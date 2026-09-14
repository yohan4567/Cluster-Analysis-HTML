#!/usr/bin/env python3
"""
Omega Workbench — local server (import-only mode)
=======================================================================
This server does NOT run any membership-filtering pipeline. You filter
stars yourself (in VS Code, or however you like) and produce a CSV.
This server's only job is to:

  1. Look for that CSV in its own folder, by the cluster's compact ID
     appearing in the filename (e.g. "ngc5139_final.csv" matches
     cluster "NGC 5139").
  2. Check the file actually has the columns the app needs.
  3. Save a canonical copy + a small stats sidecar, so next time you
     pick this cluster it loads instantly — no re-uploading.

FILE CONVENTIONS  (all files live in this folder, next to server.py)
  candidate file  : any *.csv whose name contains the cluster's compact
                     ID, e.g. "ngc5139_final.csv" for "NGC 5139"
  saved copy      : "{CompactID}_filtered.csv"
  stats sidecar   : "{CompactID}_filtered.meta.json"

ENDPOINTS
  GET  /find_data?id=NGC5139
       -> { filtered: {...} | null, candidate: {...} | null }
  POST /save_filtered
       body: { cluster_id, ra, dec, source: {type:'file', filename}
                                           | {type:'upload', csv} }
       -> validates columns, saves to disk, returns summary stats
  GET  /load_filtered?id=NGC5139
       -> instantly re-loads a previously saved file, no re-validation
"""

import http.server
import socketserver
import json
import io
import os
import re
import sys
import datetime
import webbrowser
import threading
import traceback
from urllib.parse import urlsplit, parse_qs

HOST = "127.0.0.1"          # localhost only — safe, not exposed to the network
PORT = 8000
HTML_FILE = "cluster_analysis.html"
HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, "data")   # per-cluster outputs live under here

# Columns the app actually needs to work. REQUIRED ones block import if
# missing (nothing downstream would function without them). RECOMMENDED
# ones are used by specific views (CMD, identification) — missing them
# just narrows what you can look at later, so it's a warning, not a wall.
REQUIRED_COLUMNS = ["ra", "dec", "pmra", "pmdec", "membership_prob", "is_member"]
RECOMMENDED_COLUMNS = ["source_id", "parallax", "phot_g_mean_mag", "bp_rp", "ruwe"]


# ----------------------------------------------------------------------
# filename helpers
# ----------------------------------------------------------------------
def sanitize_id(s):
    """Keep only safe characters; this becomes part of a filename."""
    return re.sub(r"[^A-Za-z0-9_-]", "", s or "")[:80]


def safe_filename(fn):
    """Strip any path components — never allow escaping HERE."""
    return os.path.basename(fn or "")


def cluster_dir(cid):
    """The per-cluster folder: data/NGC5139/ . Created on demand.
    All of a cluster's saved outputs (the filtered CSV, its stats sidecar,
    and any future graph data) live together in here."""
    d = os.path.join(DATA_DIR, cid)
    os.makedirs(d, exist_ok=True)
    return d


def filtered_path(cid):
    return os.path.join(cluster_dir(cid), "filtered.csv")

def catalog_path(cid, variant="dr3"):
    if cid == "NGC5139" and variant == "dr3_fpr":
        return os.path.join(cluster_dir(cid), "filtered_dr3_fpr.csv")
    if cid == "NGC5139" and variant == "dr3":
        p = os.path.join(cluster_dir(cid), "filtered_dr3.csv")
        return p if os.path.exists(p) else filtered_path(cid)
    return filtered_path(cid)

def catalog_variants(cid):
    if cid != "NGC5139": return []
    out = []
    for key, label in (("dr3", "Gaia DR3"), ("dr3_fpr", "Gaia DR3 + Gaia FPR")):
        p = catalog_path(cid, key)
        if os.path.exists(p): out.append({"id": key, "label": label, "filename": os.path.basename(p), "rows": count_rows(p)})
    return out


def meta_path(cid):
    return os.path.join(cluster_dir(cid), "filtered.meta.json")


def read_header(path):
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        first = f.readline()
    return [c.strip() for c in first.strip().split(",")]


def header_from_text(text):
    first = text.split("\n", 1)[0]
    return [c.strip() for c in first.strip().rstrip("\r").split(",")]


def validate_columns(header):
    missing_required = [c for c in REQUIRED_COLUMNS if c not in header]
    missing_recommended = [c for c in RECOMMENDED_COLUMNS if c not in header]
    return missing_required, missing_recommended


def count_rows(path):
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        return max(0, sum(1 for _ in f) - 1)


def find_candidate(cid):
    """Best-guess file in this folder that looks like it belongs to this
    cluster, excluding files we ourselves produced (*_filtered.csv)."""
    cid_l = cid.lower()
    if not cid_l:
        return None
    matches = []
    for fn in os.listdir(HERE):
        low = fn.lower()
        if not low.endswith(".csv") or low.endswith("_filtered.csv"):
            continue
        if cid_l in low.replace(" ", "").replace("_", ""):
            matches.append((os.path.getmtime(os.path.join(HERE, fn)), fn))
    if not matches:
        return None
    matches.sort(reverse=True)
    fn = matches[0][1]
    full = os.path.join(HERE, fn)
    try:
        rows = count_rows(full)
        header = read_header(full)
        missing_required, missing_recommended = validate_columns(header)
    except Exception as e:
        return {"filename": fn, "rows": None, "valid": False,
                "missing_required": ["(could not read file: " + str(e) + ")"],
                "missing_recommended": []}
    return {
        "filename": fn, "rows": rows,
        "valid": len(missing_required) == 0,
        "missing_required": missing_required,
        "missing_recommended": missing_recommended,
    }


def summarize(csv_text):
    """Compute the small stats block for a validated, already-filtered CSV."""
    import pandas as pd
    df = pd.read_csv(io.StringIO(csv_text), low_memory=False)
    # is_member may come in as real booleans or as "True"/"False" text
    if df["is_member"].dtype != bool:
        df["is_member"] = (
            df["is_member"].astype(str).str.strip().str.lower()
            .isin(["true", "1", "yes"])
        )
    n_kept = int(len(df))
    n_members = int(df["is_member"].sum())
    mem = df[df["is_member"]]

    def num(x):
        x = float(x)
        return None if (x != x) else x  # NaN check without importing math

    return dict(
        n_kept=n_kept,
        n_members=n_members,
        member_frac=(n_members / n_kept) if n_kept else 0.0,
        pmra=num(mem["pmra"].mean()) if len(mem) else None,
        pmdec=num(mem["pmdec"].mean()) if len(mem) else None,
    )


def summarize_path(path):
    """Same stats as summarize(), but read straight from a file on disk."""
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        return summarize(f.read())


# ----------------------------------------------------------------------
# automatic plot generation (server-side, matplotlib)
# ----------------------------------------------------------------------
# Whenever a cluster's filtered.csv is discovered, loaded, or saved, the
# server renders its standard plots (2D sky map, CMD, 3D map) as PNGs into
# the same folder. Generation runs on a background thread so responses
# stay fast, and a sidecar (plots.meta.json) remembers which version of
# the CSV the plots came from — so they regenerate only when the data
# actually changes.
PLOT_LOCK = threading.Lock()

INK_DIM = "#93a0b8"
STAR_COLOR = "#fff3d6"


def plots_meta_path(cid):
    return os.path.join(cluster_dir(cid), "plots.meta.json")


def plots_stale(cid):
    fp = filtered_path(cid)
    if not os.path.exists(fp):
        return False
    try:
        with open(plots_meta_path(cid), encoding="utf-8") as f:
            m = json.load(f)
    except Exception:
        return True
    return m.get("source_mtime") != os.path.getmtime(fp)


def ensure_plots_async(cid):
    if plots_stale(cid):
        threading.Thread(target=generate_plots, args=(cid,), daemon=True).start()


def _dark_axes(fig, ax):
    fig.patch.set_facecolor("black")
    ax.set_facecolor("black")
    for spine in getattr(ax, "spines", {}).values():
        spine.set_color("#1c2540")
    ax.tick_params(colors=INK_DIM, labelsize=8)
    ax.xaxis.label.set_color(INK_DIM)
    ax.yaxis.label.set_color(INK_DIM)
    if hasattr(ax, "zaxis"):
        ax.zaxis.label.set_color(INK_DIM)


def generate_plots(cid):
    """Render sky_map_2d.png, cmd.png (if colour data exists) and
    map_3d.png into the cluster's folder. Safe to call repeatedly."""
    with PLOT_LOCK:
        if not plots_stale(cid):
            return
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            import numpy as np
            import pandas as pd

            fp = filtered_path(cid)
            df = pd.read_csv(fp)
            if df["is_member"].dtype != bool:
                df["is_member"] = (df["is_member"].astype(str).str.strip()
                                   .str.lower().isin(["true", "1", "yes"]))
            mem = df[df["is_member"]]
            if not len(mem):
                return

            # cluster center: stored meta if available, else member median
            meta = load_or_build_meta(cid) or {}
            ra0 = meta.get("cluster_ra") or float(mem["ra"].median())
            dec0 = meta.get("cluster_dec") or float(mem["dec"].median())
            cos0 = np.cos(np.radians(dec0))
            dx = (mem["ra"].to_numpy() - ra0) * cos0 * 60.0    # arcmin
            dy = (mem["dec"].to_numpy() - dec0) * 60.0

            out = cluster_dir(cid)

            # ---- 1. 2D sky map --------------------------------------
            fig, ax = plt.subplots(figsize=(7, 7), dpi=170)
            _dark_axes(fig, ax)
            ax.scatter(dx, dy, s=1.2, c=STAR_COLOR, alpha=0.16, linewidths=0)
            ax.set_xlabel("Δ RA · cos(dec)  (arcmin)")
            ax.set_ylabel("Δ Dec  (arcmin)")
            ax.invert_xaxis()          # RA increases to the left on the sky
            ax.set_aspect("equal")
            ax.set_title(f"{cid} — {len(mem):,} member stars", color=INK_DIM, fontsize=10)
            fig.savefig(os.path.join(out, "sky_map_2d.png"),
                        facecolor="black", bbox_inches="tight")
            plt.close(fig)

            # ---- 2. CMD (needs colour + magnitude columns) ----------
            if "bp_rp" in mem.columns and "phot_g_mean_mag" in mem.columns:
                sel = mem.dropna(subset=["bp_rp", "phot_g_mean_mag"])
                if len(sel):
                    fig, ax = plt.subplots(figsize=(5.6, 7), dpi=170)
                    _dark_axes(fig, ax)
                    ax.scatter(sel["bp_rp"], sel["phot_g_mean_mag"],
                               s=1.6, c=STAR_COLOR, alpha=0.35, linewidths=0)
                    ax.set_xlabel("BP − RP  (colour)")
                    ax.set_ylabel("Gaia G magnitude")
                    ax.invert_yaxis()  # brighter is up
                    ax.set_title(f"{cid} — colour–magnitude diagram", color=INK_DIM, fontsize=10)
                    fig.savefig(os.path.join(out, "cmd.png"),
                                facecolor="black", bbox_inches="tight")
                    plt.close(fig)

            # ---- 3. 3D map ------------------------------------------
            # Per-star depth is NOT measurable (Gaia parallax errors dwarf
            # the cluster's true depth), so the Z axis is *statistically
            # reconstructed*: onion-peel the observed projected profile into
            # a 3D density nu(r), then draw each star's z from
            # p(z|R) ~ nu(sqrt(R^2+z^2)) at its own projected radius R.
            # (The earlier shuffled-offsets depth had the right one-axis
            # spread but the wrong 3D shape: disk+pole / box artifacts.)
            R = np.hypot(dx, dy)
            Rmax = float(R.max()) * 1.001 + 1e-9
            NG = 48
            edges3 = np.linspace(0, Rmax, NG + 1)
            counts3, _ = np.histogram(R, bins=edges3)
            sigma_p = counts3 / (np.pi * (edges3[1:]**2 - edges3[:-1]**2))
            sigma_p = np.convolve(sigma_p, [0.25, 0.5, 0.25], mode="same")
            # completeness correction: crowding deletes stars from dense
            # centers, carving a fake central hole -> hollow-shell artifact.
            # True profiles never rise outward: enforce monotone envelope.
            sigma_p = np.maximum.accumulate(sigma_p[::-1])[::-1]
            Rmid3 = 0.5 * (edges3[1:] + edges3[:-1])

            def _chord(j, kk):
                lo = max(edges3[kk], Rmid3[j]); hi = max(edges3[kk + 1], Rmid3[j])
                return 2 * (np.sqrt(max(hi * hi - Rmid3[j]**2, 0.0))
                            - np.sqrt(max(lo * lo - Rmid3[j]**2, 0.0)))

            nu = np.zeros(NG)
            for j in range(NG - 1, -1, -1):        # peel outside-in
                s = sigma_p[j] - sum(nu[kk] * _chord(j, kk) for kk in range(j + 1, NG))
                L = _chord(j, j)
                nu[j] = max(s / L, 0.0) if L > 0 else 0.0

            rng = np.random.default_rng(sum(ord(c) for c in cid) * 7919)  # stable seed
            NZ = 48
            dz = np.zeros(len(R))
            rbin = np.clip((R / Rmax * NG).astype(int), 0, NG - 1)
            for b in np.unique(rbin):
                sel = rbin == b
                zmax = np.sqrt(max(Rmax * Rmax - Rmid3[b]**2, 0.0))
                if zmax <= 0:
                    continue
                zg = (np.arange(NZ) + 0.5) / NZ * zmax
                w = nu[np.clip((np.hypot(Rmid3[b], zg) / Rmax * NG).astype(int), 0, NG - 1)]
                cw = np.cumsum(w)
                if cw[-1] <= 0:
                    continue
                pick = np.searchsorted(cw, rng.random(sel.sum()) * cw[-1]).clip(0, NZ - 1)
                dz[sel] = np.where(rng.random(sel.sum()) < 0.5, -1, 1) * zg[pick]
            # physical scale, if a mean parallax is usable (pc per arcmin)
            unit = "arcmin"
            sx, sy, sz = dx, dy, dz
            # physical scale: prefer the literature distance (stored by the
            # V/sigma pipeline in rotation.json) - the old median-parallax
            # estimate was noise-biased ~3x low for distant clusters
            dist_pc = None
            rj = os.path.join(cluster_dir(cid), "rotation.json")
            if os.path.exists(rj):
                try:
                    dist_pc = json.load(open(rj, encoding="utf-8"))["dist_kpc"] * 1000.0
                except Exception:
                    dist_pc = None
            if dist_pc is None and "parallax" in mem.columns:
                plx = mem["parallax"].to_numpy()
                plx = plx[np.isfinite(plx) & (plx > 0)]
                if len(plx) > 50:
                    dist_pc = 1000.0 / float(np.median(plx))
            if dist_pc:
                k = dist_pc * np.pi / (180 * 60)   # pc per arcmin
                sx, sy, sz = dx * k, dy * k, dz * k
                unit = "pc"
            fig = plt.figure(figsize=(7.5, 7), dpi=150)
            ax = fig.add_subplot(projection="3d")
            _dark_axes(fig, ax)
            ax.set_box_aspect((1, 1, 1))
            ax.scatter(sx, sy, sz, s=1.0, c=STAR_COLOR, alpha=0.25, linewidths=0)
            ax.set_xlabel(f"X ({unit})")
            ax.set_ylabel(f"Y ({unit})")
            ax.set_zlabel(f"Z modeled ({unit})")
            lim = max(np.nanmax(np.abs(v)) for v in (sx, sy, sz))
            for setter in (ax.set_xlim, ax.set_ylim, ax.set_zlim):
                setter(-lim, lim)
            ax.set_title(f"{cid} — 3D view (depth modeled, spherical symmetry)",
                         color=INK_DIM, fontsize=10)
            fig.savefig(os.path.join(out, "map_3d.png"),
                        facecolor="black", bbox_inches="tight")
            plt.close(fig)

            with open(plots_meta_path(cid), "w", encoding="utf-8") as f:
                json.dump({"source_mtime": os.path.getmtime(fp),
                           "created": datetime.datetime.now().isoformat(timespec="seconds")}, f)
        except Exception:
            traceback.print_exc()


def load_or_build_meta(cid):
    """Return the stats dict for a cluster's filtered.csv, computing and
    persisting the sidecar if it is missing or incomplete. This covers the
    case where a filtered.csv was dropped into the folder by hand (no
    accompanying .meta.json), which otherwise showed up as blank counts."""
    fp, mp = filtered_path(cid), meta_path(cid)
    if not os.path.exists(fp):
        return None
    meta = {}
    if os.path.exists(mp):
        try:
            meta = json.load(open(mp, encoding="utf-8"))
        except Exception:
            meta = {}
    if meta.get("n_members") is None:
        try:
            meta.update(summarize_path(fp))
            meta.setdefault("source_raw", "placed in folder")
            meta.setdefault("created", datetime.datetime.now().isoformat(timespec="seconds"))
            with open(mp, "w", encoding="utf-8") as f:
                json.dump(meta, f)
        except Exception:
            pass  # a malformed CSV just leaves the counts blank, as before
    return meta


# ----------------------------------------------------------------------
# HTTP handler
# ----------------------------------------------------------------------
class Handler(http.server.BaseHTTPRequestHandler):
    def _send(self, code, body, ctype="application/json"):
        data = body.encode("utf-8") if isinstance(body, str) else body
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(data)))
        self.send_header("Cache-Control", "no-store")
        # Allows the app's own HTML to reach this server even if the page
        # was opened directly as a file:// document instead of through
        # http://127.0.0.1:PORT/. Harmless: this server only listens on
        # localhost, so only this machine can ever reach it.
        self.send_header("Access-Control-Allow-Origin", "*")
        self.end_headers()
        self.wfile.write(data)

    def do_OPTIONS(self):
        self.send_response(204)
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "Content-Type")
        self.send_header("Content-Length", "0")
        self.end_headers()

    def _send_json(self, obj, code=200):
        self._send(code, json.dumps(obj))

    def do_GET(self):
        parts = urlsplit(self.path)
        path = parts.path
        qs = parse_qs(parts.query)

        if path in ("/", "/index.html", "/" + HTML_FILE):
            fpath = os.path.join(HERE, HTML_FILE)
            if not os.path.exists(fpath):
                self._send(500, f"Cannot find {HTML_FILE} next to server.py", "text/plain")
                return
            with open(fpath, "rb") as f:
                self._send(200, f.read(), "text/html; charset=utf-8")
            return

        if path == "/ping":
            self._send_json({"ok": True})
            return

        if path == "/find_data":
            cid = sanitize_id((qs.get("id") or [""])[0])
            if not cid:
                self._send_json({"error": "missing id"}, 400)
                return
            resp = {"filtered": None, "candidate": None}
            variant = (qs.get("variant") or ["dr3"])[0]
            fp = catalog_path(cid, variant)
            meta = summarize_path(fp) if os.path.exists(fp) else None
            if meta is not None:
                resp["filtered"] = {
                    "filename": os.path.basename(fp),
                    "mtime": os.path.getmtime(fp),
                    **meta,
                }
                resp["variants"] = catalog_variants(cid)
                resp["selected_variant"] = variant
                if variant == "dr3": ensure_plots_async(cid)
            else:
                resp["candidate"] = find_candidate(cid)
            self._send_json(resp)
            return

        if path == "/list_clusters":
            out = []
            for name in sorted(os.listdir(DATA_DIR)):
                d = os.path.join(DATA_DIR, name)
                if not os.path.isdir(d) or not os.path.exists(os.path.join(d, "filtered.csv")):
                    continue
                # computes + persists filtered.meta.json from the CSV if it's
                # missing, instead of silently showing a blank member count
                meta = load_or_build_meta(name) or {}
                out.append({"id": name, "n_members": meta.get("n_members"),
                            "ra": meta.get("cluster_ra"), "dec": meta.get("cluster_dec"),
                            "has_rotation": os.path.exists(os.path.join(d, "rotation.json"))})
            self._send_json({"clusters": out})
            return

        if path == "/rotation":
            cid = sanitize_id((qs.get("id") or [""])[0])
            if "bin_mode" in qs or "nbins" in qs:
                try:
                    from vsigma_study.vsigma_pipeline import analyze_cluster, CLUSTERS
                    if cid not in CLUSTERS:
                        raise ValueError("Unknown cluster")
                    mode = qs.get("bin_mode", ["equal_number"])[0]
                    count = int(qs.get("nbins", ["16"])[0])
                    variant = qs.get("variant", ["dr3"])[0]
                    if variant not in ("dr3", "dr3_fpr"):
                        raise ValueError("Unknown catalog")
                    result = analyze_cluster(cid, want_rv=False, bin_mode=mode, nbins=count,
                                             catalog=catalog_path(cid, variant), save_figure=False)
                    if result is None:
                        raise ValueError("At least 50 usable members are required")
                    result["catalog_variant"] = variant
                    self._send_json(result)
                except ValueError as e:
                    self._send_json({"error": str(e)}, 400)
                except Exception as e:
                    self._send_json({"error": str(e)}, 500)
                return

            fp = os.path.join(cluster_dir(cid), "rotation.json")
            if not cid or not os.path.exists(fp):
                self._send_json({"error": "no rotation analysis for this cluster yet"}, 404)
                return
            self._send(200, open(fp, "r", encoding="utf-8").read())
            return

        if path == "/vsigma_summary":
            fp = os.path.join(DATA_DIR, "vsigma_summary.json")
            if not os.path.exists(fp):
                self._send_json({"error": "vsigma_summary.json not found - run the pipeline"}, 404)
                return
            self._send(200, open(fp, "r", encoding="utf-8").read())
            return

        if path == "/filter_study":
            fp = os.path.join(DATA_DIR, "filter_study.json")
            if not os.path.exists(fp):
                self._send_json({"error": "no filter study data - run "
                                          "filter_study/make_filter_figures.py"}, 404)
                return
            self._send(200, open(fp, "r", encoding="utf-8").read())
            return

        if path == "/filter_image":
            # the per-cluster filter-process figure (PNG)
            cid = sanitize_id((qs.get("id") or [""])[0])
            fp = os.path.join(cluster_dir(cid), "filter_process.png")
            if not cid or not os.path.exists(fp):
                self._send_json({"error": "no filter figure for this cluster"}, 404)
                return
            with open(fp, "rb") as f:
                self._send(200, f.read(), "image/png")
            return

        if path == "/tidal_study":
            fp = os.path.join(DATA_DIR, "tidal_study.json")
            if not os.path.exists(fp):
                self._send_json({"error": "no tidal study data - run tidal_study/export_app_data.py"}, 404)
                return
            self._send(200, open(fp, "r", encoding="utf-8").read())
            return

        if path == "/vsigma_literature":
            fp = os.path.join(DATA_DIR, "vsigma_literature.json")
            if not os.path.exists(fp):
                self._send_json({"error": "no literature comparison yet"}, 404)
                return
            self._send(200, open(fp, "r", encoding="utf-8").read())
            return

        if path == "/load_filtered":
            cid = sanitize_id((qs.get("id") or [""])[0])
            variant = (qs.get("variant") or ["dr3"])[0]
            fp = catalog_path(cid, variant)
            if not cid or not os.path.exists(fp):
                self._send_json({"error": f"no cached file for '{cid}'"}, 404)
                return
            csv_text = open(fp, "r", encoding="utf-8").read()
            meta = summarize_path(fp)
            meta["catalog_variant"] = variant
            try:
                import pandas as pd
                o = pd.read_csv(fp, usecols=["catalog_origin"])["catalog_origin"].value_counts()
                meta["catalog_origin_counts"] = {str(k): int(v) for k, v in o.items()}
            except Exception: pass
            if variant == "dr3": ensure_plots_async(cid)
            self._send_json({"summary": meta, "csv": csv_text, "cached": True, "variant": variant, "filename": os.path.basename(fp), "variants": catalog_variants(cid)})
            return

        self._send_json({"error": "not found"}, 404)

    def do_POST(self):
        if self.path == "/save_plot":
            self._handle_save_plot()
            return
        if self.path != "/save_filtered":
            self._send_json({"error": "not found"}, 404)
            return
        try:
            length = int(self.headers.get("Content-Length", 0))
            payload = json.loads(self.rfile.read(length))
            cid = sanitize_id(payload.get("cluster_id", ""))
            if not cid:
                self._send_json({"error": "missing cluster_id"})
                return

            source = payload.get("source", {})
            if source.get("type") == "file":
                fn = safe_filename(source.get("filename"))
                fpath = os.path.join(HERE, fn)
                if not fn or not os.path.exists(fpath):
                    self._send_json({"error": f"file not found: {fn}"})
                    return
                csv_text = open(fpath, "r", encoding="utf-8", errors="replace").read()
                source_label = fn
            else:
                csv_text = source.get("csv", "")
                source_label = "uploaded"

            header = header_from_text(csv_text)
            missing_required, missing_recommended = validate_columns(header)
            if missing_required:
                self._send_json({
                    "error": "This file is missing required columns: " + ", ".join(missing_required),
                    "missing_required": missing_required,
                    "missing_recommended": missing_recommended,
                })
                return

            try:
                summary = summarize(csv_text)
            except Exception as e:
                self._send_json({"error": f"Could not read this CSV: {e}",
                                 "trace": traceback.format_exc()})
                return
            summary["missing_recommended"] = missing_recommended

            out_path = filtered_path(cid)
            with open(out_path, "w", encoding="utf-8", newline="") as f:
                f.write(csv_text)
            meta = dict(summary)
            meta["source_raw"] = source_label
            meta["cluster_ra"] = float(payload.get("ra", 0) or 0)
            meta["cluster_dec"] = float(payload.get("dec", 0) or 0)
            meta["created"] = datetime.datetime.now().isoformat(timespec="seconds")
            with open(meta_path(cid), "w", encoding="utf-8") as f:
                json.dump(meta, f)

            ensure_plots_async(cid)
            self._send_json({"summary": summary, "csv": csv_text,
                             "saved_as": os.path.basename(out_path)})
        except Exception as e:
            self._send_json({"error": str(e), "trace": traceback.format_exc()})

    def _handle_save_plot(self):
        """Save a rendered plot image (PNG, base64) into the cluster's
        folder, next to its filtered.csv. Used by the analysis views to
        auto-store the pictures they generate."""
        import base64
        try:
            length = int(self.headers.get("Content-Length", 0))
            payload = json.loads(self.rfile.read(length))
            cid = sanitize_id(payload.get("cluster_id", ""))
            name = re.sub(r"[^A-Za-z0-9_-]", "", payload.get("name", ""))[:60]
            b64 = payload.get("png_base64", "")
            if not cid or not name or not b64:
                self._send_json({"error": "missing cluster_id, name, or image data"})
                return
            raw = base64.b64decode(b64)
            out_path = os.path.join(cluster_dir(cid), name + ".png")
            with open(out_path, "wb") as f:
                f.write(raw)
            self._send_json({"saved_as": os.path.basename(out_path),
                             "path": out_path, "bytes": len(raw)})
        except Exception as e:
            self._send_json({"error": str(e), "trace": traceback.format_exc()})

    def log_message(self, *args):
        pass  # keep the console quiet


class ReusableServer(socketserver.ThreadingTCPServer):
    allow_reuse_address = True
    daemon_threads = True


def main():
    try:
        import pandas  # noqa: F401
    except ImportError:
        print("\n  Missing required package: pandas")
        print("  Install it with:\n\n     pip install pandas\n")
        sys.exit(1)

    global PORT
    httpd = None
    for attempt_port in range(PORT, PORT + 10):
        try:
            httpd = ReusableServer((HOST, attempt_port), Handler)
            PORT = attempt_port
            break
        except OSError:
            continue
    if httpd is None:
        print(f"  Could not bind to any port in {PORT}-{PORT+9}. Close other apps and retry.")
        sys.exit(1)

    url = f"http://{HOST}:{PORT}/"
    print("\n  GC Analysis is running.")
    print(f"  Open:  {url}")
    print("  Press Ctrl+C here to stop.\n")

    threading.Timer(0.8, lambda: webbrowser.open(url)).start()
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\n  Stopped.")
        httpd.shutdown()


if __name__ == "__main__":
    main()
