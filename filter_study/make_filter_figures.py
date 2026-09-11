"""
Filter-process study: one figure per cluster showing every stage of the
membership filter, plus the probability distribution it produces.

Panels (left to right, top to bottom):
  1. MOTION   proper-motion plot coloured by p_pm      (the GMM / window step)
  2. POSITION radial density profile, all vs members   (the King-profile step)
  3. COLOUR   colour-magnitude diagram, field vs members
  4. RESULT   histogram of the final membership probability

Outputs:
  data/<Cluster>/filter_process.png     one figure per cluster
  data/filter_study.json                summary numbers for the app screen
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

APP = Path(r"C:\Users\hoche\OneDrive\Desktop\Cluster Analysis HTML")
DATA = APP / "data"
HERE = Path(__file__).parent
sys.path.insert(0, str(APP / "membership_filter"))
from membership_filter import CLUSTERS, cluster_cfg   # noqa: E402

BG, PANEL, INK, DIM = "#0a0e18", "#101827", "#e9edf6", "#93a0b8"
COOL, WARM, GRID, OK = "#6ea8ff", "#f4a862", "#1c2540", "#7fd6a6"

# folder name in data/  ->  key used in membership_filter.CLUSTERS
def lit_key(cid):
    for name in CLUSTERS:
        if name.replace(" ", "") == cid:
            return name
    return None


def style(ax, title=None):
    ax.set_facecolor(PANEL)
    for s in ax.spines.values():
        s.set_color(GRID)
    ax.tick_params(colors=DIM, labelsize=9)
    ax.xaxis.label.set_color(DIM)
    ax.yaxis.label.set_color(DIM)
    ax.grid(color=GRID, lw=0.5, alpha=0.5)
    if title:
        ax.set_title(title, color=INK, fontsize=11, pad=8)


def make_figure(cid, df, cfg, out_png):
    mem = df[df.is_member]
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 10), dpi=120)
    fig.patch.set_facecolor(BG)
    (a1, a2), (a3, a4) = axes

    # ---------- 1. MOTION: proper motion coloured by p_pm ----------
    style(a1, "STEP 3 · MOTION — how each star drifts across the sky")
    s = df.sample(min(60000, len(df)), random_state=1)
    sc = a1.scatter(s.pmra, s.pmdec, s=3, c=s.p_pm, cmap="viridis",
                    vmin=0, vmax=1, alpha=0.6, lw=0)
    pm0 = cfg.get("pm")
    if pm0:
        a1.plot(pm0[0], pm0[1], "+", color=WARM, ms=15, mew=2.5)
        a1.annotate("known cluster motion", pm0, color=WARM, fontsize=9,
                    xytext=(14, 16), textcoords="offset points",
                    arrowprops=dict(color=WARM, arrowstyle="->", lw=1.1))
    # zoom to the interesting region around the cluster clump
    if pm0:
        a1.set_xlim(pm0[0] - 12, pm0[0] + 12)
        a1.set_ylim(pm0[1] - 12, pm0[1] + 12)
    a1.set_xlabel("drift left-right,  pmRA  (mas/yr)")
    a1.set_ylabel("drift up-down,  pmDec  (mas/yr)")
    cb = fig.colorbar(sc, ax=a1, pad=0.02)
    cb.set_label("motion score  p_pm", color=DIM, fontsize=9)
    cb.ax.tick_params(colors=DIM, labelsize=8)
    cb.outline.set_edgecolor(GRID)

    # ---------- 2. POSITION: radial density profile ----------
    style(a2, "STEP 4 · POSITION — do they crowd toward the centre?")
    r = df["r_proj_arcmin"].to_numpy()
    r = r[np.isfinite(r) & (r > 0)]
    if len(r) > 50:
        edges = np.geomspace(max(r.min(), 0.2), r.max(), 18)
        area = np.pi * (edges[1:] ** 2 - edges[:-1] ** 2)
        rm = np.sqrt(edges[1:] * edges[:-1])
        call, _ = np.histogram(r, bins=edges)
        cmem, _ = np.histogram(mem["r_proj_arcmin"].dropna(), bins=edges)
        a2.plot(rm, call / area, "o-", color=COOL, ms=4, lw=1.3, label="all stars")
        a2.plot(rm, np.maximum(cmem / area, 1e-6), "s-", color=WARM, ms=4, lw=1.3,
                label="members only")
        a2.set_xscale("log")
        a2.set_yscale("log")
        a2.legend(facecolor=PANEL, labelcolor=INK, edgecolor=GRID, fontsize=9)
    a2.set_xlabel("distance from centre  (arcmin)")
    a2.set_ylabel("stars per arcmin$^2$")

    # ---------- 3. COLOUR: CMD ----------
    style(a3, "STEP 5 · COLOUR — do they share one family sequence?")
    cmd = df.dropna(subset=["bp_rp", "phot_g_mean_mag"])
    fld = cmd[~cmd.is_member]
    mmm = cmd[cmd.is_member]
    if len(fld):
        a3.scatter(fld.bp_rp.sample(min(20000, len(fld)), random_state=2),
                   fld.phot_g_mean_mag.sample(min(20000, len(fld)), random_state=2),
                   s=2, c=DIM, alpha=0.28, lw=0, label="rejected (field)")
    if len(mmm):
        k = min(30000, len(mmm))
        sub = mmm.sample(k, random_state=2)
        a3.scatter(sub.bp_rp, sub.phot_g_mean_mag, s=2, c=WARM, alpha=0.5, lw=0,
                   label="members")
    a3.invert_yaxis()
    a3.set_xlabel("colour  (BP - RP)")
    a3.set_ylabel("brightness  (G mag)")
    a3.legend(facecolor=PANEL, labelcolor=INK, edgecolor=GRID, fontsize=9, markerscale=4)

    # ---------- 4. RESULT: membership probability histogram ----------
    style(a4, "STEP 6 · RESULT — final membership probability")
    p = df["membership_prob"].dropna().to_numpy()
    bins = np.linspace(0, 1, 41)
    thr = cfg.get("threshold", 0.5)
    below = p[p < thr]
    above = p[p >= thr]
    a4.hist(below, bins=bins, color=DIM, alpha=0.75, label=f"rejected  (p < {thr:g})")
    a4.hist(above, bins=bins, color=WARM, alpha=0.85, label=f"members  (p >= {thr:g})")
    a4.axvline(thr, color=OK, lw=1.4, ls="--")
    a4.set_yscale("log")
    a4.set_xlabel("membership probability")
    a4.set_ylabel("number of stars")
    a4.legend(facecolor=PANEL, labelcolor=INK, edgecolor=GRID, fontsize=9)

    frac = len(mem) / len(df) if len(df) else 0
    method = str(df["pm_method"].dropna().iloc[0]) if df["pm_method"].notna().any() else "-"
    method_txt = {"gmm2": "2-blob GMM", "gmm3": "3-blob GMM",
                  "window": "direct motion window", "too_few": "too few stars"}.get(method, method)
    fig.suptitle(
        f"{cid}  —  filtering process     "
        f"[{cfg.get('strategy','?')} strategy · motion step: {method_txt}]\n"
        f"{len(df):,} stars examined  →  {len(mem):,} kept as members ({frac:.1%})",
        color=INK, fontsize=13, y=0.985)
    fig.tight_layout(rect=[0, 0, 1, 0.945])
    fig.savefig(out_png, facecolor=BG, bbox_inches="tight")
    plt.close(fig)

    def jnum(x):
        """JSON-safe: bare NaN is invalid JSON and breaks JSON.parse in the
        browser, so anything non-finite becomes null."""
        if x is None:
            return None
        x = float(x)
        return x if np.isfinite(x) else None

    def col_mean(name):
        if name not in df.columns or not df[name].notna().any():
            return None
        return jnum(np.nanmean(df[name]))

    return dict(
        cluster=cid,
        strategy=cfg.get("strategy", "?"),
        pm_method=method, pm_method_label=method_txt,
        threshold=jnum(thr),
        cmd_weight=jnum(cfg.get("cmd_w", 1.0)),
        n_examined=int(len(df)),
        n_members=int(len(mem)),
        member_frac=jnum(frac),
        mean_p_pm=col_mean("p_pm"),
        # all-NaN for crowding-limited cores, where the spatial test stands down
        mean_p_spatial=col_mean("p_spatial"),
        spatial_used=bool("p_spatial" in df.columns and df["p_spatial"].notna().any()),
        mean_p_cmd=col_mean("p_cmd"),
        cmd_used=bool("p_cmd" in df.columns and df["p_cmd"].notna().any()
                      and float(cfg.get("cmd_w", 1.0)) > 0),
        # how decisive the filter was: fraction of stars in the ambiguous middle
        frac_ambiguous=jnum(np.mean((p > 0.2) & (p < 0.8))) if len(p) else None,
    )


def main():
    rows = []
    for d in sorted(DATA.iterdir()):
        fp = d / "filtered.csv"
        if not d.is_dir() or not fp.exists():
            continue
        cid = d.name
        key = lit_key(cid)
        cfg = cluster_cfg(key) if key else cluster_cfg("__unknown__")
        try:
            df = pd.read_csv(fp)
            if df["is_member"].dtype != bool:
                df["is_member"] = (df["is_member"].astype(str).str.strip()
                                   .str.lower().isin(["true", "1", "yes"]))
            rows.append(make_figure(cid, df, cfg, d / "filter_process.png"))
            print(f"[{cid}] figure written  ({rows[-1]['n_members']:,}/"
                  f"{rows[-1]['n_examined']:,} members, {rows[-1]['pm_method_label']})")
        except Exception as exc:
            print(f"[{cid}] FAILED: {exc!r}")

    (DATA / "filter_study.json").write_text(json.dumps(rows), encoding="utf-8")
    (HERE / "filter_study_summary.json").write_text(json.dumps(rows, indent=1), encoding="utf-8")
    print(f"\n{len(rows)} clusters -> data/filter_study.json")


if __name__ == "__main__":
    main()
