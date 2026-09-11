"""
Tidal-stripping search around uncontaminated globular clusters.
================================================================
Idea: stars stripped by the Milky Way's tides leave the cluster but keep
(nearly) the cluster's proper motion and its colour-magnitude sequence.
So we query a WIDE field (degrees, far beyond the tidal radius), keep only
stars that move like the cluster AND lie on the cluster's CMD ridge line,
and then look at where the survivors sit on the sky. A stream of them in
a line = tidal tails. Their CMD turnoff matching the cluster's = same age
and metallicity = proof they were born in the cluster.

Targets are deliberately high-latitude, contamination-free clusters with
literature-reported tails to validate against:
  Pal 5    - THE textbook tidal-tail cluster (Odenkirchen+2001, Ibata+2019)
  NGC 6341 (M92) - stream reported by Thomas et al. 2020
  NGC 288  - suspected short tails (Jordi & Grebel 2010)

Outputs in this folder: <cid>_tail_map.png, <cid>_tail_cmd.png,
tidal_results.json, plus cached wide-field queries in data/.
"""
import json
import socket
import time
import traceback
from pathlib import Path

# A dead Gaia connection must RAISE, not hang forever (the previous run sat
# 9 hours on a zombie socket with no timeout).
socket.setdefaulttimeout(600)

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

HERE = Path(__file__).parent
CACHE = HERE / "gaia_cache"
APPDATA = Path(r"C:\Users\hoche\OneDrive\Desktop\Cluster Analysis HTML\data")
K = 4.74047

# cid: (ra, dec, r_t arcmin, dist kpc, search radius deg, pm window mas/yr)
TARGETS = {
    "Pal5":    (229.0219,  -0.1116,  7.583, 23.2, 4.0, 2.0),
    "NGC6341": (259.2808,  43.1359, 12.444,  8.3, 3.5, 2.0),
    "NGC288":  ( 13.1885, -26.5826, 13.193,  8.9, 3.0, 2.0),
}

INK = "#93a0b8"


def cluster_pm(cid):
    """Systemic PM from our own high-confidence members."""
    df = pd.read_csv(APPDATA / cid / "filtered.csv",
                     usecols=["pmra", "pmdec", "membership_prob"])
    m = df[df["membership_prob"] >= 0.9]
    if len(m) < 100:
        m = df[df["membership_prob"] >= 0.5]
    return float(m["pmra"].mean()), float(m["pmdec"].mean())


def wide_query(cid):
    """Wide-field Gaia cone with PM / parallax / quality prefilters baked
    into the ADQL so the download stays small. Cached to disk."""
    out = CACHE / f"{cid}_wide.csv"
    if out.exists():
        print(f"[{cid}] cached wide query: {out.name}")
        return pd.read_csv(out)
    ra0, dec0, rt, dist, radius, pmwin = TARGETS[cid]
    pm0 = cluster_pm(cid)
    from astroquery.gaia import Gaia
    Gaia.ROW_LIMIT = -1

    def build(r_out, r_in):
        # ONLY the cheap, highly-selective cuts live in ADQL (PM box kills
        # ~95% of rows; G limit is indexed). ruwe / parallax / bp_rp cuts
        # moved client-side - they made the server-side query so expensive
        # that jobs sat queued for hours.
        ring = (f"""
      AND 0=CONTAINS(POINT('ICRS', ra, dec),
                     CIRCLE('ICRS', {ra0}, {dec0}, {r_in}))""" if r_in > 0 else "")
        return f"""SELECT source_id, ra, dec, pmra, pmdec, pmra_error, pmdec_error,
        parallax, parallax_error, phot_g_mean_mag, bp_rp, ruwe
    FROM gaiadr3.gaia_source
    WHERE 1=CONTAINS(POINT('ICRS', ra, dec),
                     CIRCLE('ICRS', {ra0}, {dec0}, {r_out})){ring}
      AND pmra BETWEEN {pm0[0] - pmwin} AND {pm0[0] + pmwin}
      AND pmdec BETWEEN {pm0[1] - pmwin} AND {pm0[1] + pmwin}
      AND phot_g_mean_mag < 20.5"""

    # equal-area annuli keep each download small and restartable
    fracs = np.sqrt(np.linspace(0.0, 1.0, 9)) * radius
    print(f"[{cid}] querying Gaia in 8 annuli to {radius} deg, "
          f"PM box +-{pmwin} around ({pm0[0]:.2f},{pm0[1]:.2f})")
    chunks = []
    for i, (r_in, r_out) in enumerate(zip(fracs[:-1], fracs[1:])):
        part = CACHE / f"{cid}_ann{i}.csv"
        if part.exists():                       # survives crashes/restarts
            df = pd.read_csv(part)
            print(f"[{cid}]   annulus {r_in:.2f}-{r_out:.2f}: {len(df):,} stars (cached)")
            chunks.append(df)
            continue
        last = None
        for attempt in range(5):
            try:
                t0 = time.time()
                df = Gaia.launch_job_async(build(r_out, r_in)).get_results().to_pandas()
                print(f"[{cid}]   annulus {r_in:.2f}-{r_out:.2f}: {len(df):,} stars "
                      f"in {time.time()-t0:.0f}s")
                df.to_csv(part, index=False)
                chunks.append(df)
                last = None
                break
            except Exception as e:
                last = e
                print(f"[{cid}]   annulus {r_in:.2f}-{r_out:.2f} attempt {attempt+1} failed: {e!r}")
                time.sleep(60 * (attempt + 1))   # 503s = server maintenance; back off hard
        if last is not None:
            raise last
    df = pd.concat(chunks, ignore_index=True).drop_duplicates(subset="source_id")
    # quality cuts applied client-side (cheap here, expensive in ADQL)
    df = df[(df["ruwe"] < 1.4) & df["bp_rp"].notna()
            & df["parallax"].notna() & (df["parallax"] < 0.6) & (df["parallax"] > -1)]
    df.to_csv(out, index=False)
    print(f"[{cid}] total {len(df):,} stars cached after client-side quality cuts")
    return df


def cmd_ridge(cid):
    """Median colour vs magnitude ridge of the cluster's own members,
    with a per-bin tolerance. Used as the matched filter."""
    df = pd.read_csv(APPDATA / cid / "filtered.csv",
                     usecols=["phot_g_mean_mag", "bp_rp", "membership_prob"])
    m = df[(df["membership_prob"] >= 0.9)].dropna()
    if len(m) < 200:
        m = df[(df["membership_prob"] >= 0.5)].dropna()
    edges = np.arange(m["phot_g_mean_mag"].min(), 20.6, 0.4)
    mids, med, tol = [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        sl = m[(m["phot_g_mean_mag"] >= lo) & (m["phot_g_mean_mag"] < hi)]["bp_rp"]
        if len(sl) < 15:
            continue
        mids.append(0.5 * (lo + hi))
        med.append(sl.median())
        tol.append(max(2.0 * 1.4826 * np.median(np.abs(sl - sl.median())), 0.06))
    return np.array(mids), np.array(med), np.array(tol), m


def turnoff(mids, med):
    """Turnoff = bluest point of the ridge, searched only in the FAINT 40%
    of the ridge - otherwise a blue horizontal branch (e.g. M92's) wins the
    "bluest point" and the result isn't the turnoff at all."""
    lo = np.quantile(mids, 0.6)
    sel = mids >= lo
    i = int(np.argmin(np.where(sel, med, np.inf)))
    return float(mids[i]), float(med[i])


def analyze(cid):
    ra0, dec0, rt_arcmin, dist, radius, pmwin = TARGETS[cid]
    rt_deg = rt_arcmin / 60.0
    wide = wide_query(cid)
    pm0 = cluster_pm(cid)
    mids, med, tol, members = cmd_ridge(cid)

    # candidate selection: tight PM circle + CMD matched filter
    dpm = np.hypot(wide["pmra"] - pm0[0], wide["pmdec"] - pm0[1])
    ridge_col = np.interp(wide["phot_g_mean_mag"], mids, med)
    ridge_tol = np.interp(wide["phot_g_mean_mag"], mids, tol)
    in_cmd = np.abs(wide["bp_rp"] - ridge_col) < ridge_tol
    in_pm = dpm < 1.5
    cand = wide[in_cmd & in_pm].copy()

    x = (cand["ra"] - ra0) * np.cos(np.radians(dec0))     # deg, East+
    y = cand["dec"] - dec0
    r = np.hypot(x, y)
    cand["x"], cand["y"], cand["r"] = x, y, r
    outside = cand[cand["r"] > 1.5 * rt_deg]

    # density map (smoothed), background from the field edge
    B = 60
    H, xe, ye = np.histogram2d(cand["x"], cand["y"],
                               bins=B, range=[[-radius, radius], [-radius, radius]])
    S = gaussian_filter(H.astype(float), 2.0)
    xc = 0.5 * (xe[:-1] + xe[1:]); yc = 0.5 * (ye[:-1] + ye[1:])
    XX, YY = np.meshgrid(xc, yc, indexing="ij")
    RR = np.hypot(XX, YY)
    edge = (RR > 0.75 * radius) & (RR < radius)
    bg, bg_sd = float(np.median(S[edge])), float(np.std(S[edge]))
    sig = (S - bg) / max(bg_sd, 1e-9)

    # tail direction: PA of the 2nd-moment principal axis of significant
    # cells OUTSIDE the tidal radius (the cluster core is masked out)
    mask = (sig > 2.0) & (RR > 1.5 * rt_deg) & (RR < 0.75 * radius)
    n_sig_cells = int(mask.sum())
    pa_tail = None
    if n_sig_cells >= 4:
        w = np.clip(sig[mask], 0, None)
        mx, my = XX[mask], YY[mask]
        cxx = np.sum(w * mx * mx); cyy = np.sum(w * my * my); cxy = np.sum(w * mx * my)
        # principal axis angle in (E,N); convert to PA East of North
        ang = 0.5 * np.arctan2(2 * cxy, cxx - cyy)
        pa_tail = float(np.degrees(np.arctan2(np.sin(ang), np.cos(ang))))
        pa_tail = (90.0 - pa_tail) % 180.0          # to PA (E of N), axis so mod 180

    # direction of cluster motion on the sky, PA East of North (mod 180
    # for comparison with the tail AXIS, tails are leading+trailing)
    pa_pm = float(np.degrees(np.arctan2(pm0[0], pm0[1])) % 180.0)

    # membership excess outside the tidal radius vs background expectation
    area_out = np.pi * ((0.75 * radius)**2 - (1.5 * rt_deg)**2)
    cell_area = (2 * radius / B)**2
    n_out = int(((outside["r"] < 0.75 * radius)).sum())
    n_bg_expect = bg / cell_area / (gaussian_filter(np.ones((B, B)), 2.0).mean()) * area_out
    # simpler robust estimate: median smoothed density is per-cell counts
    n_bg_expect = float(bg * area_out / cell_area)

    # CMD turnoff comparison: cluster members vs stars INSIDE the detected
    # tail overdensity (comparing against the whole outer field would dilute
    # the sample with background stars and shift the turnoff colour)
    to_g, to_col = turnoff(mids, med)
    xi = np.clip(np.searchsorted(xe, outside["x"], side="right") - 1, 0, B - 1)
    yi = np.clip(np.searchsorted(ye, outside["y"], side="right") - 1, 0, B - 1)
    in_tail_cell = mask[xi, yi] if n_sig_cells else np.zeros(len(outside), bool)
    tail_stars = outside[in_tail_cell]
    tail_near_to = tail_stars[(np.abs(tail_stars["phot_g_mean_mag"] - to_g) < 1.0)]
    to_col_tail = float(tail_near_to["bp_rp"].median()) if len(tail_near_to) >= 10 else None

    # all-magnitudes CMD consistency: normalized colour offset from the
    # cluster ridge, in units of the local tolerance. True siblings center
    # near 0; pure background spreads uniformly across [-1, 1].
    # (More robust than the turnoff for distant clusters like Pal 5, whose
    # real turnoff sits at Gaia's magnitude floor.)
    cmd_off = None
    if len(tail_stars) >= 30:
        offs = ((tail_stars["bp_rp"] - np.interp(tail_stars["phot_g_mean_mag"], mids, med))
                / np.interp(tail_stars["phot_g_mean_mag"], mids, tol))
        cmd_off = float(np.median(offs))

    res = dict(cluster=cid, dist_kpc=dist, r_t_deg=round(rt_deg, 3),
               search_radius_deg=radius, n_wide=int(len(wide)),
               n_candidates=int(len(cand)), n_outside_rt=int(len(outside)),
               n_outside_expected_bg=round(n_bg_expect, 1),
               excess_outside=round(len(outside[outside["r"] < 0.75 * radius]) - n_bg_expect, 1),
               n_sig_cells=n_sig_cells,
               tail_pa_deg=None if pa_tail is None else round(pa_tail, 1),
               pm_pa_deg=round(pa_pm, 1),
               pa_agreement_deg=None if pa_tail is None else
                   round(min(abs(pa_tail - pa_pm), 180 - abs(pa_tail - pa_pm)), 1),
               turnoff_G=round(to_g, 2), turnoff_color_cluster=round(to_col, 3),
               turnoff_color_tail=to_col_tail if to_col_tail is None else round(to_col_tail, 3),
               cmd_offset_tail_norm=cmd_off if cmd_off is None else round(cmd_off, 3),
               n_tail_stars=int(len(tail_stars)))

    # ---- figures -----------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 7.2), dpi=150)
    fig.patch.set_facecolor("black"); ax.set_facecolor("#05070d")
    pc = ax.pcolormesh(xe, ye, sig.T, cmap="magma", vmin=-1, vmax=max(4, sig.max() * 0.8))
    cb = fig.colorbar(pc, ax=ax); cb.set_label("significance (sigma)", color=INK)
    cb.ax.tick_params(colors=INK)
    th = np.linspace(0, 2 * np.pi, 100)
    ax.plot(1.5 * rt_deg * np.cos(th), 1.5 * rt_deg * np.sin(th), color="#6ea8ff", lw=1.2, ls="--")
    L = radius * 0.55
    ax.plot([L * np.sin(np.radians(pa_pm)), -L * np.sin(np.radians(pa_pm))],
            [L * np.cos(np.radians(pa_pm)), -L * np.cos(np.radians(pa_pm))],
            color="#7fd6a6", lw=1.4, ls=":", label=f"PM direction PA={pa_pm:.0f}°")
    if pa_tail is not None:
        ax.plot([L * np.sin(np.radians(pa_tail)), -L * np.sin(np.radians(pa_tail))],
                [L * np.cos(np.radians(pa_tail)), -L * np.cos(np.radians(pa_tail))],
                color="#ffc890", lw=1.6, label=f"tail axis PA={pa_tail:.0f}°")
    ax.invert_xaxis()      # sky convention: East to the left
    ax.set_xlabel("ΔRA·cos(dec) (deg)  [E left]"); ax.set_ylabel("ΔDec (deg)")
    for sp in ax.spines.values(): sp.set_color("#1c2540")
    ax.tick_params(colors=INK); ax.xaxis.label.set_color(INK); ax.yaxis.label.set_color(INK)
    ax.legend(facecolor="black", labelcolor=INK, edgecolor="#1c2540", fontsize=8)
    ax.set_title(f"{cid}: PM+CMD-selected star density beyond the cluster "
                 f"(dashed circle = 1.5 tidal radii)", color=INK, fontsize=10)
    fig.tight_layout(); fig.savefig(HERE / f"{cid}_tail_map.png", facecolor="black",
                                    bbox_inches="tight"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 7), dpi=150)
    fig.patch.set_facecolor("black"); ax.set_facecolor("#05070d")
    cmd_sample = tail_stars if len(tail_stars) >= 30 else outside
    cmd_label = (f"stars in detected tail (n={len(tail_stars)})" if len(tail_stars) >= 30
                 else f"candidates beyond 1.5 r_t (n={len(outside)})")
    ax.scatter(members["bp_rp"], members["phot_g_mean_mag"], s=2, c="#5f6b86",
               alpha=0.25, linewidths=0, label="cluster members")
    ax.scatter(cmd_sample["bp_rp"], cmd_sample["phot_g_mean_mag"], s=8, c="#ffc890",
               alpha=0.85, linewidths=0, label=cmd_label)
    ax.plot(med, mids, color="#6ea8ff", lw=1.2, label="CMD ridge")
    ax.axhline(to_g, color="#7fd6a6", lw=0.9, ls="--", label=f"turnoff G={to_g:.1f}")
    ax.invert_yaxis(); ax.set_xlabel("BP − RP"); ax.set_ylabel("G")
    for sp in ax.spines.values(): sp.set_color("#1c2540")
    ax.tick_params(colors=INK); ax.xaxis.label.set_color(INK); ax.yaxis.label.set_color(INK)
    ax.legend(facecolor="black", labelcolor=INK, edgecolor="#1c2540", fontsize=8)
    ax.set_title(f"{cid}: tail candidates share the cluster CMD", color=INK, fontsize=10)
    fig.tight_layout(); fig.savefig(HERE / f"{cid}_tail_cmd.png", facecolor="black",
                                    bbox_inches="tight"); plt.close(fig)
    return res


def main():
    CACHE.mkdir(exist_ok=True)
    # stage 1: all queries first, so every download is cached even if a
    # later analysis step hits a bug
    for cid in TARGETS:
        try:
            wide_query(cid)
        except Exception:
            traceback.print_exc()
    # stage 2: analysis
    results = []
    for cid in TARGETS:
        try:
            r = analyze(cid)
            results.append(r)
            print(f"[{cid}] candidates={r['n_candidates']:,} outside_rt={r['n_outside_rt']}"
                  f" (bg expect {r['n_outside_expected_bg']})  tail_PA={r['tail_pa_deg']}"
                  f" PM_PA={r['pm_pa_deg']} agree={r['pa_agreement_deg']}deg"
                  f"  TO_col cluster={r['turnoff_color_cluster']} tail={r['turnoff_color_tail']}")
        except Exception:
            traceback.print_exc()
    (HERE / "tidal_results.json").write_text(json.dumps(results, indent=1), encoding="utf-8")
    print(f"\n{len(results)} clusters analyzed -> tidal_results.json")


if __name__ == "__main__":
    main()
