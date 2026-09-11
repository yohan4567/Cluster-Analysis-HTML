"""Export the tidal-stripping study into one compact JSON the app can plot
interactively (data/tidal_study.json). Reuses the cached wide-field queries
and the same selection logic as tidal_pipeline.py."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter

import tidal_pipeline as tp

HERE = Path(__file__).parent
OUT = Path(r"C:\Users\hoche\OneDrive\Desktop\Cluster Analysis HTML\data\tidal_study.json")

results = {r["cluster"]: r for r in json.loads((HERE / "tidal_results.json").read_text(encoding="utf-8"))}

payload = []
for cid in tp.TARGETS:
    ra0, dec0, rt_arcmin, dist, radius, pmwin = tp.TARGETS[cid]
    rt_deg = rt_arcmin / 60.0
    wide = tp.wide_query(cid)                      # cached
    pm0 = tp.cluster_pm(cid)
    mids, med, tol, members = tp.cmd_ridge(cid)

    dpm = np.hypot(wide["pmra"] - pm0[0], wide["pmdec"] - pm0[1])
    ridge_col = np.interp(wide["phot_g_mean_mag"], mids, med)
    ridge_tol = np.interp(wide["phot_g_mean_mag"], mids, tol)
    cand = wide[(np.abs(wide["bp_rp"] - ridge_col) < ridge_tol) & (dpm < 1.5)].copy()

    x = ((cand["ra"] - ra0) * np.cos(np.radians(dec0))).to_numpy()
    y = (cand["dec"] - dec0).to_numpy()

    # same density/significance grid as the pipeline, to tag in-tail stars
    B = 60
    H, xe, ye = np.histogram2d(x, y, bins=B, range=[[-radius, radius], [-radius, radius]])
    S = gaussian_filter(H.astype(float), 2.0)
    xc = 0.5 * (xe[:-1] + xe[1:]); yc = 0.5 * (ye[:-1] + ye[1:])
    XX, YY = np.meshgrid(xc, yc, indexing="ij")
    RR = np.hypot(XX, YY)
    edge = (RR > 0.75 * radius) & (RR < radius)
    sig = (S - np.median(S[edge])) / max(np.std(S[edge]), 1e-9)
    mask = (sig > 2.0) & (RR > 1.5 * rt_deg) & (RR < 0.75 * radius)
    xi = np.clip(np.searchsorted(xe, x, side="right") - 1, 0, B - 1)
    yi = np.clip(np.searchsorted(ye, y, side="right") - 1, 0, B - 1)
    in_tail = mask[xi, yi]

    # a light CMD sample of cluster members for the background layer
    msamp = members.sample(min(4000, len(members)), random_state=1)

    r = results.get(cid, {})
    payload.append(dict(
        cluster=cid, dist_kpc=dist, rt_deg=round(rt_deg, 4),
        radius_deg=radius, stats=r,
        stars=dict(x=[round(v, 4) for v in x], y=[round(v, 4) for v in y],
                   g=[round(v, 3) for v in cand["phot_g_mean_mag"]],
                   bprp=[round(v, 3) for v in cand["bp_rp"]],
                   in_tail=[bool(v) for v in in_tail]),
        ridge=dict(g=[round(v, 2) for v in mids], color=[round(v, 3) for v in med]),
        members_cmd=dict(g=[round(v, 3) for v in msamp["phot_g_mean_mag"]],
                         bprp=[round(v, 3) for v in msamp["bp_rp"]]),
    ))

OUT.write_text(json.dumps(payload), encoding="utf-8")
print(f"wrote {OUT}  ({OUT.stat().st_size/1e6:.1f} MB, {len(payload)} clusters)")
