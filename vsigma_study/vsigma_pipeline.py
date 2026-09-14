"""
V/sigma pipeline for the GC rotation study.
============================================
Computes, for every cluster with a filtered.csv in the app's data folder:

  * error-deconvolved velocity-dispersion profile sigma(R)
  * PM rotation profile  V_rot(R) = <v_tan> per radial bin
    (Bianchini et al. 2018 method -- NOT the sinusoid fit; see METHODS.md 2.2)
  * V/sigma (peak and global) with bootstrap uncertainties
  * 2D rotation map (annular sectors, mean v_tan)
  * counter-rotation / sign-flip detection            (多旋轉軸 part 1)
  * 3D spin-axis fit per radial shell from Gaia RVs    (多旋轉軸 part 2)
    with perspective-rotation correction

Outputs per cluster -> <data>/<Cluster>/rotation.json + rotation.png
Cross-cluster       -> <data>/vsigma_summary.json, vsigma_summary.csv,
                       correlations in the same JSON, plots in this folder.

Conventions (fixed throughout):
  x = dRA*cos(dec0)  [East +],  y = dDec  [North +]   (arcmin)
  v_tan = (x*v_y - y*v_x)/R    + = counterclockwise with East right, North up
  v_rad = (x*v_x + y*v_y)/R    + = outward
  z     = line of sight, away from observer
  Spin vector Omega:  <v_x> = -Oz*y, <v_y> = +Oz*x, v_z = Ox*y - Oy*x
"""
import json
import shutil
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA = Path(__file__).resolve().parents[1] / "data"
OUT = Path(__file__).resolve().parent
K = 4.74047                      # (mas/yr)*(kpc) -> km/s
P_CUT = 0.90
N_BOOT = 200
RNG = np.random.default_rng(42)

# ----------------------------------------------------------------------
# Cluster parameters from the team's target list (gc_target_list.xlsx)
# name: (ra, dec, r_c', r_h', conc, r_t', [Fe/H], d_kpc, lit_rotation)
# V_LOS (km/s, heliocentric, Harris 2010) added for perspective corrections.
# ----------------------------------------------------------------------
CLUSTERS = {
    "NGC5139":  (201.6968, -47.4796, 2.37, 5.00, 1.31, 48.389, -1.53,  5.2, "Confirmed rotator",        232.1),
    "NGC104":   (  6.0236, -72.0813, 0.36, 3.17, 2.07, 42.296, -0.72,  4.5, "Confirmed rotator",        -18.0),
    "NGC7078":  (322.4930,  12.1670, 0.14, 1.00, 2.29, 27.298, -2.37, 10.4, "Confirmed rotator",       -107.0),
    "NGC7089":  (323.3626,  -0.8233, 0.32, 1.06, 1.59, 12.449, -1.65, 11.5, "Confirmed rotator",         -5.3),
    "NGC5904":  (229.6384,   2.0810, 0.44, 1.77, 1.73, 23.629, -1.29,  7.5, "Confirmed rotator",         53.2),
    "NGC6656":  (279.0998, -23.9048, 1.33, 3.36, 1.38, 31.904, -1.70,  3.2, "Confirmed rotator",       -146.3),
    "NGC6273":  (255.6575, -26.2680, 0.43, 1.32, 1.53, 14.570, -1.74,  8.8, "Confirmed rotator",        135.0),
    "NGC5272":  (205.5484,  28.3773, 0.37, 2.31, 1.89, 28.721, -1.50, 10.2, "Borderline rotator",      -147.6),
    "NGC6752":  (287.7171, -59.9846, 0.17, 1.91, 2.50, 53.759, -1.54,  4.0, "Borderline rotator",       -26.7),
    "NGC6809":  (294.9988, -30.9648, 1.80, 2.83, 0.93, 15.320, -1.94,  5.4, "Borderline / non-rotator", 174.7),
    "NGC288":   ( 13.1885, -26.5826, 1.35, 2.23, 0.99, 13.193, -1.32,  8.9, "Non-rotator",              -45.4),
    "NGC362":   ( 15.8094, -70.8488, 0.18, 0.82, 1.76, 10.358, -1.26,  8.6, "Not tested / unclear",     223.5),
    "NGC6397":  (265.1754, -53.6743, 0.05, 2.90, 2.50, 15.811, -2.02,  2.3, "Not tested / unclear",      18.8),
    "NGC6341":  (259.2808,  43.1359, 0.26, 1.02, 1.68, 12.444, -2.31,  8.3, "Not tested / unclear",    -120.0),
    "NGC5024":  (198.2302,  18.1682, 0.35, 1.31, 1.72, 18.368, -2.10, 17.9, "Not tested / unclear",     -62.9),
    "NGC2419":  (114.5353,  38.8824, 0.32, 0.89, 1.37,  7.502, -2.15, 82.6, "Not tested / unclear",     -20.2),
    "NGC1851":  ( 78.5282, -40.0466, 0.09, 0.51, 1.86,  6.520, -1.18, 12.1, "Not tested / unclear",     320.5),
    "NGC6266":  (255.3033, -30.1137, 0.22, 0.92, 1.71, 11.283, -1.18,  6.8, "Not tested / unclear",     -70.1),
    "Terzan5":  (267.0200, -24.7792, 0.16, 0.72, 1.62,  6.670, -0.23,  6.9, "Not tested / unclear",     -82.0),
    "NGC6440":  (267.2196, -20.3603, 0.14, 0.48, 1.62,  5.836, -0.36,  8.5, "Not tested / unclear",     -76.6),
    "NGC5286":  (206.6117, -51.3743, 0.28, 0.73, 1.41,  7.197, -1.69, 11.7, "Not tested / unclear",      57.4),
    "NGC6121":  (245.8968, -26.5258, 1.16, 4.33, 1.65, 51.815, -1.16,  2.2, "Not tested / unclear",      70.7),
    "Pal5":     (229.0219,  -0.1116, 2.29, 2.73, 0.52,  7.583, -1.41, 23.2, "Not tested / unclear",     -58.7),
    "NGC6544":  (271.8358, -24.9973, 0.05, 1.21, 1.63,  2.133, -1.40,  3.0, "Not tested / unclear",     -38.0),
}

# ----------------------------------------------------------------------
# Dynamical age = Age / t_rh: how many half-mass relaxation times the
# cluster has lived through (Bianchini+2018 predict V/sigma decays with it).
# t_rh: Baumgardt & Hilker catalog, log T_rh column
#       (people.smp.uq.edu.au/HolgerBaumgardt/globular/parameter.html, 2026-07-18).
# Age:  Forbes & Bridges 2010 compilation, except Terzan5 (Ferraro+2016,
#       dominant old population) and NGC6440 (Pallanca+2021).
# name: (age_gyr, trh_gyr)
DYN_AGE = {
    "NGC104":  (13.06,  5.012), "NGC288":  (10.62,  3.090),
    "NGC362":  (10.37,  1.023), "NGC1851": ( 9.98,  1.000),
    "NGC2419": (12.30, 42.658), "NGC5024": (12.67,  8.511),
    "NGC5139": (11.52, 20.893), "NGC5272": (11.39,  3.162),
    "NGC5286": (12.54,  1.413), "NGC5904": (10.62,  3.162),
    "NGC6121": (12.54,  0.871), "NGC6266": (11.78,  1.175),
    "NGC6273": (11.90,  2.951), "NGC6341": (13.18,  1.259),
    "NGC6397": (12.67,  0.513), "NGC6440": (13.00,  1.023),
    "NGC6544": (10.37,  0.245), "NGC6656": (12.67,  3.090),
    "NGC6752": (11.78,  1.995), "NGC6809": (12.29,  3.162),
    "NGC7078": (12.93,  1.622), "NGC7089": (11.78,  2.818),
    "Pal5":    ( 9.80,  7.943), "Terzan5": (12.00,  2.692),
}


def jnum(x):
    """JSON-safe number (None for NaN/inf)."""
    x = float(x)
    return x if np.isfinite(x) else None


def load_members(cid, catalog=None):
    """Returns (members, p_cut_used). Falling back to 0.5 means the
    membership model separated poorly -> sample flagged as contaminated."""
    fp = Path(catalog) if catalog else DATA / cid / "filtered.csv"
    if not fp.exists():
        return None, None
    df = pd.read_csv(fp, low_memory=False)
    for c in ["ra", "dec", "pmra", "pmdec", "pmra_error", "pmdec_error", "membership_prob"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["ra", "dec", "pmra", "pmdec", "pmra_error", "pmdec_error"])
    cols = ["ra", "dec", "pmra", "pmdec", "pmra_error", "pmdec_error"]
    df = df[np.isfinite(df[cols]).all(axis=1) & (df["pmra_error"] > 0) & (df["pmdec_error"] > 0)]
    m = df[df["membership_prob"] >= P_CUT]
    cut = P_CUT
    if len(m) < 100:
        m = df[df["membership_prob"] >= 0.5]
        cut = 0.5
    return m.reset_index(drop=True), cut


def deconvolved_sigma(vx, vy, ex2):
    """1D dispersion from two components, minus mean measurement error^2.
    vx, vy: residual velocities (km/s); ex2: per-star mean err^2 (km/s)^2."""
    obs2 = 0.5 * (np.var(vx) + np.var(vy))
    true2 = obs2 - np.mean(ex2)
    return (np.sqrt(true2) if true2 > 0 else np.nan), obs2


def analyze_cluster(cid, want_rv=True, bin_mode="equal_number", nbins=None, catalog=None, save_figure=True):
    if bin_mode not in ("equal_number", "equal_radius"):
        raise ValueError("Unknown bin mode")
    if nbins is not None and (not isinstance(nbins, int) or not 2 <= nbins <= 100):
        raise ValueError("Bin count must be an integer from 2 to 100")
    rng = np.random.default_rng(42)
    ra0, dec0, rc, rh, conc, rt, feh, dist, litrot, vlos = CLUSTERS[cid]
    m, p_cut_used = load_members(cid, catalog)
    if m is None or len(m) < 50:
        return None
    n = len(m)
    if nbins is not None and nbins > n // 2:
        raise ValueError(f"Choose at most {n // 2} bins for {n} usable stars")
    kd = K * dist                                     # (mas/yr) -> km/s
    pc_per_arcmin = dist * 1000 * np.pi / (180 * 60)

    x = (m["ra"].to_numpy() - ra0) * np.cos(np.radians(dec0)) * 60.0   # arcmin E+
    y = (m["dec"].to_numpy() - dec0) * 60.0                            # arcmin N+
    R = np.hypot(x, y)

    # systemic PM: 3-sigma clipped mean
    pmx, pmy = m["pmra"].to_numpy(), m["pmdec"].to_numpy()
    keep = np.ones(n, bool)
    for _ in range(3):
        mx, my = pmx[keep].mean(), pmy[keep].mean()
        sx, sy = pmx[keep].std(), pmy[keep].std()
        keep = (np.abs(pmx - mx) < 3 * sx) & (np.abs(pmy - my) < 3 * sy)
    pm0 = (pmx[keep].mean(), pmy[keep].mean())

    # residual PM after removing the systemic (bulk) motion
    dpmx = pmx - pm0[0]                                # mas/yr, East
    dpmy = pmy - pm0[1]                                # mas/yr, North

    # perspective-rotation correction (van de Ven+2006; Bianchini+2018):
    # the cluster's own bulk PM, viewed across a wide field, shears the
    # residual PM field and masquerades as a tangential (rotation) signal.
    # First-order term scales as mu_sys * (angular offset in radians), so it
    # is negligible for compact/slow clusters but large for omega Cen -- the
    # biggest, closest, fastest-moving-on-sky GC -- which is exactly what
    # distorted its rotation profile. (This is the transverse-motion shear;
    # the separate line-of-sight perspective on <v_rad> is vrad_persp below.)
    x_rad = np.radians(x / 60.0)                       # East offset, radians
    y_rad = np.radians(y / 60.0)                       # North offset, radians
    dpmx = dpmx + x_rad * pm0[0]
    dpmy = dpmy + y_rad * pm0[1]

    vx = kd * dpmx                                     # km/s, East
    vy = kd * dpmy                                     # km/s, North
    err2 = (kd**2) * 0.5 * (m["pmra_error"].to_numpy()**2 + m["pmdec_error"].to_numpy()**2)

    Rsafe = np.where(R > 1e-6, R, 1e-6)
    v_tan = (x * vy - y * vx) / Rsafe
    v_rad = (x * vx + y * vy) / Rsafe

    # ---- radial bins (equal count) ----------------------------------
    nbins = nbins or int(np.clip(n // 400, 4, 16))
    order = np.argsort(R, kind="stable")
    if bin_mode == "equal_number":
        bin_idx = np.array_split(order, nbins)
        edges = [float(R.min())] + [float((R[a[-1]] + R[b[0]]) / 2)
                 for a, b in zip(bin_idx[:-1], bin_idx[1:])] + [float(R.max())]
    else:
        edges = np.linspace(0, R.max(), nbins + 1)
        labels = np.clip(np.searchsorted(edges, R, side="right") - 1, 0, nbins - 1)
        bin_idx = [np.flatnonzero(labels == i) for i in range(nbins)]
    # Resample original assignments, preserving equal-number ties and empty rings.
    labels = np.empty(n, dtype=int)
    for i, members in enumerate(bin_idx):
        labels[members] = i

    prof = dict(r_mid=[], n=[], vrot=[], vrot_err=[], sig=[], sig_err=[],
                vrad=[], vrad_err=[])
    for b in bin_idx:
        if len(b) < 2:
            for key in prof:
                prof[key].append(len(b) if key == "n" else np.nan)
            continue
        rb = R[b]
        prof["r_mid"].append(np.median(rb))
        prof["n"].append(len(b))
        prof["vrot"].append(v_tan[b].mean())
        prof["vrot_err"].append(v_tan[b].std() / np.sqrt(len(b)))
        s, obs2 = deconvolved_sigma(vx[b], vy[b], err2[b])
        prof["sig"].append(s)
        # SE of the DECONVOLVED sigma: Var(sigma_obs^2-hat) ~ sigma_obs^4/N
        # (2N samples from two components), then /(2*sigma_true) -- much
        # larger than s/sqrt(2N) when errors dominate the observed spread
        prof["sig_err"].append(obs2 / (2 * s * np.sqrt(len(b))) if s > 0 else np.nan)
        prof["vrad"].append(v_rad[b].mean())
        prof["vrad_err"].append(v_rad[b].std() / np.sqrt(len(b)))
    for k in prof:
        prof[k] = np.array(prof[k], float)

    # predicted perspective signal in <v_rad>(R): receding cluster -> apparent contraction
    prof["vrad_persp"] = -vlos * (prof["r_mid"] / 60.0) * np.pi / 180.0

    # ---- headline stats ----------------------------------------------
    elig = prof["n"] >= min(200, max(30, n / 10))
    if not elig.any():
        elig = prof["n"] >= 2
    i_pk = int(np.nanargmax(np.where(elig, np.abs(prof["vrot"]), np.nan)))
    v_peak = abs(prof["vrot"][i_pk])
    r_peak = prof["r_mid"][i_pk]

    inner = R <= max(rh, np.quantile(R, 0.25))
    sigma0, _ = deconvolved_sigma(vx[inner], vy[inner], err2[inner])
    sig_glob, _ = deconvolved_sigma(vx, vy, err2)
    vsig_peak = v_peak / sigma0 if sigma0 > 0 else np.nan
    vsig_glob = abs(v_tan.mean()) / sig_glob if sig_glob > 0 else np.nan

    # sinusoid diagnostic (colleague's fit == residual systemic-PM dipole;
    # also a sensitive contamination indicator: clean clusters give < ~0.5 km/s)
    th = np.arctan2(x, y)
    A = np.column_stack([np.sin(th), np.cos(th)])
    (a_s, b_s), *_ = np.linalg.lstsq(A, v_tan, rcond=None)
    dipole = float(np.hypot(a_s, b_s))

    # rotation-detection significance: chi^2 of the bin means against zero.
    # Guards against the max-statistic bias of V_peak - a non-rotator can
    # still show V_peak ~ 2x bin error just by picking the noisiest bin.
    from scipy.stats import chi2 as chi2_dist
    ok_b = prof["vrot_err"] > 0
    chisq = float(np.sum((prof["vrot"][ok_b] / prof["vrot_err"][ok_b])**2))
    p_rot = float(chi2_dist.sf(chisq, int(ok_b.sum()))) if ok_b.sum() else 1.0

    # ---- measurement-quality cascade -----------------------------------
    # unmeasurable   sigma deconvolves to zero (errors >> intrinsic spread)
    # contaminated   membership fallback cut used, OR dispersion profile
    #                RISES outward (field stars), OR dipole > 1 km/s
    #                (all three are empirical contamination signatures;
    #                clean clusters in this sample: dipole 0.07-0.6 km/s)
    # error-dominated  deconvolved variance < 5x its own sampling noise
    med_err_vel = float(np.sqrt(np.median(err2)))
    n_in = max(int(inner.sum()), 1)
    obs2_in = 0.5 * (np.var(vx[inner]) + np.var(vy[inner]))
    snr_sigma = (sigma0**2) / (obs2_in * np.sqrt(2.0 / n_in)) if obs2_in > 0 else 0.0
    third = max(len(prof["sig"]) // 3, 1)
    sig_in3, sig_out3 = np.median(prof["sig"][:third]), np.median(prof["sig"][-third:])
    err_io = np.hypot(np.median(prof["sig_err"][:third]), np.median(prof["sig_err"][-third:]))
    rising = bool(np.isfinite(err_io) and sig_out3 > sig_in3 + 2 * err_io)
    if not np.isfinite(sigma0) or sigma0 <= 0:
        quality = "unmeasurable"
    elif p_cut_used < P_CUT or rising or dipole > 1.0:
        quality = "contaminated"
    elif snr_sigma < 5:
        quality = "error-dominated"
    else:
        quality = "good"

    # ---- bootstrap V/sigma uncertainty --------------------------------
    boots = []
    edges_arr = np.array(edges)
    for _ in range(N_BOOT):
        s = rng.integers(0, n, n)
        vt_b, vx_b, vy_b, e2_b, R_b = v_tan[s], vx[s], vy[s], err2[s], R[s]
        which = labels[s]
        vb = np.full(nbins, np.nan)
        for i in range(nbins):
            sel = which == i
            if sel.sum() > 5:
                vb[i] = vt_b[sel].mean()
        vpk = np.nanmax(np.abs(np.where(elig, vb, np.nan)))
        inn = R_b <= max(rh, np.quantile(R, 0.25))
        s0, _ = deconvolved_sigma(vx_b[inn], vy_b[inn], e2_b[inn])
        if s0 > 0 and np.isfinite(vpk):
            boots.append(vpk / s0)
    vsig_err = float(np.std(boots)) if len(boots) > 20 else np.nan
    v_peak_err = float(prof["vrot_err"][i_pk])

    # ---- counter-rotation / sign structure (多旋轉軸 part 1) ----------
    sig_bins = np.abs(prof["vrot"]) > 2 * prof["vrot_err"]
    signs = np.sign(prof["vrot"]) * sig_bins
    counter = bool((signs > 0).any() and (signs < 0).any())

    # ---- rotation map (annular sectors) -------------------------------
    map_nr = min(nbins, 8)
    r_edges = (np.linspace(0, R.max(), map_nr + 1) if bin_mode == "equal_radius"
               else np.quantile(R, np.linspace(0, 1, map_nr + 1)))
    t_edges = np.linspace(-np.pi, np.pi, 13)
    vmap = np.full((map_nr, 12), np.nan)
    ri = np.clip(np.searchsorted(r_edges, R, side="right") - 1, 0, map_nr - 1)
    # sector angle = position angle East of North (atan2(E, N)) so the polar
    # displays (matplotlib zero at N; HTML rotation 90/clockwise) are honest
    ti = np.clip(np.searchsorted(t_edges, np.arctan2(x, y), side="right") - 1, 0, 11)
    for i in range(map_nr):
        for j in range(12):
            sel = (ri == i) & (ti == j)
            if sel.sum() >= 10:
                vmap[i, j] = v_tan[sel].mean()

    res = dict(
        cluster=cid, n_members=int(n), dist_kpc=dist, feh=feh, conc=conc,
        r_c_arcmin=rc, r_h_arcmin=rh, r_t_arcmin=rt, lit_rotation=litrot,
        vlos_lit=vlos, pm_sys=dict(pmra=jnum(pm0[0]), pmdec=jnum(pm0[1])),
        diagnostics=dict(
            estimator="moment subtraction; not a selection-corrected likelihood fit",
            kinematic_catalog=("Gaia DR3 + FPR" if "catalog_origin" in m and (m["catalog_origin"] == "gaia_fpr").any() else "Gaia DR3"),
            origin_counts=({str(k): int(v) for k, v in m["catalog_origin"].value_counts().items()}
                           if "catalog_origin" in m else {"gaia_dr3": n}),
            unresolved_bins=int(np.sum(~np.isfinite(prof["sig"]))),
            warning="Exploratory estimate: per-star reported error variances are averaged within each bin. Membership selection can bias dispersion; FPR scores are not calibrated DR3 probabilities. Unresolved values are not zero.",
        ),
        binning=dict(mode=bin_mode, count=nbins, bootstrap_samples=N_BOOT,
                     sparse_bins=int((prof["n"] < 30).sum()), map_rings=map_nr),
        bins=dict(
            r_edges_arcmin=[jnum(v) for v in edges],
            r_mid_arcmin=[jnum(v) for v in prof["r_mid"]],
            r_mid_pc=[jnum(v * pc_per_arcmin) for v in prof["r_mid"]],
            n=[int(v) for v in prof["n"]],
            v_rot=[jnum(v) for v in prof["vrot"]],
            v_rot_err=[jnum(v) for v in prof["vrot_err"]],
            sigma=[jnum(v) for v in prof["sig"]],
            sigma_err=[jnum(v) for v in prof["sig_err"]],
            v_rad=[jnum(v) for v in prof["vrad"]],
            v_rad_err=[jnum(v) for v in prof["vrad_err"]],
            v_rad_perspective=[jnum(v) for v in prof["vrad_persp"]],
        ),
        map=dict(r_edges_arcmin=[jnum(v) for v in r_edges],
                 theta_edges_deg=[jnum(np.degrees(v)) for v in t_edges],
                 v_tan_mean=[[jnum(v) for v in row] for row in vmap]),
        stats=dict(v_peak_kms=jnum(v_peak), v_peak_err_kms=jnum(v_peak_err),
                   r_peak_arcmin=jnum(r_peak), sigma0_kms=jnum(sigma0),
                   vsig_peak=jnum(vsig_peak) if quality == "good" else None,
                   vsig_peak_err=(jnum(vsig_err) if quality == "good" and jnum(vsig_peak) is not None else None),
                   median_err_vel_kms=jnum(med_err_vel),
                   sigma0_snr=jnum(snr_sigma),
                   p_cut_used=p_cut_used,
                   rotation_detection_p=jnum(p_rot),
                   rotation_detected=bool(p_rot < 0.01),
                   quality=quality,
                   dipole_diag_kms=jnum(dipole),
                   counter_rotation=counter,
                   rotation_sense="CCW (E-right convention)" if prof["vrot"][i_pk] > 0 else "CW (E-right convention)"),
    )

    # ---- 3D spin axis from RVs (多旋轉軸 part 2) ----------------------
    res["axis3d"] = dict(available=False, reason="no RV data")
    if want_rv:
        try:
            res["axis3d"] = fit_axis3d(cid, m, x, y, R, v_tan, pm0, dist,
                                       pc_per_arcmin, rh)
        except Exception as e:
            res["axis3d"] = dict(available=False, reason=f"error: {e}")
    # a twist between shells is only meaningful if rotation itself is
    # detected - otherwise it is the angle between two noise vectors
    if res["axis3d"].get("twist_significant") and not res["stats"]["rotation_detected"]:
        res["axis3d"]["twist_significant"] = False
        res["axis3d"]["twist_note"] = "rotation not detected; twist is noise"

    if save_figure:
        make_cluster_figure(cid, prof, r_edges, t_edges, vmap, res, rh)
    return res


# ----------------------------------------------------------------------
# Gaia RVs: query (cached), crossmatch, 3D spin-vector fit per shell
# ----------------------------------------------------------------------
def get_rv_table(cid):
    cache = DATA / cid / "gaia_rv.csv"
    if cache.exists():
        return pd.read_csv(cache)
    ra0, dec0 = CLUSTERS[cid][0], CLUSTERS[cid][1]
    rt = CLUSTERS[cid][5]
    from astroquery.gaia import Gaia
    Gaia.ROW_LIMIT = -1
    q = f"""SELECT source_id, ra, dec, radial_velocity, radial_velocity_error
    FROM gaiadr3.gaia_source
    WHERE 1=CONTAINS(POINT('ICRS', ra, dec),
                     CIRCLE('ICRS', {ra0}, {dec0}, {min(rt / 60.0, 1.0)}))
      AND radial_velocity IS NOT NULL"""
    last = None
    for attempt in range(2):
        try:
            df = Gaia.launch_job_async(q).get_results().to_pandas()
            df.to_csv(cache, index=False)
            return df
        except Exception as e:
            last = e
            time.sleep(20)
    raise last


def _fit_omega_xy(xp, yp, rv):
    """Fit rv = c + Ox*yp - Oy*xp  (xp, yp in pc; rv km/s) with 3-sigma clip."""
    keep = np.ones(len(rv), bool)
    for _ in range(3):
        A = np.column_stack([np.ones(keep.sum()), yp[keep], -xp[keep]])
        coef, *_ = np.linalg.lstsq(A, rv[keep], rcond=None)
        resid = rv - (coef[0] + coef[1] * yp - coef[2] * xp)
        s = resid[keep].std()
        keep = np.abs(resid) < 3 * max(s, 1e-9)
    # final fit on clipped sample + parameter errors from residual scatter
    A = np.column_stack([np.ones(keep.sum()), yp[keep], -xp[keep]])
    coef, *_ = np.linalg.lstsq(A, rv[keep], rcond=None)
    resid = rv[keep] - A @ coef
    cov = np.linalg.inv(A.T @ A) * resid.var()
    return coef, np.sqrt(np.diag(cov)), int(keep.sum())


def fit_axis3d(cid, members, x, y, R, v_tan, pm0, dist, pc_per_arcmin, rh):
    rv = get_rv_table(cid)
    if "source_id" not in members.columns:
        return dict(available=False, reason="no source_id column")
    merged = members.merge(rv, on="source_id", suffixes=("", "_rv"))
    merged = merged.dropna(subset=["radial_velocity"])
    if len(merged) < 50:
        return dict(available=False, reason=f"only {len(merged)} RV members")

    ra0, dec0 = CLUSTERS[cid][0], CLUSTERS[cid][1]
    vlos = CLUSTERS[cid][9]
    xm = (merged["ra"].to_numpy() - ra0) * np.cos(np.radians(dec0)) * 60.0
    ym = (merged["dec"].to_numpy() - dec0) * 60.0
    xp, yp = xm * pc_per_arcmin, ym * pc_per_arcmin
    rvv = merged["radial_velocity"].to_numpy().astype(float)

    # perspective correction: apparent RV gradient from systemic transverse motion
    vtE, vtN = K * dist * pm0[0], K * dist * pm0[1]          # km/s
    rv_corr = rvv - (vtE * np.radians(xm / 60.0) + vtN * np.radians(ym / 60.0))

    # Omega_z from PM rotation curve, solid-body fit through origin (inner half)
    Rp = R * pc_per_arcmin
    inner = Rp <= np.quantile(Rp, 0.5)
    Oz = float((v_tan[inner] * Rp[inner]).sum() / (Rp[inner]**2).sum())

    def one_fit(sel_x, sel_y, sel_rv, oz):
        coef, errs, nk = _fit_omega_xy(sel_x, sel_y, sel_rv)
        Ox, Oy = float(coef[1]), float(coef[2])
        Om = np.array([Ox, Oy, oz])
        mag = float(np.linalg.norm(Om))
        return dict(omega_kms_per_pc=[jnum(Ox), jnum(Oy), jnum(oz)],
                    omega_err=[jnum(errs[1]), jnum(errs[2]), None],
                    n_used=nk, mag=jnum(mag),
                    incl_deg=jnum(np.degrees(np.arccos(np.clip(oz / mag, -1, 1)))) if mag > 0 else None,
                    pa_deg=jnum(np.degrees(np.arctan2(Ox, Oy)) % 360.0))

    out = dict(available=True, n_rv=int(len(merged)))
    out["global"] = one_fit(xp, yp, rv_corr, Oz)

    # radial shells -> axis twist (多旋轉軸 test)
    if len(merged) >= 150:
        Rm = np.hypot(xp, yp)
        med = np.median(Rm)
        sh = []
        for name, sel in (("inner", Rm <= med), ("outer", Rm > med)):
            # shell-specific Omega_z from PM curve in that radial range
            selR = (Rp <= med) if name == "inner" else (Rp > med)
            oz_s = float((v_tan[selR] * Rp[selR]).sum() / max((Rp[selR]**2).sum(), 1e-9))
            f = one_fit(xp[sel], yp[sel], rv_corr[sel], oz_s)
            f["shell"] = name
            f["r_range_pc"] = [jnum(Rm[sel].min()), jnum(Rm[sel].max())]
            sh.append(f)
        out["shells"] = sh
        a = np.array(sh[0]["omega_kms_per_pc"], float)
        b = np.array(sh[1]["omega_kms_per_pc"], float)
        na, nb = np.linalg.norm(a), np.linalg.norm(b)
        if na > 0 and nb > 0:
            out["twist_deg"] = jnum(np.degrees(np.arccos(np.clip(a @ b / (na * nb), -1, 1))))
            # Monte-Carlo error on the twist: perturb each shell's spin vector
            # by its fit errors (Oz err approximated at 25%) - a large twist
            # between two barely-detected vectors is NOT a detection.
            ea = np.array([sh[0]["omega_err"][0] or 0, sh[0]["omega_err"][1] or 0, 0.25 * abs(a[2])])
            eb = np.array([sh[1]["omega_err"][0] or 0, sh[1]["omega_err"][1] or 0, 0.25 * abs(b[2])])
            angs = []
            for _ in range(1000):
                pa_v = a + RNG.normal(0, 1, 3) * ea
                pb_v = b + RNG.normal(0, 1, 3) * eb
                npa, npb = np.linalg.norm(pa_v), np.linalg.norm(pb_v)
                if npa > 0 and npb > 0:
                    angs.append(np.degrees(np.arccos(np.clip(pa_v @ pb_v / (npa * npb), -1, 1))))
            out["twist_err_deg"] = jnum(np.std(angs))
            out["twist_significant"] = bool(out["twist_deg"] is not None and
                                            out["twist_err_deg"] is not None and
                                            out["twist_deg"] > 2 * out["twist_err_deg"])
    return out


# ----------------------------------------------------------------------
# figures
# ----------------------------------------------------------------------
INK = "#93a0b8"

def _dark(fig, axes):
    fig.patch.set_facecolor("black")
    for ax in np.atleast_1d(axes).ravel():
        ax.set_facecolor("#05070d")
        for sp in ax.spines.values():
            sp.set_color("#1c2540")
        ax.tick_params(colors=INK, labelsize=8)
        ax.xaxis.label.set_color(INK)
        ax.yaxis.label.set_color(INK)
        ax.title.set_color(INK)
        ax.grid(color="#1c2540", lw=0.5, alpha=0.6)


def make_cluster_figure(cid, prof, r_edges, t_edges, vmap, res, rh):
    fig = plt.figure(figsize=(12, 10), dpi=140)
    ax1 = fig.add_subplot(2, 2, 1)
    ax2 = fig.add_subplot(2, 2, 2)
    ax3 = fig.add_subplot(2, 2, 3)
    ax4 = fig.add_subplot(2, 2, 4, projection="polar")
    _dark(fig, [ax1, ax2, ax3])
    ax4.set_facecolor("#05070d")
    ax4.tick_params(colors=INK, labelsize=7)
    ax4.grid(color="#1c2540", lw=0.5, alpha=0.6)

    r = prof["r_mid"]
    kd_fig = K * res["dist_kpc"]                      # km/s -> mas/yr divisor
    ax1.axhline(0, color="#5f6b86", lw=1)
    ax1.errorbar(r, np.abs(prof["vrot"]) / kd_fig, prof["vrot_err"] / kd_fig,
                 fmt="o-", ms=4, lw=1.2,
                 color="#6ea8ff", ecolor="#6ea8ff", capsize=2)
    ax1.axvline(rh, color="#f4a862", lw=0.8, ls="--")
    ax1.set_xlabel("R (arcmin)"); ax1.set_ylabel("|mu_t|  (mas/yr)")
    ax1.set_title("Rotation profile (mean tangential PM)")

    ax2.errorbar(r, prof["sig"], prof["sig_err"], fmt="s-", ms=4, lw=1.2,
                 color="#7fd6a6", ecolor="#7fd6a6", capsize=2)
    ax2.axvline(rh, color="#f4a862", lw=0.8, ls="--", label="r_h")
    ax2.set_xlabel("R (arcmin)"); ax2.set_ylabel("sigma (km/s, error-deconvolved)")
    ax2.set_title("Velocity-dispersion profile")
    ax2.legend(facecolor="black", labelcolor=INK, edgecolor="#1c2540", fontsize=8)

    vs = np.abs(prof["vrot"]) / np.where(prof["sig"] > 0, prof["sig"], np.nan)
    vs_err = vs * np.sqrt((prof["vrot_err"] / np.where(np.abs(prof["vrot"]) > 0, np.abs(prof["vrot"]), np.nan))**2
                          + (prof["sig_err"] / np.where(prof["sig"] > 0, prof["sig"], np.nan))**2)
    ax3.errorbar(r, vs, vs_err, fmt="^-", ms=4, lw=1.2, color="#ffc890",
                 ecolor="#ffc890", capsize=2)
    ax3.set_xlabel("R (arcmin)"); ax3.set_ylabel("|V_rot| / sigma")
    st = res["stats"]
    fmt3 = lambda v: "n/a" if v is None else f"{v:.3f}"
    ax3.set_title(f"V/sigma profile   (peak {fmt3(st['vsig_peak'])} ± {fmt3(st['vsig_peak_err'])})"
                  + ("" if st["quality"] == "good" else f"  [{st['quality']}]"))

    T, Rr = np.meshgrid(t_edges, r_edges)
    vmax = np.nanmax(np.abs(vmap)) if np.isfinite(vmap).any() else 1
    pc = ax4.pcolormesh(T, Rr, vmap, cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax4.set_theta_zero_location("N")   # position angle: 0 = North...
    ax4.set_theta_direction(-1)        # ...increasing toward East (clockwise here)
    ax4.set_title("Rotation map: mean v_tan (km/s), angle = PA E of N", color=INK, fontsize=9)
    cb = fig.colorbar(pc, ax=ax4, pad=0.1, shrink=0.8)
    cb.ax.tick_params(colors=INK, labelsize=7)

    fig.suptitle(f"{cid} — {res['n_members']:,} members — V/sigma(peak) = {fmt3(st['vsig_peak'])}"
                 + f" — literature: {res['lit_rotation']}",
                 color=INK, fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out1 = DATA / cid / "rotation.png"
    fig.savefig(out1, facecolor="black", bbox_inches="tight")
    fig.savefig(OUT / "per_cluster_plots" / f"{cid}_rotation.png", facecolor="black", bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------
# cross-cluster summary + correlations
# ----------------------------------------------------------------------
def correlations(rows):
    from scipy.stats import spearmanr
    df = pd.DataFrame(rows)
    df["r_h_pc"] = df["r_h_arcmin"] * df["dist_kpc"] * 1000 * np.pi / (180 * 60)
    df["log_n"] = np.log10(df["n_members"])
    params = [("feh", "[Fe/H]"), ("conc", "concentration c"),
              ("r_h_pc", "half-light radius (pc)"), ("sigma0", "sigma0 (km/s)"),
              ("log_n", "log10 N members"), ("dist_kpc", "distance (kpc)"),
              ("dyn_age", "dynamical age (Age/t_rh)")]
    corr = []
    reliable = df[df["quality"] == "good"]
    for key, label in params:
        ok = df[[key, "vsig_peak"]].dropna()
        rho, p = spearmanr(ok[key], ok["vsig_peak"])
        okr = reliable[[key, "vsig_peak"]].dropna()
        rho_r, p_r = (spearmanr(okr[key], okr["vsig_peak"]) if len(okr) > 4
                      else (np.nan, np.nan))
        corr.append(dict(param=key, label=label, spearman_rho=jnum(rho),
                         p_value=jnum(p), n=len(ok),
                         spearman_rho_reliable=jnum(rho_r),
                         p_value_reliable=jnum(p_r), n_reliable=len(okr)))
    return df, corr


def correlation_figure(df, corr):
    ncol = 3
    nrow = int(np.ceil(len(corr) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(15, 4.5 * nrow), dpi=140)
    _dark(fig, axes)
    for ax in np.atleast_1d(axes).ravel()[len(corr):]:
        ax.set_visible(False)
    color_map = {"Confirmed rotator": "#6ea8ff", "Borderline rotator": "#ffc890",
                 "Borderline / non-rotator": "#ffc890", "Non-rotator": "#ff7a7a",
                 "Not tested / unclear": "#93a0b8"}
    keys = [c["param"] for c in corr]
    labels = {c["param"]: c["label"] for c in corr}
    stats = {c["param"]: c for c in corr}
    for ax, key in zip(np.atleast_1d(axes).ravel(), keys):
        for _, row in df.iterrows():
            if row["vsig_peak"] is None or not np.isfinite(row["vsig_peak"]):
                continue
            ax.errorbar(row[key], row["vsig_peak"], yerr=row["vsig_peak_err"] or 0,
                        fmt="o", ms=5, color=color_map.get(row["lit_rotation"], "#93a0b8"),
                        ecolor="#5f6b86", capsize=2)
            ax.annotate(row["cluster"].replace("NGC", "N"), (row[key], row["vsig_peak"]),
                        fontsize=6, color="#5f6b86", xytext=(3, 3),
                        textcoords="offset points")
        c = stats[key]
        ax.set_xlabel(labels[key]); ax.set_ylabel("V/sigma (peak)")
        ax.set_title(f"Spearman rho = {c['spearman_rho']:.2f}  (p = {c['p_value']:.3f})")
    fig.suptitle("V/sigma vs cluster parameters — colored by literature rotation status "
                 "(blue=confirmed, orange=borderline, red=non, grey=untested)",
                 color=INK, fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(OUT / "vsigma_correlations.png", facecolor="black", bbox_inches="tight")
    plt.close(fig)


def ranking_figure(df):
    d = df.dropna(subset=["vsig_peak"]).sort_values("vsig_peak", ascending=True)
    fig, ax = plt.subplots(figsize=(9, 10), dpi=140)
    _dark(fig, ax)
    color_map = {"Confirmed rotator": "#6ea8ff", "Borderline rotator": "#ffc890",
                 "Borderline / non-rotator": "#ffc890", "Non-rotator": "#ff7a7a",
                 "Not tested / unclear": "#93a0b8"}
    colors = [color_map.get(s, "#93a0b8") for s in d["lit_rotation"]]
    ax.barh(d["cluster"], d["vsig_peak"], xerr=d["vsig_peak_err"].fillna(0),
            color=colors, ecolor="#5f6b86", capsize=2)
    ax.set_xlabel("V/sigma (peak)")
    ax.set_title("Clusters ranked by V/sigma — blue = literature-confirmed rotators",
                 color=INK)
    fig.tight_layout()
    fig.savefig(OUT / "vsigma_ranking.png", facecolor="black", bbox_inches="tight")
    plt.close(fig)


def multi_axis_report(results):
    lines = ["# 多旋轉軸 (Multiple / Nested Rotation Axes) — Findings", ""]
    lines += [
        "Two testable versions of the idea (see METHODS.md §5):",
        "1. **Sign/amplitude structure of V_rot(R)** — does the line-of-sight spin",
        "   component change with radius (e.g. a counter-rotating core)?",
        "2. **3D spin-axis tilt between radial shells** — with Gaia radial velocities,",
        "   fit the full spin vector per shell; the angle between shells is the",
        "   'nested axis' tilt.", ""]
    lines.append("## Per-cluster results\n")
    lines.append("| Cluster | rotation detected (PM) | counter-rotation | N_RV | 3D axis incl (deg) | axis PA (deg) | shell twist ± err (deg) | twist significant? |")
    lines.append("|---|---|---|---|---|---|---|---|")
    for r in results:
        ax = r.get("axis3d", {})
        g = ax.get("global", {})
        tw, twe = ax.get("twist_deg"), ax.get("twist_err_deg")
        lines.append("| {} | {} | {} | {} | {} | {} | {} | {} |".format(
            r["cluster"],
            "YES" if r["stats"]["rotation_detected"] else "no",
            "YES" if r["stats"]["counter_rotation"] else "no",
            ax.get("n_rv", "—") if ax.get("available") else "—",
            f"{g.get('incl_deg'):.0f}" if g.get("incl_deg") is not None else "—",
            f"{g.get('pa_deg'):.0f}" if g.get("pa_deg") is not None else "—",
            (f"{tw:.0f} ± {twe:.0f}" if tw is not None and twe is not None else "—"),
            ("**YES**" if ax.get("twist_significant") else "no") if tw is not None else "—"))
    lines += ["",
        "**How to read this:** a large shell twist with small errors means the inner and",
        "outer parts of the cluster do NOT share one rotation axis — the honest version",
        "of the team's nested-axes idea. A 'YES' in counter-rotation means the",
        "line-of-sight spin actually reverses sign with radius.", "",
        "**Caveats:** RV samples are small (only bright stars have Gaia RVs); the",
        "perspective-rotation correction has been applied (it would otherwise fake a",
        "rotation gradient of up to ~1.6 km/s for nearby fast-moving clusters).",
        "Stars do not literally orbit 'sub-poles' (the moon-earth analogy) — relaxation",
        "destroys such hierarchies — but a radius-dependent spin vector is real physics",
        "and is exactly what these two tests measure."]
    (OUT / "multi_axis_report.md").write_text("\n".join(lines), encoding="utf-8")


# ----------------------------------------------------------------------
def main():
    (OUT / "per_cluster_plots").mkdir(parents=True, exist_ok=True)
    results, failures = [], []
    for cid in CLUSTERS:
        if not (DATA / cid / "filtered.csv").exists():
            print(f"[{cid}] no filtered.csv — skipped")
            continue
        t0 = time.time()
        try:
            r = analyze_cluster(cid)
            if r is None:
                failures.append((cid, "too few members"))
                continue
            (DATA / cid / "rotation.json").write_text(json.dumps(r), encoding="utf-8")
            results.append(r)
            st = r["stats"]
            ax = r.get("axis3d", {})
            vs_txt = f"{st['vsig_peak']:.3f}" if st["vsig_peak"] is not None else "n/a"
            print(f"[{cid}] N={r['n_members']:,}  V/sig={vs_txt}"
                  f"±{st['vsig_peak_err'] or 0:.3f}  q={st['quality']}"
                  f"  det_p={st['rotation_detection_p']:.1e}  sigma0={st['sigma0_kms']:.2f}"
                  f"  dipole={st['dipole_diag_kms']:.3f}  counter={st['counter_rotation']}"
                  f"  RV={'n=%d' % ax.get('n_rv', 0) if ax.get('available') else ax.get('reason', '')}"
                  f"  [{time.time()-t0:.1f}s]  lit={r['lit_rotation']}")
        except Exception as e:
            failures.append((cid, repr(e)))
            traceback.print_exc()

    rows = [dict(cluster=r["cluster"], n_members=r["n_members"],
                 dist_kpc=r["dist_kpc"], feh=r["feh"], conc=r["conc"],
                 r_h_arcmin=r["r_h_arcmin"], lit_rotation=r["lit_rotation"],
                 age_gyr=DYN_AGE.get(r["cluster"], (None, None))[0],
                 trh_gyr=DYN_AGE.get(r["cluster"], (None, None))[1],
                 dyn_age=(round(DYN_AGE[r["cluster"]][0] / DYN_AGE[r["cluster"]][1], 3)
                          if r["cluster"] in DYN_AGE else None),
                 v_peak=r["stats"]["v_peak_kms"], sigma0=r["stats"]["sigma0_kms"],
                 vsig_peak=r["stats"]["vsig_peak"], vsig_peak_err=r["stats"]["vsig_peak_err"],
                 quality=r["stats"]["quality"],
                 rotation_detected=r["stats"]["rotation_detected"],
                 rotation_detection_p=r["stats"]["rotation_detection_p"],
                 counter_rotation=r["stats"]["counter_rotation"],
                 axis_incl_deg=(r["axis3d"].get("global", {}) or {}).get("incl_deg"),
                 axis_twist_deg=r["axis3d"].get("twist_deg"))
            for r in results]
    df, corr = correlations(rows)
    df.to_csv(OUT / "vsigma_summary.csv", index=False)
    summary = dict(generated="pipeline", clusters=rows, correlations=corr)
    (DATA / "vsigma_summary.json").write_text(json.dumps(summary), encoding="utf-8")
    (OUT / "vsigma_summary.json").write_text(json.dumps(summary, indent=1), encoding="utf-8")
    correlation_figure(df, corr)
    ranking_figure(df)
    multi_axis_report(results)

    print(f"\n{len(results)} clusters done, {len(failures)} failures: {failures}")
    print("Correlations (Spearman):")
    for c in corr:
        print(f"  V/sigma vs {c['label']:26s} rho={c['spearman_rho']:+.2f}  p={c['p_value']:.3f}")


if __name__ == "__main__":
    main()
