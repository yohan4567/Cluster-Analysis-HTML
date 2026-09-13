# =====================================================================
# ADAPTIVE GAIA DR3 GLOBULAR CLUSTER MEMBERSHIP FILTER (v2)
# =====================================================================
# Each cluster is filtered with a strategy suited to its environment:
#
#   clean     - high Galactic latitude, simple field -> anchored 2-comp GMM
#   crowded   - bulge/disk field -> 3-comp GMM + hard PM gate, CMD off/down
#   companion - SMC in the field (47 Tuc, NGC 362) -> 3-comp GMM
#   sparse    - too few stars for a GMM (Pal 5) -> direct PM window
#   distant   - PM errors dominate (NGC 2419) -> error-aware PM window
#
# Every cluster is anchored to its literature mean proper motion
# (Vasiliev & Baumgardt 2021), so the fit can never lock onto the
# wrong clump. If the fitted cluster component lands too far from the
# literature value, the code falls back to the direct window method.
#
# Reuses the raw Gaia CSVs cached in gc_filter_outputs/ - only queries
# Gaia for clusters whose raw CSV is missing (e.g. NGC 6121).
#
# Usage:
#   python gc_gaia_dr3_adaptive_filter.py
#   python gc_gaia_dr3_adaptive_filter.py --only "NGC 6121" "Terzan 5"
# =====================================================================

import argparse
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent

# ======================= GLOBAL SETTINGS =======================

# the target list now lives next to this script (self-contained folder)
TARGET_LIST_DEFAULT = str(Path(__file__).parent / "gc_target_list.xlsx")
RAW_CACHE_FOLDERS = [SCRIPT_DIR / "gc_filter_outputs"]   # where v1 raw CSVs live
OUTPUT_FOLDER_DEFAULT = SCRIPT_DIR / "gc_filter_outputs_v2"

RUWE_MAX = 1.4
PM_ERR_MAX = 1.5
PARALLAX_NSIGMA = 3.0
GAIA_G_MAG_MAX = 20.5
GMM_N_INIT = 3
CMD_TRAIN_PROB = 0.50

# How far (mas/yr) a fitted GMM component may sit from the literature PM
# before we distrust the fit and fall back to the window method.
ANCHOR_TOLERANCE = 1.5


# ======================= PER-CLUSTER CONFIG =======================
# pm       : literature mean proper motion (pmra, pmdec) in mas/yr
#            from Vasiliev & Baumgardt (2021)
# pm_disp  : rough central PM dispersion of the cluster, mas/yr
# strategy : which filtering method to use
# cmd_w    : weight of the colour-magnitude test, 0 = off (heavy
#            differential reddening), 1 = full
# Optional overrides: n_comp, pm_k (gate size in sigma), spatial_w,
#                     threshold

STRATEGY_DEFAULTS = {
    "clean":     dict(method="gmm",    n_comp=2, pm_k=4.0, cmd_w=1.0, spatial_w=1.0, threshold=0.5),
    "crowded":   dict(method="gmm",    n_comp=3, pm_k=3.0, cmd_w=0.0, spatial_w=1.0, threshold=0.5),
    "companion": dict(method="gmm",    n_comp=3, pm_k=4.0, cmd_w=1.0, spatial_w=1.0, threshold=0.5),
    "sparse":    dict(method="window",           pm_k=3.0, cmd_w=0.5, spatial_w=0.3, threshold=0.3),
    # distant: PM errors dominate the window score, so typical true members
    # only reach p ~ 0.3-0.5 -> lower acceptance threshold
    "distant":   dict(method="window",           pm_k=3.0, cmd_w=1.0, spatial_w=1.0, threshold=0.25),
}

CLUSTERS = {
    # ---- clean halo fields ----
    "NGC 288":  dict(strategy="clean", pm=(4.164, -5.705),  pm_disp=0.15),
    "NGC 1851": dict(strategy="clean", pm=(2.145, -0.650),  pm_disp=0.20),
    "NGC 5024": dict(strategy="clean", pm=(-0.141, -1.331), pm_disp=0.15),
    "NGC 5139": dict(strategy="clean", pm=(-3.250, -6.746), pm_disp=0.75),  # Omega Cen
    "NGC 5272": dict(strategy="clean", pm=(-0.152, -2.670), pm_disp=0.20),  # M3
    "NGC 5904": dict(strategy="clean", pm=(4.086, -9.870),  pm_disp=0.25),  # M5
    "NGC 6341": dict(strategy="clean", pm=(-4.935, -0.625), pm_disp=0.20),  # M92
    "NGC 6752": dict(strategy="clean", pm=(-3.161, -4.027), pm_disp=0.35),
    "NGC 6809": dict(strategy="clean", pm=(-3.432, -9.311), pm_disp=0.25),  # M55
    "NGC 7078": dict(strategy="clean", pm=(-0.659, -3.803), pm_disp=0.25),  # M15
    "NGC 7089": dict(strategy="clean", pm=(3.435, -2.159),  pm_disp=0.20),  # M2

    # ---- bulge / disk crowded fields ----
    "Terzan 5": dict(strategy="crowded", pm=(-1.989, -5.243),   pm_disp=0.30),
    "NGC 6440": dict(strategy="crowded", pm=(-1.187, -4.020),   pm_disp=0.25),
    "NGC 6544": dict(strategy="crowded", pm=(-2.304, -18.604),  pm_disp=0.45),
    "NGC 6266": dict(strategy="crowded", pm=(-4.978, -2.947),   pm_disp=0.40),  # M62
    "NGC 6273": dict(strategy="crowded", pm=(-3.249, 1.660),    pm_disp=0.30),  # M19
    "NGC 6656": dict(strategy="crowded", pm=(9.851, -5.617),    pm_disp=0.55),  # M22
    "NGC 5286": dict(strategy="crowded", pm=(0.197, -0.153),    pm_disp=0.25, cmd_w=0.5),
    "NGC 6397": dict(strategy="crowded", pm=(3.260, -17.664),   pm_disp=0.45, cmd_w=0.5),
    "NGC 6121": dict(strategy="crowded", pm=(-12.514, -19.022), pm_disp=0.50, cmd_w=0.3),  # M4

    # ---- SMC companion in the field ----
    "NGC 104": dict(strategy="companion", pm=(5.252, -2.551), pm_disp=0.60),  # 47 Tuc
    "NGC 362": dict(strategy="companion", pm=(6.694, -2.535), pm_disp=0.25),

    # ---- special cases ----
    "Pal 5":    dict(strategy="sparse",  pm=(-2.730, -2.654), pm_disp=0.10),
    "NGC 2419": dict(strategy="distant", pm=(-0.008, -0.559), pm_disp=0.10),
}


def cluster_cfg(name):
    entry = CLUSTERS.get(name)
    if entry is None:
        # unknown cluster: behave like v1 (clean, no literature anchor)
        cfg = dict(STRATEGY_DEFAULTS["clean"])
        cfg.update(strategy="clean", pm=None, pm_disp=0.30)
        return cfg
    cfg = dict(STRATEGY_DEFAULTS[entry["strategy"]])
    cfg.update(entry)
    return cfg


@dataclass(frozen=True)
class ClusterConfig:
    name: str
    common_name: str
    ra_deg: float
    dec_deg: float
    distance_kpc: float
    core_radius_arcmin: float
    tidal_radius_arcmin: float
    search_radius_deg: float

    @property
    def safe_name(self):
        return "".join(ch for ch in self.name.replace(" ", "_") if ch.isalnum() or ch in "_-")


def read_targets(target_list):
    df = pd.read_excel(target_list, sheet_name="GC Target List", header=3)

    numeric_cols = [
        "RA (deg)",
        "Dec (deg)",
        "Core radius r_c (arcmin)",
        "Tidal radius r_t (arcmin)",
        "Heliocentric distance (kpc)",
    ]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=["Cluster ID", "RA (deg)", "Dec (deg)"]).reset_index(drop=True)

    # The spreadsheet's r_t = r_c * 10^c breaks for core-collapsed clusters
    # (r_c is a placeholder there, as the sheet's own notes warn). Two came
    # out physically impossible (r_t < 3x half-light radius), truncating the
    # search circle inside the cluster itself. Literature values instead:
    #   NGC 6544: ~18' (Baumgardt & Hilker 2018 structural catalogue)
    #   Pal 5:    ~16' (Odenkirchen et al. 2002)
    RT_OVERRIDES_ARCMIN = {"NGC 6544": 18.0, "Pal 5": 16.0}

    targets = []
    for _, row in df.iterrows():
        rt = float(row["Tidal radius r_t (arcmin)"])
        rt = RT_OVERRIDES_ARCMIN.get(str(row["Cluster ID"]).strip(), rt)
        targets.append(
            ClusterConfig(
                name=str(row["Cluster ID"]).strip(),
                common_name="" if pd.isna(row.get("Common Name")) else str(row["Common Name"]).strip(),
                ra_deg=float(row["RA (deg)"]),
                dec_deg=float(row["Dec (deg)"]),
                distance_kpc=float(row["Heliocentric distance (kpc)"]),
                core_radius_arcmin=float(row["Core radius r_c (arcmin)"]),
                tidal_radius_arcmin=rt,
                search_radius_deg=rt / 60.0,
            )
        )
    return targets


# ======================= GAIA QUERY (with cache reuse) =======================

def find_cached_raw(cluster, output_folder):
    candidates = [Path(output_folder) / f"{cluster.safe_name}_gaia_dr3_raw.csv"]
    candidates += [folder / f"{cluster.safe_name}_gaia_dr3_raw.csv" for folder in RAW_CACHE_FOLDERS]
    for path in candidates:
        if path.exists():
            return path
    return None


def query_gaia_dr3(cluster, output_folder):
    cached = find_cached_raw(cluster, output_folder)
    if cached is not None:
        print(f"[{cluster.name}] using cached Gaia query: {cached}")
        return pd.read_csv(cached)

    from astroquery.gaia import Gaia

    print(f"[{cluster.name}] querying Gaia DR3...")
    Gaia.ROW_LIMIT = -1

    def build_query(r_out, r_in=None):
        ring = ""
        if r_in is not None and r_in > 0:
            ring = f"""
      AND 0 = CONTAINS(
        POINT('ICRS', ra, dec),
        CIRCLE('ICRS', {cluster.ra_deg}, {cluster.dec_deg}, {r_in})
      )"""
        return f"""
    SELECT
        source_id, ra, dec, ra_error, dec_error,
        parallax, parallax_error,
        pmra, pmra_error, pmdec, pmdec_error,
        phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag,
        bp_rp, ruwe
    FROM gaiadr3.gaia_source
    WHERE 1 = CONTAINS(
        POINT('ICRS', ra, dec),
        CIRCLE('ICRS', {cluster.ra_deg}, {cluster.dec_deg}, {r_out})
    ){ring}
      AND pmra IS NOT NULL AND pmdec IS NOT NULL
      AND pmra_error IS NOT NULL AND pmdec_error IS NOT NULL
      AND parallax IS NOT NULL AND parallax_error IS NOT NULL
      AND ruwe IS NOT NULL AND phot_g_mean_mag IS NOT NULL AND bp_rp IS NOT NULL
      AND ruwe < {RUWE_MAX}
      AND pmra_error < {PM_ERR_MAX} AND pmdec_error < {PM_ERR_MAX}
      AND phot_g_mean_mag < {GAIA_G_MAG_MAX}
    """

    def run_query(adql, label, attempts=3):
        last_exc = None
        for attempt in range(1, attempts + 1):
            try:
                job = Gaia.launch_job_async(adql)
                return job.get_results().to_pandas()
            except Exception as exc:
                last_exc = exc
                print(f"[{cluster.name}] {label} attempt {attempt} failed: {exc!r}")
                time.sleep(15 * attempt)
        raise last_exc

    t0 = time.time()
    R = cluster.search_radius_deg
    try:
        df = run_query(build_query(R), "query")
    except Exception:
        # Large fields in dense regions can keep dropping mid-download.
        # Fall back to four annuli of equal area (smaller result chunks).
        print(f"[{cluster.name}] full-field query failed -> retrying in 4 annulus chunks")
        fracs = np.sqrt([0.0, 0.25, 0.5, 0.75, 1.0]) * R
        chunks = []
        for r_in, r_out in zip(fracs[:-1], fracs[1:]):
            chunk = run_query(
                build_query(r_out, r_in),
                f"annulus {r_in:.3f}-{r_out:.3f} deg",
            )
            print(f"[{cluster.name}]   annulus {r_in:.3f}-{r_out:.3f} deg: {len(chunk):,} stars")
            chunks.append(chunk)
        df = pd.concat(chunks, ignore_index=True).drop_duplicates(subset="source_id")
    print(f"[{cluster.name}] Gaia query returned {len(df):,} stars in {time.time() - t0:.1f} s")

    raw_csv = Path(output_folder) / f"{cluster.safe_name}_gaia_dr3_raw.csv"
    raw_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(raw_csv, index=False)
    print(f"[{cluster.name}] saved raw Gaia query: {raw_csv}")
    return df


# ======================= FILTER STAGES =======================

def quality_cuts(df, cluster):
    df = df.copy()
    for col in ["ruwe", "pmra_error", "pmdec_error", "parallax", "parallax_error"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    cluster_plx = 1.0 / cluster.distance_kpc
    good_ruwe = df["ruwe"].isna() | (df["ruwe"] < RUWE_MAX)
    good_pm = (df["pmra_error"].isna() | (df["pmra_error"] < PM_ERR_MAX)) & (
        df["pmdec_error"].isna() | (df["pmdec_error"] < PM_ERR_MAX)
    )
    plx_err = df["parallax_error"].replace(0, np.nan)
    plx_excess_sig = (df["parallax"] - cluster_plx) / plx_err
    not_fg = df["parallax"].isna() | plx_err.isna() | (plx_excess_sig < PARALLAX_NSIGMA)

    return df[good_ruwe & good_pm & not_fg].reset_index(drop=True)


def pm_window_score(df, valid, pm0, pm_disp, pm_k):
    """Per-star Gaussian score around the literature PM, folding in each
    star's own measurement errors. Also returns the hard-gate mask."""
    dx = df.loc[valid, "pmra"].to_numpy() - pm0[0]
    dy = df.loc[valid, "pmdec"].to_numpy() - pm0[1]
    sx = np.sqrt(pm_disp**2 + df.loc[valid, "pmra_error"].to_numpy() ** 2)
    sy = np.sqrt(pm_disp**2 + df.loc[valid, "pmdec_error"].to_numpy() ** 2)
    z2 = (dx / sx) ** 2 + (dy / sy) ** 2
    score = np.exp(-0.5 * z2)
    inside = z2 <= pm_k**2
    return score, inside


def pm_membership(df, cluster, cfg):
    df = df.copy()
    for col in ["ra", "dec", "pmra", "pmdec", "pmra_error", "pmdec_error"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    valid = df[["pmra", "pmdec", "pmra_error", "pmdec_error"]].notna().all(axis=1)
    df["p_pm"] = 0.0
    df["pm_method"] = ""
    if valid.sum() < 5:
        df.loc[valid, "p_pm"] = 1.0
        df.loc[valid, "pm_method"] = "too_few"
        return df

    pm0 = cfg.get("pm")
    if pm0 is None:
        # no literature anchor: estimate from stars near the centre (v1 style)
        dx = (df["ra"] - cluster.ra_deg) * np.cos(np.radians(cluster.dec_deg)) * 60.0
        dy = (df["dec"] - cluster.dec_deg) * 60.0
        r = np.sqrt(dx**2 + dy**2)
        core = df[valid & (r < cluster.core_radius_arcmin)]
        if len(core) < 10:
            core = df[valid & (r <= np.nanpercentile(r[valid], 20))]
        pm0 = (core["pmra"].median(), core["pmdec"].median())
    pm0 = np.asarray(pm0, dtype=float)

    pm_disp = cfg["pm_disp"]
    pm_k = cfg["pm_k"]

    method = cfg["method"]
    if method == "gmm" and valid.sum() >= 20 * cfg.get("n_comp", 2):
        from sklearn.mixture import GaussianMixture

        X = df.loc[valid, ["pmra", "pmdec"]].to_numpy()
        gmm = GaussianMixture(
            n_components=cfg["n_comp"],
            covariance_type="full",
            random_state=0,
            n_init=GMM_N_INIT,
            reg_covar=1e-4,
        ).fit(X)

        comp = int(np.argmin(np.linalg.norm(gmm.means_ - pm0, axis=1)))
        anchor_dist = float(np.linalg.norm(gmm.means_[comp] - pm0))

        if anchor_dist <= ANCHOR_TOLERANCE:
            p = gmm.predict_proba(X)[:, comp]
            # hard gate: kill stars far from the literature PM even if the
            # broad field component happens to give them high posterior
            _, inside = pm_window_score(df, valid, pm0, pm_disp, pm_k)
            p = np.where(inside, p, 0.0)
            df.loc[valid, "p_pm"] = p
            df.loc[valid, "pm_method"] = f"gmm{cfg['n_comp']}"
            print(
                f"[{cluster.name}] PM: GMM({cfg['n_comp']}) anchored, "
                f"component at {anchor_dist:.2f} mas/yr from literature PM"
            )
            return df

        print(
            f"[{cluster.name}] PM: GMM component {anchor_dist:.2f} mas/yr from "
            f"literature PM (> {ANCHOR_TOLERANCE}) -> window fallback"
        )

    score, inside = pm_window_score(df, valid, pm0, pm_disp, pm_k)
    df.loc[valid, "p_pm"] = np.where(inside, np.maximum(score, 0.01), 0.0)
    df.loc[valid, "pm_method"] = "window"
    if method != "window":
        pass  # message already printed above
    else:
        print(f"[{cluster.name}] PM: error-aware window around literature PM")
    return df


def _king(r, amp, rc, bkg, rt):
    inner = 1.0 / np.sqrt(1.0 + (r / rc) ** 2)
    edge = 1.0 / np.sqrt(1.0 + (rt / rc) ** 2)
    return amp * np.clip(inner - edge, 0, None) ** 2 + bkg


def spatial_membership(df, cluster, n_bins=25):
    from scipy.optimize import curve_fit

    df = df.copy()
    dx = (df["ra"] - cluster.ra_deg) * np.cos(np.radians(cluster.dec_deg)) * 60.0
    dy = (df["dec"] - cluster.dec_deg) * 60.0
    r = np.sqrt(dx**2 + dy**2)
    df["r_proj_arcmin"] = r

    finite = np.isfinite(r)
    df["p_spatial"] = 0.0
    if finite.sum() < 5:
        df.loc[finite, "p_spatial"] = 1.0
        return df

    r_fit = r[finite]
    rt = min(cluster.tidal_radius_arcmin, max(np.nanmax(r_fit), 1.0))
    n_bins = max(5, min(n_bins, int(np.sqrt(len(r_fit)))))

    # Log-spaced bins so the core is resolved even when it is much smaller
    # than the field of view (linear bins lumped compact cores into the
    # first bin, and the King fit then pinned rc at its lower bound).
    rc_known = max(cluster.core_radius_arcmin, 0.05)
    r_lo = min(max(0.05, rc_known / 4.0), np.nanmax(r_fit) / 10.0)
    edges = np.geomspace(r_lo, np.nanmax(r_fit), n_bins + 1)
    edges = np.concatenate([[0.0], edges])
    counts, _ = np.histogram(r, bins=edges)
    areas = np.pi * (edges[1:] ** 2 - edges[:-1] ** 2)
    density = counts / areas
    # Poisson uncertainty per bin: without this the fit is dominated by the
    # tiny, ultra-dense inner bins and the King profile collapses.
    sigma = np.sqrt(np.maximum(counts, 1)) / areas
    r_mid = 0.5 * (edges[1:] + edges[:-1])
    keep = counts > 0
    density, sigma, r_mid = density[keep], sigma[keep], r_mid[keep]

    try:
        popt, _ = curve_fit(
            lambda rr, amp, rc, bkg: _king(rr, amp, rc, bkg, rt),
            r_mid,
            density,
            sigma=sigma,
            p0=[
                max(density.max() - np.median(density), 1e-6),
                np.clip(rc_known, 0.1, 0.5 * rt),
                max(np.median(density), 1e-6),
            ],
            # rc may exceed the catalogue core radius: for compact clusters
            # Gaia only resolves an inflated "effective" core (crowding),
            # so the upper bound is set by the tidal radius instead
            bounds=([0, 0.05, 0], [np.inf, max(0.5 * rt, 0.2), np.inf]),
            maxfev=10000,
        )
        amp, rc, bkg = popt
        edge = 1.0 / np.sqrt(1.0 + (rt / rc) ** 2)

        # Sanity check: in crowding-limited cores (e.g. M62) Gaia loses the
        # central stars and the density profile is flat or even inverted.
        # The King fit then reports "no cluster" - in that case the spatial
        # term carries no information, so leave it neutral (NaN -> PM-only).
        central_frac = amp * (1.0 - edge) ** 2
        central_frac = central_frac / max(central_frac + bkg, 1e-12)
        if central_frac < 0.3:
            print(f"[{cluster.name}] spatial: no central concentration in "
                  f"Gaia counts (crowding-limited core?) -> PM-only")
            df["p_spatial"] = np.nan
            return df

        inner = 1.0 / np.sqrt(1.0 + (r / rc) ** 2)
        cl = amp * np.clip(inner - edge, 0, None) ** 2
        # p_spatial = LOCAL cluster fraction at this radius (raw, not
        # normalised): what fraction of stars at radius r are cluster stars
        p = cl / np.maximum(cl + bkg, 1e-12)
        df.loc[finite, "p_spatial"] = np.clip(p[finite], 0.01, 0.99)
    except Exception:
        scale = max(np.nanpercentile(r_fit, 50), 1.0)
        df.loc[finite, "p_spatial"] = np.clip(
            np.exp(-0.5 * (r_fit / scale) ** 2), 0.01, 0.99
        )

    return df


def cmd_membership(df, mag_bin=0.5):
    df = df.copy()
    df["p_cmd"] = np.nan
    if "bp_rp" not in df.columns or "phot_g_mean_mag" not in df.columns:
        return df

    df["bp_rp"] = pd.to_numeric(df["bp_rp"], errors="coerce")
    df["phot_g_mean_mag"] = pd.to_numeric(df["phot_g_mean_mag"], errors="coerce")

    p_sp = df["p_spatial"].fillna(1.0)  # NaN spatial = uninformative, don't veto
    train = df[(df["p_pm"] * p_sp) > CMD_TRAIN_PROB].dropna(
        subset=["phot_g_mean_mag", "bp_rp"]
    )
    if len(train) < 50:
        train = df[(df["p_pm"] > 0.50) & (p_sp > 0.25)].dropna(
            subset=["phot_g_mean_mag", "bp_rp"]
        )
    if len(train) < 50:
        return df

    g0, g1 = train["phot_g_mean_mag"].min(), train["phot_g_mean_mag"].max()
    edges = np.arange(g0, g1 + mag_bin, mag_bin)
    mids, cols, sigs = [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        sl = train[(train["phot_g_mean_mag"] >= lo) & (train["phot_g_mean_mag"] < hi)]
        if len(sl) < 10:
            continue
        med = sl["bp_rp"].median()
        sig = max(1.4826 * np.median(np.abs(sl["bp_rp"] - med)), 0.03)
        mids.append(0.5 * (lo + hi))
        cols.append(med)
        sigs.append(sig)

    if len(mids) < 2:
        return df

    mids, cols, sigs = map(np.array, (mids, cols, sigs))
    has = df["phot_g_mean_mag"].notna() & df["bp_rp"].notna()
    g = df.loc[has, "phot_g_mean_mag"].to_numpy()
    c = df.loc[has, "bp_rp"].to_numpy()
    z = (c - np.interp(g, mids, cols)) / np.interp(g, mids, sigs)
    df.loc[has, "p_cmd"] = np.exp(-0.5 * z**2)
    return df


def combine_and_flag(df, cfg):
    """Bayesian combination: the PM posterior provides the base odds, and
    the spatial term updates those odds by how cluster-dominated the field
    is at each star's radius, relative to the field as a whole. This way a
    star with decisive proper motion survives even in the outskirts, while
    a PM-ambiguous star needs the spatial prior on its side."""
    df = df.copy()

    p_pm = df["p_pm"].fillna(0.0).clip(1e-9, 1 - 1e-9)
    f_glob = float(np.clip(p_pm.mean(), 0.01, 0.99))
    # NaN p_spatial = "spatial term uninformative": fill with the global
    # fraction so the prior ratio is exactly 1 (pure PM decision)
    p_loc = df["p_spatial"].fillna(f_glob).clip(0.01, 0.99)

    prior_ratio = (p_loc / (1 - p_loc)) / (f_glob / (1 - f_glob))
    odds = (p_pm / (1 - p_pm)) * prior_ratio ** cfg["spatial_w"]
    p = odds / (1.0 + odds)

    # stars hard-gated out in PM space stay out
    p = p.where(df["p_pm"] > 0, 0.0)

    cmd_w = cfg["cmd_w"]
    if cmd_w > 0:
        has_cmd = df["p_cmd"].notna()
        cmd_factor = (0.5 + 0.5 * df["p_cmd"]) ** cmd_w
        p = p.where(~has_cmd, p * cmd_factor)

    df["membership_prob"] = p.fillna(0.0).clip(0.0, 1.0)
    df["is_member"] = df["membership_prob"] >= cfg["threshold"]
    return df


def add_distance_columns(df, min_snr=5.0):
    """Per-star distance from parallax — honest version.

    1/parallax is only meaningful when the parallax is measured well
    (positive, and at least min_snr times its error, i.e. <= 20% error).
    At globular-cluster distances that is true for only a few percent of
    stars (mostly foreground); everything else gets NaN on purpose.
    For cluster members, adopt the cluster's literature distance instead.
    """
    df = df.copy()
    if "parallax" in df.columns and "parallax_error" in df.columns:
        plx = pd.to_numeric(df["parallax"], errors="coerce")
        err = pd.to_numeric(df["parallax_error"], errors="coerce")
        ok = (plx > 0) & (err > 0) & (plx / err >= min_snr)
        df["dist_pc"] = np.where(ok, 1000.0 / plx, np.nan)
        df["dist_pc_err"] = np.where(ok, 1000.0 * err / plx**2, np.nan)
    return df


def filter_cluster(df, cluster, cfg):
    df = quality_cuts(df, cluster)
    if len(df) == 0:
        for col in ["p_pm", "pm_method", "r_proj_arcmin", "p_spatial", "p_cmd",
                    "membership_prob", "is_member"]:
            df[col] = []
        return df

    df = pm_membership(df, cluster, cfg)
    df = spatial_membership(df, cluster)
    df = cmd_membership(df)
    df = combine_and_flag(df, cfg)
    df = add_distance_columns(df)
    return df


# ======================= DRIVER =======================

def run_one_cluster(cluster, output_folder):
    cfg = cluster_cfg(cluster.name)
    print(f"[{cluster.name}] strategy: {cfg['strategy']}"
          + (f" ({cluster.common_name})" if cluster.common_name else ""))

    raw = query_gaia_dr3(cluster, output_folder)
    print(f"[{cluster.name}] raw stars: {len(raw):,}")

    t0 = time.time()
    filtered = filter_cluster(raw, cluster, cfg)

    out_dir = Path(output_folder)
    out_dir.mkdir(parents=True, exist_ok=True)
    filtered_csv = out_dir / f"{cluster.safe_name}_filtered.csv"
    filtered.to_csv(filtered_csv, index=False)

    n_members = int(filtered["is_member"].sum()) if len(filtered) else 0
    frac = n_members / len(filtered) if len(filtered) else 0.0
    print(f"[{cluster.name}] stars after quality cuts: {len(filtered):,}")
    print(f"[{cluster.name}] members: {n_members:,} ({frac:.1%})  "
          f"[{time.time() - t0:.1f} s]")
    return {
        "cluster": cluster.name,
        "common_name": cluster.common_name,
        "strategy": cfg["strategy"],
        "raw_stars": len(raw),
        "filtered_stars": len(filtered),
        "members": n_members,
        "member_fraction": frac,
        "filtered_csv": str(filtered_csv),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--target-list", default=TARGET_LIST_DEFAULT)
    parser.add_argument("--output-folder", default=str(OUTPUT_FOLDER_DEFAULT))
    parser.add_argument("--only", nargs="*", help="Optional cluster IDs, e.g. NGC 5139")
    parser.add_argument("--include-fpr", action="store_true",
                        help="Build a versioned NGC 5139 DR3+FPR CSV from the existing DR3 baseline")
    parser.add_argument("--dr3-baseline", type=Path,
                        default=SCRIPT_DIR.parent / "data/NGC5139/filtered.csv")
    parser.add_argument("--fpr-mirror", choices=["esa", "aip"], default="esa")
    args = parser.parse_args()

    if args.include_fpr:
        if args.only and args.only != ["NGC 5139"]:
            parser.error("The published FPR crowded-field extension supports NGC 5139 only.")
        from fpr_extension import build_catalog
        build_catalog(args.dr3_baseline, Path(args.output_folder), args.target_list, args.fpr_mirror)
        return

    targets = read_targets(args.target_list)
    if args.only:
        wanted = {name.strip().lower() for name in args.only}
        targets = [t for t in targets if t.name.lower() in wanted]

    print(f"Loaded {len(targets)} cluster target(s).")
    if not targets:
        raise SystemExit("No matching cluster targets found.")

    summaries, failures = [], []
    for cluster in targets:
        try:
            summaries.append(run_one_cluster(cluster, args.output_folder))
        except Exception as exc:
            failures.append({"cluster": cluster.name, "error": repr(exc)})
            print(f"[{cluster.name}] ERROR: {exc!r}")

    out_dir = Path(args.output_folder)
    out_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(summaries).to_csv(out_dir / "batch_summary.csv", index=False)
    if failures:
        pd.DataFrame(failures).to_csv(out_dir / "batch_failures.csv", index=False)

    print(f"\nBatch summary saved: {out_dir / 'batch_summary.csv'}")
    if failures:
        print(f"Failures saved: {out_dir / 'batch_failures.csv'}")


if __name__ == "__main__":
    main()
