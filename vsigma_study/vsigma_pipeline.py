"""Versioned Gaia PM rotation/dispersion entry point.

Interactive and batch analysis use pm_kinematics.py (joint Gaussian MCMC).
Batch runs save per-cluster rotation.json/png only; legacy cross-cluster
rankings are not regenerated from a different V/sigma definition.
"""
import json
import traceback
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

DATA = Path(__file__).resolve().parents[1] / "data"
OUT = Path(__file__).resolve().parent
K = 4.74047                      # (mas/yr)*(kpc) -> km/s
P_CUT = 0.90

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


def analyze_cluster(cid, want_rv=False, bin_mode="auto", nbins=None, catalog=None, save_figure=True, max_steps=4000):
    """Versioned joint PM inference; want_rv is retained for API compatibility."""
    try:
        from .pm_kinematics import cached_analyze, save_profile_figure
    except ImportError:
        from pm_kinematics import cached_analyze, save_profile_figure
    if cid not in CLUSTERS:
        raise ValueError("Unknown cluster")
    result = cached_analyze(cid, catalog or DATA / cid / "filtered.csv", CLUSTERS[cid], bin_mode, nbins, max_steps)
    if save_figure:
        save_profile_figure(result, DATA / cid / "rotation.png")
    return result


def main():
    for cid in CLUSTERS:
        if not (DATA / cid / "filtered.csv").exists():
            continue
        try:
            result = analyze_cluster(cid)
            (DATA / cid / "rotation.json").write_text(json.dumps(result, allow_nan=False), encoding="utf-8")
            print(f"{cid}: {result['binning']['count']} bins; {result['diagnostics']['unresolved_bins']} unavailable dispersion bins")
        except Exception:
            traceback.print_exc()
    print("New PM profiles saved. Legacy peak/central cross-cluster studies were not overwritten.")


if __name__ == "__main__":
    main()
