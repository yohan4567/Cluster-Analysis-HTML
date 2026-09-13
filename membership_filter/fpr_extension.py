"""Versioned Omega Cen FPR extension; never overwrites the DR3 baseline.

Schema: https://gea.esac.esa.int/archive/documentation/FPR/chap_datamodel/
sec_dm_focused_product_release/ssec_dm_crowded_field_source.html
FPR membership_prob is a compatibility SCORE, not a calibrated posterior.
"""
import io
import json
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd
import requests

TAP = "https://gea.esac.esa.int/tap-server/tap/sync"
VERSION = "dr3_fpr_v1"
# FPR crowded_field_source does not publish DR3's RUWE column.
COLUMNS = "source_id,ref_epoch,ra,dec,ra_error,dec_error,parallax,parallax_error,pmra,pmra_error,pmdec,pmdec_error,pmra_pmdec_corr,phot_g_mean_mag"


def tap_csv(query, endpoint=TAP):
    for attempt in range(3):
        try:
            response = requests.get(endpoint, params=dict(REQUEST="doQuery", LANG="ADQL",
                FORMAT="votable", MAXREC=1000000, QUERY=query), timeout=(30, 240))
            response.raise_for_status()
            root = ET.fromstring(response.content)
            for element in root.iter():
                if element.tag.split("}")[-1] == "INFO" and element.get("name") == "QUERY_STATUS":
                    if element.get("value") != "OK":
                        raise ValueError("Gaia query error/truncation: " + str(element.attrib) + str(element.text))
            from astropy.io.votable import parse_single_table
            df = parse_single_table(io.BytesIO(response.content)).to_table(use_names_over_ids=True).to_pandas()
            if "source_id" in df:
                df["source_id"] = df.source_id.astype("string")
            return df
        except (requests.RequestException, ValueError):
            if attempt == 2:
                raise
            time.sleep(2 * (attempt + 1))


def download_fpr(cluster, folder, endpoint=TAP):
    cache = folder / "NGC_5139_gaia_fpr_raw.csv"
    metadata = cache.with_suffix(".json")
    where = ("region_name = 'NGC5139' AND 1=CONTAINS(POINT('ICRS',ra,dec),"
             f"CIRCLE('ICRS',{cluster.ra_deg},{cluster.dec_deg},{cluster.search_radius_deg}))")
    query = f"SELECT {COLUMNS} FROM gaiafpr.crowded_field_source WHERE {where}"
    if cache.exists() and metadata.exists():
        info = json.loads(metadata.read_text(encoding="utf-8"))
        df = pd.read_csv(cache, dtype={"source_id": "string"})
        if info.get("query") == query and len(df) == info.get("rows"):
            return df
    print(f"Downloading official FPR sources from {endpoint} (before quality cuts)...", flush=True)
    expected = int(tap_csv("SELECT COUNT(*) AS n FROM gaiafpr.crowded_field_source WHERE " + where, endpoint).iloc[0, 0])
    df = tap_csv(query, endpoint)
    if len(df) != expected or df.source_id.isna().any() or df.source_id.duplicated().any():
        raise ValueError(f"Incomplete or duplicate FPR download: expected {expected}, received {len(df)}")
    if set(COLUMNS.split(",")) - set(df):
        raise ValueError("FPR response is missing required columns")
    df.to_csv(cache, index=False)
    metadata.write_text(json.dumps(dict(query=query, rows=len(df), endpoint=endpoint), indent=2), encoding="utf-8")
    return df


def filter_fpr(raw, cluster, cfg):
    from membership_filter import RUWE_MAX, PM_ERR_MAX, GAIA_G_MAG_MAX, PARALLAX_NSIGMA, add_distance_columns
    df = raw.copy()
    for col in COLUMNS.split(","):
        if col != "source_id":
            df[col] = pd.to_numeric(df[col], errors="coerce")
    required = [c for c in COLUMNS.split(",") if c != "source_id"]
    good = np.isfinite(df[required]).all(axis=1)
    good &= (df.ra >= 0) & (df.ra < 360) & df.dec.between(-90, 90)
    for col in ["pmra_error", "pmdec_error"]:
        good &= (df[col] > 0) & (df[col] < PM_ERR_MAX)
    good &= (df.parallax_error > 0) & (df.phot_g_mean_mag < GAIA_G_MAG_MAX)
    good &= df.pmra_pmdec_corr.between(-1, 1)
    good &= (df.parallax - 1 / cluster.distance_kpc) / df.parallax_error < PARALLAX_NSIGMA
    df = df.loc[good].copy().reset_index(drop=True)
    # Keep a schema-compatible empty column; RUWE is unavailable for FPR.
    df["ruwe"] = np.nan
    # First-order epoch alignment is sufficient over 1.5 yr here. Native
    # coordinates/epoch are retained; PM components are not re-zeroed.
    df["catalog_ra"], df["catalog_dec"] = df.ra, df.dec
    df["catalog_ref_epoch"] = df.ref_epoch
    dt = 2016.0 - df.ref_epoch
    df["ra"] = (df.ra + dt * df.pmra / (3.6e6 * np.cos(np.radians(df.dec)))) % 360
    df["dec"] = df.dec + dt * df.pmdec / 3.6e6
    df["ref_epoch"] = 2016.0
    x = (df.ra - cluster.ra_deg) * np.cos(np.radians(cluster.dec_deg)) * 60
    y = (df.dec - cluster.dec_deg) * 60
    df["r_proj_arcmin"] = np.hypot(x, y)
    df = df.loc[df.r_proj_arcmin <= cluster.tidal_radius_arcmin].copy()
    # Mahalanobis distance: intrinsic isotropic dispersion + per-star PM
    # covariance. Do not fit a naive GMM to heteroscedastic FPR measurements.
    a = cfg["pm_disp"] ** 2 + df.pmra_error ** 2
    b = cfg["pm_disp"] ** 2 + df.pmdec_error ** 2
    c = df.pmra_pmdec_corr * df.pmra_error * df.pmdec_error
    dx, dy = df.pmra - cfg["pm"][0], df.pmdec - cfg["pm"][1]
    z2 = np.maximum((b * dx**2 + a * dy**2 - 2*c*dx*dy) / (a*b-c*c), 0)
    df["p_pm"] = np.where(z2 <= cfg["pm_k"]**2, np.exp(-0.5*z2), 0)
    df["pm_method"] = "fpr_covariance_window"
    for col in ["bp_rp", "phot_bp_mean_mag", "phot_rp_mean_mag", "p_cmd", "p_spatial"]:
        df[col] = np.nan
    df["membership_prob"] = df.p_pm
    df["is_member"] = df.membership_prob >= cfg["threshold"]
    df["catalog_origin"] = "gaia_fpr"
    df["membership_method"] = "pm_score_only_uncalibrated"
    df["cmd_used"] = False
    df["spatial_used"] = False
    df["kinematics_validated"] = False
    return add_distance_columns(df)


def build_catalog(baseline, folder, target_list, mirror="esa"):
    from membership_filter import read_targets, cluster_cfg
    cluster = next(c for c in read_targets(target_list) if c.name == "NGC 5139")
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    dr3 = pd.read_csv(baseline, dtype={"source_id": "string"})
    required = {"source_id", "ra", "dec", "pmra", "pmdec", "membership_prob", "is_member"}
    if required - set(dr3):
        raise ValueError("Baseline must be an existing filtered DR3 catalog")
    if "catalog_origin" in dr3 and not dr3.catalog_origin.eq("gaia_dr3").all():
        raise ValueError("Baseline contains non-DR3 rows; use the original catalog")
    if dr3.source_id.isna().any() or dr3.source_id.duplicated().any():
        raise ValueError("Baseline has missing/duplicate source IDs")
    flags = dr3.is_member.astype(str).str.lower()
    if not flags.isin(["true", "false", "1", "0"]).all():
        raise ValueError("Invalid baseline is_member flags")
    dr3["is_member"] = flags.isin(["true", "1"])
    dr3["catalog_origin"] = "gaia_dr3"
    dr3["ref_epoch"] = dr3.get("ref_epoch", 2016.0)
    dr3["membership_method"] = "existing_dr3_filter"
    dr3["cmd_used"] = dr3.get("p_cmd", pd.Series(np.nan, index=dr3.index)).notna()
    dr3["spatial_used"] = dr3.get("p_spatial", pd.Series(np.nan, index=dr3.index)).notna()
    dr3["kinematics_validated"] = False
    endpoint = {"esa": TAP, "aip": "https://gaia.aip.de/tap/sync"}[mirror]
    raw = download_fpr(cluster, folder, endpoint)
    fpr = filter_fpr(raw, cluster, cluster_cfg(cluster.name))
    # ESA publishes FPR as an add-on excluding nominal sources. Preserve
    # close resolved companions; never blindly delete nearest neighbours.
    from scipy.spatial import cKDTree
    def unit_vectors(frame):
        ra, dec = np.radians(frame.ra.to_numpy()), np.radians(frame.dec.to_numpy())
        return np.column_stack([np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)])
    if len(fpr):
        chord, _ = cKDTree(unit_vectors(dr3)).query(unit_vectors(fpr))
        fpr["nearest_dr3_arcsec"] = np.degrees(2*np.arcsin(np.clip(chord/2, 0, 1))) * 3600
        fpr["close_dr3_review"] = fpr.nearest_dr3_arcsec < 0.1
    combined = pd.concat([dr3, fpr], ignore_index=True)
    combined["source_key"] = combined.catalog_origin + ":" + combined.source_id
    if combined.source_key.duplicated().any():
        raise ValueError("Duplicate catalog-qualified IDs")
    combined["filter_version"] = VERSION
    output = folder / "NGC_5139_filtered_dr3_fpr_v1.csv"
    combined.to_csv(output, index=False)
    counts = []
    for origin, part in combined.groupby("catalog_origin"):
        x = (part.ra-cluster.ra_deg)*np.cos(np.radians(cluster.dec_deg))*60
        y = (part.dec-cluster.dec_deg)*60
        core = np.hypot(x, y) <= cluster.core_radius_arcmin
        counts.append(dict(catalog=origin, rows=len(part), members=int(part.is_member.sum()),
                           core_rows=int(core.sum()), core_members=int((core & part.is_member).sum())))
    report = dict(version=VERSION, baseline=str(Path(baseline).resolve()), output=str(output.resolve()),
        fpr_downloaded=len(raw), counts=counts, core_radius_arcmin=cluster.core_radius_arcmin,
        fpr_cuts=dict(ruwe_max=1.4, pm_error_max=1.5, g_max=20.5, foreground_sigma=3,
                      membership_score_min=0.5, intrinsic_pm_dispersion=0.75),
        notes=["DR3 selections and scores preserved.",
               "FPR score is not a calibrated membership posterior; no colour or spatial prior used.",
               "FPR coordinates propagated to J2016.0; original coordinates and epoch retained.",
               "ESA FPR add-on excludes nominal sources. Close neighbours are flagged, not deleted.",
               "Not validated for combined kinematic analysis; old rotation products are unchanged."])
    output.with_suffix(".json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2), flush=True)
    return combined
