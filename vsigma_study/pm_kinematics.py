"""Projected PM kinematics. Conditional on catalog membership and reported errors.

See METHODS_MCMC.md for priors, bin selection, and limitations. No catalog edits.
"""
from pathlib import Path
import hashlib
import json
import threading

import emcee
import numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2

VERSION = "pm-mcmc-2"
K = 4.74047
LOCK = threading.Lock()


def basis(ra, dec):
    a, d = np.radians(ra), np.radians(dec)
    n = np.stack([np.cos(d)*np.cos(a), np.cos(d)*np.sin(a), np.sin(d)], axis=-1)
    e = np.stack([-np.sin(a), np.cos(a), np.zeros_like(a)], axis=-1)
    north = np.stack([-np.sin(d)*np.cos(a), -np.sin(d)*np.sin(a), np.cos(d)], axis=-1)
    return n, e, north


def transform_pm(ra, dec, pmra, pmdec, errors, corr, center, distance, vlos):
    """Exact spherical bulk-velocity projection, assuming a common distance.

    Positive tangential direction: E to N on an E-right sky map.
    errors is (N,2), in mas/yr. Returns PMs and covariance in local R,t basis.
    """
    n, e, north = basis(ra, dec)
    n0, e0, north0 = basis(*center)
    cosine = np.clip(n @ n0, -1, 1)
    sine = np.linalg.norm(np.cross(n0, n), axis=1)
    radius = np.arctan2(sine, cosine) * 180 / np.pi * 60
    if np.any(sine < 1e-12):
        raise ValueError("A star lies exactly at the center; its tangential direction is undefined")
    radial = (cosine[:, None]*n - n0) / sine[:, None]
    tangent = np.cross(n, radial)
    rotation = np.stack([np.stack([np.sum(radial*e, axis=1), np.sum(radial*north, axis=1)], axis=1),
                         np.stack([np.sum(tangent*e, axis=1), np.sum(tangent*north, axis=1)], axis=1)], axis=1)
    # Solve the projected bulk PM, rather than subtracting coordinate means.
    design = np.stack([np.stack([e @ e0, e @ north0], axis=1),
                       np.stack([north @ e0, north @ north0], axis=1)], axis=1)
    perspective = vlos / (K*distance) * np.stack([e @ n0, north @ n0], axis=1)
    observed = np.column_stack([pmra, pmdec])
    target = observed - perspective
    keep = np.ones(len(ra), dtype=bool)
    for _ in range(4):
        bulk = np.linalg.lstsq(design[keep].reshape(-1, 2), target[keep].ravel(), rcond=None)[0]
        residual = target - np.einsum("nij,j->ni", design, bulk)
        scale = np.maximum(1.4826*np.median(np.abs(residual[keep]-np.median(residual[keep], axis=0)), axis=0), 1e-6)
        candidate = np.all(np.abs(residual) < 4*scale, axis=1)
        if candidate.sum() < 10:
            break
        keep = candidate
    cov = np.zeros((len(ra), 2, 2))
    cov[:, 0, 0], cov[:, 1, 1] = errors[:, 0]**2, errors[:, 1]**2
    cov[:, 0, 1] = cov[:, 1, 0] = corr*errors[:, 0]*errors[:, 1]
    return (radius, np.einsum("nij,nj->ni", rotation, residual),
            rotation @ cov @ rotation.transpose(0, 2, 1), bulk)


def select_bins(radius, covariance, rh, mode="auto", count=None):
    """Pilot information uses only errors and a fixed 0.2 mas/yr reference.

    No fitting, rotation significance, or positive-dispersion test selects bins.
    Auto joins adjacent equal-count seeds, subject to a 0.3 dex width cap
    in log10(1+R/Rh). Low-information bins at the cap remain flagged.
    """
    n = len(radius)
    if mode not in ("auto", "equal_number", "equal_radius"):
        raise ValueError("Unknown bin mode")
    if mode != "auto" and (not isinstance(count, int) or not 2 <= count <= min(100, n//2)):
        raise ValueError("Manual bin count must be 2–100 and at most half the usable stars")
    order = np.argsort(radius, kind="stable")
    information = np.mean((0.2**2 / (0.2**2 + np.diagonal(covariance, axis1=1, axis2=2)))**2, axis=1)
    if mode == "equal_radius":
        edges = np.linspace(0, radius.max(), count+1)
        labels = np.clip(np.searchsorted(edges, radius, side="right")-1, 0, count-1)
        groups = [np.flatnonzero(labels == i) for i in range(count)]
    else:
        seeds = np.array_split(order, count if mode != "auto" else max(1, min(24, n//200)))
        groups = []
        pending = seeds[0]
        for nxt in seeds[1:]:
            enough = len(pending) >= 200 and information[pending].sum() >= 100
            width = np.log10((rh + radius[nxt[-1]]) / (rh + radius[pending[0]]))
            if mode != "auto" or enough or width > 0.3:
                groups.append(pending)
                pending = nxt
            else:
                pending = np.concatenate([pending, nxt])
        groups.append(pending)
        if mode == "auto" and len(groups) > 1 and information[groups[-1]].sum() < 100:
            joined = np.concatenate(groups[-2:])
            if np.log10((rh+radius[joined[-1]])/(rh+radius[joined[0]])) <= 0.3:
                groups[-2:] = [joined]
        edges = [radius.min()] + [(radius[a[-1]]+radius[b[0]])/2 for a, b in zip(groups[:-1], groups[1:])] + [radius.max()]
    diagnostics = [dict(n=len(g), effective_information=float(information[g].sum()),
                        low_information=bool(information[g].sum() < 100 or len(g) < 200),
                        broad=bool(len(g) and np.log10((rh+radius[g].max())/(rh+radius[g].min())) > 0.3)) for g in groups]
    return groups, np.asarray(edges), diagnostics


def fit_bin(u, covariance, seed=42, max_steps=4000, prior_scale=1.0):
    """Gaussian R,t fit using emcee. Means N(0,5²), s half-normal(scale).

    Flat-space sampling in s>=0 includes the zero boundary without a log floor.
    A conservative likelihood-ratio gate controls reporting of dispersion.
    Autocorrelation checks are diagnostics, not a proof of convergence.
    """
    if len(u) < 30:
        return dict(status="too_few_stars", converged=False)
    a, b, c = covariance[:, 0, 0], covariance[:, 1, 1], covariance[:, 0, 1]
    def likelihood(p):
        mr, mt, sr, st = p
        if sr < 0 or st < 0:
            return -np.inf
        aa, bb = a+sr*sr, b+st*st
        det = aa*bb-c*c
        dr, dt = u[:, 0]-mr, u[:, 1]-mt
        return -0.5*np.sum(np.log(det)+(bb*dr*dr+aa*dt*dt-2*c*dr*dt)/det)
    def posterior(p):
        ll = likelihood(p)
        return ll - .5*np.sum((p[:2]/5)**2) - .5*np.sum((p[2:]/prior_scale)**2) if np.isfinite(ll) else -np.inf
    initial = np.r_[np.median(u, axis=0), np.sqrt(np.maximum(np.var(u, axis=0)-[a.mean(), b.mean()], 1e-5))]
    bounds = [(None, None), (None, None), (0, None), (0, None)]
    opt = minimize(lambda p: -likelihood(p), initial, bounds=bounds, method="L-BFGS-B")
    null = minimize(lambda m: -likelihood(np.r_[m, 0., 0.]), opt.x[:2], method="BFGS")
    lr = max(0., 2*(likelihood(opt.x)-likelihood(np.r_[null.x, 0., 0.])))
    dr, dt = u[:, 0]-null.x[0], u[:, 1]-null.x[1]
    noise_chi = np.sum((b*dr*dr+a*dt*dt-2*c*dr*dt)/(a*b-c*c))
    inconsistent = bool(noise_chi < chi2.ppf(.001, max(2*len(u)-2, 1)))
    rng = np.random.default_rng(seed)
    start = opt.x + rng.normal(size=(20, 4))*np.maximum(np.abs(opt.x)*.025, .002)
    start[:, 2:] = np.abs(start[:, 2:]) + 1e-6
    sampler = emcee.EnsembleSampler(20, 4, posterior)
    sampler.random_state = np.random.RandomState(seed).get_state()
    state = start
    previous_tau = None
    converged = False
    tau = np.full(4, np.nan)
    for offset in range(0, max_steps, 500):
        state = sampler.run_mcmc(state, min(500, max_steps-offset), progress=False)
        if sampler.iteration < 1000:
            continue
        tau = sampler.get_autocorr_time(tol=0, discard=sampler.iteration//3)
        length = sampler.iteration - sampler.iteration//3
        stable = previous_tau is not None and np.all(np.abs(tau-previous_tau)/tau < .1)
        if np.all(length > 50*tau) and stable:
            converged = True
            break
        previous_tau = tau
    samples = sampler.get_chain(discard=sampler.iteration//3, flat=True)
    combined = np.sqrt(np.mean(samples[:, 2:]**2, axis=1))
    # Weak-prior sensitivity check by importance reweighting to scales /2 and *2.
    sensitivity = 0.
    for scale in (prior_scale/2, prior_scale*2):
        logw = -.5*np.sum(samples[:, 2:]**2, axis=1)*(1/scale**2-1/prior_scale**2)
        weights = np.exp(logw-logw.max())
        ix = np.argsort(combined)
        alt = np.interp(.5, np.cumsum(weights[ix])/weights.sum(), combined[ix])
        sensitivity = max(sensitivity, abs(alt-np.median(combined))/max(np.std(combined), 1e-8))
    status = ("not_converged" if not converged else "error_model_conflict" if inconsistent
              else "prior_sensitive" if sensitivity > .5 else "resolved" if lr > 9 else "unresolved")
    quantiles = lambda values: np.quantile(values, [.16, .5, .84]).tolist()
    return dict(status=status, converged=converged, steps=sampler.iteration,
                tau=tau.tolist(), effective_samples=(len(samples)/tau).tolist(),
                acceptance=float(np.mean(sampler.acceptance_fraction)),
                likelihood_ratio=lr, prior_sensitivity=sensitivity,
                mean_r=quantiles(samples[:, 0]), mean_t=quantiles(samples[:, 1]),
                sigma_r=quantiles(samples[:, 2]), sigma_t=quantiles(samples[:, 3]),
                sigma=quantiles(combined), sigma_upper95=float(np.quantile(combined, .95)),
                ratio=quantiles(np.abs(samples[:, 1])/combined),
                positive_probability=float(np.mean(samples[:, 1] > 0)))


def analyze(cid, catalog, parameters, mode="auto", count=None, max_steps=4000):
    import pandas as pd
    ra0, dec0, rc, rh, conc, rt, feh, distance, litrot, vlos = parameters
    df = pd.read_csv(catalog, low_memory=False)
    required = ["ra", "dec", "pmra", "pmdec", "pmra_error", "pmdec_error", "membership_prob"]
    if any(k not in df for k in required):
        raise ValueError("Catalog requires PM errors and membership_prob")
    for k in required:
        df[k] = pd.to_numeric(df[k], errors="coerce")
    if "pmra_pmdec_corr" not in df:
        df["pmra_pmdec_corr"] = np.nan
    df["pmra_pmdec_corr"] = pd.to_numeric(df.pmra_pmdec_corr, errors="coerce")
    missing_correlation = df.pmra_pmdec_corr.isna()
    # Legacy DR3 exports did not retain rho. This is an explicit approximation,
    # reported in the API/UI; nonmissing invalid correlations remain excluded.
    df["pmra_pmdec_corr"] = df.pmra_pmdec_corr.fillna(0.)
    members = df.membership_prob >= .9
    valid = np.isfinite(df[required]).all(axis=1) & (df.pmra_error > 0) & (df.pmdec_error > 0) & (df.pmra_pmdec_corr.abs() < 1)
    valid &= (df.ra.between(0, 360) & df.dec.between(-90, 90))
    valid &= ((df.ra-ra0).abs() + (df.dec-dec0).abs()) > 1e-10
    m = df[members & valid]
    if len(m) < 50:
        raise ValueError("At least 50 members with valid PM covariance are required; no relaxed membership fallback is used")
    radius, u, cov, bulk = transform_pm(m.ra.to_numpy(), m.dec.to_numpy(), m.pmra.to_numpy(), m.pmdec.to_numpy(),
                                      m[["pmra_error", "pmdec_error"]].to_numpy(), m.pmra_pmdec_corr.to_numpy(),
                                      (ra0, dec0), distance, vlos)
    groups, edges, bin_diagnostics = select_bins(radius, cov, rh, mode, count)
    kd = K*distance
    bins = {key: [] for key in ("r_mid_arcmin", "r_mid_pc", "r_over_rh", "n", "v_rot", "v_rot_err", "v_rot_lo", "v_rot_hi",
                               "v_rad", "v_rad_err", "sigma", "sigma_err", "sigma_lo", "sigma_hi", "sigma_r", "sigma_t",
                               "sigma_r_lo", "sigma_r_hi", "sigma_t_lo", "sigma_t_hi", "sigma_upper95", "ratio", "ratio_lo", "ratio_hi",
                               "observed_sigma", "rms_error", "status", "rotation_detected")}
    fits = []
    for i, g in enumerate(groups):
        fit = fit_bin(u[g], cov[g], seed=42+i, max_steps=max_steps)
        fits.append({**bin_diagnostics[i], **fit})
        good_mean = fit["converged"] and fit["status"] not in ("error_model_conflict", "prior_sensitive")
        resolved = fit["status"] == "resolved"
        radius_mid = float(np.median(radius[g])) if len(g) else None
        row = dict(r_mid_arcmin=radius_mid, r_mid_pc=radius_mid*distance*1000*np.pi/10800 if len(g) else None,
                   r_over_rh=radius_mid/rh if len(g) else None, n=len(g), status=fit["status"],
                   observed_sigma=float(np.sqrt(np.mean(np.var(u[g], axis=0)))*kd) if len(g) else None,
                   rms_error=float(np.sqrt(np.mean(np.trace(cov[g], axis1=1, axis2=2)/2))*kd) if len(g) else None,
                   sigma_upper95=fit.get("sigma_upper95", 0)*kd if fit["status"] == "unresolved" else None,
                   rotation_detected=bool(good_mean and min(fit["positive_probability"], 1-fit["positive_probability"]) < .003/(2*len(groups))))
        for key, source, usable in (("v_rot", "mean_t", good_mean), ("v_rad", "mean_r", good_mean),
                                    ("sigma", "sigma", resolved), ("sigma_r", "sigma_r", resolved), ("sigma_t", "sigma_t", resolved),
                                    ("ratio", "ratio", resolved)):
            lo, mid, hi = np.array(fit[source])*(1 if key == "ratio" else kd) if usable else (None, None, None)
            row[key], row[key+"_lo"], row[key+"_hi"] = mid, lo, hi
            row[key+"_err"] = (hi-lo)/2 if usable else None
        for key in bins:
            bins[key].append(row.get(key))
    bins["r_edges_arcmin"] = edges.tolist()
    # Descriptive map only; this is not another dispersion estimator.
    nr = min(len(groups), 8)
    re = np.unique(np.quantile(radius, np.linspace(0, 1, nr+1)))
    theta = np.arctan2(np.sin(np.radians(m.ra.to_numpy()-ra0))*np.cos(np.radians(m.dec.to_numpy())),
                       np.cos(np.radians(dec0))*np.sin(np.radians(m.dec.to_numpy()))-
                       np.sin(np.radians(dec0))*np.cos(np.radians(m.dec.to_numpy()))*np.cos(np.radians(m.ra.to_numpy()-ra0)))
    te = np.linspace(-np.pi, np.pi, 13)
    ri = np.clip(np.searchsorted(re, radius, side="right")-1, 0, len(re)-2)
    ti = np.clip(np.searchsorted(te, theta, side="right")-1, 0, 11)
    vmap = [[float(np.mean(u[(ri == i) & (ti == j), 1])*kd) if np.sum((ri == i) & (ti == j)) >= 10 else None for j in range(12)] for i in range(len(re)-1)]
    origin = {str(k): int(v) for k, v in m.catalog_origin.value_counts().items()} if "catalog_origin" in m else {"gaia_dr3": len(m)}
    finite = [i for i, v in enumerate(bins["v_rot"]) if v is not None]
    peak = max(finite, key=lambda i: abs(bins["v_rot"][i])) if finite else None
    signs = [np.sign(v) for v, detected in zip(bins["v_rot"], bins["rotation_detected"]) if detected]
    return dict(model_version=VERSION, cluster=cid, n_members=len(m), dist_kpc=distance, feh=feh, conc=conc,
                r_c_arcmin=rc, r_h_arcmin=rh, r_t_arcmin=rt, lit_rotation=litrot, vlos_lit=vlos,
                pm_sys=dict(pmra=float(bulk[0]), pmdec=float(bulk[1])),
                binning=dict(mode=mode, count=len(groups), sparse_bins=sum(len(g)<30 for g in groups), map_rings=len(re)-1,
                             rule="Errors-only pilot: reference 0.2 mas/yr, target information 100, minimum 200 stars; merge cap 0.3 dex in 1+R/Rh. Provisional, not an optimized physical bin count."),
                diagnostics=dict(estimator="Joint radial/tangential Gaussian likelihood with emcee",
                                 kinematic_catalog="Gaia DR3 + FPR" if origin.get("gaia_fpr") else "Gaia DR3", origin_counts=origin,
                                 rejected_covariance=int((members & ~valid).sum()),
                                 diagonal_covariance_stars=int((members & valid & missing_correlation).sum()),
                                 unresolved_bins=sum(s != "resolved" for s in bins["status"]), bin_fits=fits,
                                 warning="Conditional on membership and reported errors. Missing PM correlations use an explicitly approximate diagonal input covariance. Spatial systematics, selection truncation, distance and bulk-motion uncertainty are not marginalized. FPR membership scores are not calibrated DR3 probabilities. Auto-bin thresholds are provisional; no aperture or central V/sigma is inferred."),
                bins=bins, map=dict(r_edges_arcmin=re.tolist(), theta_edges_deg=np.degrees(te).tolist(), v_tan_mean=vmap),
                stats=dict(v_peak_kms=abs(bins["v_rot"][peak]) if peak is not None else None,
                           v_peak_err_kms=bins["v_rot_err"][peak] if peak is not None else None,
                           r_peak_arcmin=bins["r_mid_arcmin"][peak] if peak is not None else None,
                           sigma0_kms=None, vsig_peak=None, vsig_peak_err=None, sigma0_snr=None,
                           median_err_vel_kms=float(np.sqrt(np.median(np.trace(cov, axis1=1, axis2=2)/2))*kd),
                           p_cut_used=.9, rotation_detection_p=None, rotation_detected=any(bins["rotation_detected"]),
                           quality="conditional", dipole_diag_kms=None, counter_rotation=bool(1 in signs and -1 in signs),
                           rotation_sense="See signed rotation profile"),
                axis3d=dict(available=False, reason="Legacy 3D estimator is not run with the new PM model"))


def cached_analyze(cid, catalog, parameters, mode="auto", count=None, max_steps=4000):
    if max_steps not in (4000, 8000):
        raise ValueError("Sampling limit must be 4000 or 8000 steps")
    catalog = Path(catalog)
    stat = catalog.stat()
    identity = [VERSION, str(catalog.resolve()), stat.st_size, stat.st_mtime_ns, parameters, mode, count if mode != "auto" else None]
    if max_steps != 4000:
        identity.append(max_steps)
    key = hashlib.sha256(json.dumps(identity).encode()).hexdigest()
    directory = Path(__file__).resolve().parent / ".kinematics_cache"
    destination = directory / (key+".json")
    with LOCK:
        if destination.exists():
            return json.loads(destination.read_text(encoding="utf-8"))
        result = analyze(cid, catalog, parameters, mode, count, max_steps=max_steps)
        directory.mkdir(exist_ok=True)
        temp = destination.with_suffix(".tmp")
        temp.write_text(json.dumps(result, allow_nan=False), encoding="utf-8")
        temp.replace(destination)
        return result


def save_profile_figure(result, path):
    import matplotlib.pyplot as plt
    b = result["bins"]
    fig, axes = plt.subplots(3, 1, figsize=(8, 10), sharex=True)
    x = b["r_over_rh"]
    for ax, keys, label in zip(axes, [("v_rot",), ("sigma_r", "sigma_t"), ("ratio",)],
                               ["Signed Vt (km/s)", "Intrinsic dispersion (km/s)", "PM |Vt| / sigma_1D"]):
        for key in keys:
            values = np.array([np.nan if v is None else v for v in b[key]])
            low = np.array([np.nan if v is None else v for v in b[key+"_lo"]])
            high = np.array([np.nan if v is None else v for v in b[key+"_hi"]])
            ax.errorbar(x, values, yerr=[values-low, high-values], fmt="o-", label=key)
        ax.set_ylabel(label)
        ax.legend()
    axes[-1].set_xlabel("R / projected half-light radius")
    fig.suptitle(result["cluster"]+" — conditional PM MCMC; gaps are unavailable estimates")
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)
