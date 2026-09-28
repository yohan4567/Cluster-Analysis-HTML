# Conditional proper-motion kinematics

The interactive estimator is versioned in `pm_kinematics.py`. It supersedes the
moment-subtraction profile; old cross-cluster summary files are not comparable.
This is an exploratory model conditional on the selected sample, not a fully
selection-corrected or spatial-systematics-marginalized inference.

## Coordinates and data

Use membership_prob >= 0.9 without a lower-threshold fallback. Keep finite
positions, proper motions and positive PM errors. Valid reported correlations
must be strictly between -1 and 1. Legacy exports without correlations use zero
input cross-covariance **as an explicitly reported approximation**; they are not
claimed to have complete measured covariance. Future DR3 downloads retain
pmra_pmdec_corr. Existing raw caches are not automatically redownloaded.

Fit a constant 3D bulk velocity's projection into each star's local east/north
basis (catalog distance and LOS speed fixed, two transverse components fitted
by iteratively clipped least squares). Subtract that projected velocity. This
includes perspective contraction/expansion and the changing local sky basis.
Transform the residual vector and covariance into the great-circle radial and
perpendicular tangential directions. Tangential positive means E to N on an
E-right sky map. This assumes common distance; per-star depth is not inferred.
Bulk-motion uncertainty is not propagated and can couple radial bins.

## Likelihood and prior

Within each ring, assume a constant mean (mu_R, mu_t) and diagonal intrinsic
covariance diag(s_R²,s_t²). Add each star's transformed full measurement
covariance. Stars are assumed independent. Fit in mas/yr, converting to km/s
using 4.74047 D_kpc. Mean priors are independent N(0,5²) mas/yr; dispersion
priors are independent half-normal(scale=1 mas/yr), with support at zero.
Do not subtract measurement variance after fitting.

Use L-BFGS-B initialization and 20 emcee walkers, up to 4000 steps (8000 with the Longer sampling option). After at
least 1000 steps, check every 500: discard the first third, require retained
length >50 times every parameter's autocorrelation time, and <10% relative
change of all autocorrelation estimates. Report steps, autocorrelation times,
effective samples and acceptance fraction. These checks are not a proof of
convergence or protection from model misspecification. A bounded run can fail.

Importance-reweight posterior samples to half/double the dispersion prior scale;
flag a median combined-dispersion shift >0.5 original posterior standard
deviations. This is a local sensitivity diagnostic, not a replacement for
independent chains or a comprehensive prior study.

Report posterior 16th/50th/84th percentiles. A conservative **heuristic** gate
requires 2 Delta log likelihood >9 relative to both intrinsic dispersions zero
before reporting resolved combined dispersion. This boundary test is not
advertised as an exact chi-square significance or a universal detection limit.
If the observed residual quadratic form under measurement errors alone is
below the 0.1% chi-square quantile (2N-2 degrees of freedom), flag an error-model
conflict. It can indicate truncation, error overestimation or other misspecification;
it does not identify a unique cause. Compatible unresolved fits show a conditional
95% upper limit. Nonconverged/conflicting/prior-sensitive fits do not give ratios.

The component intervals remain posterior estimates: resolution of their average
does not establish a separate detection of each component. Mean-sign rotation
flags use posterior tail probability <0.003/(2 times bin count), conditional on
this model, not a calibrated global frequentist detection with systematics.

## Automatic binning

Construct min(24, floor(N/200)) equal-count seeds, at least one. A per-star pilot
information contribution is the average of [s_ref²/(s_ref²+error_component²)]²,
with fixed s_ref=0.2 mas/yr. Join adjacent seeds until at least 200 stars and
summed information >=100, unless the prospective width exceeds 0.3 dex in
log10(1+R/Rh). A final underfilled group can merge backward only within this cap.
Broad initial seeds and low-information groups remain flagged, not forced to
yield positive dispersion. Stars are never discarded by bin selection.

The pilot follows the information scaling of Gaussian variance estimation,
but ignores cross-covariance and does not know the true dispersion. Thresholds
are **provisional defaults**, with regression checks, not a calibrated optimum
for all 24 clusters. Rules use errors/positions, never the fitted rotation peak
or dispersion detection. Different catalogs may produce different boundaries;
compare their reported boundaries and matched samples before attributing any
change to FPR. A shared-edge comparison control is not implemented yet.

## Outputs and remaining scientific checks

Local PM ratio = abs(Vt)/sqrt((sigma_R²+sigma_t²)/2), computed sample by sample.
This differs from Bianchini et al.'s V_peak,PM/sigma_0,LOS. No central or
aperture statistic is inferred. The peak is descriptive and depends on binning.
The sector map is a descriptive mean, not a separate likelihood fit.

Remaining work before publication: selection-function/field-mixture modeling;
FPR-specific error and membership calibration; spatial systematic covariance;
bulk/distance nuisance parameters; within-ring gradients and position-angle
structure; matched brightness/radius and membership-threshold sensitivity;
independent-chain convergence and repeated-simulation interval coverage;
shared DR3/FPR bin edges; broader prior and bin-threshold validation.

Regression tests cover exact spherical bulk removal (including perspective),
pure rotation, unequal intrinsic dispersions with correlated heteroscedastic
errors, artificially truncated samples, and conservation of stars in auto/manual
bins. They are necessary implementation checks, not proof of astrophysical validity.

References: [Bianchini et al. 2018](https://arxiv.org/abs/1806.02580),
[Vasiliev 2019](https://arxiv.org/abs/1811.05345),
[Vasiliev & Baumgardt 2021](https://arxiv.org/abs/2102.09568).
