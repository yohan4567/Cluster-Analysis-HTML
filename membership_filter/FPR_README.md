# Omega Centauri: DR3 + FPR extension

Run from the `Cluster Analysis HTML` directory:

```powershell
python -m pip install -r requirements.txt
python membership_filter/membership_filter.py --include-fpr --only "NGC 5139" --fpr-mirror aip --output-folder membership_filter/fpr_run
```

The default mirror is ESA (`--fpr-mirror esa`); AIP hosts the same
`gaiafpr.crowded_field_source` table. Network access is required for the
first run. A validated local raw cache is reused on subsequent runs.
Interrupted/incomplete downloads are not accepted as valid caches.

Outputs in `membership_filter/fpr_run`:

- `NGC_5139_filtered_dr3_fpr_v1.csv`: quality-passing DR3 and FPR rows,
  including nonmembers, with `membership_prob` and `is_member`.
- `NGC_5139_filtered_dr3_fpr_v1.json`: source-specific member/core counts
  and processing notes.
- `NGC_5139_gaia_fpr_raw.csv` and `.json`: downloaded observations and query provenance.

The original `data/NGC5139/filtered.csv` is the DR3 baseline and is not
overwritten. Existing DR3 memberships are retained, not refitted. This is
an opt-in extension mode; ordinary DR3 filtering is unchanged. Other
clusters are unsupported by this FPR release. The app's active data and
previous rotation outputs are not replaced automatically.

## FPR selection

Finite astrometry and positive errors are required; each PM error < 1.5
mas/yr, G < 20.5, and foreground parallax excess < 3 sigma. FPR does not
publish RUWE, so no RUWE cut is applied to FPR rows (the output column is
empty).
FPR coordinates are linearly propagated from their catalog epoch to
J2016.0, retaining `catalog_ra`, `catalog_dec`, `catalog_ref_epoch`.

FPR has no reliable BP/RP colours. Those columns and `p_cmd` are left
empty. No King spatial prior is fitted to the FPR footprint. Instead an
error-aware PM window around (-3.250, -6.746) mas/yr, with intrinsic
dispersion 0.75 mas/yr and full supplied PM measurement covariance, gives
`exp(-z_squared/2)` (zero outside the existing 4-sigma gate). Scores >=
0.5 are flagged as candidate members. These conservative defaults keep
the project's magnitude/error limits rather than admitting all faint FPR
sources. Adjustments require a contamination/completeness assessment.

`membership_prob` is a compatibility column: for FPR it is an uncalibrated
PM consistency score, NOT a posterior probability or a quantity directly
comparable to the DR3 Bayesian-style membership estimate. `catalog_origin`
and `membership_method` identify the distinction. No kinematic validation
is claimed; `kinematics_validated` is false. Do not use the combined
catalog for precision rotation/dispersion without validating source-specific
selection and measurement errors.

Source IDs are preserved as strings; `source_key` includes the catalog.
ESA's FPR is an add-on excluding nominal sources. FPR stars within 0.1
arcsec of a retained DR3 source at the common epoch are flagged for review,
not automatically deleted (they could be resolved close companions).

## Verification

```powershell
cd membership_filter
python -m unittest test_fpr_extension -q
```

Tests cover missing colours, bad/missing errors, foreground rejection,
PM covariance, epoch handling, source IDs, and baseline preservation.
Real full-catalog download/end-to-end results still need verification.

Reference: https://doi.org/10.1051/0004-6361/202347203
