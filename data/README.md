# Data directory

Each cluster directory contains the canonical `filtered.csv` catalog plus analysis outputs and plots used by the app.

## Omega Centauri (NGC5139)

`NGC5139` has two selectable catalogs:

- `filtered_dr3.csv`: Gaia DR3 filtered catalog.
- `filtered_dr3_fpr.csv`: Gaia DR3 combined with the Gaia FPR crowded-field extension (tracked with Git LFS).
- `filtered.csv`: backward-compatible canonical DR3 copy used by existing analysis scripts.

The FPR extension has astrometry and proper motions but no calibrated BP/RP colours, so CMD analysis uses the DR3 catalog.

## Project-wide outputs

The JSON files directly under this folder are generated summaries consumed by the app:

- `vsigma_summary.json` and `vsigma_literature.json`
- `tidal_study.json`
- `filter_study.json`

Generated caches and plots remain next to the cluster catalog that produced them so the server can locate them automatically.
