# Globular Clusters Analysis System

A Gaia DR3-based system for studying the physical and dynamical properties of globular clusters (GCs), including membership filtering, rotation analysis, tidal-stripping studies, visualization, and cross-cluster comparisons.

---

## Globular Cluster Analysis

This project focuses on the analysis of key parameters and dynamical properties of **globular clusters (GCs)** using data from **ESA Gaia Data Release 3 (Gaia DR3)**.

The project currently includes **24 pre-filtered globular clusters**, with their corresponding datasets stored in `./data`. The analysis pipeline is designed to be flexible, allowing both input data and filtering methods to be modified directly.

### Features

#### Data Filtering Pipeline

The project includes a membership-filtering pipeline for processing Gaia DR3 data and identifying suitable stellar members of globular clusters.

The filtering criteria can be modified according to the requirements of different analyses. After modifying the data or filtering methods, the project's cache should automatically load the updated results.

#### Cluster Motion Analysis

One of the main focuses of this project is the **kinematic properties of globular clusters**.

The primary parameter currently used is the **\(v/\sigma\) ratio**, which investigates the relative importance of ordered rotational motion and random stellar motion within a cluster.

Further kinematic analyses are planned for future development.

#### 2D and 3D Visualization

The project provides both **2D graphs and 3D models** for visualizing globular clusters.

Various parameters can be adjusted directly, allowing users to explore how different assumptions and parameters affect the resulting representations.

#### Globular Cluster Comparisons

The project provides comparisons between different globular clusters, allowing their physical and kinematic properties to be examined side by side.

Detailed comparisons and related analyses can be found in the project documentation.

---

## Data Source

All astronomical data used in this project are based on **ESA Gaia Data Release 3 (Gaia DR3)**.

Gaia DR3 is the third major public data release from the **European Space Agency's Gaia mission**, released on **13 June 2022**.

Gaia's primary goal is to create a highly precise three-dimensional map of the Milky Way by measuring properties including:

* Stellar positions
* Parallax and distances
* Proper motions
* Radial velocities
* Brightness
* Other astrophysical properties

These measurements provide the foundation for studying the spatial distribution and kinematics of stars within globular clusters.

---

# Getting Started

A little note here, if you are trying to edit the code, make sure you run:
```text
pip install -r requirements.txt
```

## Starting the App

### Windows

Double-click:

```text
Start App - Windows.bat
```

### macOS

Run:

```text
Start App - Mac.command
```

The application will open in a browser tab.

Select a cluster, load its data, and explore the available analyses and visualizations.

### Main Views

The top navigation bar provides access to:

* **V/sigma Study** — Cross-cluster rotation and kinematic analysis
* **Tidal Study** — Tidal-stripping analysis for Pal 5, M92, and NGC 288
* **Filter Process** — Visualization and explanation of how cluster membership was determined

---

# Project Structure

```text
.
├── data/
│   ├── NGC5139/
│   │   ├── filtered.csv
│   │   ├── rotation.json
│   │   ├── rotation.png
│   │   ├── sky_map_2d.png
│   │   ├── cmd.png
│   │   ├── map_3d.png
│   │   ├── gaia_rv.csv
│   │   └── ...
│   │
│   ├── <other clusters>/
│   ├── vsigma_summary.json
│   └── vsigma_literature.json
│
├── membership_filter/
│   ├── membership_filter.py
│   ├── gc_target_list.xlsx
│   ├── run_summary.csv
│   └── archive/
│
├── vsigma_study/
│   ├── VSigma_Study_Explained.docx
│   ├── METHODS.md
│   ├── vsigma_pipeline.py
│   ├── vsigma_summary.csv
│   ├── vsigma_summary.json
│   ├── vsigma_ranking.png
│   ├── vsigma_correlations.png
│   ├── multi_axis_report.md
│   ├── literature_claims.json
│   └── per_cluster_plots/
│
├── tidal_study/
│   ├── Tidal_Stripping_Explained.docx
│   ├── tidal_pipeline.py
│   ├── tidal_results.json
│   ├── literature_refs.json
│   ├── gaia_cache/
│   ├── make_explainer_docx.py
│   └── export_app_data.py
│
├── filter_study/
│   ├── make_filter_figures.py
│   └── filter_study_summary.json
│
├── server.py
└── cluster_analysis.html
```

## `data/`

Contains the data used by the application.

Each globular cluster has its own directory. For example:

```text
data/NGC5139/
```

`NGC5139` corresponds to **Omega Centauri**.

Typical files include:

| File             | Description                                              |
| ---------------- | -------------------------------------------------------- |
| `filtered.csv`   | Filtered stellar catalogue with membership probabilities |
| `rotation.json`  | Rotation and \(v/\sigma\) analysis results               |
| `rotation.png`   | Rotation analysis plot                                   |
| `sky_map_2d.png` | 2D sky map                                               |
| `cmd.png`        | Colour-magnitude diagram                                 |
| `map_3d.png`     | 3D cluster visualization                                 |
| `gaia_rv.csv`    | Gaia radial velocities used for 3D spin-axis analysis    |

`NGC5139` provides two selectable catalog versions:

| File | Contents |
|------|----------|
| `filtered_dr3.csv` | Gaia DR3 filtered catalog with photometry and CMD support |
| `filtered_dr3_fpr.csv` | Gaia DR3 plus Gaia FPR crowded-field extension |

The FPR extension improves coverage in the crowded core, but its rows contain astrometry and proper motions without calibrated BP/RP colours. Therefore the CMD view is unavailable when this version is selected; use the DR3 version for colour–magnitude analysis. The large FPR file is tracked with Git LFS.

`NGC5139` also contains special **MODELED** and **RAW** 3D datasets and their corresponding documentation.

The `data/` directory also contains project-wide files:

```text
vsigma_summary.json
vsigma_literature.json
```

---

# Analysis Modules

## Membership Filter

Located in:

```text
membership_filter/
```

The membership-filtering program identifies likely cluster members from Gaia data using cluster-specific strategies based on:

* Position
* Proper motion
* Colour-magnitude diagrams (CMD)

### Main Files

```text
membership_filter.py
```

Main membership-filtering program.

```text
gc_target_list.xlsx
```

The target list containing the 24 globular clusters.

```text
run_summary.csv
```

Summary of the number of selected members for each cluster.

```text
archive/
```

Older outputs retained for reference.

---

## V/Sigma Study

Located in:

```text
vsigma_study/
```

This module investigates the rotational properties and kinematics of the 24 globular clusters.

### Documentation

```text
VSigma_Study_Explained.docx
```

A plain-language introduction to the analysis.

```text
METHODS.md
```

Detailed technical methodology.

### Analysis

```text
vsigma_pipeline.py
```

Main \(v/\sigma\) analysis pipeline.

```text
vsigma_summary.csv
vsigma_summary.json
```

Results for all 24 clusters.

```text
vsigma_ranking.png
```

Ranking of clusters by \(v/\sigma\).

```text
vsigma_correlations.png
```

Correlations between \(v/\sigma\) and other cluster parameters.

```text
multi_axis_report.md
```

Investigation into nested or multiple rotation axes.

```text
literature_claims.json
```

Comparison between claims from previous literature and the results obtained in this project.

```text
per_cluster_plots/
```

Contains individual rotation plots for each cluster.

---

## Tidal Study

Located in:

```text
tidal_study/
```

This module investigates **tidal stripping** in:

* Palomar 5
* M92
* NGC 288

### Main Files

```text
Tidal_Stripping_Explained.docx
```

Plain-language explanation of the tidal-stripping study.

```text
tidal_pipeline.py
```

Main tidal-feature search and analysis program.

```text
tidal_results.json
```

Results for the three clusters.

```text
literature_refs.json
```

Published references used for comparison.

```text
gaia_cache/
```

Cached wide-field Gaia data.

```text
make_explainer_docx.py
```

Regenerates the explanatory Word document from the analysis results.

```text
export_app_data.py
```

Rebuilds the interactive tidal-study data used by the application.

---

## Filter Study

Located in:

```text
filter_study/
```

This module generates figures and summaries describing the membership-filtering process.

```text
make_filter_figures.py
```

Rebuilds the per-cluster filtering figures and generates:

```text
data/filter_study.json
```

```text
filter_study_summary.json
```

Contains the filtering results in a readable format.

Each cluster's filtering visualization is stored as:

```text
data/<cluster>/filter_process.png
```

---

# Running the Analysis

All major analysis modules can be rerun independently.

### Membership Filtering

```bash
python membership_filter/membership_filter.py
```

### V/Sigma Analysis

```bash
python vsigma_study/vsigma_pipeline.py
```

### Tidal Study

```bash
python tidal_study/tidal_pipeline.py
```

### Filter Figures

```bash
python filter_study/make_filter_figures.py
```

The pipelines reuse cached downloads whenever possible.

Deleting a relevant cache file will force the program to perform a fresh Gaia query.

---

# Future Development

Future work will include additional investigations into the physical and dynamical properties of globular clusters, including further kinematic analyses and additional cluster parameters.

---

# License

This project is licensed under the **MIT License**. See the [`LICENSE`](LICENSE) file for details.

Gaia data are distributed under the **CC BY-NC 3.0 IGO license**. Please refer to the [Gaia Data Use and License](https://www.cosmos.esa.int/web/gaia-users/license) page for the applicable terms.

---

# Contact

For questions, suggestions, or inquiries regarding this project, please contact:

**Discord:** `zz.oao`

## PM rotation and dispersion (MCMC)

The Rotation view defaults to **Auto (per cluster)**. Initial equal-count rings are merged using reported error information and a radial-width cap. The applied count is displayed. These are provisional, reproducible settings, not an empirically optimal bin count. Manual **Equal number** and **Equal radius width** modes remain available (2–100 bins).

The new estimator fits mean radial and signed tangential motion and their intrinsic dispersions together, using each star's PM error covariance. Geometry uses spherical local bases and a projected bulk space velocity, including perspective contraction/expansion. Fits run in mas/yr with `emcee`; plots convert to km/s and show radius in units of the catalog projected half-light radius.

The displayed local ratio is `abs(Vt) / sqrt((sigma_R² + sigma_t²)/2)`, with posterior intervals. It is **not** Bianchini's peak-PM/central-LOS statistic. Central/aperture headline ratios and the old 3D fit are not produced by the new model. Existing cross-cluster studies remain explicitly labeled legacy; interactive requests do not overwrite them or catalogs.

Allow several minutes for a first calculation. If convergence gaps remain, enable **Longer sampling** and update the profile (up to 8000 instead of 4000 steps; this does not guarantee resolution). Results are cached under ignored `vsigma_study/.kinematics_cache/`, keyed by model version, catalog path/size/modification time, parameters and bin settings. Changed input data invalidates this cache. Restart `server.py` after updating the code. Install dependencies with `python -m pip install -r requirements.txt`.

### Scientific limitations

DR3 and FPR are included when selected, with finite positive uncertainties. Valid reported PM correlation is used; missing correlation in legacy CSVs uses a **diagonal input covariance approximation**, whose star count is shown. Invalid nonmissing correlations are excluded. Future DR3 downloads retain the correlation column; existing CSVs/caches are not silently rewritten. No lower-membership fallback is applied.

The fit is **conditional on the selected sample and reported errors**: selection truncation, star-to-star spatial systematics, distance and bulk-motion uncertainty are not marginalized. FPR membership scores are not calibrated DR3 probabilities. A likelihood cannot undo biased selection.

Unconverged fits, error-model conflicts and unresolved dispersions are distinguished. Unresolved, converged, model-compatible bins show a conditional 95% upper limit; no ordinary ratio is shown there. Error bars are 16th–84th posterior percentiles, not observational-error bars. Rotation flags are conditional posterior-sign checks, not systematics-corrected detections.

See [the method specification](vsigma_study/METHODS_MCMC.md) for priors, bin rules and validation. Run scientific regression checks with `python -m unittest vsigma_study.test_pm_kinematics -v`.
