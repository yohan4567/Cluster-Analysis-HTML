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

`NGC5139` additionally contains special **MODELED** and **RAW** 3D datasets and their corresponding documentation.

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
