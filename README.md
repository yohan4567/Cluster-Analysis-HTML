# Globular Clusters Analysis System

## Globular Cluster Analysis

This project focuses on the analysis of key parameters and dynamical properties of **globular clusters (GCs)** using data from **ESA Gaia Data Release 3 (Gaia DR3)**.

The project currently includes **24 pre-filtered globular clusters**, with their corresponding datasets stored in `./data`. The analysis pipeline is designed to be flexible, allowing both the input data and filtering methods to be modified directly.

### Features

#### Data Filtering Pipeline

The project includes a filtering pipeline for processing Gaia DR3 data and selecting suitable stellar members of globular clusters.

The filtering criteria can be modified directly according to the requirements of different analyses. After modifying the data or filtering methods, the project's cache should automatically update when the program is run.

#### Cluster Motion Analysis

One of the main focuses of this project is studying the **kinematic properties of globular clusters**.

Currently, the primary parameter used in our analysis is the **\(v/\sigma\) ratio**, which can be used to investigate the relative importance of ordered rotation and random stellar motions within a cluster.

Further kinematic analyses are planned for future development.

#### 2D and 3D Visualization

The project provides both **2D graphs and 3D models** for visualizing globular clusters.

Various parameters can be adjusted directly, allowing users to explore how different assumptions and parameters affect the resulting representations of the clusters.

#### Globular Cluster Comparisons

The project also provides comparisons between different globular clusters, allowing their physical and kinematic properties to be examined side by side.

Detailed comparisons and related analyses can be found in the project documentation.

### Data Source

All data used in this project are based **solely on ESA Gaia Data Release 3 (Gaia DR3)**.

**Gaia Data Release 3** is the third major public data release from the **European Space Agency's Gaia mission**, released on **13 June 2022**.

Gaia's primary goal is to create a highly precise three-dimensional map of the Milky Way by measuring properties including:

* Stellar positions
* Parallax and distances
* Proper motions
* Radial velocities
* Brightness
* Other astrophysical properties

These measurements provide the foundation for studying the spatial distribution and kinematics of stars within globular clusters.

### Data

The project currently contains **24 pre-filtered globular clusters**.

The corresponding data files are located in:

```text
./data
```

The data can be modified directly if you wish to:

* Add additional globular clusters
* Remove existing clusters
* Modify or replace cluster data
* Apply different filtering criteria

The cache is designed to load automatically after changes are made.

### Future Development

Several additional analyses are planned for future development, including further investigations into the physical and dynamical properties of globular clusters.

### License

This project is licensed under the **MIT License**. See the [LICENSE](LICENSE) file for details.
Gaia data are distributed under the **CC BY-NC 3.0 IGO license**, please refer to [LICENSE](https://www.cosmos.esa.int/web/gaia-users/license) for more information.

### Contact

If you have any questions, suggestions, or inquiries regarding this project, please contact:

**Discord:** `zz.oao`


=====================================

TO START THE APP:  double-click "Start App - Windows.bat"
                   (on a Mac: "Start App - Mac.command")
It opens a browser tab. Pick a cluster, load its data, explore the views.
Top-bar buttons: "V/sigma study" (cross-cluster rotation study),
"Tidal study" (tidal-stripping results for Pal 5, M92, NGC 288), and
"Filter process" (how membership was decided, cluster by cluster).


FOLDERS
-------

data\
    What the app reads. One folder per globular cluster (NGC5139 = omega
    Centauri, etc.), each containing:
      filtered.csv        the star catalog with membership probabilities
      rotation.json/png   rotation & V/sigma analysis results and plot
      sky_map_2d.png, cmd.png, map_3d.png   auto-generated plots
      gaia_rv.csv         Gaia radial velocities (for the 3D spin axis)
    Plus two app-wide files: vsigma_summary.json, vsigma_literature.json.
    NGC5139 also holds the two special 3D files (MODELED vs RAW) + README.

membership_filter\
    membership_filter.py    the program that finds cluster members in Gaia
                            data (per-cluster strategies, PM+position+CMD)
    gc_target_list.xlsx     the 24-cluster target list it reads
    run_summary.csv         how many members each cluster got
    archive\                older outputs kept for reference

vsigma_study\
    VSigma_Study_Explained.docx   START HERE - plain-language write-up
    METHODS.md                    full technical methodology
    vsigma_pipeline.py            the analysis program (rerun anytime)
    vsigma_summary.csv/.json      results table, all 24 clusters
    vsigma_ranking.png            clusters ranked by V/sigma
    vsigma_correlations.png       V/sigma vs cluster parameters
    multi_axis_report.md          the nested-rotation-axis investigation
    literature_claims.json        what past papers claimed vs what we found
    per_cluster_plots\            one rotation figure per cluster

tidal_study\
    Tidal_Stripping_Explained.docx   START HERE - plain-language write-up
    tidal_pipeline.py                the search program (rerun anytime)
    tidal_results.json               numbers for all three clusters
    literature_refs.json             the published papers we compared against
    gaia_cache\                      downloaded wide-field Gaia data
    make_explainer_docx.py           regenerates the Word file from results
    export_app_data.py               rebuilds the app's interactive tidal

filter_study\
    make_filter_figures.py       rebuilds the per-cluster filter-process
                                 figures + data\filter_study.json
    filter_study_summary.json    the same numbers, in readable form
    (each cluster's figure is saved as filter_process.png inside its own
     folder under data\)

server.py / cluster_analysis.html
    The app itself (local web server + the page it serves).


RERUNNING THINGS
----------------
Membership filter:  python membership_filter\membership_filter.py
V/sigma study:      python vsigma_study\vsigma_pipeline.py
Tidal search:       python tidal_study\tidal_pipeline.py
Filter figures:    python filter_study\make_filter_figures.py
(All reuse cached downloads where possible; deleting a cache file forces
a fresh Gaia query.)
