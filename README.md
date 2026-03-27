[README.md](https://github.com/user-attachments/files/26305468/README.md)
[README.md](https://github.com/user-attachments/files/26305468/README.md)
# Expanding Freshwater Biodiversity Models for Species Abundance Target Delivery

This repository contains the code and data underlying the project report *Expanding Freshwater Biodiversity Models for Species Abundance Target Delivery*, produced for the Department for Environment, Food & Rural Affairs (Defra).

The project models the abundance of 266 freshwater taxa (235 invertebrates, 31 fish) that contribute to England's legally binding species abundance targets under the Environment Act 2021. It establishes workflows to project how Defra's statutory and non-statutory water environment targets could contribute to halting and reversing species declines by 2030 and 2042.

---

## Repository structure

```
.
├── paths.R                        # Central path definitions — source this in every script
├── bio_data_prep.R                # (1) Biological data preparation
├── wq_data_prep.R                 # (2) Water quality data preparation
├── site_matching.R                # (3) Site matching (biological, WQ, discharge sites)
├── prep_mod_args.R                # (4) Prepare model argument table (water quality GAMMs)
├── prep_inla_mod_args.R           # (5) Prepare model argument table (INLA species models)
├── fit_mod.R                      # (5) Fit monthly water quality GAMMs
├── fit_mod_batch.sh               # Batch script for fit_mod.R on HPC
├── problem_mods.R                 # (5b) Re-fit any failed water quality models
├── env_data_prep.R                # (6) Environmental covariate preparation
├── historical_pred_data_prep.R    # (7) Data preparation for hindcasting
├── future_rainfall.R              # (8) Derive future precipitation anomaly rasters
├── future_pred_data_prep.R        # (9) Data preparation for future projections
├── prep_future_pred_args.R        # (10) Build future scenario prediction argument grids
├── inla_mod_final_v3.R            # (11) Fit INLA species abundance models
├── inla_mod_batch_v3.sh           # Batch script for inla_mod_final_v3.R on HPC
├── inla_posteriors.R              # (12) Draw posterior samples from fitted models
├── inla_mod_diagnostics.R         # (13) Model diagnostics and summary plots
├── inla_historical_pred_v2.R      # (14) Hindcast predictions (2003–2023)
├── inla_future_pred.R             # (15) National-scale future projections
├── inla_future_pred_batch.sh      # Batch script for inla_future_pred.R on HPC
├── inla_spatial_pred.R            # (16) Spatially explicit projections
├── analysis.R                     # (17) Aggregate results, compute probabilities, figures
└── data/
    ├── bundled/                   # Small open datasets committed to this repository
    ├── external/                  # Large or licensed datasets — download separately (see below)
    │   └── README.md              # Download instructions for each external file
    └── processed/                 # Contains three QGIS-generated snapped site shapefiles
                                   # (meta_snap, ea_wq_snap, stws_snap) and the RICT PCA
                                   # raster (rict_pca_axes.tif) — all committed to the
                                   # repository. Large intermediate outputs written here
                                   # by the pipeline scripts are gitignored.
```

---

## Pipeline overview

The scripts should be run in the numbered order above. A brief description of each stage is given below.

| Step | Script | Description |
|------|--------|-------------|
| 1 | `bio_data_prep.R` | Processes EA Fish & Ecology Data Explorer records into harmonised invertebrate and fish abundance matrices and metadata |
| 2 | `wq_data_prep.R` | Cleans and subsets EA Water Quality Archive records for eight physicochemical determinands |
| 3 | `site_matching.R` | Snaps biological, water quality, and consented discharge sites to the OS Open Rivers network and matches them to river sub-segments |
| 4 | `prep_mod_args.R` | Builds a model argument table assigning a task ID to each water quality monitoring area for parallelised GAMM fitting |
| 5 | `fit_mod.R` | Fits generalised additive mixed models (GAMMs) to estimate monthly sub-segment-level water quality at biological survey locations |
| 6 | `env_data_prep.R` | Assembles all environmental covariates (RICT PCA scores, rainfall anomalies, water quality predictions, sewage/mining/habitat/barrier metrics) and assigns cross-validation folds |
| 7 | `historical_pred_data_prep.R` | Prepares scaled predictor tables for hindcasting (2003–2023) |
| 8 | `future_rainfall.R` | Derives future precipitation anomaly rasters for 2030 and 2042 from CHESS-SCAPE projections |
| 9 | `future_pred_data_prep.R` | Prepares predictor tables for future projections (2030, 2042) under three RCPs |
| 10 | `prep_future_pred_args.R` | Constructs scenario argument grids covering all combinations of pressure interventions and RCPs |
| 11 | `prep_inla_mod_args.R` | Builds a model argument table assigning a task ID to each species for parallelised INLA model fitting; creates required output directories |
| 12 | `inla_mod_final_v3.R` | Fits INLA species abundance models for each taxon under "no metals" and "with metals" covariate configurations, with five-fold cross-validation |
| 13 | `inla_posteriors.R` | Draws 200 posterior samples from each fitted model for use in projections |
| 14 | `inla_mod_diagnostics.R` | Produces diagnostic plots and summary statistics for all fitted models |
| 15 | `inla_historical_pred_v2.R` | Generates hindcast predictions and computes geometric mean relative abundance indices (2003–2023) |
| 16 | `inla_future_pred.R` | Generates national-scale future projections across all species and scenarios |
| 17 | `inla_spatial_pred.R` | Generates spatially explicit projections per biological monitoring site |
| 18 | `analysis.R` | Aggregates projections, computes four scenario probabilities, and produces figures |

Steps 5, 12 and 16 are computationally intensive and are designed to run as parallelised batch jobs on an HPC cluster using the accompanying `.sh` scripts.

---

## Getting started

### 1. Clone the repository

```bash
git clone https://github.com/wilkesma/DefraFSAP/tree/main
cd freshwater-biodiversity-models
```

### 2. Install R packages

This project uses `renv` to manage package versions. After opening the project in R or RStudio:

```r
install.packages("renv")
renv::restore()
```

R-INLA is not on CRAN and must be installed separately before running `renv::restore()`:

```r
install.packages("INLA",
  repos = c(getOption("repos"),
            INLA = "https://inla.r-inla-download.org/R/stable"),
  dep = TRUE)
```

### 3. Download external data

Several datasets are too large or subject to licensing restrictions that prevent them being committed to this repository. Instructions for downloading each file are in [`data/external/README.md`](data/external/README.md). Place all downloaded files in `data/external/` using the exact filenames specified there.

### 4. Run the scripts

`paths.R` is sourced at the top of each script — no manual path editing should be required. Run scripts in the numbered order in the table above. All outputs generated by the scripts are written to `data/processed/`.

---

## Data availability

Datasets committed to this repository (`data/bundled/`) are reproduced under open licences (Open Government Licence v3.0 or equivalent). Full download and licensing details for external datasets are provided in [`data/external/README.md`](data/external/README.md).

The `data/processed/` folder contains four small files committed to the repository: the three QGIS-generated snapped site shapefiles (`meta_snap`, `ea_wq_snap`, `stws_snap`) and the RICT PCA raster (`rict_pca_axes.tif`). All other intermediate outputs generated by running the pipeline (model objects, prediction arrays, etc.) are written to `data/processed/` but are gitignored due to their size; they will be recreated by running the scripts.

---

## Citation

If you use this code or data, please cite the project report:

> Wilkes, M.A. et al. (2025). Expanding Freshwater Biodiversity Models for Species Abundance Target Delivery. Report to Defra.

---

## Licence

The code in this repository is released under the MIT Licence (see `LICENSE`).

Data included in `data/bundled/` are subject to their original licences (e.g. Open Government Licence v3.0). External datasets are not redistributed and must be obtained from their original providers under their respective terms.
