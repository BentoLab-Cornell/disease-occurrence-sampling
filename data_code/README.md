# Population exposure and background sampling: replication package

This package reproduces the sensitivity analyses accompanying the report on Gibb et al. (2026). It contains prepared data, R code and saved coefficient estimates. Reproduction starts from the supplied covariate matrices and spatial meshes; downloading rasters and repeating GIS extraction are outside this package.

## Run

Unzip the package, open a terminal in `Gibb_sensitivity_replication`, and run:

```sh
Rscript code/install.R
Rscript run.R
```

This writes nine figures (PDF and PNG), coefficient tables and comparison summaries to `output/saved`. It uses the saved estimates and does not fit models. Figure numbers match the report: synthetic null (1), dengue background comparison (2), sampling choices (3), reconstruction (4), predictor omission (5), country omissions (6), additional draws (7), all dengue slopes (8), and prior/mesh sensitivity (9).

To refit the models from the supplied inputs:

```sh
Rscript code/install.R refit
Rscript run.R refit all
Rscript run.R results refitted
```

Refits are saved individually under `output/refits`; regenerated figures and tables go to `output/refitted`. A full refit takes many hours. Jobs run sequentially, and repeating the command resumes completed jobs. The final command requires all refits and never substitutes saved estimates for missing fits. Delete a job's output CSV before intentionally rerunning it after a code or input change.

For a single analysis, replace `all` with `dengue`, `prior`, `synthetic`, `west_nile` or `background`. To list jobs and run one example:

```sh
Rscript run.R jobs dengue
Rscript run.R refit dengue sampling__980019_area_log_positive
```

Tested on Linux with R 4.3.3, INLA 23.03.26, Matrix 1.6-5 and ggplot2 3.4.4. The installer obtains the specified INLA release from its official archive. Model refitting requires that INLA version; figures require ggplot2 only. Small numerical differences between machines can occur. Session information is written with the outputs.

## Files and analyses

| Input | Contents | Code |
|---|---|---|
| `data/dengue.rds` | Archived and reconstructed covariates, 48 sampled backgrounds, country exclusion flags, standardization and fixed mesh | `code/dengue.R` |
| `data/dengue_jobs.csv` | The 120 dengue fits, with sampling seed, arm and intervention | `code/dengue.R` |
| `data/synthetic.rds` | Simulated landscape and generating/model settings | `code/synthetic.R` |
| `data/west_nile.rds` | County data, archived background covariates, population, calibrated rates and occurrence/incidence model settings | `code/west_nile.R` |
| `data/west_nile_background.rds` | Archived occurrence data, finite background frame and simulated null cases | `code/background.R` |
| `data/results/*.csv` | Saved posterior means, standard deviations and marginal 95% credible intervals | `code/summarize.R`, `code/figures.R` |

RDS files are named R lists, readable with `readRDS()`. For example, `str(readRDS("data/dengue.rds"), max.level=1)` lists the data components. Dengue coordinates are longitude/latitude; `presence` distinguishes cases (1) from backgrounds (0). Covariate names and the archived means/standard deviations are in `model$vars` and `model$scaling`. Healthcare travel time is `log(1 + minutes)` before standardization. Missing covariates retain the original model's handling.

The dengue sampling arms have weights proportional to:

| Arm | Weight | Eligible cells |
|---|---|---|
| A | Cell area × log(1 + population), replacing zero log weights with 1.06840784199454 | Original sampling frame |
| B | Cell area × log(1 + population) | Inhabited cells |
| C | Cell area × population | Inhabited cells |
| D | Population | Inhabited cells |

The original four seeds are 980019, 980027, 980035 and 980043; the 16 additional A/D pairs use 980001–980016. Each global fit uses 11,313 cases and 12,108 backgrounds. Country omissions remove the flagged cases and backgrounds without redrawing locations. C minus B isolates the log-to-linear change; D minus A compares the complete recipes. Paired changes compare posterior means. Across-draw ranges and Monte Carlo standard errors describe background-sampling variability, separately from each model's credible interval.

## Sources

The underlying study is [Gibb et al., *The anthropogenic fingerprint on emerging infectious diseases*](https://doi.org/10.1038/s41586-026-11058-6). Archived inputs and model specifications derive from the [authors' release](https://doi.org/10.5281/zenodo.20738794) and [code version 1.0.0](https://github.com/rorygibb/fingerprint_ms/tree/v1.0.0).

The dengue matrices are a reconstruction. Five predictors remain unmatched to the authors' original layers: forest loss, precipitation change, protected areas, social vulnerability and temperature change. Substitutions use Hansen 2021 v1.9 forest loss, ERA5-Land climate, October 2026 WDPCA protected areas and a GRDI social-vulnerability layer whose historical edition is unconfirmed. Archived values are supplied for the original-location comparisons. Country flags derive from GADM 4.1 boundaries; boundary source files are not redistributed here. The report gives the scientific interpretation and source comparisons.
