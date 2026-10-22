# Population exposure and sampling in disease occurrence models

R code, prepared data, and saved model estimates for studying how population exposure, background sampling, and outcome definition affect associations in spatial models of infectious disease occurrence.

The analyses combine simulations with dengue and West Nile virus examples. They examine associations under known generating assumptions and compare estimates across sampling and modeling choices in empirical data. The main questions concern the role of population size, the distribution used to select background locations, and the distinction between reported occurrence and incidence.

## Analyses

| Component | What it examines |
| --- | --- |
| Synthetic experiments | Logarithmic and population-proportional background sampling under rare and saturating occurrence. Null scenarios hold individual infection risk and reporting constant; positive controls introduce an environmental effect. Null scenarios also include a background distribution based on the generating occurrence probabilities. |
| Dengue sensitivity analyses | Changes in estimated associations across background sampling recipes, reconstructed covariates, predictor omissions, country exclusions, and repeated background draws. Additional fits examine prior and spatial mesh choices. |
| West Nile model comparisons | Occurrence and incidence models applied to simulated county counts under Poisson and negative-binomial generating processes. Comparisons include population offsets, link functions, background distributions, and equal-population controls. An incidence model is also fitted to observed counts. |
| West Nile background resampling | Logarithmic and linear population weighting within a finite set of archived background locations, using observed occurrence records and simulated null records. Models use the original predictor set or an extension including healthcare travel time. |

The simulations identify behavior under specified assumptions. Empirical sensitivity analyses assess how estimates change with analytical choices. Population-proportional sampling is a comparator whose interpretation depends on the generating process, particularly when occurrence probabilities saturate. These analyses do not directly estimate the effects of spatial or temporal exposure mismatch.

## Reproduce figures and tables

Run the following commands from the repository root, the directory containing `run.R`:

```sh
Rscript code/install.R
Rscript run.R
```

This uses the supplied coefficient estimates without refitting models. It writes:

- Nine figures in PDF and PNG format to `output/saved/figures/`.
- Coefficient tables and comparison summaries to `output/saved/tables/`.
- R session information to `output/saved/session.txt`.

The figures cover the synthetic experiments, dengue comparisons, and prior/mesh sensitivity. Tables include all analysis components, including both West Nile modules.

The explicit form of the default command is:

```sh
Rscript run.R results saved
```

## Prepare the dengue input

The dengue input is distributed as two binary parts: `dengue.rds.part1` and `dengue.rds.part2`. Place both files in the repository's `data/` directory and combine them to recreate `data/dengue.rds` before refitting models that use the dengue input.

On macOS or Linux, run this command from the repository root:

```sh
cat data/dengue.rds.part1 data/dengue.rds.part2 > data/dengue.rds
```

Alternatively, run the following in R from the repository root (works on Windows, macOS, and Linux):

```r
parts <- file.path("data", c("dengue.rds.part1", "dengue.rds.part2"))
stopifnot(all(file.exists(parts)))
local({
  output <- file("data/dengue.rds", open = "wb")
  on.exit(close(output))
  for (part in parts) {
    writeBin(readBin(part, what = "raw", n = file.info(part)$size), output)
  }
})
```

The order of the parts matters. Each part contains a portion of the original file and cannot be opened separately with `readRDS()`. After combining them, use `readRDS("data/dengue.rds")` as usual. These commands overwrite any existing `data/dengue.rds`.

## Refit models

To fit all models from the prepared inputs and then summarize the new estimates:

```sh
Rscript code/install.R refit
Rscript run.R refit all
Rscript run.R results refitted
```

There are 1,341 fitting jobs across five modules: `dengue`, `prior`, `synthetic`, `west_nile`, and `background`. Jobs run sequentially, with one INLA thread. Individual coefficient files are written to `output/refits/<module>/<job-id>.csv`.

To list jobs or run one dengue fit:

```sh
Rscript run.R jobs dengue
Rscript run.R refit dengue sampling__980019_area_log_positive
```

To refit a single module, replace `all` with its name. Repeating a refit command resumes work by validating and skipping existing job outputs. After changing code or inputs, remove the affected output CSVs before rerunning those jobs. The resume mechanism does not check whether inputs have changed.

`Rscript run.R results refitted` requires every expected refit across all modules. It stops if a fit is missing and does not substitute saved estimates. Keep `data/results/` in place when refitting: these files also define the expected jobs and coefficient structure.

### Software

The recorded analysis environment is Linux with R 4.3.3, INLA 23.03.26, Matrix 1.6-5, and ggplot2 3.4.4.

Generating figures and tables requires ggplot2. Refitting additionally requires INLA and its dependencies; the installer explicitly checks for Matrix, foreach, sp, and lifecycle. The refitting code requires INLA 23.03.26, which the installer obtains from the official archive.

Only the INLA version is enforced. Other missing packages are installed from CRAN, so the installer does not recreate a fully fixed package environment. Numerical results may vary slightly across environments. The summary workflow records the active R session with its outputs.

## Repository contents

| Path | Contents |
| --- | --- |
| `run.R` | Entry point for listing jobs, refitting models, and generating results. |
| `code/` | Installation, model fitting, summaries, and figure scripts. |
| `data/dengue.rds.part1`, `data/dengue.rds.part2` | Binary parts of the dengue input; combine them as described above. |
| `data/dengue.rds` | Recreated from the two parts. Contains archived and reconstructed covariates, sampled backgrounds, country exclusion flags, scaling information, and spatial model inputs. |
| `data/dengue_jobs.csv` | The 120 dengue fitting jobs, including sampling seeds and interventions. |
| `data/synthetic.rds` | Simulated landscape and generating/model settings. |
| `data/west_nile.rds` | County data, archived background covariates, population, calibrated rates, and model settings. |
| `data/west_nile_background.rds` | Archived occurrence records, the finite background frame, and simulated null cases. |
| `data/results/` | Saved posterior means, standard deviations, and marginal 95% credible intervals. |

The RDS files contain named R lists. Inspect their components with, for example:

```r
x <- readRDS("data/dengue.rds")
str(x, max.level = 1)
```

Reproduction starts from the supplied matrices, spatial meshes, and simulation settings. Raw raster acquisition and the GIS steps used to prepare these inputs are outside this repository.

## Dengue background comparisons

The four sampling arms use weights proportional to:

| Arm | Weight | Eligible cells |
| --- | --- | --- |
| A | Cell area × log(1 + population), with zero log weights replaced by 1.06840784199454 | Original sampling frame |
| B | Cell area × log(1 + population) | Inhabited cells |
| C | Cell area × population | Inhabited cells |
| D | Population | Inhabited cells |

The C minus B comparison isolates the log-to-linear weighting change while retaining the cell-area factor and inhabited sampling frame. D minus A compares the complete sampling recipes. B minus A examines the handling of uninhabited cells, and D minus C removes the cell-area factor.

The four main seeds are 980019, 980027, 980035, and 980043. Sixteen additional A/D pairs use seeds 980001–980016. Each global fit uses 11,313 occurrence records and 12,108 background locations. Country exclusions remove flagged records and backgrounds without drawing new locations.

Dengue covariates use the archived scaling parameters, stored in `model$scaling`; predictor names are in `model$vars`. Healthcare travel time is transformed as `log(1 + minutes)` before standardization. Coordinates are longitude and latitude, and `presence` distinguishes occurrence records (1) from background locations (0).

Paired changes compare posterior means. Across-draw ranges and Monte Carlo standard errors describe variation from background sampling, separately from each model's marginal credible interval. Coefficient units depend on the module: the West Nile model-comparison script converts internally standardized slopes back to input predictor units.

## Data sources and reconstruction

The empirical examples use archived inputs and model specifications from Gibb et al. (2026), [*The anthropogenic fingerprint on emerging infectious diseases*](https://doi.org/10.1038/s41586-026-11058-6), including the [data and model release](https://doi.org/10.5281/zenodo.20738794) and [source code, version 1.0.0](https://github.com/rorygibb/fingerprint_ms/tree/v1.0.0).

The dengue sensitivity analyses use reconstructed covariates alongside archived values for comparisons at the original locations. Five predictors remain unmatched to the original source layers: forest loss, precipitation change, protected areas, social vulnerability, and temperature change. The reconstruction uses Hansen 2021 v1.9 forest loss, ERA5-Land climate data, October 2026 WDPCA protected areas, and a GRDI social-vulnerability layer whose historical edition is unconfirmed. These differences should be considered when interpreting the empirical comparisons.

Country exclusion flags derive from GADM 4.1 boundaries. Boundary source files are not included. The West Nile background exercise resamples the supplied finite background frame; it does not reconstruct a new global background surface.
