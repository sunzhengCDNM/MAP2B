# Pairwise classification of cervical cancer from MAP2B profiles

This repository contains the R code and processed MAP2B profiles used to
evaluate species-level abundance and coverage classifiers for cervical cancer.
The reported figure compares normal cervical samples with cervical cancer and
cervical intraepithelial neoplasia with cervical cancer.

## Analysis

Species abundance and coverage are analysed separately. Each comparison
contains 15 samples per class, with cervical cancer defined as the positive
class. Model performance is evaluated by outer leave-one-out cross-validation
(LOOCV). Preprocessing, Boruta feature selection and random-forest tuning are
performed within each training fold.

Uncertainty is estimated from 2,000 class-stratified bootstrap resamples of the
fixed out-of-fold predictions. This bootstrap quantifies uncertainty conditional
on the observed out-of-fold predictions; it does not refit the models.

## Files

```text
R/01_fit_pairwise_models.R       LOOCV, Boruta selection and random forests
R/02_bootstrap_performance.R     Bootstrap confidence intervals and mean curves
R/03_plot_pairwise_performance.R Final AUROC and AUPRC figure
run_all.R                        Complete analysis
plot_figure.R                    Figure-only reproduction
data/input/                      MAP2B abundance and coverage tables
data/derived/                    Predictions, metrics and curve summaries
results/                         Final PDF
```

## Requirements

The analysis uses R and the following packages:

```r
install.packages(c(
  "data.table",
  "randomForest",
  "Boruta",
  "matrixTests",
  "pROC",
  "PRROC",
  "ggplot2",
  "patchwork"
))
```

Tested package versions are recorded in
`environment/package_versions.tsv` and `environment/sessionInfo.txt`.
The complete dependency graph is recorded in `renv.lock`. Restore it before
running the analysis in a new R installation:

```bash
Rscript -e 'install.packages("renv")'
Rscript -e 'renv::restore()'
```

## Reproduce the figure

The public-safe workflow uses aggregate bootstrap summaries and does not
require sample-level data:

```bash
Rscript plot_figure.R
```


## Reproduce the complete analysis

After obtaining the required data-sharing approval, place `abundance.tsv` and
`coverage.tsv` in `data/input/` and run:

```bash
Rscript run_all.R
```

## Data

The input files are tab-separated MAP2B output tables despite the `.xls`
extension used by the original pipeline. They are stored here with `.tsv`
extensions. See `data/README.md` for file structure, checksums and the
data-sharing notice.

The sample-level matrices and derived tables containing sample identifiers are
retained locally but excluded from Git by default. Remove the corresponding
entries from `.gitignore` only after institutional data-sharing approval has
been confirmed.

## Computational environment

The reference analysis was run on a MacBook Pro with an Apple M4 Pro chip
(14-core CPU), 48 GB unified memory and macOS 15.6. MAP2B processing of the 45
samples required approximately 2 h 25 min. The R analysis starts from the
MAP2B abundance and coverage tables.

