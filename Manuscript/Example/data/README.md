# Data files

## MAP2B input

`input/abundance.tsv` and `input/coverage.tsv` are the species-level tables
produced by MAP2B v1.8. The original files used an `.xls` extension but contain
tab-separated text.

Each table contains:

- 1,132 species rows;
- seven taxonomy columns from kingdom to species;
- 45 pseudonymised sample columns;
- 15 normal cervical samples (`Nor`);
- 15 cervical intraepithelial neoplasia samples (`CIN`);
- 15 cervical cancer samples (`Ca`).

MD5 checksums:

```text
c966f93ac406949bfaf3ac1415dc7e4f  abundance.tsv
2fbd35a262e2ea69a47154eee1b642d1  coverage.tsv
```

These sample-level matrices are included locally for complete reproducibility
and excluded from Git by default. Confirm that the study consent, institutional
policy and repository terms permit public release before changing this rule.

## Derived data

The analysis writes the following files to `derived/`:

```text
sample_metadata.tsv
loocv_predictions.tsv
fold_parameters.tsv
fold_selected_features.tsv
feature_selection_frequency.tsv
performance_metrics.tsv
analysis_parameters.tsv
bootstrap_metrics.tsv
bootstrap_curve_points.tsv
bootstrap_parameters.tsv
```

`loocv_predictions.tsv` contains one prediction for every held-out sample in
each data type and pairwise comparison. Cervical cancer is the positive class.
The bootstrap files are derived from these fixed out-of-fold predictions.

The aggregate bootstrap metrics and curve points contain no sample identifiers
and are sufficient to reproduce the PDF. Derived files containing sample
identifiers are excluded from Git until their public release is approved.
