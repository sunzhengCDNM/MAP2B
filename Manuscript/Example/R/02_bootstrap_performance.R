# Estimate uncertainty from class-stratified resamples of the fixed
# out-of-fold predictions. The models are not refitted in this step.

required_packages <- c("data.table", "pROC", "PRROC")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(data.table)
  library(pROC)
  library(PRROC)
})

args <- commandArgs(trailingOnly = TRUE)
project_dir <- normalizePath(if (length(args) == 0L) "." else args[[1L]])
data_dir <- file.path(project_dir, "data", "derived")
prediction_file <- file.path(data_dir, "loocv_predictions.tsv")

bootstrap_seed <- 20260722L
n_bootstrap <- 2000L
curve_grid <- seq(0, 1, length.out = 201L)

calculate_curves <- function(true_binary, prob_case) {
  true_binary <- factor(true_binary, levels = c("control", "case"))
  roc_fit <- pROC::roc(
    response = true_binary,
    predictor = prob_case,
    levels = c("control", "case"),
    direction = "<",
    quiet = TRUE
  )
  pr_fit <- PRROC::pr.curve(
    scores.class0 = prob_case[true_binary == "case"],
    scores.class1 = prob_case[true_binary == "control"],
    curve = TRUE
  )

  list(
    auroc = as.numeric(pROC::auc(roc_fit)),
    auprc = as.numeric(pr_fit$auc.integral),
    roc = data.frame(
      x = 1 - roc_fit$specificities,
      y = roc_fit$sensitivities
    ),
    pr = data.frame(
      x = pr_fit$curve[, 1L],
      y = pr_fit$curve[, 2L]
    )
  )
}

interpolate_roc <- function(curve) {
  collapsed <- aggregate(y ~ x, data = curve, FUN = max)
  collapsed <- collapsed[order(collapsed$x), , drop = FALSE]
  collapsed$y <- cummax(collapsed$y)
  approx(
    x = collapsed$x,
    y = collapsed$y,
    xout = curve_grid,
    method = "linear",
    rule = 2,
    ties = "ordered"
  )$y
}

interpolate_pr <- function(curve) {
  collapsed <- aggregate(y ~ x, data = curve, FUN = mean)
  collapsed <- collapsed[order(collapsed$x), , drop = FALSE]
  approx(
    x = collapsed$x,
    y = collapsed$y,
    xout = curve_grid,
    method = "linear",
    rule = 2,
    ties = "ordered"
  )$y
}

bootstrap_predictions <- function(group_data, seed) {
  observed <- calculate_curves(group_data$true_binary, group_data$prob_case)
  case_index <- which(group_data$true_binary == "case")
  control_index <- which(group_data$true_binary == "control")

  auroc <- numeric(n_bootstrap)
  auprc <- numeric(n_bootstrap)
  roc_values <- matrix(
    NA_real_,
    nrow = n_bootstrap,
    ncol = length(curve_grid)
  )
  pr_values <- matrix(
    NA_real_,
    nrow = n_bootstrap,
    ncol = length(curve_grid)
  )

  set.seed(seed)
  for (i in seq_len(n_bootstrap)) {
    index <- c(
      sample(case_index, length(case_index), replace = TRUE),
      sample(control_index, length(control_index), replace = TRUE)
    )
    result <- calculate_curves(
      group_data$true_binary[index],
      group_data$prob_case[index]
    )
    auroc[[i]] <- result$auroc
    auprc[[i]] <- result$auprc
    roc_values[i, ] <- interpolate_roc(result$roc)
    pr_values[i, ] <- interpolate_pr(result$pr)
  }

  summarise_curve <- function(values, curve_type) {
    data.table(
      curve_type = curve_type,
      x = curve_grid,
      mean = colMeans(values),
      lower = apply(values, 2L, quantile, probs = 0.025, names = FALSE),
      upper = apply(values, 2L, quantile, probs = 0.975, names = FALSE)
    )
  }

  list(
    metrics = data.table(
      auroc = observed$auroc,
      auroc_lower = unname(quantile(auroc, 0.025)),
      auroc_upper = unname(quantile(auroc, 0.975)),
      auprc = observed$auprc,
      auprc_lower = unname(quantile(auprc, 0.025)),
      auprc_upper = unname(quantile(auprc, 0.975))
    ),
    curves = rbind(
      summarise_curve(roc_values, "roc"),
      summarise_curve(pr_values, "pr")
    )
  )
}

predictions <- fread(prediction_file)
required_columns <- c(
  "data_type", "comparison", "model", "sample_id", "true_binary", "prob_case"
)
if (!all(required_columns %in% names(predictions))) {
  stop("The out-of-fold prediction table is incomplete")
}

group_counts <- predictions[
  ,
  .(
    n_samples = uniqueN(sample_id),
    n_case = sum(true_binary == "case"),
    n_control = sum(true_binary == "control")
  ),
  by = .(data_type, comparison, model)
]
if (
  nrow(group_counts) != 4L ||
    any(group_counts$n_samples != 30L) ||
    any(group_counts$n_case != 15L) ||
    any(group_counts$n_control != 15L)
) {
  stop("Each analysis must contain 15 cases and 15 controls")
}

# Offsets preserve the random streams used for the reported figure.
analysis_groups <- data.table(
  data_type = c("Abundance", "Coverage", "Abundance", "Coverage"),
  comparison = c("CIN_vs_Ca", "CIN_vs_Ca", "Nor_vs_Ca", "Nor_vs_Ca"),
  model = "Boruta_RF",
  seed_offset = c(1000L, 2000L, 9000L, 10000L)
)

metric_records <- vector("list", nrow(analysis_groups))
curve_records <- vector("list", nrow(analysis_groups))

for (i in seq_len(nrow(analysis_groups))) {
  group <- analysis_groups[i]
  group_data <- predictions[
    data_type == group$data_type &
      comparison == group$comparison &
      model == group$model
  ]
  result <- bootstrap_predictions(
    group_data = group_data,
    seed = bootstrap_seed + group$seed_offset
  )

  identifiers <- group[, .(data_type, comparison, model)]
  metric_records[[i]] <- cbind(identifiers, result$metrics)
  curve_records[[i]] <- cbind(identifiers, result$curves)

  cat(group$comparison, group$data_type, "complete\n")
}

bootstrap_metrics <- rbindlist(metric_records)
bootstrap_curves <- rbindlist(curve_records)
bootstrap_parameters <- data.table(
  parameter = c(
    "bootstrap_seed", "n_bootstrap", "sampling",
    "curve_grid_points", "curve_interval"
  ),
  value = c(
    bootstrap_seed,
    n_bootstrap,
    "class-stratified resampling of fixed out-of-fold predictions",
    length(curve_grid),
    "pointwise 95% percentile interval"
  )
)

fwrite(
  bootstrap_metrics,
  file.path(data_dir, "bootstrap_metrics.tsv"),
  sep = "\t"
)
fwrite(
  bootstrap_curves,
  file.path(data_dir, "bootstrap_curve_points.tsv"),
  sep = "\t"
)
fwrite(
  bootstrap_parameters,
  file.path(data_dir, "bootstrap_parameters.tsv"),
  sep = "\t"
)

cat("Bootstrap analysis complete\n")
print(bootstrap_metrics)
