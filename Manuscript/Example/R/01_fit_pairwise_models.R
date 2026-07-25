# Fit Boruta-selected random forests for the Nor-versus-Ca and
# CIN-versus-Ca comparisons. All data-dependent steps are estimated within
# each outer leave-one-out training fold.

required_packages <- c(
  "data.table", "randomForest", "Boruta", "matrixTests", "pROC", "PRROC"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(data.table)
  library(randomForest)
  library(Boruta)
  library(matrixTests)
  library(pROC)
  library(PRROC)
})

args <- commandArgs(trailingOnly = TRUE)
project_dir <- normalizePath(if (length(args) == 0L) "." else args[[1L]])
input_dir <- file.path(project_dir, "data", "input")
output_dir <- file.path(project_dir, "data", "derived")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

analysis_seed <- 20260722L
prevalence_cutoff <- 0.10
boruta_max_runs <- 200L
tuning_trees <- 100L
final_trees <- 500L
terminal_node_size <- 5L

read_profile_table <- function(path, data_type) {
  taxonomy_columns <- c(
    "#Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"
  )
  profile <- fread(path, check.names = FALSE)
  if (!identical(names(profile)[seq_along(taxonomy_columns)], taxonomy_columns)) {
    stop(data_type, " table does not contain the expected taxonomy columns")
  }

  species <- as.character(profile$Species)
  sample_ids <- names(profile)[-(seq_along(taxonomy_columns))]
  status <- sub("-[0-9]+$", "", sample_ids)

  if (anyNA(species) || any(species == "") || anyDuplicated(species)) {
    stop(data_type, " table contains invalid or duplicated species names")
  }
  if (anyDuplicated(sample_ids) || !all(status %in% c("Nor", "CIN", "Ca"))) {
    stop(data_type, " table contains invalid sample identifiers")
  }

  feature_by_sample <- as.matrix(profile[, ..sample_ids])
  storage.mode(feature_by_sample) <- "numeric"
  if (any(!is.finite(feature_by_sample)) || any(feature_by_sample < 0)) {
    stop(data_type, " table contains invalid feature values")
  }

  feature_id <- make.names(species, unique = TRUE)
  x <- t(feature_by_sample)
  rownames(x) <- sample_ids
  colnames(x) <- feature_id

  list(
    x = x,
    metadata = data.table(
      sample_id = sample_ids,
      status = factor(status, levels = c("Nor", "CIN", "Ca"))
    ),
    taxonomy = profile[, ..taxonomy_columns],
    feature_map = data.table(feature_id = feature_id, species = species),
    path = normalizePath(path),
    data_type = data_type
  )
}

check_profile_alignment <- function(abundance, coverage) {
  if (!identical(abundance$metadata$sample_id, coverage$metadata$sample_id)) {
    stop("Abundance and coverage sample columns are not aligned")
  }
  if (!identical(abundance$feature_map$species, coverage$feature_map$species)) {
    stop("Abundance and coverage species rows are not aligned")
  }
  if (!identical(abundance$taxonomy, coverage$taxonomy)) {
    stop("Abundance and coverage taxonomy columns are not aligned")
  }

  counts <- table(abundance$metadata$status)
  if (!identical(as.integer(counts[c("Nor", "CIN", "Ca")]), c(15L, 15L, 15L))) {
    stop("Expected 15 Nor, 15 CIN and 15 Ca samples")
  }
}

preprocess_fold <- function(x_train_raw, x_test_raw) {
  prevalence <- colMeans(x_train_raw > 0)
  retained <- names(prevalence)[prevalence > prevalence_cutoff]
  if (length(retained) == 0L) {
    stop("No features passed the training-fold prevalence filter")
  }

  x_train <- x_train_raw[, retained, drop = FALSE]
  x_test <- x_test_raw[, retained, drop = FALSE]
  minimum_nonzero <- min(x_train[x_train > 0])
  pseudocount <- minimum_nonzero / 10

  x_train <- log10(x_train + pseudocount)
  x_test <- log10(x_test + pseudocount)

  feature_variance <- apply(x_train, 2L, var)
  variable <- names(feature_variance)[
    is.finite(feature_variance) & feature_variance > 0
  ]
  if (length(variable) == 0L) {
    stop("No variable features remained after transformation")
  }

  list(
    x_train = x_train[, variable, drop = FALSE],
    x_test = x_test[, variable, drop = FALSE],
    pseudocount = pseudocount,
    n_prevalence = length(retained),
    n_variable = length(variable)
  )
}

rank_features_by_wilcoxon <- function(x_train, y_train) {
  result <- suppressWarnings(
    matrixTests::col_wilcoxon_twosample(
      x_train[y_train == "control", , drop = FALSE],
      x_train[y_train == "case", , drop = FALSE]
    )
  )
  p_value <- result$pvalue
  p_value[!is.finite(p_value)] <- 1
  names(p_value) <- colnames(x_train)
  names(sort(p_value))
}

select_boruta_features <- function(x_train, y_train, seed) {
  set.seed(seed)
  boruta_error <- NA_character_
  fit <- tryCatch(
    Boruta::Boruta(
      x = as.data.frame(x_train, check.names = FALSE),
      y = y_train,
      pValue = 0.01,
      mcAdj = TRUE,
      maxRuns = boruta_max_runs,
      doTrace = 0
    ),
    error = function(e) {
      boruta_error <<- conditionMessage(e)
      NULL
    }
  )

  selected <- if (is.null(fit)) {
    character()
  } else {
    Boruta::getSelectedAttributes(fit, withTentative = FALSE)
  }

  if (length(selected) < 5L) {
    ranked <- rank_features_by_wilcoxon(x_train, y_train)
    selected <- head(ranked, min(10L, length(ranked)))
    method <- "wilcoxon_fallback"
  } else {
    method <- "boruta"
  }

  list(features = selected, method = method, error = boruta_error)
}

make_mtry_grid <- function(n_features) {
  candidates <- floor(c(
    log2(n_features),
    sqrt(n_features) / 2,
    sqrt(n_features),
    sqrt(n_features) * 2,
    n_features * 0.10,
    n_features * 0.20,
    n_features / 3
  ))
  as.integer(sort(unique(candidates[candidates >= 1L & candidates <= n_features])))
}

fit_random_forest <- function(x_train, y_train, x_test, seed) {
  mtry_grid <- make_mtry_grid(ncol(x_train))
  oob_error <- rep(Inf, length(mtry_grid))

  for (i in seq_along(mtry_grid)) {
    set.seed(seed + i)
    tuning_fit <- randomForest(
      x = x_train,
      y = y_train,
      ntree = tuning_trees,
      mtry = mtry_grid[[i]],
      nodesize = terminal_node_size,
      importance = FALSE
    )
    candidate_error <- tuning_fit$err.rate[tuning_trees, "OOB"]
    if (is.finite(candidate_error)) {
      oob_error[[i]] <- candidate_error
    }
  }

  best_index <- which.min(oob_error)
  best_mtry <- mtry_grid[[best_index]]
  set.seed(seed + 1000L)
  final_fit <- randomForest(
    x = x_train,
    y = y_train,
    ntree = final_trees,
    mtry = best_mtry,
    nodesize = terminal_node_size,
    importance = FALSE
  )
  probability <- predict(final_fit, newdata = x_test, type = "prob")

  list(
    prob_case = as.numeric(probability[, "case"]),
    best_mtry = best_mtry,
    best_oob_error = oob_error[[best_index]],
    mtry_grid = paste(mtry_grid, collapse = ",")
  )
}

calculate_performance <- function(true_binary, prob_case) {
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
    curve = FALSE
  )
  list(
    auroc = as.numeric(pROC::auc(roc_fit)),
    auprc = as.numeric(pr_fit$auc.integral)
  )
}

# Read and align the two MAP2B feature matrices.
abundance <- read_profile_table(
  file.path(input_dir, "abundance.tsv"), "Abundance"
)
coverage <- read_profile_table(
  file.path(input_dir, "coverage.tsv"), "Coverage"
)
check_profile_alignment(abundance, coverage)

profiles <- list(Abundance = abundance, Coverage = coverage)
data_seed_offset <- c(Abundance = 100000L, Coverage = 200000L)
comparisons <- list(
  Nor_vs_Ca = list(
    groups = c("Nor", "Ca"),
    positive = "Ca",
    seed_offset = 20000L
  ),
  CIN_vs_Ca = list(
    groups = c("CIN", "Ca"),
    positive = "Ca",
    seed_offset = 30000L
  )
)

prediction_records <- list()
parameter_records <- list()
feature_records <- list()
record_index <- 0L

# Outer leave-one-out cross-validation.
for (data_type in names(profiles)) {
  profile <- profiles[[data_type]]

  for (comparison_name in names(comparisons)) {
    comparison <- comparisons[[comparison_name]]
    pair_index <- which(profile$metadata$status %in% comparison$groups)
    pair_metadata <- profile$metadata[pair_index]
    x_pair <- profile$x[pair_index, , drop = FALSE]
    y_pair <- factor(
      ifelse(pair_metadata$status == comparison$positive, "case", "control"),
      levels = c("control", "case")
    )

    if (!identical(as.integer(table(y_pair)), c(15L, 15L))) {
      stop(comparison_name, " must contain 15 controls and 15 cases")
    }

    for (test_index in seq_len(nrow(x_pair))) {
      train_index <- setdiff(seq_len(nrow(x_pair)), test_index)
      y_train <- factor(y_pair[train_index], levels = c("control", "case"))
      fold <- preprocess_fold(
        x_train_raw = x_pair[train_index, , drop = FALSE],
        x_test_raw = x_pair[test_index, , drop = FALSE]
      )

      fold_seed <- analysis_seed +
        data_seed_offset[[data_type]] +
        comparison$seed_offset +
        test_index * 100L +
        2L

      selection <- select_boruta_features(
        x_train = fold$x_train,
        y_train = y_train,
        seed = fold_seed
      )
      selected <- selection$features
      rf <- fit_random_forest(
        x_train = fold$x_train[, selected, drop = FALSE],
        y_train = y_train,
        x_test = fold$x_test[, selected, drop = FALSE],
        seed = fold_seed
      )

      record_index <- record_index + 1L
      prediction_records[[record_index]] <- data.table(
        data_type = data_type,
        comparison = comparison_name,
        model = "Boruta_RF",
        fold = test_index,
        sample_id = pair_metadata$sample_id[[test_index]],
        true_status = as.character(pair_metadata$status[[test_index]]),
        positive_status = comparison$positive,
        true_binary = as.character(y_pair[[test_index]]),
        prob_case = rf$prob_case
      )

      parameter_records[[record_index]] <- data.table(
        data_type = data_type,
        comparison = comparison_name,
        model = "Boruta_RF",
        fold = test_index,
        held_out_sample = pair_metadata$sample_id[[test_index]],
        seed = fold_seed,
        prevalence_cutoff = prevalence_cutoff,
        pseudocount = fold$pseudocount,
        n_prevalence_features = fold$n_prevalence,
        n_variable_features = fold$n_variable,
        selection_method = selection$method,
        n_selected_features = length(selected),
        boruta_error = selection$error,
        best_mtry = rf$best_mtry,
        best_oob_error = rf$best_oob_error,
        mtry_grid = rf$mtry_grid
      )

      feature_map <- profile$feature_map[
        match(selected, profile$feature_map$feature_id)
      ]
      feature_records[[record_index]] <- data.table(
        data_type = data_type,
        comparison = comparison_name,
        model = "Boruta_RF",
        fold = test_index,
        held_out_sample = pair_metadata$sample_id[[test_index]],
        selection_method = selection$method,
        feature_id = feature_map$feature_id,
        species = feature_map$species
      )
    }
  }
}

predictions <- rbindlist(prediction_records)
fold_parameters <- rbindlist(parameter_records)
fold_selected_features <- rbindlist(feature_records)

performance_metrics <- predictions[
  ,
  {
    performance <- calculate_performance(true_binary, prob_case)
    .(
      positive_status = unique(positive_status),
      n_samples = .N,
      n_positive = sum(true_binary == "case"),
      n_negative = sum(true_binary == "control"),
      auroc = performance$auroc,
      auprc = performance$auprc
    )
  },
  by = .(data_type, comparison, model)
]

feature_selection_frequency <- fold_selected_features[
  ,
  .(n_selected = .N),
  by = .(data_type, comparison, model, feature_id, species)
][
  ,
  selection_frequency := n_selected / 30
][
  order(data_type, comparison, -selection_frequency, species)
]

sample_metadata <- copy(abundance$metadata)
sample_metadata[, status := as.character(status)]

analysis_parameters <- data.table(
  parameter = c(
    "analysis_seed", "prevalence_cutoff", "boruta_p_value",
    "boruta_multiple_testing_adjustment", "boruta_max_runs",
    "boruta_tentative_features", "wilcoxon_fallback_threshold",
    "wilcoxon_fallback_features", "tuning_trees", "final_trees",
    "terminal_node_size", "abundance_md5", "coverage_md5"
  ),
  value = c(
    analysis_seed, prevalence_cutoff, 0.01, TRUE, boruta_max_runs,
    "excluded", 5, 10, tuning_trees, final_trees, terminal_node_size,
    unname(tools::md5sum(abundance$path)),
    unname(tools::md5sum(coverage$path))
  )
)

fwrite(sample_metadata, file.path(output_dir, "sample_metadata.tsv"), sep = "\t")
fwrite(predictions, file.path(output_dir, "loocv_predictions.tsv"), sep = "\t")
fwrite(fold_parameters, file.path(output_dir, "fold_parameters.tsv"), sep = "\t")
fwrite(
  fold_selected_features,
  file.path(output_dir, "fold_selected_features.tsv"),
  sep = "\t"
)
fwrite(
  feature_selection_frequency,
  file.path(output_dir, "feature_selection_frequency.tsv"),
  sep = "\t"
)
fwrite(
  performance_metrics,
  file.path(output_dir, "performance_metrics.tsv"),
  sep = "\t"
)
fwrite(
  analysis_parameters,
  file.path(output_dir, "analysis_parameters.tsv"),
  sep = "\t"
)

cat("Pairwise LOOCV complete\n")
print(performance_metrics)
