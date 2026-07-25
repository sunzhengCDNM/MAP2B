# Draw the publication figure from the saved bootstrap summaries.

required_packages <- c("data.table", "ggplot2", "patchwork")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
project_dir <- normalizePath(if (length(args) == 0L) "." else args[[1L]])
data_dir <- file.path(project_dir, "data", "derived")
result_dir <- file.path(project_dir, "results")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

metrics <- fread(file.path(data_dir, "bootstrap_metrics.tsv"))
curves <- fread(file.path(data_dir, "bootstrap_curve_points.tsv"))

series_levels <- c(
  "Nor_vs_Ca__Abundance",
  "Nor_vs_Ca__Coverage",
  "CIN_vs_Ca__Abundance",
  "CIN_vs_Ca__Coverage"
)
series_names <- c(
  Nor_vs_Ca__Abundance = "Nor vs Ca | Abundance",
  Nor_vs_Ca__Coverage = "Nor vs Ca | Coverage",
  CIN_vs_Ca__Abundance = "CIN vs Ca | Abundance",
  CIN_vs_Ca__Coverage = "CIN vs Ca | Coverage"
)
series_colours <- c(
  Nor_vs_Ca__Abundance = "#08519C",
  Nor_vs_Ca__Coverage = "#6BAED6",
  CIN_vs_Ca__Abundance = "#D94801",
  CIN_vs_Ca__Coverage = "#F6A21A"
)

add_series <- function(x) {
  x[, series := factor(
    paste(comparison, data_type, sep = "__"),
    levels = series_levels
  )]
  x
}

metrics <- add_series(metrics)
curves <- add_series(curves)

if (
  nrow(metrics) != 4L ||
    anyNA(metrics$series) ||
    anyNA(curves$series)
) {
  stop("Expected four Boruta-RF analysis series")
}

make_panel <- function(curve_type) {
  requested_curve_type <- curve_type
  metric_name <- if (curve_type == "roc") "AUROC" else "AUPRC"
  estimate_column <- tolower(metric_name)
  lower_column <- paste0(estimate_column, "_lower")
  upper_column <- paste0(estimate_column, "_upper")

  legend_data <- copy(metrics)
  legend_data[, label := sprintf(
    "%s (%s %.2f, 95%% CI %.2f-%.2f)",
    series_names[as.character(series)],
    metric_name,
    get(estimate_column),
    get(lower_column),
    get(upper_column)
  )]
  legend_labels <- setNames(
    legend_data$label,
    as.character(legend_data$series)
  )
  panel_data <- curves[curve_type == requested_curve_type]

  panel <- ggplot(
    panel_data,
    aes(x = x, colour = series, fill = series, group = series)
  ) +
    geom_ribbon(
      aes(ymin = lower, ymax = upper),
      alpha = 0.11,
      colour = NA
    ) +
    geom_line(aes(y = mean), linewidth = 1) +
    scale_colour_manual(
      values = series_colours,
      breaks = series_levels,
      labels = legend_labels[series_levels]
    ) +
    scale_fill_manual(values = series_colours) +
    guides(
      fill = "none",
      colour = guide_legend(nrow = 4L, byrow = TRUE)
    ) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    labs(title = metric_name, colour = NULL) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom",
      legend.text = element_text(size = 8.2),
      legend.key.width = grid::unit(10, "mm"),
      plot.margin = margin(6, 8, 6, 8)
    )

  if (curve_type == "roc") {
    panel +
      geom_abline(
        intercept = 0,
        slope = 1,
        linetype = "dashed",
        colour = "grey60"
      ) +
      labs(x = "1 - Specificity", y = "Sensitivity")
  } else {
    panel +
      geom_hline(
        yintercept = 0.5,
        linetype = "dashed",
        colour = "grey60"
      ) +
      labs(x = "Recall", y = "Precision")
  }
}

figure <- (make_panel("roc") | make_panel("pr")) +
  plot_annotation(
    title = "Pairwise Classification of Cervical Cancer",
    subtitle = paste(
      "LOOCV performance of Boruta-selected random forests using",
      "species abundance and coverage"
    )
  )

output_file <- file.path(result_dir, "Figure_pairwise_classification.pdf")
ggsave(
  output_file,
  figure,
  width = 14,
  height = 7.2,
  units = "in",
  device = grDevices::pdf
)

cat("Figure written to", output_file, "\n")
