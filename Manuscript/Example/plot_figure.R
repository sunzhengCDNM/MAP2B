#!/usr/bin/env Rscript

# Recreate the PDF from the committed bootstrap summaries.

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_argument) != 1L) {
  stop("Run this file with Rscript")
}
project_dir <- dirname(normalizePath(sub("^--file=", "", script_argument)))
rscript <- file.path(R.home("bin"), "Rscript")
plot_script <- file.path(project_dir, "R", "03_plot_pairwise_performance.R")

status <- system2(
  rscript,
  args = c(shQuote(plot_script), shQuote(project_dir))
)
if (!identical(status, 0L)) {
  stop("Figure generation failed with exit status ", status)
}
