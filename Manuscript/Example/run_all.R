# Run the complete analysis from MAP2B profiles to the final PDF.

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_argument) != 1L) {
  stop("Run this file with Rscript")
}
project_dir <- dirname(normalizePath(sub("^--file=", "", script_argument)))
rscript <- file.path(R.home("bin"), "Rscript")

scripts <- c(
  "01_fit_pairwise_models.R",
  "02_bootstrap_performance.R",
  "03_plot_pairwise_performance.R"
)

for (script in scripts) {
  path <- file.path(project_dir, "R", script)
  cat("\nRunning", script, "\n")
  status <- system2(
    rscript,
    args = c(shQuote(path), shQuote(project_dir))
  )
  if (!identical(status, 0L)) {
    stop(script, " failed with exit status ", status)
  }
}

cat("\nAnalysis complete\n")
