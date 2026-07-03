# =============================================================================
# run_all.R  --  end-to-end synthetic M2FPCA mHealth analysis
#
# Reproduces (on fully synthetic data) the workflow of Dey, Ghosal, Merikangas &
# Zipunnikov (2026), "Multivariate Functional Principal Component Analysis for
# Mixed-Type mHealth Data: An Application to Mood Disorders" (arXiv:2603.11385):
# estimate the joint covariance of mixed-type
# EMA + actigraphy functional data, extract shared principal components,
# compute subject-level scores as digital biomarkers, and use them to predict
# diagnostic group.
#
# Usage (from this directory):
#     Rscript run_all.R
#
# Outputs land in ./output (RDS artifacts) and ./figures (PNGs).
# NO REAL PARTICIPANT DATA IS USED — see 01_generate_synthetic_data.R.
# =============================================================================

here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) getwd())
Sys.setenv(SDA_OUT = file.path(here, "output"),
           SDA_FIG = file.path(here, "figures"))

steps <- c("01_generate_synthetic_data.R",
           "02_run_m2fpca.R",
           "03_scores_and_prediction.R",
           "04_figures.R")

for (s in steps) {
  cat("\n############################################################\n")
  cat("##", s, "\n")
  cat("############################################################\n")
  source(file.path(here, s), echo = FALSE)
}

cat("\nDONE. See ./output and ./figures.\n")
