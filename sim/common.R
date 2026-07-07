# common.R -- shared prelude for the Phase-4.B simulation drivers
#
# Run the drivers from the repo root:
#   RENV_CONFIG_AUTOLOADER_ENABLED=false Rscript sim/sim_multivariate.R
#   RENV_CONFIG_AUTOLOADER_ENABLED=false Rscript sim/sim_contamination.R
# They need only the package functions (sourced below, so no installed
# gtmeta is required) plus metafor and mvmeta on the library path.
# Replication counts are overridable via GTMETA_NREP* env vars; every
# cell sets its own seed, so results are reproducible per cell.

for (f in list.files("R", full.names = TRUE)) source(f)
dir.create(file.path("sim", "results"), showWarnings = FALSE,
           recursive = TRUE)

# k x m estimates + list of covariances -> list of igmi_gaussian objects
studies_from_sim <- function(y, Sigma) {
  lapply(seq_len(nrow(y)), function(i) {
    igmi_gaussian(y[i, ], Sigma = Sigma[[i]])
  })
}

print_section <- function(title) {
  cat("\n==", title, "==\n")
}
