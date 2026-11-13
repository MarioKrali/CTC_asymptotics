# Common setup, sourced at the top of every script in scripts/.
#
# Environment variables:
#   CTC_QUICK=1   run strongly reduced settings (smoke test; NOT the paper results).
#                 Outputs then go to results/quick/ and figures/quick/.
#   CTC_CORES=n   number of cores for the simulation scripts (default 1).
#                 Each configuration re-seeds itself, so results are identical
#                 for any number of cores. Forking is not available on Windows,
#                 where it always uses one core.

if (!file.exists(file.path("R", "setup.R"))) {
  stop("Run the scripts from the repository root (the folder containing R/ and scripts/).",
       call. = FALSE)
}

QUICK   <- identical(Sys.getenv("CTC_QUICK"), "1")
N_CORES <- max(1L, as.integer(Sys.getenv("CTC_CORES", "1")))

FIG_DIR  <- if (QUICK) file.path("figures", "quick") else "figures"
RES_DIR  <- if (QUICK) file.path("results", "quick") else "results"
DATA_DIR <- "data"
for (dir in c(FIG_DIR, RES_DIR, file.path(RES_DIR, "tables"))) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
}

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) {
  if (basename(f) != "setup.R") source(f)
}

if (QUICK) message("CTC_QUICK=1: reduced settings for testing only (not the paper results).")
