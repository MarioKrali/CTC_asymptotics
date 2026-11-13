# Run every analysis script in order (from the repository root):
#   Rscript run_all.R
# Set CTC_CORES to parallelise and CTC_QUICK=1 for a fast smoke test (see README).

scripts <- sort(list.files("scripts", pattern = "^[0-9]+_.*\\.R$", full.names = TRUE))
for (s in scripts) {
  message("\n==== ", s, " ====")
  t0 <- Sys.time()
  # each script runs in a fresh R process so no objects leak between scripts
  status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(s))
  if (status != 0) stop(s, " failed (exit status ", status, ")")
  message("Finished in ", format(round(Sys.time() - t0, 1)))
}
