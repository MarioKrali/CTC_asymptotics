# Small helpers shared by the scripts.

#' lapply() over independent simulation configurations, forked over N_CORES
#' cores on Unix. Every configuration calls set.seed() itself, so results do
#' not depend on the number of cores.
par_lapply <- function(X, FUN, cores = N_CORES) {
  if (cores > 1L && .Platform$OS.type == "unix") {
    out <- parallel::mclapply(X, FUN, mc.cores = cores, mc.preschedule = FALSE)
    failed <- vapply(out, inherits, logical(1), what = "try-error")
    if (any(failed)) stop("A parallel job failed:\n", out[[which(failed)[1]]])
    out
  } else {
    lapply(X, FUN)
  }
}

#' Save a ggplot to FIG_DIR (width/height in inches).
save_fig <- function(plot, file, width, height) {
  path <- file.path(FIG_DIR, file)
  ggplot2::ggsave(path, plot, width = width, height = height, units = "in")
  message("Saved ", path)
  invisible(path)
}

#' Save an R object to RES_DIR as .rds.
save_result <- function(object, file) {
  path <- file.path(RES_DIR, file)
  saveRDS(object, path)
  message("Saved ", path)
  invisible(path)
}

#' Read an .rds result written by an earlier script, with a helpful error.
load_result <- function(file, produced_by) {
  path <- file.path(RES_DIR, file)
  if (!file.exists(path)) {
    stop(path, " not found. Run ", produced_by, " first.", call. = FALSE)
  }
  readRDS(path)
}

#' Format numbers in (0, 1] as the array labels used throughout (0.1 -> ".1", 1 -> "1").
p_label <- function(p) sub("^0", "", as.character(p))
