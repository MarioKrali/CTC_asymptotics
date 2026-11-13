# Simulation study: recovery of the ancestral graph of heavy-tailed LSEMs with
# the bootstrap lower/upper CTC graphs, for d in {5, 10, 20}, p in {0.1, 0.2},
# n in {500, 5000, 50000}. k = floor(0.5 * n^0.4).
#
# Output: results/lsem_graph_coverage.rds, an array
#   [replicate, metric, d, p, n]   (metric names: LSEM_METRICS in R/lsem.R)
# Figures and tables: scripts/02_lsem_figures_tables.R.

source("R/setup.R")

d_vec <- c(5, 10, 20)
p_vec <- c(0.1, 0.2)
n_vec <- c(500, 5000, 50000)
reps  <- if (QUICK) 5 else 1000   # Monte Carlo replications
B     <- if (QUICK) 50 else 1000  # bootstrap replicates
gamma_threshold <- 0.95
c_var <- 0.1
alpha <- 0.1
k_of_n <- function(n) floor(0.5 * n^0.4)

run_config <- function(d, p, n) {
  set.seed(112233)  # same seed for every configuration (common random numbers)
  counts <- matrix(NA_real_, reps, 10)
  nshd   <- numeric(reps)
  for (g in seq_len(reps)) {
    A      <- simulate_lsem_A(d, p)
    true_A <- 1 * (ctc_population(A) >= gamma_threshold)
    nshd[g] <- nshd_ctc_vs_ancestral(true_A, A)

    X   <- simulate_lsem_X(A, n)
    fit <- tail_moments_dboot(reverse_ranks(X), k_of_n(n), B, gamma_threshold, c_var, alpha)
    cnt <- graph_recovery_counts(true_A,
                                 low = t(fit$ancestral_graph_lower),
                                 up  = t(fit$ancestral_graph_upper))
    counts[g, ] <- cnt
  }
  colnames(counts) <- names(cnt)
  recovery_metrics(counts, nshd)
}

configs <- expand.grid(n = n_vec, p = p_vec, d = d_vec)
t0 <- Sys.time()
out <- par_lapply(seq_len(nrow(configs)), function(i) {
  with(configs[i, ], run_config(d, p, n))
})
message("Simulation time: ", format(round(Sys.time() - t0, 1)))

agg_res <- array(NA_real_, dim = c(reps, length(LSEM_METRICS), length(d_vec), length(p_vec), length(n_vec)),
                 dimnames = list(NULL, LSEM_METRICS, as.character(d_vec), p_label(p_vec), as.character(n_vec)))
for (i in seq_len(nrow(configs))) {
  cf <- configs[i, ]
  agg_res[, , as.character(cf$d), p_label(cf$p), as.character(cf$n)] <- out[[i]]
}

save_result(agg_res, "lsem_graph_coverage.rds")
