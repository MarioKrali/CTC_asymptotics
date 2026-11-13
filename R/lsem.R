# Heavy-tailed linear structural equation models (LSEMs) and graph-recovery metrics.

#' Random LSEM coefficient matrix A = (I - C)^{-1} diag(scale).
#' C is strictly upper triangular (a DAG): each edge is present with
#' probability p and has weight Unif(0.2, 1.2). Innovation scales are Unif(0.5, 1.5).
simulate_lsem_A <- function(d, p) {
  C <- matrix(sample(c(0, 1), d^2, replace = TRUE, prob = c(1 - p, p)) *
                runif(d^2, 0.2, 1.2), nrow = d)
  C[lower.tri(C, diag = TRUE)] <- 0
  marginal_scale <- runif(d, 0.5, 1.5)
  solve(diag(1, d) - C) %*% diag(marginal_scale)
}

#' n observations of X = A Z with i.i.d. half-Cauchy innovations Z (n x d matrix).
simulate_lsem_X <- function(A, n) {
  d <- nrow(A)
  t(A %*% matrix(abs(rt(n * d, 1)), nrow = d))
}

#' Population CTC matrix of the LSEM:
#' Gamma[i, j] = 1/2 + 1/2 * sum_{s: A[i, s] > 0} A[j, s] / sum_s A[j, s],
#' i.e. conditioning on node j, outcome node i (the transpose of hat_Hn).
#' The diagonal is 0.
ctc_population <- function(A) {
  d <- nrow(A)
  Gamma <- matrix(0, d, d)
  for (i in seq_len(d)) {
    for (j in seq_len(d)[-i]) {
      Gamma[i, j] <- 1/2 + 1/2 * sum(A[j, A[i, ] > 0]) / sum(A[j, ])
    }
  }
  Gamma
}

#' Normalised SHD between the CTC-thresholded graph and the true ancestral graph of A.
nshd_ctc_vs_ancestral <- function(true_A, A) {
  anc <- (A - diag(diag(A))) > 0
  sum(abs(true_A - anc)) / max(sum(true_A) + sum(anc), 1)
}

#' Recovery statistics of a lower/upper graph pair against the true graph.
#' All inputs are d x d 0/1 (or logical) matrices with zero diagonal.
graph_recovery_counts <- function(true_A, low, up) {
  d <- nrow(true_A)
  n_pairs <- d * (d - 1)
  TP_L <- sum(low * true_A);          TP_U <- sum(up * true_A)
  FP_L <- sum(pmax(low - true_A, 0)); FP_U <- sum(pmax(up - true_A, 0))
  FN_L <- sum(pmax(true_A - low, 0)); FN_U <- sum(pmax(true_A - up, 0))
  c(
    covered  = (max(low - true_A) <= 0) && (min(up - true_A) >= 0),  # low ⊆ true ⊆ up
    pair_cov = (sum((true_A >= low) * (true_A <= up)) - d) / n_pairs,
    TP_L = TP_L, TP_U = TP_U, FP_L = FP_L, FP_U = FP_U, FN_L = FN_L, FN_U = FN_U,
    TN_L = n_pairs - TP_L - FP_L - FN_L,
    TN_U = n_pairs - TP_U - FP_U - FN_U
  )
}

#' Names of the metrics stored in the simulation arrays (in this order).
LSEM_METRICS <- c("nSHD", "Cov_joint", "Cov_pair", "Precision_L", "Precision_U",
                  "TPR_L", "TPR_U", "FPR_L", "FPR_U", "FNR_L", "FNR_U")

#' Turn a reps x stats matrix (rows from graph_recovery_counts(), plus nSHD)
#' into a reps x LSEM_METRICS matrix. Rates are NaN when undefined (0/0).
#' Cov_joint is the coverage frequency over all replicates, repeated in every row.
recovery_metrics <- function(counts, nshd) {
  with(as.data.frame(counts), cbind(
    nSHD        = nshd,
    Cov_joint   = mean(covered),
    Cov_pair    = pair_cov,
    Precision_L = TP_L / (TP_L + FP_L), Precision_U = TP_U / (TP_U + FP_U),
    TPR_L       = TP_L / (TP_L + FN_L), TPR_U       = TP_U / (TP_U + FN_U),
    FPR_L       = FP_L / (FP_L + TN_L), FPR_U       = FP_U / (FP_U + TN_U),
    FNR_L       = FN_L / (FN_L + TP_L), FNR_U       = FN_U / (FN_U + TP_U)
  ))
}
