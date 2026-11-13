# Empirical CTC matrix with multiplier-bootstrap confidence graphs.

#' Estimate the CTC matrix and lower/upper ancestral graphs
#'
#' @param ranks_margs n x d matrix of ranks. Small ranks are the extreme
#'   observations: pass ranks of -X to work with the survival copula.
#' @param k number of upper order statistics conditioned on.
#' @param B_rep number of multiplier-bootstrap replicates.
#' @param gamma_threshold threshold on the CTC used to declare an ancestral relation.
#' @param c_var studentisation floor: each standard deviation is floored at
#'   c_var times the largest one, which avoids division by zero.
#' @param alpha level; the (1 - alpha) quantile of the bootstrapped maximum is used.
#'
#' @return A list with
#'   hat_Hn                d x d empirical CTC matrix; entry [j, l] conditions on node j.
#'   var_Hn                d x d empirical variances of the conditional exceedances.
#'   boot_vars             d x d bootstrap variances.
#'   ancestral_graph_lower logical d x d adjacency matrix of the lower graph.
#'   ancestral_graph_upper logical d x d adjacency matrix of the upper graph.
#'
#' The multipliers are drawn at once as a k x B_rep matrix whose column r
#' holds the multipliers of bootstrap replicate r.
tail_moments_dboot <- function(ranks_margs, k, B_rep, gamma_threshold, c_var, alpha) {
  n <- nrow(ranks_margs)
  d <- ncol(ranks_margs)
  u_margs <- ranks_margs / (n + 1)

  hat_Hn    <- matrix(0, d, d)
  var_Hn    <- matrix(0, d, d)
  boot_vars <- matrix(0, d, d)
  boot_mat  <- array(0, dim = c(d, d, B_rep))  # bootstrap process; diagonal stays 0

  xi <- matrix(rnorm(k * B_rep), nrow = k, ncol = B_rep)  # column r = multipliers of replicate r

  for (j in seq_len(d)) {
    # conditional exceedances W: the other margins on the k most extreme observations of margin j
    idx <- order(u_margs[, j])[seq_len(k)]
    W_j <- 1 - u_margs[idx, -j, drop = FALSE]                # k x (d - 1)

    hat_Hn[j, -j] <- matrixStats::colMeans2(W_j)
    var_Hn[j, -j] <- matrixStats::colVars(W_j)

    centred <- sweep(W_j, 2, hat_Hn[j, -j])
    boot_j  <- k^(-0.5) * crossprod(centred, xi)             # (d - 1) x B_rep
    boot_mat[j, -j, ] <- boot_j
    boot_vars[j, -j]  <- matrixStats::rowVars(boot_j)
  }

  # studentised maximum of the bootstrap process, one value per replicate
  sd_floor    <- pmax(boot_vars^0.5, c_var * max(boot_vars)^0.5)
  max_process <- apply(abs(boot_mat) / as.vector(sd_floor), 3, max)
  quant_alpha <- quantile(max_process, 1 - alpha)

  half_width <- quant_alpha * boot_vars^0.5 / sqrt(k)
  list(
    hat_Hn                = hat_Hn,
    var_Hn                = var_Hn,
    boot_vars             = boot_vars,
    ancestral_graph_lower = (hat_Hn - half_width > gamma_threshold),
    ancestral_graph_upper = (hat_Hn + half_width > gamma_threshold)
  )
}

#' Ranks of -X column-wise, so that rank 1 is the largest observation.
reverse_ranks <- function(X) {
  apply(-X, 2, function(x) rank(x, ties.method = "first"))
}
