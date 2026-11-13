# Bivariate max-linear model: estimators of Gamma_12 (or Gamma_21) and its variance.
#
# Model: X1 = max(a11 Z1, (1 - a11) Z2), X2 = Z2, Z1 and Z2 i.i.d. unit Frechet.
# "direction 12" estimates Gamma_12 (condition on U extreme), "21" swaps the
# components and conditions on V.

#' n draws of the survival copula (U, V) of the max-linear model.
rMaxLinearCopula <- function(n, a11) {
  a12 <- 1 - a11
  Z1 <- 1 / rexp(n)
  Z2 <- 1 / rexp(n)
  X1 <- pmax(a11 * Z1, a12 * Z2)
  X2 <- Z2
  cbind(U = 1 - exp(-1 / X1), V = 1 - exp(-1 / X2))
}

# Population quantities for direction 12 (for direction 21: Gamma = 1, variance 0)
m_theory     <- function(a11) 1 - a11 / 2                         # Gamma_12
C_oracle_fun <- function(a11) 1 - 2 * a11 / 3                     # E[(1 - V)^2 | U extreme]
C_rank_fun   <- function(a11) C_oracle_fun(a11) - m_theory(a11)^2 # Var(W)

theory_values <- function(a11, direction) {
  if (direction == "21") {
    c(m = 1, C_oracle = 1, C_rank = 0)
  } else {
    c(m = m_theory(a11), C_oracle = C_oracle_fun(a11), C_rank = C_rank_fun(a11))
  }
}

#' Oracle (known margins) plug-in using the threshold U <= k/n.
oracle_plugin_indicator <- function(U, V, k) {
  n <- length(U)
  W <- (1 - V)[U <= k / n]
  list(m_hat = sum(W) / k, C_hat = sum(W^2) / k)
}

#' Rank plug-in using exactly the k smallest ranks of U.
#' `ord` (= order of the ranks of U) and `Vh` (= normalised ranks of V) can be
#' passed in precomputed when several k are evaluated on the same sample.
rank_plugin_exactk <- function(U, V, k,
                               ord = order(rank(U, ties.method = "first")),
                               Vh  = rank(V, ties.method = "first") / (length(V) + 1)) {
  W <- 1 - Vh[ord[seq_len(k)]]
  m_hat <- mean(W)
  list(m_hat = m_hat, C_hat = mean(W^2) - m_hat^2, W = W)
}

#' Multiplier-bootstrap variance of k^{-1/2} sum (W - mean W) xi.
#' Same multiplier draws as calling rnorm(k) B times.
boot_var_rank <- function(W, B = 1000L) {
  k <- length(W)
  xi <- matrix(rnorm(k * B), nrow = k, ncol = B)
  var(colSums((W - mean(W)) * xi) / sqrt(k))
}

#' Oracle and rank estimates (incl. bootstrap variance) for one sample and one k.
#' Returns a 2-row data frame (oracle, rank).
estimate_both <- function(U, V, k, B, ord, Vh) {
  oracle <- oracle_plugin_indicator(U, V, k)
  rk     <- rank_plugin_exactk(U, V, k, ord = ord, Vh = Vh)
  C_boot <- boot_var_rank(rk$W, B = B)
  est <- c(oracle$m_hat, oracle$C_hat, rk$m_hat, rk$C_hat, C_boot)
  if (!all(is.finite(est))) stop("Non-finite estimate (k = ", k, ")")
  data.frame(method = c("oracle", "rank"),
             m_hat  = c(oracle$m_hat, rk$m_hat),
             C_hat  = c(oracle$C_hat, rk$C_hat),
             C_boot = c(NA_real_, C_boot))
}

#' One sample from the max-linear model with the requested orientation.
draw_uv <- function(n, a11, direction) {
  UV <- rMaxLinearCopula(n, a11)
  if (direction == "21") UV <- UV[, c(2, 1)]
  list(U = as.numeric(UV[, 1]), V = as.numeric(UV[, 2]))
}

#' Simulation with k = max(floor(c * n^0.4), 10), a fresh sample for every (a11, c, rep).
#' Returns one row per (a11, c, rep, method).
run_grid_c <- function(n, a11_vec, c_vec, reps, B, seed = 123, direction = "12") {
  set.seed(seed)
  out <- vector("list", length(a11_vec) * length(c_vec) * reps)
  i <- 0L
  for (a11 in a11_vec) {
    for (c_val in c_vec) {
      k <- max(floor(c_val * n^0.4), 10L)
      for (r in seq_len(reps)) {
        s  <- draw_uv(n, a11, direction)
        Vh <- rank(s$V, ties.method = "first") / (n + 1)
        ord <- order(rank(s$U, ties.method = "first"))
        i <- i + 1L
        out[[i]] <- cbind(n = n, a11 = a11, c = c_val, k = k, rep = r,
                          estimate_both(s$U, s$V, k, B, ord, Vh))
      }
    }
  }
  res <- do.call(rbind, out)
  add_theory_columns(res, direction)
}

#' Simulation on a grid k = max(floor(sc * n^alpha), 10). All (alpha, sc) pairs
#' are evaluated on the same sample, so ranks are computed once per sample.
run_grid_alpha_sc <- function(n, a11_vec, alpha_vec, sc_vec, reps, B, seed = 123,
                              direction = "12") {
  set.seed(seed)
  grid <- expand.grid(sc = sc_vec, alpha = alpha_vec)[, c("alpha", "sc")]
  out <- vector("list", length(a11_vec) * reps)
  i <- 0L
  for (a11 in a11_vec) {
    for (r in seq_len(reps)) {
      s   <- draw_uv(n, a11, direction)
      Vh  <- rank(s$V, ties.method = "first") / (n + 1)
      ord <- order(rank(s$U, ties.method = "first"))
      rows <- lapply(seq_len(nrow(grid)), function(g) {
        k <- max(floor(grid$sc[g] * n^grid$alpha[g]), 10L)
        cbind(n = n, a11 = a11, alpha = grid$alpha[g], sc = grid$sc[g], k = k, rep = r,
              estimate_both(s$U, s$V, k, B, ord, Vh))
      })
      i <- i + 1L
      out[[i]] <- do.call(rbind, rows)
    }
  }
  res <- do.call(rbind, out)
  add_theory_columns(res, direction)
}

add_theory_columns <- function(res, direction) {
  th <- t(vapply(res$a11, theory_values, numeric(3), direction = direction))
  res$m_theory <- th[, "m"]
  res$C_theory <- ifelse(res$method == "oracle", th[, "C_oracle"], th[, "C_rank"])
  res
}
