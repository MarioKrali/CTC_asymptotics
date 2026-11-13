# LSEM simulation over a grid of tuning parameters k = max(floor(c * n^kappa), 10),
# kappa in {0.3, 0.4, 0.5, 0.6}, c in {0.25, 0.5, 0.75, 1}, with B = 100 bootstrap
# replicates. All (kappa, c) pairs are evaluated on the same simulated data set.
#
# Output: results/lsem_tuning_grid.rds, an array [replicate, metric, d, p, n, kappa, c]
# Figure: plot_covgrid.pdf (joint coverage across all (d, p, n), per (kappa, c))

source("R/setup.R")
library(dplyr)
library(ggplot2)

d_vec     <- c(5, 10, 20)
p_vec     <- c(0.1, 0.2)
n_vec     <- c(500, 5000, 50000)
kappa_vec <- c(0.3, 0.4, 0.5, 0.6)
c_vec     <- c(0.25, 0.5, 0.75, 1)
reps <- if (QUICK) 4 else 100
B    <- if (QUICK) 20 else 100
gamma_threshold <- 0.95

run_config <- function(d, p, n) {
  set.seed(112233)
  n_k <- length(kappa_vec) * length(c_vec)
  counts <- array(NA_real_, dim = c(reps, n_k, 10))
  nshd   <- numeric(reps)
  for (g in seq_len(reps)) {
    A      <- simulate_lsem_A(d, p)
    true_A <- 1 * (ctc_population(A) >= gamma_threshold)
    nshd[g] <- nshd_ctc_vs_ancestral(true_A, A)
    ranks  <- reverse_ranks(simulate_lsem_X(A, n))
    j <- 0L
    for (kappa in kappa_vec) {
      for (c_val in c_vec) {
        j <- j + 1L
        fit <- tail_moments_dboot(ranks, max(floor(c_val * n^kappa), 10), B,
                                  gamma_threshold, 0.1, 0.1)
        counts[g, j, ] <- graph_recovery_counts(true_A,
                                                low = t(fit$ancestral_graph_lower),
                                                up  = t(fit$ancestral_graph_upper))
      }
    }
  }
  cnt_names <- names(graph_recovery_counts(diag(0, 2), diag(0, 2), diag(0, 2)))
  # [reps, metric, kappa, c]; c varies fastest in the loop above
  res <- vapply(seq_len(n_k), function(j) {
    m <- counts[, j, ]; colnames(m) <- cnt_names
    recovery_metrics(m, nshd)
  }, matrix(0, reps, length(LSEM_METRICS)))
  array(res, dim = c(reps, length(LSEM_METRICS), length(c_vec), length(kappa_vec))) |>
    aperm(c(1, 2, 4, 3))
}

configs <- expand.grid(n = n_vec, p = p_vec, d = d_vec)
t0 <- Sys.time()
out <- par_lapply(seq_len(nrow(configs)), function(i) {
  with(configs[i, ], run_config(d, p, n))
})
message("Simulation time: ", format(round(Sys.time() - t0, 1)))

agg_res <- array(NA_real_,
                 dim = c(reps, length(LSEM_METRICS), length(d_vec), length(p_vec), length(n_vec),
                         length(kappa_vec), length(c_vec)),
                 dimnames = list(NULL, LSEM_METRICS, as.character(d_vec), p_label(p_vec),
                                 as.character(n_vec), p_label(kappa_vec), p_label(c_vec)))
for (i in seq_len(nrow(configs))) {
  cf <- configs[i, ]
  agg_res[, , as.character(cf$d), p_label(cf$p), as.character(cf$n), , ] <- out[[i]]
}
save_result(agg_res, "lsem_tuning_grid.rds")

# ---------------------------------------------------------------------------
# Figure: joint coverage (in %) per (kappa, c); boxplot over the 18 (d, p, n)
# configurations, highlighted point = d = 20, p = 0.2, n = 50000.
# ---------------------------------------------------------------------------
cov_df <- as.data.frame.table(agg_res[1, "Cov_joint", , , , , ], responseName = "value",
                              stringsAsFactors = FALSE)
names(cov_df) <- c("d", "p", "n", "kappa", "c", "value")
cov_df <- mutate(cov_df, value = 100 * value, metric = "Joint coverage")

# As in the figure of the paper, each configuration's value is repeated once
# per Monte Carlo replication (reps times) before computing the boxplot statistics.
box_df <- cov_df[rep(seq_len(nrow(cov_df)), each = reps), ]

highlight <- filter(cov_df, n == "50000", p == ".2", d == "20")

plot_grid_coverage <- ggplot(box_df, aes(x = c, y = value, group = c)) +
  geom_boxplot(outlier.shape = NA, width = .4, fill = "gray75", alpha = .2, color = "black") +
  geom_point(data = highlight, aes(x = c, y = value),
             color = "gray50", shape = 18, size = 2, position = position_dodge(width = 0.8)) +
  facet_grid(rows = vars(metric), cols = vars(kappa),
             labeller = labeller(kappa = as_labeller(function(x) paste0("kappa == 0", x),
                                                     default = label_parsed))) +
  scale_x_discrete(breaks = p_label(c_vec), labels = as.character(c_vec)) +
  labs(y = NULL, x = "c") +
  theme_bw(base_size = 17, base_line_size = 0.5) +
  theme(
    panel.spacing = unit(.5, "lines"),
    strip.background = element_rect(fill = "gray100", color = "gray50"),
    strip.text = element_text(size = 17, margin = margin(10, 0, 10, 0)),
    strip.text.y = element_text(size = 17, margin = margin(0, 10, 0, 10)),
    axis.text = element_text(size = 17),
    axis.text.x = element_text(angle = 0, hjust = .5, vjust = 1),
    panel.background = element_rect(fill = "white"),
    plot.margin = margin(1, 1, 2, 1, "cm")
  )

save_fig(plot_grid_coverage, "plot_covgrid.pdf", width = 11, height = 4.5)
