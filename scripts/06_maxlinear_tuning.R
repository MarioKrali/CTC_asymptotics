# Choice of the tuning parameters in k = max(floor(c * n^kappa), 10) for the
# bivariate max-linear model: MSE of the rank estimator of Gamma_12 and Gamma_21
# at n = 50000 over a grid of (kappa, c).
#
# Output: results/maxlinear_tuning.rds  list(dir12 = data.frame, dir21 = data.frame)
# Figure: ML_bivariate_choice_c_alpha.pdf

source("R/setup.R")
library(dplyr)
library(ggplot2)
library(patchwork)

a11_vec   <- c(0.25, 0.4, 0.6, 0.85)
kappa_vec <- c(0.3, 0.4, 0.5, 0.6)
c_vec     <- c(0.25, 0.5, 0.75, 1)
# The figure only uses n = 50000. Each n is seeded separately, so adding
# further sample sizes here does not change the n = 50000 results.
n_vec <- if (QUICK) 5000 else 50000
reps  <- if (QUICK) 5 else 1000
B     <- if (QUICK) 20 else 100

jobs <- expand.grid(n = n_vec, direction = c("12", "21"), stringsAsFactors = FALSE)
t0 <- Sys.time()
out <- par_lapply(seq_len(nrow(jobs)), function(i) {
  run_grid_alpha_sc(jobs$n[i], a11_vec, kappa_vec, c_vec, reps = reps, B = B,
                    seed = 123, direction = jobs$direction[i])
})
message("Simulation time: ", format(round(Sys.time() - t0, 1)))

res <- list(dir12 = do.call(rbind, out[jobs$direction == "12"]),
            dir21 = do.call(rbind, out[jobs$direction == "21"]))
save_result(res, "maxlinear_tuning.rds")

# ---------------------------------------------------------------------------
# Figure: MSE of the rank estimator vs c, one line per kappa, one column per theta
# ---------------------------------------------------------------------------
col_kappa <- c("0.3" = "black", "0.4" = "gray35", "0.5" = "gray60", "0.6" = "gray75")

mse_plot <- function(res, strip_label) {
  mse_df <- res |>
    filter(n == max(n_vec), method == "rank") |>
    group_by(a11, alpha, sc) |>
    summarise(mse_m = mean((m_hat - m_theory)^2, na.rm = TRUE), .groups = "drop") |>
    mutate(mse_strip = strip_label)

  a11_values <- as.character(sort(unique(mse_df$a11)))
  a11_labels <- setNames(paste0("theta == ", a11_values), a11_values)

  ggplot(mse_df, aes(x = sc, y = mse_m, colour = factor(alpha), group = factor(alpha))) +
    geom_line(linewidth = .75) +
    geom_point(size = 1) +
    facet_grid(rows = vars(mse_strip), cols = vars(a11),
               labeller = labeller(mse_strip = label_parsed,
                                   a11 = as_labeller(a11_labels, label_parsed))) +
    scale_colour_manual(values = col_kappa, name = expression(kappa)) +
    labs(title = NULL, x = expression(c), y = NULL) +
    theme_bw(base_size = 17) +
    theme(strip.background = element_rect(fill = "gray100"),
          strip.text = element_text(face = "bold"),
          strip.placement = "outside",
          strip.text.y.right = element_text(angle = -90, margin = margin(l = 6, r = 6)),
          legend.position = "right")
}

p_top <- mse_plot(res$dir12, "MSE(hat(Gamma)[12])") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
p_bottom <- mse_plot(res$dir21, "MSE(hat(Gamma)[21])") +
  theme(strip.text.x = element_blank(), strip.background.x = element_blank())

combined_plot <- wrap_plots(p_top, p_bottom, ncol = 1, guides = "collect") &
  theme(legend.position = "right",
        axis.text    = element_text(size = 17),
        axis.title   = element_text(size = 17),
        strip.text   = element_text(size = 17, face = "bold"),
        legend.text  = element_text(size = 17),
        legend.title = element_text(size = 17))

save_fig(combined_plot, "ML_bivariate_choice_c_alpha.pdf", width = 12, height = 6)
