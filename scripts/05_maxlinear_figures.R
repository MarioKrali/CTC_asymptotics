# Figures for the bivariate max-linear simulation (scripts/04_maxlinear_simulation.R).
#
#   ML_emp_boot_10x10.pdf, ML_emp21_boot_10x10.pdf  variance estimators (directions 12, 21)
#   ML_Gamma_7x10.pdf, ML_Gamma21_7x10.pdf          Gamma estimators (directions 12, 21)
#   A_04Gamma_7x5.pdf                               sampling density of Gamma_12 (a11 = 0.4)

source("R/setup.R")
library(dplyr)
library(tidyr)
library(ggplot2)

res <- load_result("maxlinear_c_grid.rds", "scripts/04_maxlinear_simulation.R")
density_res <- load_result("maxlinear_density.rds", "scripts/04_maxlinear_simulation.R")

col_c <- c("0.25" = "black", "0.5" = "gray35", "0.75" = "gray60", "1" = "gray75")

#' Boxplots of `diff` by theta (= a11), coloured by c, faceted method x n.
#' Every row (method) gets a symmetric y-range shared by all its columns.
facet_boxplot <- function(df, method_levels, method_lab) {
  df <- df |>
    mutate(a11_f  = factor(a11),
           c_f    = factor(c),
           method = factor(method, levels = method_levels))

  # invisible points forcing a symmetric, row-wise identical y-range (max |diff| + 5 %)
  limits_df <- df |>
    group_by(method) |>
    summarise(lim = max(abs(diff), na.rm = TRUE) * 1.05, .groups = "drop") |>
    tidyr::crossing(n = unique(df$n)) |>
    reframe(diff = c(-lim, lim), .by = c(method, n)) |>
    mutate(a11_f = factor(levels(df$a11_f)[1], levels = levels(df$a11_f)))

  n_lab <- setNames(paste0("n=", sort(unique(df$n))), sort(unique(df$n)))

  ggplot(df, aes(x = a11_f, y = diff, colour = c_f)) +
    geom_blank(data = limits_df, aes(x = a11_f, y = diff), inherit.aes = FALSE) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
    geom_boxplot(position = position_dodge(width = 0.6), width = 0.45,
                 outlier.shape = 1, outlier.size = 1.3, fill = NA, linewidth = 0.5) +
    scale_colour_manual(values = col_c, name = "   c") +
    labs(x = expression(theta), y = NULL) +
    facet_grid(rows = vars(method), cols = vars(n), scales = "free_y",
               labeller = labeller(method = as_labeller(method_lab, label_parsed),
                                   n = as_labeller(n_lab))) +
    theme_bw(base_size = 17) +
    theme(
      panel.grid.major = element_line(colour = scales::alpha("grey40", 0.25), linewidth = 0.3),
      panel.grid.minor = element_line(colour = scales::alpha("grey40", 0.15), linewidth = 0.2),
      strip.placement  = "outside",
      strip.text.y     = element_text(size = 17, margin = margin(l = 6, r = 6)),
      strip.background = element_rect(fill = "grey100", colour = "black"),
      strip.text.x     = element_text(size = 17),
      axis.title.x     = element_text(margin = margin(t = 5)),
      legend.position  = "right",
      legend.key.size  = unit(.8, "cm"),
      panel.spacing    = unit(.5, "lines")
    ) +
    scale_y_continuous(expand = c(0, 0))
}

# ---------------------------------------------------------------------------
# 1. Variance: rank estimator, bootstrap estimator and their difference
# ---------------------------------------------------------------------------
variance_plot <- function(res) {
  rank_rows <- filter(res, method == "rank")
  df <- bind_rows(
    mutate(rank_rows, diff = C_hat  - C_theory, method = "rank"),
    mutate(rank_rows, diff = C_boot - C_theory, method = "boot"),
    mutate(rank_rows, diff = C_hat  - C_boot,   method = "difference")
  )
  facet_boxplot(df, c("rank", "boot", "difference"), c(
    rank       = "Rank:   ~ hat(sigma)^2 - sigma^2",
    boot       = "Bootstrap: ~ hat(sigma)[boot]^2 - sigma^2",
    difference = "Difference: ~ hat(sigma)^2 - hat(sigma)[boot]^2"
  ))
}

save_fig(variance_plot(res$dir12), "ML_emp_boot_10x10.pdf",   width = 10, height = 10)
save_fig(variance_plot(res$dir21), "ML_emp21_boot_10x10.pdf", width = 10, height = 10)

# ---------------------------------------------------------------------------
# 2. Gamma: oracle vs rank estimator
# ---------------------------------------------------------------------------
gamma_plot <- function(res) {
  facet_boxplot(mutate(res, diff = m_hat - m_theory), c("oracle", "rank"), c(
    oracle = "Oracle: ~ hat(Gamma)[known]-Gamma",
    rank   = "Rank:   ~ hat(Gamma)[unknown] - Gamma"
  ))
}

save_fig(gamma_plot(res$dir12), "ML_Gamma_7x10.pdf",   width = 10, height = 7)
save_fig(gamma_plot(res$dir21), "ML_Gamma21_7x10.pdf", width = 10, height = 7)

# ---------------------------------------------------------------------------
# 3. Sampling density of the Gamma_12 estimators with normal approximations
# ---------------------------------------------------------------------------
hist_breaks <- seq(0.1, 1.7, by = 0.01)
gamma_true  <- density_res$m_theory[1]

est <- density_res |>
  mutate(margins = ifelse(method == "oracle", "Known margins", "Unknown margins"))
normal_curves <- est |>
  group_by(margins) |>
  summarise(mu = mean(m_hat), sigma = sd(m_hat), .groups = "drop") |>
  tidyr::crossing(x = hist_breaks) |>
  mutate(y = dnorm(x, mu, sigma))

dens_plot <- ggplot() +
  geom_histogram(data = filter(est, margins == "Known margins"),
                 aes(x = m_hat, y = after_stat(density), fill = margins),
                 breaks = hist_breaks, color = "gray90", alpha = 0.9, linewidth = 0.05) +
  geom_line(data = filter(normal_curves, margins == "Known margins"),
            aes(x = x, y = y), color = "gray70", linewidth = 0.25) +
  geom_histogram(data = filter(est, margins == "Unknown margins"),
                 aes(x = m_hat, y = after_stat(density), fill = margins),
                 breaks = hist_breaks, color = "gray70", alpha = 0.4, linewidth = 0.025) +
  geom_line(data = filter(normal_curves, margins == "Unknown margins"),
            aes(x = x, y = y), color = "gray20", linewidth = 0.25) +
  geom_vline(xintercept = gamma_true, linetype = "dashed", color = "black", linewidth = 0.5) +
  annotate("text", x = gamma_true, y = -Inf, label = "Gamma", vjust = 3, hjust = 0.5,
           color = "black", size = 5.5, parse = TRUE) +
  scale_x_continuous(limits = c(0, 1.7)) +
  scale_fill_manual(name = NULL,
                    values = c("Known margins" = "gray96", "Unknown margins" = "gray70")) +
  labs(x = NULL, y = "Density") +
  theme_minimal(base_size = 17) +
  theme(
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.2),
    panel.grid.major = element_line(color = "gray80", linewidth = 0.2),
    panel.grid.minor = element_line(color = scales::alpha("gray80", 0.5), linewidth = 0.2),
    plot.margin = margin(0.5, 1.5, 1.5, 0.5, "cm"),
    legend.position = "right",
    legend.justification = "center",
    legend.direction = "vertical",
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(size = 17)
  ) +
  coord_cartesian(clip = "off")

save_fig(dens_plot, "A_04Gamma_7x5.pdf", width = 10, height = 5)
