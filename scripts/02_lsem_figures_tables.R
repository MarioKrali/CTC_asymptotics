# Figures and LaTeX tables for the LSEM simulation (scripts/01_lsem_simulation.R).
#
# Figures: LSEM_metrics_p01_11x10.pdf, LSEM_metrics_p02_11x10.pdf,
#          LSEM_coverage_10x6.pdf, LSEM_nshd_12x4.pdf
# Tables:  results/tables/lsem_table_p01.tex, lsem_table_p02.tex

source("R/setup.R")
library(dplyr)
library(ggplot2)

agg_res <- load_result("lsem_graph_coverage.rds", "scripts/01_lsem_simulation.R")

d_levels <- dimnames(agg_res)[[3]]
p_levels <- dimnames(agg_res)[[4]]
n_levels <- dimnames(agg_res)[[5]]
grey2    <- setNames(c("gray25", "gray70"), p_levels)
p_names  <- paste0("0", p_levels)  # ".1" -> "0.1" for legends

# Long format: one row per (replicate, metric, d, p, n)
long_df <- as.data.frame(as.table(agg_res), stringsAsFactors = FALSE) |>
  rename(replicate = Var1, metric = Var2, d = Var3, p = Var4, sample_size = Var5, value = Freq) |>
  mutate(d = factor(d, levels = d_levels),
         p = factor(p, levels = p_levels),
         sample_size = factor(sample_size, levels = n_levels))

theme_lsem <- theme_bw(base_size = 17) +
  theme(panel.spacing = unit(.5, "lines"),
        strip.background = element_rect(fill = "white"),
        strip.text.y = element_text(angle = -90, size = 17),
        legend.position = "right")

# ---------------------------------------------------------------------------
# 1. Precision / TPR / FPR / FNR of the lower and upper graphs (one figure per p)
# ---------------------------------------------------------------------------
rate_df <- long_df |>
  filter(metric %in% c("Precision_L", "Precision_U", "TPR_L", "TPR_U",
                       "FPR_L", "FPR_U", "FNR_L", "FNR_U")) |>
  mutate(metric_group = sub("_[LU]$", "", metric),
         bound_type   = ifelse(endsWith(metric, "_L"), "Lower", "Upper")) |>
  filter(is.finite(value), value >= 0, value <= 1)   # drops undefined (0/0) rates

rate_means <- rate_df |>
  group_by(metric_group, d, p, sample_size, bound_type) |>
  summarise(mean_value = mean(value), .groups = "drop")

bound_cols <- c("Lower" = "gray25", "Upper" = "gray70")
dodge <- position_dodge(width = 0.5)

for (p_val in p_levels) {
  pl <- ggplot() +
    geom_boxplot(data = filter(rate_df, p == p_val),
                 aes(x = sample_size, y = value, fill = bound_type),
                 width = 0.4, outlier.shape = NA, alpha = 0.5, position = dodge) +
    geom_line(data = filter(rate_means, p == p_val),
              aes(x = sample_size, y = mean_value, color = bound_type, group = bound_type),
              linewidth = 0.25, position = dodge) +
    geom_point(data = filter(rate_means, p == p_val),
               aes(x = sample_size, y = mean_value, color = bound_type, group = bound_type),
               size = 2, shape = 18, position = dodge) +
    scale_x_discrete(name = "Sample size (n)") +
    scale_y_continuous(name = "Metric Value", limits = c(0, 1)) +
    scale_fill_manual(values = bound_cols, name = "Graph type") +
    scale_color_manual(values = bound_cols, name = "Graph type") +
    facet_grid(metric_group ~ d, labeller = labeller(d = function(d) paste0("d=", d))) +
    theme_lsem +
    theme(axis.text = element_text(size = 17), strip.text.x = element_text(size = 17))

  save_fig(pl, sprintf("LSEM_metrics_p%s_11x10.pdf", sub(".", "0", p_val, fixed = TRUE)),
           width = 10, height = 11)
}

# ---------------------------------------------------------------------------
# 2. Joint and pairwise coverage
# ---------------------------------------------------------------------------
cov_df <- long_df |>
  filter(metric %in% c("Cov_joint", "Cov_pair")) |>
  mutate(metric = factor(metric, levels = c("Cov_joint", "Cov_pair")))

cov_means <- cov_df |>
  group_by(metric, d, p, sample_size) |>
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop")

p_coverage <- ggplot() +
  # joint coverage: one value per configuration
  geom_line(data = filter(cov_means, metric == "Cov_joint"),
            aes(x = sample_size, y = mean_value, group = interaction(d, p), color = p),
            linewidth = 0.75) +
  geom_point(data = filter(cov_means, metric == "Cov_joint"),
             aes(x = sample_size, y = mean_value, group = interaction(d, p), color = p),
             size = 2, shape = 18) +
  # pairwise coverage: distribution over replicates plus mean
  geom_boxplot(data = filter(cov_df, metric == "Cov_pair"),
               aes(x = sample_size, y = value, fill = p),
               width = 0.4, position = dodge, outlier.shape = NA, alpha = 0.5) +
  geom_line(data = filter(cov_means, metric == "Cov_pair"),
            aes(x = sample_size, y = mean_value, group = interaction(d, p), color = p),
            linewidth = 0.25, position = dodge) +
  geom_point(data = filter(cov_means, metric == "Cov_pair"),
             aes(x = sample_size, y = mean_value, group = interaction(d, p), color = p),
             size = 2, shape = 18, position = dodge) +
  facet_grid(metric ~ d, scales = "free_y",
             labeller = labeller(
               metric = c(Cov_joint = "Joint coverage", Cov_pair = "Pairwise coverage"),
               d = function(x) paste0("d = ", x))) +
  scale_x_discrete(name = "Sample size (n)") +
  scale_y_continuous(name = NULL) +
  scale_color_manual(name = "p", values = grey2, labels = p_names) +
  scale_fill_manual(name = "p", values = grey2, labels = p_names) +
  theme_lsem +
  theme(axis.text = element_text(size = 17),
        strip.text.x = element_text(size = 17),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.key.size = unit(1, "lines"),
        legend.text = element_text(size = 17))

save_fig(p_coverage, "LSEM_coverage_10x6.pdf", width = 10, height = 7)

# ---------------------------------------------------------------------------
# 3. nSHD between the CTC graph (Gamma >= 0.95) and the true ancestral graph.
#    It depends only on the data-generating DAG, not on n.
# ---------------------------------------------------------------------------
nshd_df <- filter(long_df, metric == "nSHD")
nshd_means <- nshd_df |>
  group_by(metric, d, p) |>
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop")

plot_nSHD <- ggplot() +
  geom_boxplot(data = nshd_df, aes(x = p, y = value, fill = p),
               width = 0.5, position = position_dodge(width = 0.7),
               outlier.shape = NA, alpha = 0.6) +
  geom_point(data = nshd_means, aes(x = p, y = mean_value),
             size = 2, shape = 18, position = position_dodge(width = 0.7)) +
  facet_grid(metric ~ d, scales = "free_y",
             labeller = labeller(d = function(x) paste0("d = ", x))) +
  scale_x_discrete(name = "p", labels = p_names) +
  scale_y_continuous(name = NULL) +
  scale_fill_manual(name = "p", values = grey2, labels = p_names) +
  theme_lsem +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

save_fig(plot_nSHD, "LSEM_nshd_12x4.pdf", width = 10, height = 4)

# ---------------------------------------------------------------------------
# 4. LaTeX tables (means over replicates, in %)
# ---------------------------------------------------------------------------
latex_table <- function(agg_res, p_value) {
  agg_mean <- apply(agg_res, c(2, 3, 4, 5), mean, na.rm = TRUE)
  metrics <- c("Cov_joint", "Cov_pair", "Precision_L", "Precision_U",
               "TPR_L", "TPR_U", "FPR_L", "FPR_U", "FNR_L", "FNR_U", "nSHD")
  panel <- if (p_value == p_levels[1]) "A" else "B"

  header <- c(
    "\\begin{table}[t]", "\\centering",
    "\\caption{Performance metrics (in \\%) for different model configurations.}",
    "\\label{tab:graph-recovery}", "\\small", "\\setlength{\\tabcolsep}{4pt}",
    "\\begin{tabular}{@{}cc*{11}{r}@{}}", "\\toprule",
    " & & \\multicolumn{2}{c}{Coverage}", "   & \\multicolumn{2}{c}{Precision}",
    "   & \\multicolumn{2}{c}{TPR}", "   & \\multicolumn{2}{c}{FPR}",
    "   & \\multicolumn{2}{c}{FNR}", "   & \\multicolumn{1}{c}{nSHD} \\\\",
    "\\cmidrule(lr){3-4}", "\\cmidrule(lr){5-6}", "\\cmidrule(lr){7-8}",
    "\\cmidrule(lr){9-10}", "\\cmidrule(lr){11-12}", "\\cmidrule(l){13-13}",
    "$d$ & $n$",
    "& ${\\mathrm{Cov}}_{\\mathrm{joint}}$", "& ${\\mathrm{Cov}}_{\\mathrm{pair}}$",
    "& $\\mathrm{Prec}_L$", "& $\\mathrm{Prec}_U$", "& $\\mathrm{TPR}_L$", "& $\\mathrm{TPR}_U$",
    "& $\\mathrm{FPR}_L$", "& $\\mathrm{FPR}_U$", "& $\\mathrm{FNR}_L$", "& $\\mathrm{FNR}_U$",
    "& $\\mathrm{nSHD}$ \\\\", "\\midrule",
    paste0("\\multicolumn{13}{@{}l}{\\textit{Panel ", panel, ": $p=", p_value, "$}} \\\\\\midrule")
  )

  body <- unlist(lapply(d_levels, function(d) {
    rows <- vapply(n_levels, function(n) {
      vals <- agg_mean[metrics, d, p_value, n]
      cells <- ifelse(is.na(vals), "---", sprintf("%.1f", 100 * vals))
      paste0(" & $", n, "$ & ", paste(cells, collapse = " & "), " \\\\")
    }, character(1))
    c(paste0("\\multirow{3}{*}{$", d, "$}"), rows, "\\addlinespace")
  }))

  c(header, body, "\\bottomrule", "\\end{tabular}", "\\end{table}")
}

for (p_val in p_levels) {
  tex  <- latex_table(agg_res, p_val)
  path <- file.path(RES_DIR, "tables", sprintf("lsem_table_p%s.tex", sub(".", "0", p_val, fixed = TRUE)))
  writeLines(tex, path)
  message("Saved ", path)
}
