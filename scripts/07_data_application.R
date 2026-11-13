# Data application: extremal ancestral graph of six NHANES dietary nutrients
# (data/nhanes_dr1tot.csv, columns RET, AC, LZ, VA, BC, VK; see data/README.md).
#
# The CTC graph is estimated for k = 10, 20, ..., 210. Its stability across k is
# measured with the nSHD between graphs at neighbouring k.
#
# Figures: combined_mean_sds.pdf         rolling mean / SD of neighbouring nSHD
#          nhanes_graphs_k80.pdf         estimated, lower and upper graph at k = 80
#          nhanes_nshd_heatmap.pdf       nSHD between all pairs of k
#          nhanes_graph_centrality.pdf   mean nSHD of each k to all others
#          nhanes_edge_frequency.pdf     how often each edge is selected over the k grid

source("R/setup.R")
library(dplyr)
library(ggplot2)

food <- read.csv(file.path(DATA_DIR, "nhanes_dr1tot.csv"), header = TRUE)
# 1,217 respondents have no dietary recall (all six values missing); keep complete
# cases only (8,327 rows). Otherwise rank() would place the NAs last and they
# would still count towards n.
n_raw <- nrow(food)
food  <- food[complete.cases(food), ]
message("Using ", nrow(food), " complete rows (", n_raw - nrow(food), " incomplete rows dropped).")
node_names <- colnames(food)
ranks <- reverse_ranks(as.matrix(food))

gamma_threshold <- 0.95
B <- if (QUICK) 100 else 1000
k_choices <- 10 * (1:21)
n_k <- length(k_choices)

set.seed(12123)
fits <- lapply(k_choices, function(k) tail_moments_dboot(ranks, k, B, gamma_threshold, 0.1, 0.1))
names(fits) <- k_choices
save_result(fits, "nhanes_fits.rds")

# Point-estimate graphs; entry [j, l] = 1 means j -> l (l is extreme when j is)
graphs <- lapply(fits, function(f) 1 * (f$hat_Hn >= gamma_threshold))

# ---------------------------------------------------------------------------
# 1. Stability across k: nSHD between neighbouring k, rolling windows of 5
# ---------------------------------------------------------------------------
d_next <- vapply(seq_len(n_k - 1), function(i) nSHD(graphs[[i]], graphs[[i + 1]]), numeric(1))
n_windows <- n_k - 5   # number of complete windows of 5 neighbouring distances
window_stats <- data.frame(
  i    = seq_len(n_windows),
  mean = vapply(seq_len(n_windows), function(i) mean(d_next[i:(i + 4)]), numeric(1)),
  sd   = vapply(seq_len(n_windows), function(i) sd(d_next[i:(i + 4)]), numeric(1))
)

stability_plot <- function(y, ylab) {
  ggplot(window_stats, aes(x = i, y = .data[[y]])) +
    geom_line(color = "black", linewidth = 0.7) +
    geom_point(color = "black", size = 2) +
    scale_x_continuous(breaks = window_stats$i) +
    labs(x = expression(i), y = ylab, title = NULL) +
    theme_bw(base_size = 17) +
    theme(panel.grid.minor = element_blank())
}

combined_mean_sds <- cowplot::plot_grid(stability_plot("mean", "Mean nSHD"),
                                        stability_plot("sd", "SD of nSHD"))
save_fig(combined_mean_sds, "combined_mean_sds.pdf", width = 10, height = 3.5)

# ---------------------------------------------------------------------------
# 2. Estimated graph and its lower / upper confidence graphs at k = 80
# ---------------------------------------------------------------------------
fit80 <- fits[["80"]]
pdf(file.path(FIG_DIR, "nhanes_graphs_k80.pdf"), width = 6, height = 5)
label_graph <- function(G) { G <- 1 * G; dimnames(G) <- list(node_names, node_names); G }
print(lattice::levelplot(label_graph(graphs[["80"]]), main = "Estimated graph (k = 80)"))
print(lattice::levelplot(label_graph(fit80$ancestral_graph_lower), main = "Lower graph (k = 80)"))
print(lattice::levelplot(label_graph(fit80$ancestral_graph_upper), main = "Upper graph (k = 80)"))
invisible(dev.off())

# ---------------------------------------------------------------------------
# 3. nSHD between all pairs of k, and mean distance to all other graphs
# ---------------------------------------------------------------------------
D <- outer(seq_len(n_k), seq_len(n_k),
           Vectorize(function(i, j) nSHD(graphs[[i]], graphs[[j]])))

d_long <- data.frame(k_i = factor(rep(k_choices, times = n_k), levels = k_choices),
                     k_j = factor(rep(k_choices, each = n_k), levels = k_choices),
                     nSHD = as.vector(D))

p_heat <- ggplot(d_long, aes(x = k_i, y = k_j, fill = nSHD)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient(low = "white", high = "black", name = "nSHD") +
  coord_equal() +
  labs(x = "Threshold k", y = "Threshold k") +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p_heat, "nhanes_nshd_heatmap.pdf", width = 7, height = 6)

score_df <- data.frame(k = factor(k_choices, levels = k_choices),
                       mean_nshd = rowSums(D) / (n_k - 1))
p_score <- ggplot(score_df, aes(x = k, y = mean_nshd, group = 1)) +
  geom_line(color = "black", linewidth = 0.5) +
  geom_point(color = "black", size = 2) +
  geom_point(data = filter(score_df, mean_nshd == min(mean_nshd)),
             color = "white", fill = "black", shape = 21, size = 3) +
  labs(x = "Threshold k", y = "Mean nSHD to other graphs") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p_score, "nhanes_graph_centrality.pdf", width = 7, height = 4)

# ---------------------------------------------------------------------------
# 4. Edge frequencies over the k grid (out of 21)
# ---------------------------------------------------------------------------
freq <- Reduce(`+`, graphs)   # rows = sources, columns = destinations
dimnames(freq) <- list(node_names, node_names)
print(freq)

# Grey level / line width per frequency band (darker and wider = more frequent)
bands <- data.frame(
  band  = c("21", "18–19", "15", "5", "1"),
  lo    = c(21, 18, 15, 5, 1),
  hi    = c(21, 19, 15, 5, 1),
  grey  = c("black", "grey20", "grey50", "grey75", "grey85"),
  width = c(1.15, 1.00, 0.90, 0.75, 0.65)
)

nodes <- data.frame(name = node_names,
                    x = c(-2, 0, 2, -2, 0, 2),
                    y = c(0, 0, 0, -2, -2, -2))

ij <- which(freq > 0, arr.ind = TRUE)
edges <- data.frame(from_i = ij[, 1], to_i = ij[, 2], frequency = freq[ij]) |>
  filter(from_i != to_i) |>
  mutate(band_idx = vapply(frequency, function(f) {
    hit <- which(f >= bands$lo & f <= bands$hi)
    if (length(hit) == 1) hit else NA_integer_
  }, integer(1)))
if (anyNA(edges$band_idx)) {
  stop("Edge frequencies ", paste(unique(edges$frequency[is.na(edges$band_idx)]), collapse = ", "),
       " are not covered by `bands`; please extend it.")
}

offset   <- 0.09  # sideways shift of reciprocal arrows
node_gap <- 0.13  # gap between arrow ends and node centres
edges <- edges |>
  mutate(
    band = factor(bands$band[band_idx], levels = bands$band),
    x0 = nodes$x[from_i], y0 = nodes$y[from_i],
    x1 = nodes$x[to_i],   y1 = nodes$y[to_i],
    reciprocal = freq[cbind(to_i, from_i)] > 0,
    dx = x1 - x0, dy = y1 - y0, len = sqrt(dx^2 + dy^2),
    sign = if_else(from_i < to_i, 1, -1),
    nx = if_else(reciprocal, sign * (-dy / len) * offset, 0),
    ny = if_else(reciprocal, sign * (dx / len) * offset, 0),
    x_start = x0 + nx + node_gap * dx / len,
    y_start = y0 + ny + node_gap * dy / len,
    x_end   = x1 + nx - node_gap * dx / len,
    y_end   = y1 + ny - node_gap * dy / len
  )

p_edges <- ggplot() +
  geom_segment(data = edges,
               aes(x = x_start, y = y_start, xend = x_end, yend = y_end,
                   colour = band, linewidth = band),
               arrow = grid::arrow(type = "closed", length = unit(2.5, "mm"), angle = 25),
               lineend = "round") +
  geom_point(data = nodes, aes(x = x, y = y),
             shape = 21, size = 10, stroke = 0.8, fill = "white", colour = "black") +
  geom_text(data = nodes, aes(x = x, y = y, label = name), size = 3) +
  scale_colour_manual(name = sprintf("Frequency\n(out of %d)", n_k),
                      values = setNames(bands$grey, bands$band), drop = FALSE) +
  scale_linewidth_manual(values = setNames(bands$width, bands$band), guide = "none") +
  coord_equal(xlim = c(-2.7, 3.4), ylim = c(-2.6, 0.6), clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 12) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        plot.margin = margin(12, 18, 12, 12))
save_fig(p_edges, "nhanes_edge_frequency.pdf", width = 7, height = 4)
