# Bivariate max-linear model: oracle (known margins) vs rank (unknown margins)
# estimators of Gamma and of its asymptotic variance, plus the multiplier
# bootstrap variance. k = max(floor(c * n^0.4), 10).
#
# Output: results/maxlinear_c_grid.rds  list(dir12 = data.frame, dir21 = data.frame)
#         results/maxlinear_density.rds  data.frame for the density figure
# Figures: scripts/05_maxlinear_figures.R

source("R/setup.R")

a11_vec <- c(0.25, 0.4, 0.6, 0.85)
c_vec   <- c(0.25, 0.5, 0.75, 1)
n_vec   <- c(500, 5000, 50000)
reps    <- if (QUICK) 5 else 1000
B       <- if (QUICK) 50 else 1000

jobs <- expand.grid(n = n_vec, direction = c("12", "21"), stringsAsFactors = FALSE)
t0 <- Sys.time()
out <- par_lapply(seq_len(nrow(jobs)), function(i) {
  run_grid_c(jobs$n[i], a11_vec, c_vec, reps = reps, B = B, seed = 123,
             direction = jobs$direction[i])
})
message("Simulation time: ", format(round(Sys.time() - t0, 1)))

res <- list(dir12 = do.call(rbind, out[jobs$direction == "12"]),
            dir21 = do.call(rbind, out[jobs$direction == "21"]))
save_result(res, "maxlinear_c_grid.rds")

# Sampling distribution of Gamma_12 for a11 = 0.4, n = 5000, c = 0.5
# (B is irrelevant here; only the point estimates are used).
density_res <- run_grid_c(5000, 0.4, 0.5, reps = if (QUICK) 200 else 10000, B = 2,
                          seed = 123, direction = "12")
save_result(density_res, "maxlinear_density.rds")
