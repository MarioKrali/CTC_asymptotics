# Causal and conditional tail coefficients: asymptotics, inference and graph recovery

R code to reproduce the simulations, figures and data application of

> Axel Bücher and Mario Krali (2026). *Causal and conditional tail coefficients:
> asymptotics, inference and graph recovery.* arXiv preprint (link to follow).

The paper establishes functional central limit theorems for rank-based
estimators of the causal tail coefficient (CTC) and of more general conditional
tail coefficients, together with multiplier bootstrap versions. The code
illustrates the theory in a bivariate max-linear model. It also uses the
bootstrap to construct simultaneous *lower* and *upper* confidence graphs that
envelop the ancestral graph of an ancestral tail mixture model,
applied to dietary intake data from NHANES.

## Repository layout

```
R/                       functions (sourced by every script via R/setup.R)
  ctc_estimator.R        tail_moments_dboot(): CTC matrix + bootstrap lower/upper graphs
  lsem.R                 LSEM simulation, population CTC, graph-recovery metrics
  maxlinear.R            bivariate max-linear model: oracle / rank / bootstrap estimators
  graph_distance.R       normalised structural Hamming distance (nSHD)
  utils.R, setup.R       paths, saving, optional parallelisation
scripts/
  01_lsem_simulation.R        LSEM graph recovery simulation      -> results/lsem_graph_coverage.rds
  02_lsem_figures_tables.R    LSEM figures + LaTeX tables          (needs 01)
  03_lsem_tuning_grid.R       LSEM coverage over a (kappa, c) grid for k
  04_maxlinear_simulation.R   bivariate max-linear simulation      -> results/maxlinear_*.rds
  05_maxlinear_figures.R      max-linear figures                    (needs 04)
  06_maxlinear_tuning.R       MSE over a (kappa, c) grid for k
  07_data_application.R       NHANES dietary nutrients
data/                    nhanes_dr1tot.csv (see data/README.md)
figures/, results/       created when the scripts run (not tracked by git)
```

## Requirements

R (≥ 4.1) and the CRAN packages

```r
install.packages(c("matrixStats", "dplyr", "tidyr", "ggplot2", "patchwork", "cowplot", "scales"))
```

(`lattice`, `parallel` and `grid` ship with R.) The code was last run with
R 4.5.1 and ggplot2 4.0.0.

## Running

Run all scripts **from the repository root**:

```bash
Rscript scripts/01_lsem_simulation.R
Rscript scripts/02_lsem_figures_tables.R
# ... and so on, or everything at once:
Rscript run_all.R
```

Two environment variables control the runs:

| Variable      | Effect |
|---------------|--------|
| `CTC_CORES=8` | Spread the simulation configurations over 8 cores (Unix/macOS only). Every configuration seeds itself, so the results are **identical** for any number of cores. |
| `CTC_QUICK=1` | Very small settings as a smoke test (seconds). Output goes to `results/quick/` and `figures/quick/`. **Not** the paper results. |

Approximate single-core run times on a laptop: 01 ≈ 13 min, 03 ≈ 15 min,
04 ≈ 10 min, 06 ≈ 5 min, 07 < 1 min.

## Reproducibility

All simulations use fixed seeds (`112233` for the LSEM studies, `123` for the
max-linear studies, `12123` for the data application), so rerunning the
scripts reproduces the figures and tables of the paper exactly, regardless of
the number of cores used.

## Citation

```bibtex
@misc{BuecherKrali2026,
  title         = {Causal and conditional tail coefficients: asymptotics, inference and graph recovery},
  author        = {B{\"u}cher, Axel and Krali, Mario},
  year          = {2026},
  eprint        = {XXXX.XXXXX},
  archivePrefix = {arXiv},
  primaryClass  = {math.ST}
}
```