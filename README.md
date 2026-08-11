# FESV

**Fixed Effects Covariate Selection Using Shapley Values**

`FESV` implements a cooperative game-theoretic framework for covariate
selection in two-way fixed effects (FE) panel data models. The method is fully
frequentist: unit and time effects are held fixed (partialled out), and the
candidate covariates are the players of a cooperative game whose characteristic
function is the improvement in a fit statistic (adjusted R^2 by default) over
the FE-only baseline.

## Installation

You can install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("htwekocmu/FESV")
```

Or build and install the package locally from this repository:

```r
# in the FESV/ directory
devtools::install()
```

## Quick start

```r
library(FESV)

# Example data: simulated panel with active set {x1, x3, x5}
set.seed(1)
N <- 40; T <- 15
unit <- rep(1:N, times = T); time <- rep(1:T, each = N)
beta <- c(0.8, 0, 0.6, 0, -0.5, 0)
X <- matrix(rnorm(N * T * length(beta)), ncol = length(beta))
colnames(X) <- paste0("x", 1:length(beta))
y <- as.vector(X %*% beta) + rnorm(N * T) + unit * 0.5 + time * 0.3

res <- fe_shapley_covariate_selection(y, X, unit, time,
                                      exact = TRUE, gamma = 0.01)

res$shapley_values      # average marginal contribution of each covariate
res$selected_sign       # Rule A: phi_j > 0
res$selected_nested     # Rule B: argmax_k J(k) = v(M_k) - gamma*k
```

## Main functions

| Function | Description |
|---|---|
| `fe_shapley_covariate_selection()` | Master routine: Shapley values, diagnostics, Rule A and Rule B |
| `fe_shapley_exact()` | Exact Shapley values (p <= 20) |
| `fe_shapley_mc()` | Monte Carlo (permutation) Shapley values for larger p |
| `fe_nested_selection()` | Rule B nested penalised criterion |
| `bootstrap_shapley()` | Residual-bootstrap CIs and SEs for Shapley values |
| `fe_ols()` | Fixed-effects OLS on a selected coalition |
| `fe_check_superadditivity()`, `fe_check_core_membership()` | Game-theoretic diagnostics |
| `selection_metrics()` | Precision / recall / F1 against a known active set |

## Data

Three balanced panels (World Bank WDI, 2000-2023) ship with the package:
`developing_panel` (25 countries), `developed_panel` (20 countries), and
`all_panel` (pooled, 45 countries).

## License

GPL (>= 3)

## Citation

Ko, H. (2026). *Shapley-Value Covariate Selection for Panel Data Models: A
Cooperative Game-Theoretic Approach.* Working paper.

```r
citation("FESV")
```
