# FESV 1.0.0

* Initial release of FESV: fixed effects covariate selection using Shapley values.
* Exact Shapley values for p <= 20 candidate covariates and Monte Carlo
  permutation sampling for larger p.
* Two data-driven selection rules: Rule A (positive Shapley value) and
  Rule B (nested penalised criterion J(k) = v(M_k) - gamma*k).
* Residual-bootstrap confidence intervals and standard errors for Shapley values.
* Game-theoretic diagnostics: superadditivity and core-membership checks,
  efficiency identity.
* Fixed-effects OLS estimation on a selected coalition with corrected
  standard errors.
* Selection metrics (precision, recall, F1) and three worked panel datasets
  (developing, developed, all countries, 2000-2023).
