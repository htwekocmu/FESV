test_that("fe_ols matches the two-way within estimator", {
  set.seed(11)
  N <- 12; T <- 7
  unit <- rep(1:N, times = T); time <- rep(1:T, each = N)
  X <- matrix(rnorm(N * T * 3), ncol = 3)
  colnames(X) <- c("x1", "x2", "x3")
  y <- 0.8 * X[, 1] - 0.5 * X[, 3] + rnorm(N * T) + unit * 2 + time

  out <- fe_ols(y, X, unit, time, sel = c(1, 3))

  # benchmark: standard two-way FE (within) regression via dummies
  ref <- lm(y ~ factor(unit) + factor(time) + X[, 1] + X[, 3])
  expect_equal(out$coefficients$Estimate,
               unname(coef(ref)[c("X[, 1]", "X[, 3]")]),
               tolerance = 1e-10)
  # R^2 is net of the fixed effects (partialled data)
  y_fe <- fe_partial_vec(y, unit, time)
  X_fe <- fe_partial_mat(X, unit, time)[, c(1, 3)]
  r2_ref <- summary(lm(y_fe ~ X_fe - 1))$r.squared
  expect_equal(out$r2, r2_ref, tolerance = 1e-10)
  # df = NT - (N + T - 1) - k
  expect_equal(out$df_residual, N * T - (N + T - 1) - 2)
})

test_that("bootstrap_shapley returns sensible CI and se", {
  set.seed(12)
  N <- 10; T <- 6
  unit <- rep(1:N, times = T); time <- rep(1:T, each = N)
  X <- matrix(rnorm(N * T * 4), ncol = 4)
  colnames(X) <- paste0("x", 1:4)
  y <- as.vector(X %*% c(0.8, 0, 0.6, 0)) + rnorm(N * T) + unit + time

  boot <- bootstrap_shapley(y, X, unit, time, B = 50, verbose = FALSE)

  expect_equal(dim(boot$ci), c(2, 4))
  expect_length(boot$se, 4)
  expect_equal(dim(boot$phi_boot), c(50, 4))
  expect_true(all(boot$ci[1, ] <= boot$ci[2, ]))
  expect_true(all(boot$se >= 0))
})
