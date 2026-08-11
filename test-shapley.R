test_that("exact Shapley values satisfy the efficiency identity", {
  set.seed(6)
  N <- 12; T <- 6
  unit <- rep(1:N, times = T); time <- rep(1:T, each = N)
  y <- rnorm(N * T)
  X <- matrix(rnorm(N * T * 4), ncol = 4)
  colnames(X) <- paste0("x", 1:4)
  y_fe <- fe_partial_vec(y, unit, time)
  X_fe <- fe_partial_mat(X, unit, time)
  df_fe <- N * T - (N + T - 1)

  phi <- fe_shapley_exact(y_fe, X_fe, df_fe)
  v_N <- fe_model_fit(y_fe, X_fe, 1:4, df_fe, "adj_r2") -
    fe_model_fit(y_fe, X_fe, integer(0), df_fe, "adj_r2")

  expect_equal(sum(phi), v_N, tolerance = 1e-10)
  expect_length(phi, 4)
  expect_named(phi, colnames(X))
})

test_that("MC Shapley approximates exact Shapley for a small game", {
  set.seed(7)
  N <- 12; T <- 6
  unit <- rep(1:N, times = T); time <- rep(1:T, each = N)
  y <- rnorm(N * T)
  X <- matrix(rnorm(N * T * 4), ncol = 4)
  y_fe <- fe_partial_vec(y, unit, time)
  X_fe <- fe_partial_mat(X, unit, time)
  df_fe <- N * T - (N + T - 1)

  phi_exact <- fe_shapley_exact(y_fe, X_fe, df_fe)
  phi_mc <- fe_shapley_mc(y_fe, X_fe, df_fe, n_samples = 500)

  expect_equal(phi_mc, phi_exact, tolerance = 0.05)
})

test_that("null player receives zero Shapley value", {
  set.seed(8)
  N <- 20; T <- 10
  unit <- rep(1:N, times = T); time <- rep(1:T, each = N)
  X <- matrix(rnorm(N * T * 4), ncol = 4)
  colnames(X) <- c("x1", "x2", "x3", "noise")
  y <- 0.8 * X[, 1] - 0.5 * X[, 3] + rnorm(N * T) + unit + time
  y_fe <- fe_partial_vec(y, unit, time)
  X_fe <- fe_partial_mat(X, unit, time)
  df_fe <- N * T - (N + T - 1)

  phi <- fe_shapley_exact(y_fe, X_fe, df_fe)
  # the pure-noise covariate should receive the smallest value
  expect_equal(which.min(unname(phi)), 4)
})

test_that("master routine recovers the active set under Rule B", {
  set.seed(9)
  N <- 40; T <- 15
  unit <- rep(1:N, times = T); time <- rep(1:T, each = N)
  beta <- c(0.8, 0, 0.6, 0, -0.5, 0)
  X <- matrix(rnorm(N * T * 6), ncol = 6)
  colnames(X) <- paste0("x", 1:6)
  y <- as.vector(X %*% beta) + rnorm(N * T) + unit * 0.5 + time * 0.3

  res <- fe_shapley_covariate_selection(y, X, unit, time,
                                        exact = TRUE, gamma = 0.01,
                                        verbose = FALSE)

  active <- which(beta != 0)
  expect_setequal(res$selected_nested, active)
  expect_true(all(res$shapley_values[active] > 0))
  expect_lt(abs(res$efficiency_gap), 1e-10)
})
