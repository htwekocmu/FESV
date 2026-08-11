test_that("fe_model_fit returns 0 for the empty coalition relative to itself", {
  set.seed(3)
  N <- 10; T <- 5
  unit <- rep(1:N, times = T); time <- rep(1:T, each = N)
  y <- rnorm(N * T)
  X <- matrix(rnorm(N * T * 3), ncol = 3)
  y_fe <- fe_partial_vec(y, unit, time)
  X_fe <- fe_partial_mat(X, unit, time)
  df_fe <- N * T - (N + T - 1)

  v_empty <- fe_model_fit(y_fe, X_fe, integer(0), df_fe, "adj_r2")
  expect_equal(fe_characteristic_function(integer(0), y_fe, X_fe, df_fe), 0,
               tolerance = 1e-12)
  # baseline is the FE-only adjusted R^2
  expect_gt(v_empty, -Inf)
})

test_that("adj_r2 penalises degrees of freedom", {
  set.seed(4)
  N <- 15; T <- 8
  unit <- rep(1:N, times = T); time <- rep(1:T, each = N)
  y <- rnorm(N * T)
  X <- matrix(rnorm(N * T * 4), ncol = 4)
  y_fe <- fe_partial_vec(y, unit, time)
  X_fe <- fe_partial_mat(X, unit, time)
  df_fe <- N * T - (N + T - 1)

  fit0 <- fe_model_fit(y_fe, X_fe, integer(0), df_fe, "adj_r2")
  fit4 <- fe_model_fit(y_fe, X_fe, 1:4, df_fe, "adj_r2")
  # pure noise covariates must not improve adjusted R^2 on average
  expect_true(fit4 < fit0 + 0.01)

  # plain r2 never decreases when adding covariates
  r0 <- fe_model_fit(y_fe, X_fe, integer(0), df_fe, "r2")
  r4 <- fe_model_fit(y_fe, X_fe, 1:4, df_fe, "r2")
  expect_true(r4 >= r0)
})

test_that("fe_coalition_values indexes by bitmask with v(empty) = 0", {
  set.seed(5)
  N <- 10; T <- 5
  unit <- rep(1:N, times = T); time <- rep(1:T, each = N)
  y <- rnorm(N * T)
  X <- matrix(rnorm(N * T * 3), ncol = 3)
  y_fe <- fe_partial_vec(y, unit, time)
  X_fe <- fe_partial_mat(X, unit, time)
  df_fe <- N * T - (N + T - 1)

  vals <- fe_coalition_values(y_fe, X_fe, df_fe)
  expect_length(vals, 2^3)
  expect_equal(vals[1], 0, tolerance = 1e-12)
  # bitmask indexing reproduces the direct characteristic-function evaluation
  v_empty <- fe_model_fit(y_fe, X_fe, integer(0), df_fe, "adj_r2")
  for (m in 1:7) {
    S <- which(fe_intToBin(m, 3) == 1)
    v_direct <- fe_model_fit(y_fe, X_fe, S, df_fe, "adj_r2") - v_empty
    expect_equal(vals[m + 1], v_direct, tolerance = 1e-12)
  }
})
