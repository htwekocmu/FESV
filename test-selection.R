test_that("fe_nested_selection returns a valid ranking and optimum", {
  set.seed(10)
  N <- 15; T <- 8
  unit <- rep(1:N, times = T); time <- rep(1:T, each = N)
  y <- rnorm(N * T)
  X <- matrix(rnorm(N * T * 3), ncol = 3)
  colnames(X) <- c("a", "b", "c")
  y_fe <- fe_partial_vec(y, unit, time)
  X_fe <- fe_partial_mat(X, unit, time)
  df_fe <- N * T - (N + T - 1)
  phi <- fe_shapley_exact(y_fe, X_fe, df_fe)

  ns <- fe_nested_selection(phi, y_fe, X_fe, df_fe, gamma = 0.01,
                            verbose = FALSE)
  expect_true(ns$k_star >= 1 && ns$k_star <= 3)
  expect_equal(length(ns$selected), ns$k_star)
  expect_equal(ns$ranking, order(phi, decreasing = TRUE))
  expect_length(ns$J, 3)
  # selected coalition is the top-k* by Shapley value
  expect_setequal(ns$selected, ns$ranking[1:ns$k_star])
})

test_that("selection metrics are correct on hand-computed cases", {
  m1 <- selection_metrics(selected = c(1, 3, 5), active = c(1, 3, 5), p = 6)
  expect_equal(m1[["TP"]], 3)
  expect_equal(m1[["FP"]], 0)
  expect_equal(m1[["FN"]], 0)
  expect_equal(m1[["precision"]], 1)
  expect_equal(m1[["recall"]], 1)
  expect_equal(m1[["F1"]], 1)

  m2 <- selection_metrics(selected = c(1, 2, 3), active = c(1, 3, 5), p = 6)
  expect_equal(m2[["TP"]], 2)
  expect_equal(m2[["FP"]], 1)
  expect_equal(m2[["FN"]], 1)
  expect_equal(m2[["precision"]], 2/3)
  expect_equal(m2[["recall"]], 2/3)
  expect_equal(m2[["F1"]], 2/3)

  m3 <- selection_metrics(selected = integer(0), active = c(1, 2), p = 5)
  expect_equal(m3[["recall"]], 0)
  expect_equal(m3[["F1"]], 0)
})
