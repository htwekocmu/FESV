test_that("fe_partial_vec removes unit and time means", {
  set.seed(1)
  N <- 10; T <- 6
  unit <- rep(1:N, times = T)
  time <- rep(1:T, each = N)
  v <- rnorm(N * T) + unit * 3 + time * 2

  r <- fe_partial_vec(v, unit, time)

  expect_equal(length(r), N * T)
  # residuals sum to zero within every unit and every time period
  expect_equal(as.numeric(tapply(r, unit, sum)), rep(0, N), tolerance = 1e-12)
  expect_equal(as.numeric(tapply(r, time, sum)), rep(0, T), tolerance = 1e-12)
})

test_that("fe_partial_mat residualises every column", {
  set.seed(2)
  N <- 8; T <- 5
  unit <- rep(1:N, times = T)
  time <- rep(1:T, each = N)
  M <- matrix(rnorm(N * T * 3), ncol = 3)
  colnames(M) <- c("x1", "x2", "x3")

  Mr <- fe_partial_mat(M, unit, time)

  expect_equal(dim(Mr), dim(M))
  expect_equal(colnames(Mr), colnames(M))
  # equivalent to regressing each column on unit + time dummies
  D <- model.matrix(~ factor(unit) + factor(time) - 1)
  fit <- lm(M[, 1] ~ D - 1)
  expect_equal(as.numeric(Mr[, 1]),
               as.numeric(fit$residuals), tolerance = 1e-10)
})
