#' Fixed-Effects OLS on a Chosen Coalition
#'
#' Estimates the two-way FE panel model restricted to a selected coalition of
#' covariates, with unit and time effects partialled out and classical
#' (non-robust) standard errors whose degrees of freedom correct for the fixed
#' effects.
#'
#' @param y numeric response vector of length \eqn{NT}.
#' @param X numeric matrix of candidate covariates (\eqn{NT \times p}).
#' @param unit unit (individual) index, length \eqn{NT}.
#' @param time time index, length \eqn{NT}.
#' @param sel integer vector of column indices forming the coalition.
#'
#' @return List with components \code{coefficients} (data frame of estimates,
#'   standard errors, \eqn{t}-values and \eqn{p}-values), \code{r2},
#'   \code{adj_r2}, \code{sigma}, \code{df_residual}, and \code{ssr}.
#' @export
#'
#' @examples
#' set.seed(1)
#' N <- 15; T <- 8
#' unit <- rep(1:N, times = T); time <- rep(1:T, each = N)
#' X <- matrix(rnorm(N * T * 3), ncol = 3)
#' colnames(X) <- c("x1", "x2", "x3")
#' y <- 0.8 * X[, 1] - 0.5 * X[, 3] + rnorm(N * T) + unit + time
#' fe_ols(y, X, unit, time, sel = c(1, 3))
fe_ols <- function(y, X, unit, time, sel) {
  N <- length(unique(unit)); T <- length(unique(time)); nT <- length(y)
  y_fe <- fe_partial_vec(y, unit, time)
  X_fe <- fe_partial_mat(X, unit, time)
  Xs   <- X_fe[, sel, drop = FALSE]
  fit  <- stats::lm.fit(Xs, y_fe)
  resid <- fit$residuals
  ssr  <- sum(resid^2)
  k    <- length(sel)
  df   <- nT - (N + T - 1) - k
  sigma2 <- ssr / df
  XtX_inv <- solve(crossprod(Xs))
  se   <- sqrt(diag(XtX_inv) * sigma2)
  beta <- fit$coefficients
  tval <- beta / se
  pval <- 2 * stats::pt(-abs(tval), df = df)
  tss  <- sum(y_fe^2)
  r2   <- 1 - ssr / tss
  adj  <- 1 - (ssr / df) / (tss / (nT - (N + T - 1)))
  tab  <- data.frame(Estimate = beta, Std.Error = se,
                     t.value = tval, p.value = pval,
                     row.names = colnames(Xs))
  list(coefficients = tab, r2 = r2, adj_r2 = adj,
       sigma = sqrt(sigma2), df_residual = df, ssr = ssr)
}
