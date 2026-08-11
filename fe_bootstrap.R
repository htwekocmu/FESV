#' Residual Bootstrap for Shapley Values
#'
#' Produces bootstrap confidence intervals and standard errors for the Shapley
#' values using the residual bootstrap under the full FE model. Each bootstrap
#' replication resamples the residuals of the full FE fit, reconstructs the
#' response, and recomputes the exact Shapley values on the bootstrap sample.
#'
#' @param y numeric response vector of length \eqn{NT}.
#' @param X numeric matrix of candidate covariates (\eqn{NT \times p}).
#' @param unit unit (individual) index, length \eqn{NT}.
#' @param time time index, length \eqn{NT}.
#' @param B number of bootstrap replications (default \code{999}).
#' @param level confidence level (default \code{0.95}).
#' @param method characteristic-function type (see \code{\link{fe_model_fit}}).
#' @param gamma penalty per included covariate for Rule~B (passed through).
#' @param verbose logical; print progress.
#'
#' @return List with components \code{ci} (2 by \code{p} matrix of percentile
#'   confidence bounds), \code{se} (bootstrap standard errors),
#'   \code{phi_boot} (\eqn{B \times p} matrix of bootstrap Shapley values),
#'   and \code{phi_mean} (bootstrap mean).
#' @export
#'
#' @examples
#' set.seed(1)
#' N <- 10; T <- 6
#' unit <- rep(1:N, times = T); time <- rep(1:T, each = N)
#' X <- matrix(rnorm(N * T * 4), ncol = 4)
#' colnames(X) <- paste0("x", 1:4)
#' y <- as.vector(X %*% c(0.8, 0, 0.6, 0)) + rnorm(N * T) + unit + time
#' boot <- bootstrap_shapley(y, X, unit, time, B = 50, verbose = FALSE)
#' boot$ci
#' boot$se
bootstrap_shapley <- function(y, X, unit, time, B = 999,
                              level = 0.95, method = "adj_r2",
                              gamma = 0.01, verbose = TRUE) {
  N <- length(unique(unit))
  T <- length(unique(time))
  nT <- length(y)
  p <- ncol(X)

  # Estimate full model to get residuals
  y_fe <- fe_partial_vec(y, unit, time)
  X_fe <- fe_partial_mat(X, unit, time)
  fit_full <- stats::lm.fit(X_fe, y_fe)
  resid <- fit_full$residuals
  beta_hat <- fit_full$coefficients
  alpha_hat <- tapply(y - X %*% beta_hat, unit, mean)
  gamma_hat <- tapply(y - X %*% beta_hat, time, mean)
  gm <- mean(y - X %*% beta_hat)

  phi_boot <- matrix(NA, B, p)
  colnames(phi_boot) <- colnames(X)

  if (verbose) cat("Bootstrap replications:\n")

  for (b in 1:B) {
    # Residual bootstrap
    eps_b <- sample(resid, size = nT, replace = TRUE)
    y_b <- as.vector(X %*% beta_hat) +
      alpha_hat[as.character(unit)] +
      gamma_hat[as.character(time)] +
      eps_b - gm

    # Compute Shapley values on bootstrap sample
    res_b <- fe_shapley_covariate_selection(
      y_b, X, unit, time,
      method = method, exact = TRUE,
      gamma = gamma, verbose = FALSE)
    phi_boot[b, ] <- res_b$shapley_values

    if (verbose && b %% 100 == 0) cat("  ", b, "/", B, "\n")
  }

  alpha <- 1 - level
  ci <- apply(phi_boot, 2, function(x)
    stats::quantile(x, probs = c(alpha/2, 1 - alpha/2), na.rm = TRUE))

  # Bootstrap standard errors
  se <- apply(phi_boot, 2, stats::sd, na.rm = TRUE)

  list(ci = ci, se = se, phi_boot = phi_boot,
       phi_mean = colMeans(phi_boot, na.rm = TRUE))
}
