#' Goodness-of-Fit for the FE Model Using Covariates in a Coalition
#'
#' Computes a fit statistic for the two-way FE panel model restricted to the
#' covariates in coalition \code{S}. The unit and time effects are held fixed:
#' they are partialled out of \code{y_fe} and \code{X_fe} before calling this
#' function.
#'
#' @param y_fe residualised response (output of \code{\link{fe_partial_vec}}).
#' @param X_fe residualised covariate matrix (output of
#'   \code{\link{fe_partial_mat}}).
#' @param S integer vector of column indices forming the coalition (may be
#'   empty).
#' @param df_fe residual degrees of freedom of the FE-only model,
#'   \eqn{NT - (N + T - 1)}.
#' @param method characteristic-function type:
#'   \itemize{
#'     \item \code{"r2"}: plain \eqn{R^2} net of fixed effects;
#'     \item \code{"adj_r2"} (default): adjusted \eqn{R^2}, which penalises
#'           degrees of freedom so that adding an irrelevant covariate lowers
#'           the value and noise players receive negative Shapley values;
#'     \item \code{"neg_bic"}: negative Bayesian information criterion,
#'           relative to the FE-only model.
#'   }
#'
#' @return Numeric fit value (\code{-Inf} if the coalition is too large for
#'   the remaining degrees of freedom).
#' @export
#'
#' @examples
#' unit <- rep(1:10, each = 5); time <- rep(1:5, times = 10)
#' y <- rnorm(50); X <- matrix(rnorm(100), ncol = 2)
#' y_fe <- fe_partial_vec(y, unit, time)
#' X_fe <- fe_partial_mat(X, unit, time)
#' df_fe <- 50 - (10 + 5 - 1)
#' fe_model_fit(y_fe, X_fe, integer(0), df_fe)
#' fe_model_fit(y_fe, X_fe, 1:2, df_fe)
fe_model_fit <- function(y_fe, X_fe, S, df_fe, method = "adj_r2") {
  nT  <- length(y_fe)
  ns  <- length(S)
  tss <- sum(y_fe^2)

  if (ns == 0L) {
    ssr <- tss
  } else {
    Xs  <- X_fe[, S, drop = FALSE]
    fit <- tryCatch(stats::lm.fit(Xs, y_fe), error = function(e) NULL)
    if (is.null(fit)) return(-Inf)
    ssr <- sum(fit$residuals^2)
  }

  if (method == "r2") {
    return(1 - ssr / tss)
  }
  if (method == "adj_r2") {
    df_resid <- df_fe - ns
    if (df_resid <= 0) return(-Inf)
    return(1 - (ssr / df_resid) / (tss / df_fe))
  }
  if (method == "neg_bic") {
    sigma2 <- ssr / nT
    if (sigma2 <= 0) sigma2 <- .Machine$double.eps
    ll  <- -0.5 * nT * (log(2 * pi) + log(sigma2) + 1)
    bic <- -2 * ll + ns * log(nT)
    return(-bic)
  }
  stop("Unknown method: ", method)
}

#' Characteristic Function Value for a Coalition
#'
#' Evaluates \eqn{v(S) = \text{fit}(S) - \text{fit}(\emptyset)}, so that
#' \eqn{v(\emptyset) = 0} and the grand coalition value \eqn{v(N)} equals the
#' total fit gain over the FE-only baseline.
#'
#' @param S integer vector of column indices forming the coalition.
#' @param y_fe residualised response.
#' @param X_fe residualised covariate matrix.
#' @param df_fe residual degrees of freedom of the FE-only model.
#' @param method characteristic-function type (see \code{\link{fe_model_fit}}).
#' @param v_empty precomputed value of the empty coalition (optional; computed
#'   if \code{NULL}).
#'
#' @return Numeric characteristic-function value.
#' @export
#'
#' @examples
#' unit <- rep(1:10, each = 5); time <- rep(1:5, times = 10)
#' y <- rnorm(50); X <- matrix(rnorm(100), ncol = 2)
#' y_fe <- fe_partial_vec(y, unit, time)
#' X_fe <- fe_partial_mat(X, unit, time)
#' df_fe <- 50 - (10 + 5 - 1)
#' fe_characteristic_function(1:2, y_fe, X_fe, df_fe)
fe_characteristic_function <- function(S, y_fe, X_fe, df_fe,
                                       method = "adj_r2", v_empty = NULL) {
  if (is.null(v_empty)) {
    v_empty <- fe_model_fit(y_fe, X_fe, integer(0), df_fe, method)
  }
  fe_model_fit(y_fe, X_fe, S, df_fe, method) - v_empty
}
