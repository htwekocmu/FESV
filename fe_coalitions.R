#' Binary Expansion of a Coalition Bitmask
#'
#' Converts the integer bitmask \code{m} over \code{p} players into a binary
#' indicator vector. Bit \eqn{(i-1)} set corresponds to player \eqn{i} being in
#' the coalition.
#'
#' @param m non-negative integer coalition bitmask.
#' @param p number of players (covariates).
#'
#' @return Integer vector of length \code{p} with 0/1 entries.
#' @keywords internal
fe_intToBin <- function(m, p) {
  as.integer(intToBits(m))[1:p]
}

#' Precompute Coalition Values for All Subsets
#'
#' Enumerates all \eqn{2^p} coalitions and evaluates the characteristic
#' function \eqn{v(S)} for each. Values are indexed by the coalition bitmask
#' \eqn{m = 0, \dots, 2^p - 1}; bit \eqn{(i-1)} set \eqn{\Leftrightarrow}
#' player \eqn{i \in S}. Element \code{[1]} is \eqn{v(\emptyset) = 0} by
#' construction.
#'
#' @param y_fe residualised response.
#' @param X_fe residualised covariate matrix.
#' @param df_fe residual degrees of freedom of the FE-only model.
#' @param method characteristic-function type (see \code{\link{fe_model_fit}}).
#' @param verbose logical; print progress.
#'
#' @return Numeric vector of length \eqn{2^p} of characteristic-function values
#'   indexed by coalition bitmask.
#' @export
#'
#' @examples
#' unit <- rep(1:10, each = 5); time <- rep(1:5, times = 10)
#' y <- rnorm(50); X <- matrix(rnorm(50 * 3), ncol = 3)
#' colnames(X) <- c("x1", "x2", "x3")
#' y_fe <- fe_partial_vec(y, unit, time)
#' X_fe <- fe_partial_mat(X, unit, time)
#' df_fe <- 50 - (10 + 5 - 1)
#' vals <- fe_coalition_values(y_fe, X_fe, df_fe)
#' vals[1]  # v(empty) == 0
fe_coalition_values <- function(y_fe, X_fe, df_fe,
                                method = "adj_r2", verbose = FALSE) {
  p <- ncol(X_fe)
  M <- 2^p
  v_empty <- fe_model_fit(y_fe, X_fe, integer(0), df_fe, method)
  vals <- numeric(M)
  for (m in 0:(M - 1)) {
    S <- which(fe_intToBin(m, p) == 1)
    vals[m + 1] <- fe_model_fit(y_fe, X_fe, S, df_fe, method) - v_empty
    if (verbose && m %% 256 == 0) cat("  coalition", m, "of", M - 1, "\n")
  }
  vals
}

#' Print the Characteristic Function for All Coalitions
#'
#' Displays \eqn{v(S)} for every coalition of a small game, using either index
#' labels or covariate names.
#'
#' @param vals numeric vector of coalition values indexed by bitmask (output of
#'   \code{\link{fe_coalition_values}}).
#' @param p number of players.
#' @param names_p optional character vector of player (covariate) names.
#'
#' @return Invisibly, a data frame with one row per coalition.
#' @export
#'
#' @examples
#' unit <- rep(1:10, each = 5); time <- rep(1:5, times = 10)
#' y <- rnorm(50); X <- matrix(rnorm(50 * 3), ncol = 3)
#' colnames(X) <- c("x1", "x2", "x3")
#' y_fe <- fe_partial_vec(y, unit, time)
#' X_fe <- fe_partial_mat(X, unit, time)
#' df_fe <- 50 - (10 + 5 - 1)
#' vals <- fe_coalition_values(y_fe, X_fe, df_fe)
#' fe_print_characteristic_function(vals, 3, colnames(X))
fe_print_characteristic_function <- function(vals, p, names_p = NULL) {
  cat("\nCharacteristic function v(S) for all coalitions:\n")
  cat("------------------------------------------------\n")
  out <- data.frame(coalition = character(2^p), v = vals)
  for (m in 0:(2^p - 1)) {
    S <- which(fe_intToBin(m, p) == 1)
    lab <- if (length(S) == 0) "{}" else {
      inner <- if (is.null(names_p)) paste(S, collapse = ",")
               else paste(names_p[S], collapse = ",")
      paste0("{", inner, "}")
    }
    out$coalition[m + 1] <- lab
    cat(sprintf("v(%s) = %.4f\n", lab, vals[m + 1]))
  }
  invisible(out)
}
