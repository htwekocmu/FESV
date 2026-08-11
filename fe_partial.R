#' Partial Unit and Time Effects Out of a Vector
#'
#' Residualise a numeric vector on unit and time dummies, numerically
#' equivalent to estimating
#' \deqn{v_{it} = \alpha_i + \gamma_t + u_{it}}
#' and keeping the residuals (plus the grand mean, which cancels in every fit
#' statistic). The two-way fixed effects are held fixed and never selected.
#'
#' @param v numeric vector of length \eqn{NT}.
#' @param unit unit (individual) index, length \eqn{NT}.
#' @param time time index, length \eqn{NT}.
#'
#' @return Numeric vector of residuals from the two-way FE regression.
#' @export
#'
#' @examples
#' unit <- rep(1:10, each = 5)
#' time <- rep(1:5, times = 10)
#' v <- rnorm(50) + unit * 2 + time
#' r <- fe_partial_vec(v, unit, time)
#' head(r)
fe_partial_vec <- function(v, unit, time) {
  gm <- mean(v)
  um <- tapply(v, unit, mean)
  tm <- tapply(v, time, mean)
  v - um[as.character(unit)] - tm[as.character(time)] + gm
}

#' Partial Unit and Time Effects Out of Every Column of a Matrix
#'
#' Applies \code{\link{fe_partial_vec}} to each column of a matrix. The
#' resulting columns are residuals from the two-way FE regression of each
#' column on the unit and time effects.
#'
#' @param M numeric matrix (each column a covariate), \eqn{NT \times p}.
#' @param unit unit (individual) index, length \eqn{NT}.
#' @param time time index, length \eqn{NT}.
#'
#' @return Matrix of the same dimension as \code{M}.
#' @export
#'
#' @examples
#' unit <- rep(1:10, each = 5)
#' time <- rep(1:5, times = 10)
#' M <- matrix(rnorm(100), ncol = 2)
#' colnames(M) <- c("x1", "x2")
#' Xr <- fe_partial_mat(M, unit, time)
#' apply(Xr, 2, function(x) tapply(x, unit, sum))
fe_partial_mat <- function(M, unit, time) {
  out <- apply(M, 2, fe_partial_vec, unit = unit, time = time)
  matrix(out, nrow = nrow(M), dimnames = list(NULL, colnames(M)))
}
