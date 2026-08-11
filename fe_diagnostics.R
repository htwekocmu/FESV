#' Superadditivity Check
#'
#' Tests whether the cooperative game is superadditive, i.e.
#' \eqn{v(A \cup B) \ge v(A) + v(B)} for all disjoint \eqn{A, B}. Under a
#' degrees-of-freedom-penalised characteristic function (adjusted \eqn{R^2},
#' negative BIC) violations are \emph{expected} for coalitions containing noise
#' players; this mirrors the Occam-penalty violations that arise in Bayesian
#' evidence games.
#'
#' @param vals numeric vector of coalition values indexed by bitmask.
#' @param p number of players.
#' @param tol numerical tolerance for declaring a violation.
#' @param verbose logical; print a summary.
#'
#' @return List with components \code{superadditive} (logical),
#'   \code{n_pairs}, \code{n_violations}, and \code{worst_gap}.
#' @export
#'
#' @examples
#' unit <- rep(1:10, each = 5); time <- rep(1:5, times = 10)
#' y <- rnorm(50); X <- matrix(rnorm(50 * 3), ncol = 3)
#' y_fe <- fe_partial_vec(y, unit, time)
#' X_fe <- fe_partial_mat(X, unit, time)
#' df_fe <- 50 - (10 + 5 - 1)
#' vals <- fe_coalition_values(y_fe, X_fe, df_fe)
#' fe_check_superadditivity(vals, 3)
fe_check_superadditivity <- function(vals, p, tol = 1e-10, verbose = TRUE) {
  full <- 2^p - 1
  n_pairs <- 0L; n_viol <- 0L; worst_gap <- 0
  for (a in 1:(full - 1)) {
    rem <- bitwAnd(full, bitwNot(a))
    b <- rem
    while (b > 0) {
      if (b > a) {              # count each unordered disjoint pair once
        gap <- vals[bitwOr(a, b) + 1] - vals[a + 1] - vals[b + 1]
        n_pairs <- n_pairs + 1L
        if (gap < -tol) {
          n_viol <- n_viol + 1L
          if (gap < worst_gap) worst_gap <- gap
        }
      }
      b <- bitwAnd(b - 1L, rem)
    }
  }
  superadd <- (n_viol == 0L)
  if (verbose) {
    cat(sprintf("Superadditivity: %d of %d disjoint pairs violated (worst gap %.4f)\n",
                n_viol, n_pairs, worst_gap))
    if (!superadd)
      cat("  -> expected under the df-penalised characteristic function\n")
  }
  list(superadditive = superadd, n_pairs = n_pairs,
       n_violations = n_viol, worst_gap = worst_gap)
}

#' Core Membership Check
#'
#' Checks whether the Shapley vector lies in the core of the cooperative game,
#' i.e. whether \eqn{\sum_{i \in S} \phi_i \ge v(S)} for every coalition
#' \eqn{S} (with efficiency \eqn{\sum_i \phi_i = v(N)}). The core is empty in
#' general; the check is diagnostic.
#'
#' @param phi numeric vector of Shapley values.
#' @param vals numeric vector of coalition values indexed by bitmask.
#' @param p number of players.
#' @param tol numerical tolerance.
#' @param verbose logical; print a summary.
#'
#' @return List with components \code{in_core} (logical), \code{min_slack},
#'   and \code{worst_coalition} (indices of the coalition attaining the
#'   minimum slack).
#' @export
#'
#' @examples
#' unit <- rep(1:10, each = 5); time <- rep(1:5, times = 10)
#' y <- rnorm(50); X <- matrix(rnorm(50 * 3), ncol = 3)
#' y_fe <- fe_partial_vec(y, unit, time)
#' X_fe <- fe_partial_mat(X, unit, time)
#' df_fe <- 50 - (10 + 5 - 1)
#' vals <- fe_coalition_values(y_fe, X_fe, df_fe)
#' phi <- fe_shapley_exact(y_fe, X_fe, df_fe, vals = vals)
#' fe_check_core_membership(phi, vals, 3)
fe_check_core_membership <- function(phi, vals, p, tol = 1e-8, verbose = TRUE) {
  powers <- 2^(0:(p - 1))
  min_slack <- Inf; worst_S <- integer(0)
  for (m in 1:(2^p - 1)) {
    S <- which(bitwAnd(m, powers) > 0)
    slack <- sum(phi[S]) - vals[m + 1]
    if (slack < min_slack) { min_slack <- slack; worst_S <- S }
  }
  in_core <- (min_slack >= -tol)
  if (verbose) {
    cat(sprintf("Core membership: %s (min slack = %.4f at S = {%s})\n",
                ifelse(in_core, "phi IS in the core", "phi NOT in the core"),
                min_slack,
                paste(names(phi)[worst_S], collapse = ",")))
  }
  list(in_core = in_core, min_slack = min_slack, worst_coalition = worst_S)
}
