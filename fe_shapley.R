#' Exact Shapley Values
#'
#' Computes the exact Shapley value for every covariate by enumerating all
#' coalitions. The Shapley value of player \eqn{j} is
#' \deqn{\phi_j = \sum_{S \subseteq N \setminus \{j\}}
#'   \frac{|S|!\,(p - |S| - 1)!}{p!}
#'   \big[v(S \cup \{j\}) - v(S)\big].}
#'
#' @param y_fe residualised response.
#' @param X_fe residualised covariate matrix.
#' @param df_fe residual degrees of freedom of the FE-only model.
#' @param method characteristic-function type (see \code{\link{fe_model_fit}}).
#' @param verbose logical; print each covariate's value as it is computed.
#' @param vals optional precomputed vector of coalition values (output of
#'   \code{\link{fe_coalition_values}}); recomputed if \code{NULL}.
#'
#' @return Named numeric vector of Shapley values (length \code{p}). The
#'   efficiency identity \eqn{\sum_j \phi_j = v(N)} holds.
#' @export
#'
#' @examples
#' unit <- rep(1:10, each = 5); time <- rep(1:5, times = 10)
#' y <- rnorm(50); X <- matrix(rnorm(50 * 3), ncol = 3)
#' colnames(X) <- c("x1", "x2", "x3")
#' y_fe <- fe_partial_vec(y, unit, time)
#' X_fe <- fe_partial_mat(X, unit, time)
#' df_fe <- 50 - (10 + 5 - 1)
#' fe_shapley_exact(y_fe, X_fe, df_fe)
fe_shapley_exact <- function(y_fe, X_fe, df_fe, method = "adj_r2",
                             verbose = FALSE, vals = NULL) {
  p       <- ncol(X_fe)
  names_p <- colnames(X_fe)
  if (is.null(vals))
    vals <- fe_coalition_values(y_fe, X_fe, df_fe, method, verbose)

  # log combinatorial weights w(s) = s!(p-s-1)!/p!
  logw <- function(s) lfactorial(s) + lfactorial(p - s - 1) - lfactorial(p)

  phi <- numeric(p)
  powers <- 2^(0:(p - 1))
  for (i in 1:p) {
    bit_i  <- powers[i]
    others <- setdiff(0:(p - 1), i - 1)
    acc <- 0
    # iterate over all subsets S of N\{i} via bitmask over 'others'
    n_oth <- length(others)
    for (mm in 0:(2^n_oth - 1)) {
      # build coalition bitmask without i
      sel_bits <- others[which(bitwAnd(mm, 2^(0:(n_oth - 1))) > 0)]
      base_mask <- if (length(sel_bits)) sum(2^sel_bits) else 0
      s <- length(sel_bits)
      v_without <- vals[base_mask + 1]
      v_with    <- vals[bitwOr(base_mask, bit_i) + 1]
      acc <- acc + exp(logw(s)) * (v_with - v_without)
    }
    phi[i] <- acc
    if (verbose) cat("  phi_", i, " (", names_p[i], ") = ",
                     round(phi[i], 4), "\n", sep = "")
  }
  names(phi) <- names_p
  phi
}

#' Monte Carlo Shapley Values (Permutation Sampling)
#'
#' Approximates Shapley values by averaging marginal contributions over random
#' permutations of the players, a standard unbiased estimator that avoids
#' enumerating all \eqn{2^p} coalitions and is intended for larger \eqn{p}.
#' Coalition values are cached across permutations.
#'
#' @param y_fe residualised response.
#' @param X_fe residualised covariate matrix.
#' @param df_fe residual degrees of freedom of the FE-only model.
#' @param method characteristic-function type (see \code{\link{fe_model_fit}}).
#' @param n_samples number of permutations sampled.
#' @param verbose logical; print progress.
#'
#' @return Named numeric vector of estimated Shapley values.
#' @export
#'
#' @examples
#' unit <- rep(1:10, each = 5); time <- rep(1:5, times = 10)
#' y <- rnorm(50); X <- matrix(rnorm(50 * 6), ncol = 6)
#' colnames(X) <- paste0("x", 1:6)
#' y_fe <- fe_partial_vec(y, unit, time)
#' X_fe <- fe_partial_mat(X, unit, time)
#' df_fe <- 50 - (10 + 5 - 1)
#' set.seed(1)
#' fe_shapley_mc(y_fe, X_fe, df_fe, n_samples = 100)
fe_shapley_mc <- function(y_fe, X_fe, df_fe, method = "adj_r2",
                          n_samples = 2000, verbose = FALSE) {
  p       <- ncol(X_fe)
  names_p <- colnames(X_fe)
  v_empty <- fe_model_fit(y_fe, X_fe, integer(0), df_fe, method)

  vfun <- local({
    cache <- new.env(hash = TRUE, parent = emptyenv())
    function(S) {
      key <- paste(sort(S), collapse = "-")
      hit <- cache[[key]]
      if (!is.null(hit)) return(hit)
      val <- fe_model_fit(y_fe, X_fe, S, df_fe, method) - v_empty
      cache[[key]] <- val
      val
    }
  })

  phi <- numeric(p)
  for (b in 1:n_samples) {
    perm <- sample.int(p)
    prev <- 0
    S    <- integer(0)
    for (i in perm) {
      S_new <- c(S, i)
      v_new <- vfun(S_new)
      phi[i] <- phi[i] + (v_new - prev)
      prev <- v_new
      S <- S_new
    }
    if (verbose && b %% 500 == 0) cat("  MC sample", b, "of", n_samples, "\n")
  }
  phi <- phi / n_samples
  names(phi) <- names_p
  phi
}
