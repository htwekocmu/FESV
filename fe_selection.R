#' Nested Model Selection
#'
#' Orders covariates by their Shapley values and selects the model size that
#' maximises the penalised nested criterion
#' \deqn{J(k) = v(M_k) - \gamma k, \qquad k = 1, \dots, p,}
#' where \eqn{M_k} is the coalition of the top-\eqn{k} covariates and
#' \eqn{\gamma > 0} is a penalty per included covariate. This is Rule~B of the
#' FESV framework; it is a data-driven analogue of information-criterion
#' selection.
#'
#' @param phi numeric vector of Shapley values.
#' @param y_fe residualised response.
#' @param X_fe residualised covariate matrix.
#' @param df_fe residual degrees of freedom of the FE-only model.
#' @param method characteristic-function type (see \code{\link{fe_model_fit}}).
#' @param gamma penalty per included covariate (default \code{0.01}).
#' @param verbose logical; print the criterion path.
#'
#' @return List with components \code{k_star} (optimal size), \code{selected}
#'   (indices of the selected coalition), \code{selected_names},
#'   \code{J} (criterion path), and \code{ranking}.
#' @export
#'
#' @examples
#' unit <- rep(1:10, each = 5); time <- rep(1:5, times = 10)
#' y <- rnorm(50); X <- matrix(rnorm(50 * 3), ncol = 3)
#' colnames(X) <- c("x1", "x2", "x3")
#' y_fe <- fe_partial_vec(y, unit, time)
#' X_fe <- fe_partial_mat(X, unit, time)
#' df_fe <- 50 - (10 + 5 - 1)
#' phi <- fe_shapley_exact(y_fe, X_fe, df_fe)
#' fe_nested_selection(phi, y_fe, X_fe, df_fe, verbose = FALSE)
fe_nested_selection <- function(phi, y_fe, X_fe, df_fe,
                                method = "adj_r2", gamma = 0.01,
                                verbose = TRUE) {
  p       <- length(phi)
  ranking <- order(phi, decreasing = TRUE)
  names_p <- names(phi)
  v_empty <- fe_model_fit(y_fe, X_fe, integer(0), df_fe, method)

  J <- numeric(p)
  for (k in 1:p) {
    Mk  <- ranking[1:k]
    vMk <- fe_model_fit(y_fe, X_fe, Mk, df_fe, method) - v_empty
    J[k] <- vMk - gamma * k
    if (verbose) cat(sprintf("  k=%2d: J(k) = %9.4f   [%s]\n",
                             k, J[k], paste(names_p[Mk], collapse = ",")))
  }
  k_star   <- which.max(J)
  selected <- ranking[1:k_star]
  if (verbose) {
    cat("\nOptimal size k* =", k_star, "\n")
    cat("Selected:", paste(names_p[selected], collapse = ", "), "\n")
  }
  list(k_star = k_star, selected = selected,
       selected_names = names_p[selected], J = J, ranking = ranking)
}

#' Fixed Effects Covariate Selection Using the Shapley Value
#'
#' Master routine of the FESV framework. Unit and time effects are held fixed
#' (partialled out) and the candidate covariates play a cooperative game whose
#' characteristic function is the improvement in a fit statistic over the
#' FE-only baseline. Returns Shapley values, game-theoretic diagnostics, and
#' the two selection rules.
#'
#' @param y numeric response vector of length \eqn{NT}.
#' @param X numeric matrix of candidate covariates (\eqn{NT \times p}), with
#'   named columns.
#' @param unit unit (individual) index, length \eqn{NT}.
#' @param time time index, length \eqn{NT}.
#' @param method characteristic-function type (see \code{\link{fe_model_fit}}).
#' @param exact logical; if \code{TRUE} compute exact Shapley values (feasible
#'   for \eqn{p \le 20}), otherwise Monte Carlo sampling.
#' @param n_samples number of permutations when \code{exact = FALSE}.
#' @param gamma penalty per included covariate for Rule~B.
#' @param verbose logical; print progress and diagnostics.
#'
#' @return A list with components
#'   \itemize{
#'     \item \code{shapley_values}: named vector of Shapley values;
#'     \item \code{selected_covariates} / \code{selected_sign}: Rule~A
#'           (\eqn{\phi_j > 0}) indices;
#'     \item \code{selected_nested}: Rule~B indices;
#'     \item \code{nested}: output of \code{\link{fe_nested_selection}};
#'     \item \code{coalition_values}: vector of coalition values (if exact);
#'     \item \code{superadditivity}, \code{core}: diagnostics (if exact);
#'     \item \code{v_grand}, \code{efficiency_gap}: efficiency identity checks;
#'     \item \code{df_fe}, \code{method}, \code{N}, \code{T}, \code{p}.
#'   }
#' @export
#'
#' @examples
#' set.seed(1)
#' N <- 20; T <- 10
#' unit <- rep(1:N, times = T); time <- rep(1:T, each = N)
#' beta <- c(0.8, 0, 0.6, 0, -0.5, 0)
#' X <- matrix(rnorm(N * T * length(beta)), ncol = length(beta))
#' colnames(X) <- paste0("x", 1:length(beta))
#' y <- as.vector(X %*% beta) + rnorm(N * T) + unit * 0.5 + time
#' res <- fe_shapley_covariate_selection(y, X, unit, time,
#'                                       exact = TRUE, gamma = 0.01,
#'                                       verbose = FALSE)
#' res$shapley_values
#' res$selected_nested
fe_shapley_covariate_selection <- function(y, X, unit, time,
                                           method = "adj_r2",
                                           exact = TRUE,
                                           n_samples = 2000,
                                           gamma = 0.01,
                                           verbose = TRUE) {
  stopifnot(is.matrix(X))
  N <- length(unique(unit)); T <- length(unique(time))
  nT <- length(y); p <- ncol(X)
  if (is.null(colnames(X))) colnames(X) <- paste0("V", 1:p)

  df_fe <- nT - (N + T - 1)   # residual df of the FE-only model (excl. slopes)

  if (verbose) {
    cat("======================================================================\n")
    cat("  Fixed Effects Covariate Selection using the Shapley Value\n")
    cat("======================================================================\n")
    cat("  N =", N, " units, T =", T, " periods, NT =", nT, "\n")
    cat("  Players (p) =", p, ":", paste(colnames(X), collapse = ", "), "\n")
    cat("  Fixed (not players): unit effects alpha_i, time effects gamma_t\n")
    cat("  Method =", method, " | FE residual df =", df_fe, "\n")
    cat("  Computation:", ifelse(exact, "Exact", "Monte Carlo"), "\n\n")
  }

  # Hold unit and time effects fixed (partial them out once)
  y_fe <- fe_partial_vec(y, unit, time)
  X_fe <- fe_partial_mat(X, unit, time)

  # Coalition values (exact enumeration reused by Shapley + diagnostics)
  vals <- NULL
  if (exact) {
    vals <- fe_coalition_values(y_fe, X_fe, df_fe, method, verbose = FALSE)
    # Sanity: v(empty) = 0 by construction
    stopifnot(abs(vals[1]) < 1e-12)
    phi <- fe_shapley_exact(y_fe, X_fe, df_fe, method, verbose, vals = vals)
  } else {
    phi <- fe_shapley_mc(y_fe, X_fe, df_fe, method, n_samples, verbose)
  }

  # --- Cooperative game diagnostics -------------------------------------
  v_grand <- fe_characteristic_function(1:p, y_fe, X_fe, df_fe, method)
  superadd <- NULL; core <- NULL
  if (verbose) {
    cat("\nGame-theoretic verification:\n")
    cat("----------------------------\n")
    cat("v(empty) = 0 (by normalisation)\n")
    cat("Efficiency check:\n")
    cat("  v(N)                  =", round(v_grand, 6), "\n")
    cat("  Sum of Shapley values =", round(sum(phi), 6), "\n")
    cat("  Difference            =",
        format(abs(v_grand - sum(phi)), digits = 3), "\n")
  }
  if (exact) {
    superadd <- fe_check_superadditivity(vals, p, verbose = verbose)
    core     <- fe_check_core_membership(phi, vals, p, verbose = verbose)
  }

  # Rule A: positive Shapley value
  sel_sign <- which(phi > 0)

  # Rule B: nested penalised selection
  if (verbose) cat("\nNested model selection (J(k) = v(M_k) - gamma*k):\n")
  nested <- fe_nested_selection(phi, y_fe, X_fe, df_fe, method, gamma, verbose)

  list(shapley_values = phi,
       selected_covariates = sel_sign,
       selected_sign  = sel_sign,
       selected_nested = nested$selected,
       nested = nested,
       coalition_values = vals,
       superadditivity = superadd,
       core = core,
       v_grand = v_grand,
       efficiency_gap = abs(v_grand - sum(phi)),
       df_fe = df_fe,
       method = method,
       N = N, T = T, p = p)
}
