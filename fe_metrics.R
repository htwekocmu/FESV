#' Selection Metrics
#'
#' Computes classification metrics for a covariate-selection exercise against a
#' known active set:
#' \deqn{\text{Recall} = \frac{TP}{TP+FN},\qquad
#'       \text{Precision} = \frac{TP}{TP+FP},\qquad
#'       F_1 = \frac{2\,\text{Precision}\cdot\text{Recall}}
#'                  {\text{Precision}+\text{Recall}}.}
#'
#' @param selected integer vector of indices selected by the procedure.
#' @param active integer vector of indices in the true active set.
#' @param p total number of candidate covariates.
#'
#' @return Named numeric vector with \code{TP}, \code{FP}, \code{FN},
#'   \code{TN}, \code{precision}, \code{recall}, and \code{F1}.
#' @export
#'
#' @examples
#' selection_metrics(selected = c(1, 3, 5), active = c(1, 3, 5), p = 6)
#' selection_metrics(selected = c(1, 2, 3), active = c(1, 3, 5), p = 6)
selection_metrics <- function(selected, active, p) {
  sel <- rep(FALSE, p); sel[selected] <- TRUE
  act <- rep(FALSE, p); act[active]   <- TRUE
  TP <- sum(sel & act)
  FP <- sum(sel & !act)
  FN <- sum(!sel & act)
  TN <- sum(!sel & !act)
  precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  recall    <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  f1 <- if (!is.na(precision) && !is.na(recall) && (precision + recall) > 0)
    2 * precision * recall / (precision + recall) else 0
  c(TP = TP, FP = FP, FN = FN, TN = TN,
    precision = precision, recall = recall, F1 = f1)
}
