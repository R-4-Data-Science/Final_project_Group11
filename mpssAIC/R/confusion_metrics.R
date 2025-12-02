#' Confusion matrix metrics for binary classification
#'
#' @param y_true True 0/1 outcomes
#' @param p_hat Predicted probabilities
#' @param cutoff Probability cutoff (default 0.5)
#'
#' @return Data frame with accuracy, sensitivity, specificity, FDR, DOR
#' @export
confusion_metrics <- function(y_true, p_hat, cutoff = 0.5) {
  y_pred <- as.integer(p_hat >= cutoff)
  tab <- table(Predicted = y_pred, Actual = y_true)

  TP <- tab["1", "1"]
  TN <- tab["0", "0"]
  FP <- tab["1", "0"]
  FN <- tab["0", "1"]

  acc  <- (TP + TN) / sum(tab)
  sens <- TP / (TP + FN)
  spec <- TN / (TN + FP)
  fdr  <- FP / (FP + TP)
  dor  <- (TP / FN) / (FP / TN)

  data.frame(
    accuracy    = acc,
    sensitivity = sens,
    specificity = spec,
    FDR         = fdr,
    DOR         = dor
  )
}
