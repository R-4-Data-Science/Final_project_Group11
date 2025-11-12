#' Confusion metrics at cutoff 0.5 for logistic models
#'
#' @param y_true true binary labels (0/1 or FALSE/TRUE).
#' @param y_prob predicted probabilities in [0,1].
#' @param cutoff threshold, default 0.5.
#' @return A named list: prevalence, accuracy, sensitivity, specificity, FDR, DOR, TP, FP, TN, FN.
#' @export
confusion_metrics <- function(y_true, y_prob, cutoff = 0.5) {
  y_true <- as.integer(as.logical(y_true))
  y_hat <- as.integer(y_prob >= cutoff)
  TP <- sum(y_hat == 1 & y_true == 1)
  FP <- sum(y_hat == 1 & y_true == 0)
  TN <- sum(y_hat == 0 & y_true == 0)
  FN <- sum(y_hat == 0 & y_true == 1)
  prevalence <- mean(y_true == 1)
  accuracy <- (TP + TN) / length(y_true)
  sensitivity <- ifelse((TP + FN) > 0, TP / (TP + FN), NA_real_)
  specificity <- ifelse((TN + FP) > 0, TN / (TN + FP), NA_real_)
  FDR <- ifelse((TP + FP) > 0, FP / (TP + FP), NA_real_)
  DOR <- ifelse((TP * TN) > 0 && (FP * FN) > 0, (TP/FN) / (FP/TN), NA_real_)
  list(prevalence = prevalence,
       accuracy = accuracy,
       sensitivity = sensitivity,
       specificity = specificity,
       FDR = FDR,
       DOR = DOR,
       TP = TP, FP = FP, TN = TN, FN = FN)
}
