#' Compute AIC for a model with selected predictors
#'
#' @param X A data frame of predictors (n × p)
#' @param y Response vector of length n
#' @param vars Integer vector of column indices to include (empty = intercept-only)
#' @param model_type Either "gaussian" or "binomial"
#'
#' @return Numeric AIC value
#' @keywords internal
fit_model_aic <- function(X, y, vars, model_type = c("gaussian", "binomial")) {
  model_type <- match.arg(model_type)

  # Empty model (intercept only)
  if (length(vars) == 0) {
    df  <- data.frame(y = y)
    fam <- if (model_type == "gaussian") gaussian() else binomial()
    fit <- glm(y ~ 1, data = df, family = fam)
    return(AIC(fit))
  }

  # Data frame with selected predictors
  df  <- data.frame(y = y, X[, vars, drop = FALSE])
  fam <- if (model_type == "gaussian") gaussian() else binomial()
  fit <- glm(y ~ ., data = df, family = fam)

  AIC(fit)
}
