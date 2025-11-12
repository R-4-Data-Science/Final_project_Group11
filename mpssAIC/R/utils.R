#' @keywords internal
.encode_model <- function(vars) {
  if (length(vars) == 0) return("~1")
  paste(sort(vars), collapse = "+")
}

#' @keywords internal
.jaccard <- function(a, b) {
  a <- unique(a); b <- unique(b)
  if (length(a) == 0 && length(b) == 0) return(1)
  inter <- length(intersect(a, b))
  uni <- length(union(a, b))
  if (uni == 0) return(0)
  inter / uni
}

#' @keywords internal
.fit_and_aic <- function(X, y, family, terms) {
  df <- data.frame(y = y, X, check.names = FALSE)
  if (length(terms) == 0) {
    f <- as.formula("y ~ 1")
  } else {
    f <- as.formula(paste0("y ~ ", paste(terms, collapse = " + ")))
  }
  if (family == "gaussian") {
    fit <- lm(f, data = df)
  } else if (family == "binomial") {
    fit <- glm(f, data = df, family = binomial())
  } else {
    stop("Unsupported family: ", family)
  }
  aic <- AIC(fit)
  list(fit = fit, AIC = as.numeric(aic))
}

#' @keywords internal
.parents_aic <- function(X, y, family, parent_terms) {
  .fit_and_aic(X, y, family, parent_terms)$AIC
}

#' @keywords internal
.dedup_by_model <- function(df) {
  if (nrow(df) == 0) return(df)
  df <- df[order(df$AIC), , drop = FALSE]
  keep <- !duplicated(df$model_code)
  df[keep, , drop = FALSE]
}

#' @keywords internal
.vars_in_code <- function(code) {
  if (identical(code, "~1")) return(character(0))
  strsplit(code, "\\+")[[1]]
}
