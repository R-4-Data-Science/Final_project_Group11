#' Plausible models: AIC tolerance + mean stability filter
#'
#' @param forest result of `build_paths()`.
#' @param pi numeric vector of average feature stabilities named by predictor.
#' @param Delta AIC tolerance relative to the minimum AIC in the forest.
#' @param tau minimum mean stability threshold for a model to be kept.
#' @param drop_jaccard logical. If TRUE, drop near-duplicate models by Jaccard similarity.
#' @param jaccard_threshold numeric in (0,1]; models with Jaccard >= threshold are treated as duplicates (keep best AIC).
#' @param attach_fits logical. If TRUE, refits selected models on full X,y and attaches fits.
#' @param X,y (optional) required only if attach_fits = TRUE.
#' @return data.frame with columns model_code, AIC, mean_pi, terms (list).
#' @export
plausible_models <- function(forest, pi, Delta = 2, tau = 0.4,
                             drop_jaccard = TRUE, jaccard_threshold = 0.9,
                             attach_fits = FALSE, X = NULL, y = NULL) {
  aic_tab <- forest$aic_by_model
  if (!all(names(pi) %in% forest$meta$scope)) {
    warning("Names of pi do not fully match scope; proceeding with intersection.")
  }
  AIC_min <- min(aic_tab$AIC, na.rm = TRUE)
  keep <- which(aic_tab$AIC <= AIC_min + Delta)
  kept <- aic_tab[keep, , drop = FALSE]
  kept$terms <- lapply(kept$model_code, .vars_in_code)
  # mean stability per model
  kept$mean_pi <- vapply(kept$terms, function(S) {
    if (length(S) == 0) return(0)
    mean(pi[S], na.rm = TRUE)
  }, numeric(1))
  kept <- kept[kept$mean_pi >= tau, , drop = FALSE]

  # drop near-duplicates by Jaccard
  if (drop_jaccard && nrow(kept) > 1) {
    ord <- order(kept$AIC)
    kept <- kept[ord, , drop = FALSE]
    to_drop <- logical(nrow(kept))
    for (i in seq_len(nrow(kept))) {
      if (to_drop[i]) next
      Si <- kept$terms[[i]]
      for (j in seq((i+1), nrow(kept))) {
        if (to_drop[j]) next
        Sj <- kept$terms[[j]]
        if (.jaccard(Si, Sj) >= jaccard_threshold) {
          to_drop[j] <- TRUE
        }
      }
    }
    kept <- kept[!to_drop, , drop = FALSE]
  }

  rownames(kept) <- NULL

  if (attach_fits) {
    if (is.null(X) || is.null(y)) stop("Provide X and y when attach_fits = TRUE.")
    family <- forest$meta$family
    fits <- lapply(kept$terms, function(S) .fit_and_aic(as.data.frame(X), y, family, S)$fit)
    kept$fit <- fits
  }

  kept
}
