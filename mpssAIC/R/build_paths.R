#' Multi-path forward selection using AIC
#'
#' Builds a branching forest of forward-selection models by AIC.
#'
#' @param X data.frame or matrix of predictors.
#' @param y response vector.
#' @param family "gaussian" or "binomial".
#' @param K maximum number of steps (added variables).
#' @param eps numeric. Minimum AIC improvement required over the parent model
#'   for any children to be considered (prevents adding useless variables).
#' @param delta numeric. Keep all children within `delta` AIC of the *best child*
#'   for each parent.
#' @param L integer. If the frontier at a level becomes too large after merging
#'   across parents, keep only the best `L` by AIC.
#' @param scope optional character vector naming allowed columns of X.
#' @param trace logical. If TRUE, prints brief progress.
#'
#' @return A list with elements:
#'   - frontiers: list of data.frames (one per step), with columns:
#'       terms (list-column of character vectors), AIC, parent_code, model_code, fit
#'   - aic_by_model: data.frame of unique models across all steps.
#'   - meta: list of parameters and info.
#' @export
build_paths <- function(X, y, family = c("gaussian","binomial"),
                        K = 5, eps = 1e-6, delta = 2, L = 50,
                        scope = NULL, trace = FALSE) {
  family <- match.arg(family)
  X <- as.data.frame(X, check.names = FALSE)
  if (is.null(scope)) scope <- colnames(X)
  p <- length(scope)
  if (p == 0) stop("No predictors in scope.")
  # step 0 frontier: just the empty model
  step0 <- data.frame(
    terms = I(list(character(0))),
    AIC = .parents_aic(X, y, family, character(0)),
    parent_code = NA_character_,
    model_code = "~1",
    stringsAsFactors = FALSE
  )
  frontiers <- list(step0)
  all_models <- step0

  for (k in seq_len(K)) {
    if (trace) message("Step ", k, " ...")
    parents <- frontiers[[k]]
    children_list <- vector("list", nrow(parents))

    for (i in seq_len(nrow(parents))) {
      parent_terms <- parents$terms[[i]]
      parent_code <- parents$model_code[i]
      remaining <- setdiff(scope, parent_terms)

      # generate child candidates by adding one unused variable
      cand <- lapply(remaining, function(v) sort(c(parent_terms, v)))
      if (length(cand) == 0) next

      # fit all candidates, compute AIC
      fits <- lapply(cand, function(tt) .fit_and_aic(X, y, family, tt))
      aics <- vapply(fits, function(z) z$AIC, numeric(1))
      best_idx <- which.min(aics)
      best_aic <- aics[best_idx]
      parent_aic <- parents$AIC[i]

      # gate: only keep children if best child improves parent by >= eps
      if (!is.finite(parent_aic) || !is.finite(best_aic) || (parent_aic - best_aic) < eps) next

      keep <- which(aics <= best_aic + delta)
      df <- data.frame(
        terms = I(cand[keep]),
        AIC = aics[keep],
        parent_code = parent_code,
        model_code = vapply(cand[keep], .encode_model, character(1)),
        stringsAsFactors = FALSE
      )
      # attach fits (list-column)
      df$fit <- fits[keep]
      children_list[[i]] <- df
    }

    # merge children from all parents
    children <- do.call(rbind, children_list)
    if (is.null(children) || nrow(children) == 0) break

    # deduplicate by model_code; keep best AIC entry
    children <- .dedup_by_model(children)

    # size control
    if (nrow(children) > L) {
      ord <- order(children$AIC)
      children <- children[ord[seq_len(L)], , drop = FALSE]
    }

    frontiers[[k + 1]] <- children

    # collect into all_models
    all_models <- rbind(all_models,
                        children[, c("terms","AIC","model_code"), drop = FALSE])
    all_models <- .dedup_by_model(all_models)
  }

  aic_by_model <- all_models[, c("model_code","AIC"), drop = FALSE]
  rownames(aic_by_model) <- NULL

  list(
    frontiers = frontiers,
    aic_by_model = aic_by_model,
    meta = list(family = family, K = K, eps = eps, delta = delta, L = L, p = p, scope = scope)
  )
}
