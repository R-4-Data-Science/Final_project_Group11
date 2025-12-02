#' Internal: multi-path forward selection using AIC
#'
#' @param X Data frame of predictors
#' @param y Response vector
#' @param model_type "gaussian" or "binomial"
#' @param K Max steps (model size)
#' @param eps Minimum AIC improvement
#' @param delta AIC tolerance for near-ties
#' @param L Max models per step
#'
#' @return A list with variable names, models by step, and AICs by step
#' @keywords internal
multi_path_forward <- function(
    X, y,
    model_type = c("gaussian", "binomial"),
    K = NULL,
    eps = 1e-6,
    delta = 2,
    L = 50
) {
  model_type <- match.arg(model_type)
  X <- as.data.frame(X)

  p <- ncol(X)
  var_names <- colnames(X)
  if (is.null(K)) K <- min(p, 10)

  model_key <- function(idx) {
    if (length(idx) == 0) return("")
    paste(sort(idx), collapse = ",")
  }

  parent_models <- list(integer(0))
  parent_AICs   <- fit_model_aic(X, y, vars = integer(0), model_type = model_type)

  step_models <- list()
  step_AICs   <- list()

  for (k in seq_len(K)) {
    children_list <- list()
    children_AICs <- numeric(0)
    children_keys <- character(0)

    for (m in seq_along(parent_models)) {
      parent_idx <- parent_models[[m]]
      parent_aic <- parent_AICs[m]

      remaining <- setdiff(seq_len(p), parent_idx)
      if (length(remaining) == 0) next

      cand_models <- list()
      cand_AICs   <- numeric(0)

      for (j in remaining) {
        child_idx <- sort(c(parent_idx, j))
        aic_child <- fit_model_aic(X, y, vars = child_idx, model_type = model_type)
        cand_models[[length(cand_models) + 1]] <- child_idx
        cand_AICs[length(cand_AICs) + 1] <- aic_child
      }

      if (!length(cand_AICs)) next

      best_child_AIC <- min(cand_AICs)

      if ((parent_aic - best_child_AIC) < eps) {
        next
      }

      keep_idx <- which(cand_AICs <= best_child_AIC + delta)
      for (i_keep in keep_idx) {
        child_idx <- cand_models[[i_keep]]
        aic_child <- cand_AICs[i_keep]
        key       <- model_key(child_idx)

        children_list[[length(children_list) + 1]] <- child_idx
        children_AICs[length(children_AICs) + 1]   <- aic_child
        children_keys[length(children_keys) + 1]   <- key
      }
    }

    if (!length(children_list)) break

    df_children <- data.frame(
      key = children_keys,
      AIC = children_AICs,
      stringsAsFactors = FALSE
    )

    agg <- aggregate(AIC ~ key, data = df_children, FUN = min)
    agg <- agg[order(agg$AIC), ]

    if (!is.null(L) && nrow(agg) > L) {
      agg <- agg[seq_len(L), ]
    }

    new_parents <- vector("list", nrow(agg))
    new_AICs    <- numeric(nrow(agg))

    for (i in seq_len(nrow(agg))) {
      key <- agg$key[i]
      if (key == "") {
        idx <- integer(0)
      } else {
        idx <- as.integer(strsplit(key, ",")[[1]])
      }
      new_parents[[i]] <- idx
      new_AICs[i]      <- agg$AIC[i]
    }

    step_models[[k]] <- new_parents
    step_AICs[[k]]   <- new_AICs

    parent_models <- new_parents
    parent_AICs   <- new_AICs
  }

  structure(
    list(
      var_names   = var_names,
      step_models = step_models,
      step_AICs   = step_AICs
    ),
    class = "path_forest"
  )
}

#' Build multi-path AIC model selection forest
#'
#' @param X Data frame of predictors
#' @param y Response vector
#' @param family "gaussian" or "binomial"
#' @param K Max model size (steps)
#' @param eps Minimum AIC improvement
#' @param delta AIC tolerance for near-ties
#' @param L Max models to retain per step
#'
#' @return An object of class \code{"path_forest"} with models and AICs
#' @export
build_paths <- function(
    X, y,
    family = c("gaussian", "binomial"),
    K = NULL,
    eps = 1e-6,
    delta = 2,
    L = 50
) {
  family <- match.arg(family)
  multi_path_forward(
    X = X, y = y,
    model_type = family,
    K = K, eps = eps, delta = delta, L = L
  )
}
