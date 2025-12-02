## ------------------------------------------------------------
## Multi-path forward selection + stability + plausible models
## ------------------------------------------------------------

#' Multi-path forward selection using AIC
#'
#' @param X matrix or data.frame of predictors (n x p)
#' @param y response vector (length n)
#' @param family "gaussian" (lm) or "binomial" (glm with logit link)
#' @param K max number of forward steps (max model size)
#' @param eps minimum AIC improvement required to expand a parent
#' @param delta AIC tolerance around the best child per parent
#' @param L max number of models kept per level (global cap)
#' @param keep_fits logical; if TRUE, store fitted model objects
#'
#' @return object of class "path_forest" with components:
#'   - frontiers: list of data.frames of models by step
#'   - aic_by_model: data.frame of unique models and AIC
#'   - meta: list of metadata (family, K, eps, delta, L, varnames, ...)
build_paths <- function(
    X,
    y,
    family = c("gaussian", "binomial"),
    K = NULL,
    eps = 0,
    delta = 2,
    L = 50,
    keep_fits = TRUE
) {
  family <- match.arg(family)
  X <- as.data.frame(X)
  n <- NROW(X)
  p <- NCOL(X)
  
  if (length(y) != n) {
    stop("length(y) must equal nrow(X)")
  }
  
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("x", seq_len(p))
  }
  
  if (is.null(K)) {
    K <- p
  }
  K <- min(K, p)
  
  # Environment for caching fitted models by variable set
  model_cache <- new.env(parent = emptyenv())
  model_counter <- 0L
  
  model_key <- function(vars) {
    vars <- sort(as.integer(vars))
    if (length(vars) == 0L) {
      "0"  # special key for empty model
    } else {
      paste(vars, collapse = ",")
    }
  }
  
  fit_model <- function(vars) {
    key <- model_key(vars)
    if (exists(key, envir = model_cache, inherits = FALSE)) {
      return(get(key, envir = model_cache, inherits = FALSE))
    }
    
    dat <- cbind(y = y, X)
    
    if (length(vars) == 0L) {
      form <- y ~ 1
    } else {
      form <- as.formula(
        paste("y ~", paste(colnames(X)[vars], collapse = " + "))
      )
    }
    
    if (family == "gaussian") {
      fit <- stats::lm(form, data = dat)
    } else {
      fit <- stats::glm(form, data = dat, family = stats::binomial())
    }
    
    aic <- stats::AIC(fit)
    
    model_counter <<- model_counter + 1L
    model <- list(
      id   = model_counter,
      key  = key,
      vars = sort(as.integer(vars)),
      size = length(vars),
      aic  = aic,
      fit  = if (keep_fits) fit else NULL
    )
    
    assign(key, model, envir = model_cache)
    model
  }
  
  # Step 0: empty model
  empty <- fit_model(integer(0))
  
  frontiers <- vector("list", K + 1L)
  names(frontiers) <- paste0("k", 0:K)
  frontiers[[1L]] <- data.frame(
    id   = empty$id,
    key  = empty$key,
    step = 0L,
    size = 0L,
    aic  = empty$aic,
    vars = I(list(empty$vars)),
    stringsAsFactors = FALSE
  )
  
  # Forward steps 1..K
  for (k in seq_len(K)) {
    parents <- frontiers[[k]]
    if (is.null(parents) || nrow(parents) == 0L) break
    
    cand_list <- list()
    cl_idx <- 1L
    
    for (i in seq_len(nrow(parents))) {
      parent_vars <- parents$vars[[i]]
      parent_aic  <- parents$aic[i]
      
      unused <- setdiff(seq_len(p), parent_vars)
      if (length(unused) == 0L) next
      
      # Children from this parent
      child_models <- lapply(unused, function(j) {
        m <- fit_model(c(parent_vars, j))
        list(
          parent_row = i,
          id   = m$id,
          key  = m$key,
          vars = m$vars,
          size = m$size,
          aic  = m$aic
        )
      })
      
      child_df <- do.call(
        rbind,
        lapply(child_models, function(cm) {
          data.frame(
            parent_row = cm$parent_row,
            id   = cm$id,
            key  = cm$key,
            size = cm$size,
            aic  = cm$aic,
            stringsAsFactors = FALSE
          )
        })
      )
      child_df$vars <- I(lapply(child_models, `[[`, "vars"))
      
      cand_list[[cl_idx]] <- child_df
      cl_idx <- cl_idx + 1L
    }
    
    if (length(cand_list) == 0L) break
    
    cand <- do.call(rbind, cand_list)
    
    # Parent-wise filtering by eps and delta
    next_models <- list()
    nm_idx <- 1L
    
    for (i in seq_len(nrow(parents))) {
      parent_aic <- parents$aic[i]
      idx_i <- which(cand$parent_row == i)
      if (length(idx_i) == 0L) next
      
      sub_i <- cand[idx_i, ]
      best_aic <- min(sub_i$aic)
      
      # require best child improves by at least eps
      if ((parent_aic - best_aic) < eps) {
        next
      }
      
      keep_idx <- idx_i[which(sub_i$aic <= best_aic + delta)]
      if (length(keep_idx) == 0L) next
      
      next_models[[nm_idx]] <- cand[keep_idx, c("id", "key", "size", "aic", "vars")]
      nm_idx <- nm_idx + 1L
    }
    
    if (length(next_models) == 0L) break
    
    next_df <- do.call(rbind, next_models)
    
    # Deduplicate by key, keep best AIC per key
    o <- order(next_df$key, next_df$aic)
    next_df <- next_df[o, ]
    next_df <- next_df[!duplicated(next_df$key), ]
    
    # Global cap L by AIC
    if (!is.null(L) && nrow(next_df) > L) {
      o <- order(next_df$aic)
      next_df <- next_df[o[seq_len(L)], ]
    }
    
    next_df$step <- k
    next_df <- next_df[, c("id", "key", "step", "size", "aic", "vars")]
    
    frontiers[[k + 1L]] <- next_df
  }
  
  # Drop trailing NULL/empty levels
  nonempty <- which(vapply(
    frontiers,
    function(df) !is.null(df) && nrow(df) > 0L,
    logical(1)
  ))
  frontiers <- frontiers[nonempty]
  
  # Build aic_by_model from the cache
  all_models <- as.list(model_cache)
  aic_by_model <- do.call(
    rbind,
    lapply(all_models, function(m) {
      data.frame(
        id   = m$id,
        key  = m$key,
        size = m$size,
        aic  = m$aic,
        stringsAsFactors = FALSE
      )
    })
  )
  aic_by_model$vars <- I(lapply(all_models, `[[`, "vars"))
  if (keep_fits) {
    aic_by_model$fit <- I(lapply(all_models, `[[`, "fit"))
  } else {
    aic_by_model$fit <- I(vector("list", nrow(aic_by_model)))
  }
  
  # Order by id
  aic_by_model <- aic_by_model[order(aic_by_model$id), ]
  
  structure(
    list(
      frontiers     = frontiers,
      aic_by_model  = aic_by_model,
      meta          = list(
        family   = family,
        K        = K,
        eps      = eps,
        delta    = delta,
        L        = L,
        p        = p,
        varnames = colnames(X)
      )
    ),
    class = "path_forest"
  )
}


#' Resampling-based stability (model-set proportions)
#'
#' For each resample b, runs build_paths() and computes z_j^(b) =
#' proportion of models containing feature j (over unique non-empty models).
#' Then π_j = (1/B) Σ_b z_j^(b).
#'
#' @param X matrix or data.frame of predictors
#' @param y response vector
#' @param B number of resamples
#' @param family "gaussian" or "binomial"
#' @param resample "bootstrap" or "subsample"
#' @param m subsample size if resample == "subsample"
#' @param build_paths_args list of extra arguments passed to build_paths()
#'
#' @return object of class "path_stability" with components:
#'   - pi: length-p vector of stability scores
#'   - Z: B x p matrix of z_j^(b)
#'   - meta: list with resampling info
stability <- function(
    X,
    y,
    B = 50L,
    family = c("gaussian", "binomial"),
    resample = c("bootstrap", "subsample"),
    m = NULL,
    build_paths_args = list()
) {
  family   <- match.arg(family)
  resample <- match.arg(resample)
  
  X <- as.data.frame(X)
  n <- NROW(X)
  p <- NCOL(X)
  
  if (length(y) != n) {
    stop("length(y) must equal nrow(X)")
  }
  
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("x", seq_len(p))
  }
  
  if (resample == "subsample") {
    if (is.null(m)) {
      m <- floor(0.8 * n)
    }
    if (m <= 0L || m > n) stop("Invalid subsample size m.")
  }
  
  Z <- matrix(NA_real_, nrow = B, ncol = p)
  colnames(Z) <- colnames(X)
  
  for (b in seq_len(B)) {
    if (resample == "bootstrap") {
      idx <- sample.int(n, n, replace = TRUE)
    } else {
      idx <- sample.int(n, m, replace = FALSE)
    }
    
    Xb <- X[idx, , drop = FALSE]
    yb <- y[idx]
    
    forest_b <- do.call(
      build_paths,
      c(list(X = Xb, y = yb, family = family), build_paths_args)
    )
    
    models_df <- forest_b$aic_by_model
    
    if (nrow(models_df) == 0L) {
      Z[b, ] <- 0
    } else {
      # Use non-empty models for stability
      if ("size" %in% names(models_df)) {
        models_df <- models_df[models_df$size > 0L, ]
      }
      if (nrow(models_df) == 0L) {
        Z[b, ] <- 0
      } else {
        model_vars <- models_df$vars
        z_j <- numeric(p)
        for (j in seq_len(p)) {
          z_j[j] <- mean(vapply(
            model_vars,
            function(v) j %in% v,
            logical(1)
          ))
        }
        Z[b, ] <- z_j
      }
    }
  }
  
  pi_hat <- colMeans(Z)
  
  structure(
    list(
      pi   = pi_hat,
      Z    = Z,
      meta = list(
        B               = B,
        family          = family,
        resample        = resample,
        m               = m,
        p               = p,
        varnames        = colnames(X),
        build_paths_args = build_paths_args
      )
    ),
    class = "path_stability"
  )
}


#' Plausible model selection: AIC filter + average stability
#'
#' @param forest path_forest object from build_paths()
#' @param pi numeric vector of length p with stability scores
#' @param Delta AIC tolerance Δ above AIC_min
#' @param tau stability threshold τ on mean π(S)
#' @param drop_duplicates logical; if TRUE, drop near-duplicates
#' @param jaccard_threshold Jaccard similarity threshold for near-duplicates
#'
#' @return data.frame of plausible models with columns:
#'   id, key, size, aic, vars (indices), vars_char (names), mean_pi, fit (if stored)
plausible_models <- function(
    forest,
    pi,
    Delta = 2,
    tau = 0.5,
    drop_duplicates = TRUE,
    jaccard_threshold = 0.9
) {
  if (missing(pi)) {
    stop("Argument 'pi' (stability scores) must be supplied.")
  }
  
  models_df <- forest$aic_by_model
  if (!("aic" %in% names(models_df))) {
    stop("forest$aic_by_model must contain an 'aic' column.")
  }
  
  p <- length(pi)
  
  AIC_min <- min(models_df$aic)
  keep_aic <- models_df$aic <= (AIC_min + Delta)
  
  model_vars <- models_df$vars
  mean_pi <- vapply(
    model_vars,
    function(v) {
      if (length(v) == 0L) {
        NA_real_
      } else {
        mean(pi[v])
      }
    },
    numeric(1)
  )
  
  keep_stab <- (!is.na(mean_pi)) & (mean_pi >= tau)
  cand_idx <- which(keep_aic & keep_stab)
  
  if (length(cand_idx) == 0L) {
    out <- models_df[0, , drop = FALSE]
    out$mean_pi   <- numeric(0)
    out$vars_idx  <- I(list())
    out$vars_char <- I(list())
    return(out)
  }
  
  cand_df <- models_df[cand_idx, ]
  cand_df$mean_pi <- mean_pi[cand_idx]
  
  varnames <- forest$meta$varnames
  cand_df$vars_idx  <- cand_df$vars
  cand_df$vars_char <- I(lapply(
    cand_df$vars_idx,
    function(v) if (length(v) == 0L) character(0) else varnames[v]
  ))
  
  if (drop_duplicates && nrow(cand_df) > 1L) {
    o <- order(cand_df$aic)
    cand_df <- cand_df[o, ]
    
    keep <- logical(nrow(cand_df))
    kept_vars <- list()
    
    for (i in seq_len(nrow(cand_df))) {
      v <- cand_df$vars_idx[[i]]
      
      if (length(kept_vars) == 0L) {
        keep[i] <- TRUE
        kept_vars[[1L]] <- v
      } else {
        jac_high <- FALSE
        for (kv in kept_vars) {
          inter <- length(intersect(v, kv))
          uni   <- length(union(v, kv))
          jac   <- if (uni == 0L) 1 else inter / uni
          if (jac >= jaccard_threshold) {
            jac_high <- TRUE
            break
          }
        }
        if (!jac_high) {
          keep[i] <- TRUE
          kept_vars[[length(kept_vars) + 1L]] <- v
        }
      }
    }
    
    cand_df <- cand_df[keep, ]
  }
  
  row.names(cand_df) <- NULL
  cand_df
}


#' (Optional) Confusion metrics for logistic models
#'
#' Computes confusion matrix and basic metrics at cutoff 0.5:
#' prevalence, accuracy, sensitivity, specificity, FDR, DOR.
#'
#' @param fit glm object with binomial family
#' @param y optional true labels (0/1 or 2-level factor). If NULL,
#'   tries to use fit$y.
#' @param newdata optional newdata for prediction; if NULL, uses training data.
#' @param cutoff probability threshold (default 0.5)
#'
#' @return list with 'confusion_matrix' and 'metrics' (named numeric vector)
confusion_metrics <- function(
    fit,
    y = NULL,
    newdata = NULL,
    cutoff = 0.5
) {
  if (!inherits(fit, "glm")) {
    stop("'fit' must be a glm object.")
  }
  fam <- stats::family(fit)$family
  if (fam != "binomial") {
    stop("'fit' must have binomial family.")
  }
  
  if (is.null(newdata)) {
    probs <- stats::fitted(fit)
  } else {
    probs <- stats::predict(fit, newdata = newdata, type = "response")
  }
  
  if (is.null(y)) {
    if (!is.null(fit$y)) {
      y <- fit$y
    } else {
      stop("Argument 'y' must be supplied if fit$y is not available.")
    }
  }
  
  if (is.factor(y)) {
    if (nlevels(y) != 2L) stop("y must be binary.")
    y <- as.integer(y == levels(y)[2L])
  } else {
    y <- as.numeric(y)
  }
  
  if (length(y) != length(probs)) {
    stop("length(y) must equal length of predicted probabilities.")
  }
  
  pred <- as.integer(probs >= cutoff)
  
  TP <- sum(pred == 1 & y == 1, na.rm = TRUE)
  TN <- sum(pred == 0 & y == 0, na.rm = TRUE)
  FP <- sum(pred == 1 & y == 0, na.rm = TRUE)
  FN <- sum(pred == 0 & y == 1, na.rm = TRUE)
  
  total <- TP + TN + FP + FN
  
  prevalence  <- (TP + FN) / total
  accuracy    <- (TP + TN) / total
  sensitivity <- if ((TP + FN) == 0) NA_real_ else TP / (TP + FN)
  specificity <- if ((TN + FP) == 0) NA_real_ else TN / (TN + FP)
  FDR         <- if ((TP + FP) == 0) NA_real_ else FP / (TP + FP)
  DOR         <- if (TP == 0 || TN == 0 || FP == 0 || FN == 0) {
    NA_real_
  } else {
    (TP * TN) / (FP * FN)
  }
  
  cm <- matrix(
    c(TN, FP, FN, TP),
    nrow = 2, byrow = TRUE,
    dimnames = list(
      Predicted = c("0", "1"),
      True      = c("0", "1")
    )
  )
  
  list(
    confusion_matrix = cm,
    metrics = c(
      prevalence  = prevalence,
      accuracy    = accuracy,
      sensitivity = sensitivity,
      specificity = specificity,
      FDR         = FDR,
      DOR         = DOR
    )
  )
}

# Convert each frontier level into a Section-5 style table:
# first column = variable(s), second column = AIC.
frontier_tables <- function(forest) {
  frontiers <- forest$frontiers
  varnames  <- forest$meta$varnames
  
  lapply(frontiers, function(df) {
    # df$vars is a list of index vectors
    vars_char <- vapply(
      df$vars,
      function(v) {
        if (length(v) == 0L) {
          "(empty)"
        } else {
          paste(varnames[v], collapse = " + ")
        }
      },
      character(1)
    )
    
    out <- data.frame(
      variable = vars_char,
      AIC      = df$aic,
      stringsAsFactors = FALSE
    )
    
    # Sort by AIC, like the HTML examples
    out[order(out$AIC), , drop = FALSE]
  })
}


## ------------------------------------------------------------
## Minimal runnable examples (Section 5-style, base R only)
## ------------------------------------------------------------

# Example 5.1: Linear regression (Gaussian)
example_gaussian <- function(print_tables = TRUE) {
  set.seed(1)
  n <- 120; p <- 8
  X <- matrix(rnorm(n * p), n, p)
  beta <- c(2, -1.5, 0, 0, 1, rep(0, p - 5))
  y <- as.numeric(X %*% beta + rnorm(n, sd = 1))
  colnames(X) <- paste0("x", 1:p)
  
  # Build the path forest
  forest <- build_paths(
    X, y,
    family = "gaussian",
    K = 5, eps = 1, delta = 2, L = 20
  )
  
  # Stability (π_j)
  stab <- stability(
    X, y,
    B = 50,
    family = "gaussian",
    resample = "bootstrap",
    build_paths_args = list(K = 5, eps = 1, delta = 2, L = 20)
  )

  # Plausible models: AIC + stability
  plaus <- plausible_models(
    forest,
    pi = stab$pi,
    Delta = 2,
    tau = 0.6
  )
  
  if (print_tables) {
    cat("========================================\n")
    cat("Gaussian example: forward-search frontiers\n")
    cat("========================================\n\n")
    
    tabs <- frontier_tables(forest)
    for (nm in names(tabs)) {
      cat("Frontier", nm, "(step", gsub("^k", "", nm), ")\n")
      print(tabs[[nm]], row.names = FALSE)
      cat("\n")
    }
    
    cat("========================================\n")
    cat("Gaussian example: plausible models\n")
    cat("========================================\n\n")
    
    if (nrow(plaus) == 0L) {
      cat("No plausible models under the chosen (Delta, tau).\n")
    } else {
      # Show variables, AIC, and mean stability
      df_plaus <- data.frame(
        model      = vapply(
          plaus$vars_char,
          function(v) if (length(v) == 0L) "(empty)" else paste(v, collapse = " + "),
          character(1)
        ),
        AIC        = plaus$aic,
        mean_pi    = plaus$mean_pi,
        stringsAsFactors = FALSE
      )
      df_plaus <- df_plaus[order(df_plaus$AIC), , drop = FALSE]
      print(df_plaus, row.names = FALSE)
    }
  }
  
  invisible(list(forest = forest, stability = stab, plausible = plaus))
}

# Example 5.2: Logistic regression (Binomial)
example_binomial <- function(print_tables = TRUE) {
  set.seed(2)
  n <- 200; p <- 6
  X <- matrix(rnorm(n * p), n, p)
  linpred <- 1.2 * X[, 1] - 1 * X[, 2] + 0.8 * X[, 3]
  prob <- 1 / (1 + exp(-linpred))
  y <- rbinom(n, size = 1, prob = prob)
  colnames(X) <- paste0("x", 1:p)
  
  # Build the path forest
  forest <- build_paths(
    X, y,
    family = "binomial",
    K = 4, eps = 1, delta = 2, L = 20
  )
  
  # Stability (π_j)
  
  if (FALSE) {
  stab <- stability(
    X, y,
    B = 50,
    family = "binomial",
    resample = "bootstrap",
    build_paths_args = list(K = 4, eps = 1, delta = 2, L = 20)
  )
  } else {
  
  stab <- stability(
    X, y,
    B = 200,               # more resamples
    family = "binomial",
    resample = "bootstrap",
    build_paths_args = list(K = 4, eps = 1, delta = 2, L = 20)
  )
  }
  
  
  # Plausible models
  plaus <- plausible_models(
    forest,
    pi = stab$pi,
    Delta = 2,
    tau = 0.35, #0.6
    drop_duplicates   = FALSE   # make it easy to see *something*
  )
  
  if (print_tables) {
    cat("========================================\n")
    cat("Binomial example: forward-search frontiers\n")
    cat("========================================\n\n")
    
    tabs <- frontier_tables(forest)
    for (nm in names(tabs)) {
      cat("Frontier", nm, "(step", gsub("^k", "", nm), ")\n")
      print(tabs[[nm]], row.names = FALSE)
      cat("\n")
    }
    
    cat("========================================\n")
    cat("Binomial example: plausible models\n")
    cat("========================================\n\n")
    
    if (nrow(plaus) == 0L) {
      cat("No plausible models under the chosen (Delta, tau).\n")
    } else {
      df_plaus <- data.frame(
        model      = vapply(
          plaus$vars_char,
          function(v) if (length(v) == 0L) "(empty)" else paste(v, collapse = " + "),
          character(1)
        ),
        AIC        = plaus$aic,
        mean_pi    = plaus$mean_pi,
        stringsAsFactors = FALSE
      )
      df_plaus <- df_plaus[order(df_plaus$AIC), , drop = FALSE]
      print(df_plaus, row.names = FALSE)
    }
  }
  
  invisible(list(forest = forest, stability = stab, plausible = plaus))
}

# Convert each frontier level into a Section-5 style table:
# first column = variable(s), second column = AIC.
frontier_tables <- function(forest) {
  frontiers <- forest$frontiers
  varnames  <- forest$meta$varnames
  
  lapply(frontiers, function(df) {
    # df$vars is a list of index vectors
    vars_char <- vapply(
      df$vars,
      function(v) {
        if (length(v) == 0L) {
          "(empty)"
        } else {
          paste(varnames[v], collapse = " + ")
        }
      },
      character(1)
    )
    
    out <- data.frame(
      variable = vars_char,
      AIC      = df$aic,
      stringsAsFactors = FALSE
    )
    
    # Sort by AIC, like the HTML examples
    out[order(out$AIC), , drop = FALSE]
  })
}


#example_gaussian()

example_binomial()

# Assuming you captured the result like:
res <- example_binomial(print_tables = FALSE)
forest <- res$forest
stab   <- res$stability

# 1. Look at stability scores
stab$pi
summary(stab$pi)

# 2. Look at models passing just the AIC filter
models_df <- forest$aic_by_model
AIC_min   <- min(models_df$aic)
Delta     <- 2

cand_AIC <- subset(models_df, aic <= AIC_min + Delta)
nrow(cand_AIC)


pi <- stab$pi

mean_pi <- vapply(
  cand_AIC$vars,
  function(v) {
    if (length(v) == 0L) NA_real_ else mean(pi[v])
  },
  numeric(1)
)

summary(mean_pi)

plaus <- plausible_models(
  forest,
  pi = stab$pi,
  Delta = 4,   # looser AIC tolerance
  tau   = 0.4  # looser stability threshold
)
plaus