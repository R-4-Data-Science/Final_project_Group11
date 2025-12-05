# inst/shiny-examples/mpssAIC-demo/app.R

library(shiny)
library(mpssAIC)  # your package with multi_path_AIC_pipeline(), confusion_metrics(), etc.

# -------------------------------------------------------------

# Core Helper: fit_model_aic()

# -------------------------------------------------------------

#' Compute AIC for a model with selected predictors
#'
#' @param X A data frame of predictors (n × p)
#' @param y Response vector of length n
#' @param vars Integer vector of column indices to include (empty = intercept-only)
#' @param model_type Either "gaussian" or "binomial"
#'
#' @return Numeric AIC value
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


#' Algorithm 2: Compute variable stability via resampling
#'
#' @param X Data frame of predictors
#' @param y Response vector
#' @param model_type "gaussian" or "binomial"
#' @param B Number of resamples
#' @param resample_type "bootstrap" or "subsample"
#' @param m Subsample size (if using subsample)
#' @param K, eps, delta, L Algorithm 1 parameters
#'
#' @return Numeric vector of stability scores (one per predictor)
compute_stability <- function(
    X, y,
    model_type = c("gaussian", "binomial"),
    B = 50,
    resample_type = c("bootstrap", "subsample"),
    m = NULL,
    K = NULL,
    eps = 1e-6,
    delta = 2,
    L = 50
) {
  model_type <- match.arg(model_type)
  resample_type <- match.arg(resample_type)
  
  X <- as.data.frame(X)
  n <- nrow(X)
  p <- ncol(X)
  var_names <- colnames(X)
  
  if (is.null(m)) m <- ceiling(sqrt(n))
  # Matrix to store proportions: rows = resamples, cols = variables
  
  Z <- matrix(0, nrow = B, ncol = p)
  colnames(Z) <- var_names
  
  for (b in seq_len(B)) {
    if (b %% 10 == 0) cat("Resample", b, "/", B, "\n")
    
    # Generate resample
    if (resample_type == "bootstrap") {
      idx <- sample(seq_len(n), size = n, replace = TRUE)
    } else {
      idx <- sample(seq_len(n), size = m, replace = FALSE)
    }
    
    Xb <- X[idx, , drop = FALSE]
    yb <- y[idx]
    
    # Run multi-path search on resample
    mp_b <- multi_path_forward(
      X = Xb, y = yb,
      model_type = model_type,
      K = K,
      eps = eps,
      delta = delta,
      L = L
    )
    
    # Collect all models from this resample
    all_models_b <- unlist(mp_b$step_models, recursive = FALSE)
    if (!length(all_models_b)) next
    
    M <- length(all_models_b)
    
    # For each variable, count how many models contain it
    for (j in seq_len(p)) {
      count_j <- sum(vapply(
        all_models_b,
        function(idx) j %in% idx,
        logical(1)
      ))
      Z[b, j] <- count_j / M
    }
    
  }
  
  #Average over resamples
  
  pi <- colMeans(Z)
  
  cat("\n=== Algorithm 2 Complete ===\n")
  cat("Stability scores (π):\n")
  print(round(pi, 3))
  
  pi
}

#-------------------------------------------------------------
#Algorithm 1: multi_path_forward()
#-------------------------------------------------------------

#' Algorithm 1: Multiple-Path Forward Selection
#'
#' @param X Data frame of predictors
#' @param y Response vector
#' @param model_type "gaussian" or "binomial"
#' @param K Maximum number of steps (defaults to min(p, 10))
#' @param eps Minimum AIC improvement required to expand a model
#' @param delta AIC tolerance for keeping near-ties
#' @param L Maximum number of models to keep per step
#'
#' @return List with:
#' - var_names: predictor names
#' - step_models: models at each step (list of lists)
#' - step_AICs: corresponding AICs
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
  
  #Helper to create unique key for a model
  
  model_key <- function(idx) {
    if (length(idx) == 0) return("")
    paste(sort(idx), collapse = ",")
  }
  
  #Initialize: parent is empty model (intercept only)
  
  parent_models <- list(integer(0))
  parent_AICs <- fit_model_aic(X, y, vars = integer(0), model_type = model_type)
  
  step_models <- list()
  step_AICs <- list()
  
  for (k in seq_len(K)) {
    cat("Step", k, ": processing", length(parent_models), "parent model(s)...\n")
    
    children_list <- list()
    children_AICs <- numeric(0)
    children_keys <- character(0)
    
    # For each parent, generate all possible children
    for (m in seq_along(parent_models)) {
      parent_idx <- parent_models[[m]]
      parent_aic <- parent_AICs[m]
      
      remaining <- setdiff(seq_len(p), parent_idx)
      if (length(remaining) == 0) next  # no more variables to add
      
      cand_models <- list()
      cand_AICs   <- numeric(0)
      
      # Try adding each remaining variable
      for (j in remaining) {
        child_idx <- sort(c(parent_idx, j))
        aic_child <- fit_model_aic(X, y, vars = child_idx, model_type = model_type)
        cand_models[[length(cand_models) + 1]] <- child_idx
        cand_AICs[length(cand_AICs) + 1]       <- aic_child
      }
      
      if (!length(cand_AICs)) next
      
      best_child_AIC <- min(cand_AICs)
      
      # Check: is there at least eps improvement?
      if ((parent_aic - best_child_AIC) < eps) {
        cat("  Parent", m, ": no improvement ≥ eps; stopping expansion.\n")
        next
      }
      
      # Keep children within delta of best child AIC
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
    
    # If no children generated, stop
    if (!length(children_list)) {
      cat("No further improvement. Stopping.\n")
      break
    }
    
    # Deduplicate: keep best AIC per unique model key
    df_children <- data.frame(
      key = children_keys,
      AIC = children_AICs,
      stringsAsFactors = FALSE
    )
    
    agg <- aggregate(AIC ~ key, data = df_children, FUN = min)
    agg <- agg[order(agg$AIC), ]
    
    # Cap at L models
    if (!is.null(L) && nrow(agg) > L) {
      agg <- agg[seq_len(L), ]
    }
    
    cat("  Step", k, "produced", nrow(agg), "unique model(s) after dedup/cap.\n")
    
    # Convert keys back to index vectors
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
  
  cat("\n=== Algorithm 1 Complete ===\n")
  
  list(
    var_names = var_names,
    step_models = step_models,
    step_AICs = step_AICs
  )
}


#-------------------------------------------------------------
#Algorithm 4: multi_path_AIC_pipeline()
#-------------------------------------------------------------

#' Algorithm 4: Complete Multi-Path AIC Pipeline
#'
#' @param X data frame of predictors
#' @param y response vector
#' @param model_type "gaussian" or "binomial"
#' @param K, eps, delta, L Algorithm 1 parameters
#' @param B, resample_type, m Algorithm 2 parameters
#' @param Delta, tau Algorithm 3 parameters
#'
#' @return list with:
#'   - mp_full: Algorithm 1 output
#'   - stability_pi: Algorithm 2 output
#'   - plausible: Algorithm 3 output
multi_path_AIC_pipeline <- function(
    X, y,
    model_type = c("gaussian", "binomial"),
    K = NULL,
    eps = 1e-6,
    delta = 2,
    L = 50,
    B = 50,
    resample_type = c("bootstrap", "subsample"),
    m = NULL,
    Delta = 2,
    tau = 0.6
) {
  model_type    <- match.arg(model_type)
  resample_type <- match.arg(resample_type)
  
  cat("\n", strrep("=", 60), "\n")
  cat("MULTI-PATH AIC PIPELINE\n")
  cat("Model type:", model_type, "\n")
  cat("Resamples (B):", B, "\n")
  cat(strrep("=", 60), "\n\n")
  
  # Step 1: Multi-path forward selection
  cat(">>> STEP 1: Multi-Path Forward Selection <<<\n\n")
  mp_full <- multi_path_forward(
    X = X, y = y,
    model_type = model_type,
    K = K,
    eps = eps,
    delta = delta,
    L = L
  )
  
  # Step 2: Stability via resampling
  cat("\n>>> STEP 2: Stability Estimation via Resampling <<<\n\n")
  stability_pi <- compute_stability(
    X = X, y = y,
    model_type    = model_type,
    B             = B,
    resample_type = resample_type,
    m             = m,
    K = K,
    eps = eps,
    delta = delta,
    L = L
  )
  
  # Step 3: Plausible models
  cat("\n>>> STEP 3: Plausible Model Selection <<<\n\n")
  plausible <- select_plausible_models(
    mp_full      = mp_full,
    stability_pi = stability_pi,
    Delta        = Delta,
    tau          = tau
  )
  
  cat("\n", strrep("=", 60), "\n")
  cat("PIPELINE COMPLETE\n")
  cat(strrep("=", 60), "\n")
  
  list(
    mp_full      = mp_full,
    stability_pi = stability_pi,
    plausible    = plausible
  )
}

# ---------------------------------------------------------
# Algorithm 3: Plausible Model Selection
# ---------------------------------------------------------
# Combines:
#  • AIC quality (within Δ of best model)
#  • Stability (average π ≥ τ)
# ---------------------------------------------------------

#' Algorithm 3: Select plausible, stable models
#'
#' @param mp_full Output list from Algorithm 1
#' @param stability_pi Named vector of stability scores (Algorithm 2)
#' @param Delta AIC tolerance (default = 2)
#' @param tau Minimum average stability threshold (default = 0.6)
#'
#' @return Data frame of plausible models:
#'   columns = key, AIC, AIC_diff, size, vars, avg_stability
#'
select_plausible_models <- function(
    mp_full,
    stability_pi,
    Delta = 2,
    tau = 0.6
) {
  
  var_names   <- mp_full$var_names
  step_models <- mp_full$step_models
  step_AICs   <- mp_full$step_AICs
  
  # Helper: unique key for a set of variables
  model_key <- function(idx) {
    if (length(idx) == 0) return("")
    paste(sort(idx), collapse = ",")
  }
  
  # Flatten all models and AICs
  all_models <- unlist(step_models, recursive = FALSE)
  all_AICs   <- unlist(step_AICs)
  
  if (length(all_models) == 0) {
    cat("No models generated!\n")
    return(data.frame())
  }
  
  # Create model keys
  keys <- vapply(all_models, model_key, character(1))
  
  df <- data.frame(
    key = keys,
    AIC = all_AICs,
    stringsAsFactors = FALSE
  )
  
  # Aggregate by model key → keep minimum AIC per model
  agg <- aggregate(AIC ~ key, data = df, FUN = min)
  
  # Remove intercept-only model
  agg <- agg[agg$key != "", ]
  
  if (nrow(agg) == 0) {
    cat("No non-intercept models!\n")
    return(data.frame())
  }
  
  # Compute AIC differences
  best_AIC      <- min(agg$AIC)
  agg$AIC_diff  <- agg$AIC - best_AIC
  
  # Filter by AIC ≤ best + Delta
  plausible <- agg[agg$AIC_diff <= Delta, ]
  
  if (nrow(plausible) == 0) {
    cat("No models within AIC window!\n")
    return(data.frame())
  }
  
  # Compute size and stability
  size      <- integer(nrow(plausible))
  avg_stab  <- numeric(nrow(plausible))
  vars_str  <- character(nrow(plausible))
  
  for (i in seq_len(nrow(plausible))) {
    key <- plausible$key[i]
    idx <- as.integer(strsplit(key, ",")[[1]])
    
    size[i] <- length(idx)
    
    these_vars <- var_names[idx]
    vars_str[i] <- paste(these_vars, collapse = " + ")
    
    # Stability of the variables in this model
    avg_stab[i] <- mean(stability_pi[these_vars])
  }
  
  plausible$size          <- size
  plausible$vars          <- vars_str
  plausible$avg_stability <- avg_stab
  
  # Filter by stability threshold τ
  plausible_final <- plausible[plausible$avg_stability >= tau, ]
  
  if (nrow(plausible_final) == 0) {
    cat("Warning: No models meet stability threshold τ =", tau, "\n")
    cat("Returning all models in AIC window.\n")
    plausible_final <- plausible
  }
  
  # Order by AIC
  plausible_final <- plausible_final[order(plausible_final$AIC), ]
  
  cat("\n=== Algorithm 3 Complete ===\n")
  cat("Selected", nrow(plausible_final), "plausible model(s).\n")
  
  plausible_final
}


ui <- fluidPage(
  titlePanel("mpssAIC – Multi-Path AIC Pipeline Demo"),
  
  tabsetPanel(
    tabPanel(
      "Linear Regression (Gaussian)",
      br(),
      h3("Plausible Models (Gaussian)"),
      tableOutput("linear_plausible"),
      br(),
      h3("Stability Scores"),
      tableOutput("linear_stability")
    ),
    
    tabPanel(
      "Logistic Regression (Binomial)",
      br(),
      h3("Plausible Models (Binomial)"),
      tableOutput("logistic_plausible"),
      br(),
      h3("Stability Scores"),
      tableOutput("logistic_stability"),
      br(),
      h3("Confusion Matrix (Best Plausible Model)"),
      tableOutput("logistic_cm"),
      br(),
      h3("Classification Metrics"),
      tableOutput("logistic_metrics"),
      br(),
      h3("Confusion Matrix Plot"),
      plotOutput("logistic_cm_plot", height = "300px")
    )
  )
)

server <- function(input, output, session) {
  
  ## --------------------------------------------------------------
  ## 1. Simulate synthetic data
  ##    (self-contained; no need for extra exported objects)
  ## --------------------------------------------------------------
  
  set.seed(123)
  
  # ---- Linear regression (Gaussian) ----
  n_lin <- 300
  p_lin <- 6
  X_linear <- matrix(rnorm(n_lin * p_lin), nrow = n_lin, ncol = p_lin)
  colnames(X_linear) <- paste0("x", 1:p_lin)
  
  beta_lin <- c(1.0, -1.5, 0.75, 0.5, 0, 0)
  y_linear <- as.numeric(X_linear %*% beta_lin + rnorm(n_lin, sd = 1))
  
  # ---- Logistic regression (Binomial) ----
  n_log <- 300
  p_log <- 6
  X_logistic <- matrix(rnorm(n_log * p_log), nrow = n_log, ncol = p_log)
  colnames(X_logistic) <- paste0("x", 1:p_log)
  
  beta_log <- c(1.0, -1.5, 0.75, 0, 0, 0)
  eta <- X_logistic %*% beta_log
  prob <- plogis(eta)
  y_logistic <- rbinom(n_log, size = 1, prob = prob)
  
  ## --------------------------------------------------------------
  ## 2. Run multi-path AIC pipeline using mpssAIC functions
  ## --------------------------------------------------------------
  
  # Gaussian
  result_linear <- multi_path_AIC_pipeline(
    X             = X_linear,
    y             = y_linear,
    model_type    = "gaussian",
    K             = NULL,      # e.g. min(p, 10)
    eps           = 1e-6,
    delta         = 2,
    L             = 50,
    B             = 20,        # increase for more stable results
    resample_type = "bootstrap",
    m             = NULL,
    Delta         = 2,
    tau           = 0.6
  )
  
  # Binomial
  result_logistic <- multi_path_AIC_pipeline(
    X             = X_logistic,
    y             = y_logistic,
    model_type    = "binomial",
    K             = NULL,
    eps           = 1e-6,
    delta         = 2,
    L             = 50,
    B             = 20,
    resample_type = "bootstrap",
    m             = NULL,
    Delta         = 2,
    tau           = 0.6
  )
  
  ## --------------------------------------------------------------
  ## 3. Linear outputs
  ## --------------------------------------------------------------
  
  output$linear_plausible <- renderTable({
    if (is.null(result_linear$plausible) ||
        nrow(result_linear$plausible) == 0) {
      return(NULL)
    }
    tab <- result_linear$plausible[
      , c("AIC", "AIC_diff", "size", "vars", "avg_stability")
    ]
    rownames(tab) <- NULL
    tab
  })
  
  output$linear_stability <- renderTable({
    if (is.null(result_linear$stability_pi)) return(NULL)
    data.frame(
      Variable  = names(result_linear$stability_pi),
      Stability = round(as.numeric(result_linear$stability_pi), 3),
      row.names = NULL
    )
  })
  
  ## --------------------------------------------------------------
  ## 4. Logistic outputs (incl. confusion matrix)
  ## --------------------------------------------------------------
  
  output$logistic_plausible <- renderTable({
    if (is.null(result_logistic$plausible) ||
        nrow(result_logistic$plausible) == 0) {
      return(NULL)
    }
    tab <- result_logistic$plausible[
      , c("AIC", "AIC_diff", "size", "vars", "avg_stability")
    ]
    rownames(tab) <- NULL
    tab
  })
  
  output$logistic_stability <- renderTable({
    if (is.null(result_logistic$stability_pi)) return(NULL)
    data.frame(
      Variable  = names(result_logistic$stability_pi),
      Stability = round(as.numeric(result_logistic$stability_pi), 3),
      row.names = NULL
    )
  })
  
  # Confusion matrix & metrics for best plausible model
  cm_results <- reactive({
    if (is.null(result_logistic$plausible) ||
        nrow(result_logistic$plausible) == 0) {
      return(NULL)
    }
    
    best_key <- result_logistic$plausible$key[1]
    best_idx <- as.integer(strsplit(best_key, ",")[[1]])
    
    X_best <- X_logistic[, best_idx, drop = FALSE]
    df_logistic <- data.frame(y = y_logistic, X_best)
    
    fit_best <- glm(y ~ ., data = df_logistic, family = binomial())
    
    confusion_metrics(
      fit    = fit_best,
      y      = y_logistic,
      cutoff = 0.5
    )
  })
  
  output$logistic_cm <- renderTable({
    res <- cm_results()
    if (is.null(res)) return(NULL)
    res$confusion_matrix
  }, rownames = TRUE)
  
  output$logistic_metrics <- renderTable({
    res <- cm_results()
    if (is.null(res)) return(NULL)
    data.frame(
      Metric = names(res$metrics),
      Value  = round(as.numeric(res$metrics), 4),
      row.names = NULL
    )
  })
  
  output$logistic_cm_plot <- renderPlot({
    res <- cm_results()
    if (is.null(res)) return(NULL)
    cm <- res$confusion_matrix
    
    par(mar = c(5, 5, 4, 2))
    
    plot(1, 1, type = "n", xlim = c(0, 2), ylim = c(0, 2),
         xlab = "Predicted", ylab = "Actual",
         main = "Confusion Matrix (Best Plausible Model)",
         xaxt = "n", yaxt = "n", bty = "n")
    
    rect(c(0, 1), c(0, 0), c(1, 2), c(1, 1),
         col = c("lightgreen", "lightcoral"))
    rect(c(0, 1), c(1, 1), c(1, 2), c(2, 2),
         col = c("lightcoral", "lightgreen"))
    
    text(0.5, 0.5, paste0("TN\n", cm[1, 1]), cex = 1.5, font = 2)
    text(1.5, 0.5, paste0("FP\n", cm[1, 2]), cex = 1.5, font = 2)
    text(0.5, 1.5, paste0("FN\n", cm[2, 1]), cex = 1.5, font = 2)
    text(1.5, 1.5, paste0("TP\n", cm[2, 2]), cex = 1.5, font = 2)
    
    axis(1, at = c(0.5, 1.5), labels = c("0", "1"))
    axis(2, at = c(0.5, 1.5), labels = c("0", "1"))
  })
}

shinyApp(ui, server)
