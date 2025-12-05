# ============================================================================
# Wrapper Aliases for Assignment Compliance
# ============================================================================
# The assignment requires specific function names: build_paths(), stability(),
# and plausible_models(). Our implementation uses different names, so we
# create wrapper aliases here.
# ============================================================================

#' Multi-Path Forward Selection (Wrapper)
#'
#' @description
#' Wrapper for \code{\link{multi_path_forward}} to match assignment specification.
#' Builds multiple forward-selection paths by adding variables one at a time
#' and retaining all children within a ΔAIC window.
#'
#' @inheritParams multi_path_forward
#' @return Same as \code{\link{multi_path_forward}}
#' @export
#' @seealso \code{\link{multi_path_forward}}
#'
#' @examples
#' \dontrun{
#' # Gaussian example
#' X <- data.frame(x1 = rnorm(100), x2 = rnorm(100))
#' y <- 2*X$x1 - 1.5*X$x2 + rnorm(100)
#' forest <- build_paths(X, y, family = "gaussian", K = 5)
#' }
build_paths <- function(
    X, y,
    family = c("gaussian", "binomial"),
    K = NULL,
    eps = 1e-6,
    delta = 2,
    L = 50
) {
  # Call the actual implementation
  result <- multi_path_forward(
    X = X, 
    y = y,
    model_type = family,
    K = K,
    eps = eps,
    delta = delta,
    L = L
  )
  
  # Rename output to match expected format
  names(result)[names(result) == "step_models"] <- "frontiers"
  names(result)[names(result) == "step_AICs"] <- "aic_by_model"
  
  # Add meta information
  result$meta <- list(
    K = if (is.null(K)) min(ncol(X), 10) else K,
    eps = eps,
    delta = delta,
    L = L,
    family = match.arg(family)
  )
  
  # Set class for method dispatch
  class(result) <- c("path_forest", "list")
  
  result
}


#' Variable Stability via Resampling (Wrapper)
#'
#' @description
#' Wrapper for \code{\link{compute_stability}} to match assignment specification.
#' Computes stability scores via bootstrap or subsample resampling.
#'
#' @inheritParams compute_stability
#' @return Same as \code{\link{compute_stability}}
#' @export
#' @seealso \code{\link{compute_stability}}
#'
#' @examples
#' \dontrun{
#' # Gaussian example
#' X <- data.frame(x1 = rnorm(100), x2 = rnorm(100))
#' y <- 2*X$x1 - 1.5*X$x2 + rnorm(100)
#' pi_scores <- stability(X, y, family = "gaussian", B = 30)
#' }
stability <- function(
    X, y,
    family = c("gaussian", "binomial"),
    B = 50,
    resample = c("bootstrap", "subsample"),
    m = NULL,
    K = NULL,
    eps = 1e-6,
    delta = 2,
    L = 50
) {
  # Call the actual implementation
  pi_vec <- compute_stability(
    X = X,
    y = y,
    model_type = family,
    B = B,
    resample_type = resample,
    m = m,
    K = K,
    eps = eps,
    delta = delta,
    L = L
  )
  
  # Create structured output
  result <- list(
    pi = pi_vec,
    B = B,
    resample_type = match.arg(resample),
    m = m
  )
  
  class(result) <- c("path_stability", "list")
  
  result
}


#' Plausible Model Selection (Wrapper)
#'
#' @description
#' Wrapper for \code{\link{select_plausible_models}} to match assignment specification.
#' Selects models that are both high-quality (low AIC) and stable.
#'
#' @param forest Output from \code{\link{build_paths}}
#' @param pi Either output from \code{\link{stability}} or a numeric vector of stability scores
#' @param Delta AIC tolerance (default = 2)
#' @param tau Minimum average stability threshold (default = 0.6)
#'
#' @return Data frame of plausible models
#' @export
#' @seealso \code{\link{select_plausible_models}}
#'
#' @examples
#' \dontrun{
#' # Complete workflow
#' X <- data.frame(x1 = rnorm(100), x2 = rnorm(100))
#' y <- 2*X$x1 - 1.5*X$x2 + rnorm(100)
#' 
#' forest <- build_paths(X, y, family = "gaussian")
#' pi_scores <- stability(X, y, family = "gaussian", B = 30)
#' models <- plausible_models(forest, pi_scores, Delta = 2, tau = 0.6)
#' }
plausible_models <- function(
    forest,
    pi,
    Delta = 2,
    tau = 0.6
) {
  # Extract pi vector if input is path_stability object
  if (inherits(pi, "path_stability")) {
    pi_vec <- pi$pi
  } else {
    pi_vec <- pi
  }
  
  # Convert forest back to format expected by select_plausible_models
  mp_full <- list(
    var_names = forest$var_names,
    step_models = if ("frontiers" %in% names(forest)) forest$frontiers else forest$step_models,
    step_AICs = if ("aic_by_model" %in% names(forest)) forest$aic_by_model else forest$step_AICs
  )
  
  # Call the actual implementation
  result <- select_plausible_models(
    mp_full = mp_full,
    stability_pi = pi_vec,
    Delta = Delta,
    tau = tau
  )
  
  result
}