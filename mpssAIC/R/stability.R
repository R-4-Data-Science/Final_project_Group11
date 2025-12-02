#' Internal: compute variable stability via resampling
#'
#' @inheritParams build_paths
#' @param B Number of resamples
#' @param resample_type "bootstrap" or "subsample"
#' @param m Subsample size (if using subsample)
#'
#' @return Numeric vector of length p with stability scores between 0 and 1.
#' @keywords internal
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
  model_type    <- match.arg(model_type)
  resample_type <- match.arg(resample_type)

  X <- as.data.frame(X)
  n <- nrow(X)
  p <- ncol(X)
  var_names <- colnames(X)

  if (is.null(m)) m <- ceiling(sqrt(n))

  Z <- matrix(0, nrow = B, ncol = p)
  colnames(Z) <- var_names

  for (b in seq_len(B)) {
    if (resample_type == "bootstrap") {
      idx <- sample(seq_len(n), size = n, replace = TRUE)
    } else {
      idx <- sample(seq_len(n), size = m, replace = FALSE)
    }

    Xb <- X[idx, , drop = FALSE]
    yb <- y[idx]

    mp_b <- multi_path_forward(
      X = Xb, y = yb,
      model_type = model_type,
      K = K, eps = eps, delta = delta, L = L
    )

    all_models_b <- unlist(mp_b$step_models, recursive = FALSE)
    if (!length(all_models_b)) next

    M <- length(all_models_b)

    for (j in seq_len(p)) {
      count_j <- sum(vapply(all_models_b, function(idx) j %in% idx, logical(1)))
      Z[b, j] <- count_j / M
    }
  }

  pi <- colMeans(Z)
  pi
}

#' Resampling-based variable stability
#'
#' @inheritParams build_paths
#' @param B Number of resamples
#' @param resample "bootstrap" or "subsample"
#' @param m Subsample size if using subsample
#'
#' @return Named numeric vector of stability scores between 0 and 1.
#' @export
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
  family   <- match.arg(family)
  resample <- match.arg(resample)

  compute_stability(
    X = X, y = y,
    model_type    = family,
    B             = B,
    resample_type = resample,
    m             = m,
    K = K, eps = eps, delta = delta, L = L
  )
}
