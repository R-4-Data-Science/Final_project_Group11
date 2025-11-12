#' Resampling-based feature stability
#'
#' For each resample, runs `build_paths()` and computes, per feature j,
#' the proportion of models on that resample that include j.
#' Aggregates across B resamples to obtain average stability pi_j.
#'
#' @param X data.frame or matrix of predictors.
#' @param y response vector.
#' @param B number of resamples.
#' @param resample one of "bootstrap" or "subsample".
#' @param m for "subsample", the number of rows per resample (m < n).
#' @param build_args named list of arguments forwarded to `build_paths()`.
#' @return A list with:
#'   - pi: numeric vector (length p) of average stabilities in [0,1].
#'   - resample_stats: data.frame per-resample with proportions per feature (optional compact form).
#'   - meta: list with resampling settings.
#' @export
stability <- function(X, y, B = 50, resample = c("bootstrap","subsample"),
                      m = NULL, build_args = list()) {
  resample <- match.arg(resample)
  X <- as.data.frame(X, check.names = FALSE)
  n <- nrow(X)
  scope <- colnames(X)
  p <- length(scope)

  prop_mat <- matrix(NA_real_, nrow = B, ncol = p,
                     dimnames = list(paste0("b", seq_len(B)), scope))

  for (b in seq_len(B)) {
    if (resample == "bootstrap") {
      idx <- sample.int(n, replace = TRUE)
    } else {
      if (is.null(m)) stop("Provide m for subsample.")
      idx <- sample.int(n, size = m, replace = FALSE)
    }
    Xb <- X[idx, , drop = FALSE]
    yb <- y[idx]

    forest <- do.call(build_paths, c(list(X = Xb, y = yb), build_args))
    # gather unique models in this forest
    aic_tab <- forest$aic_by_model
    codes <- aic_tab$model_code
    models <- lapply(codes, .vars_in_code)
    # proportion per feature:
    # z_j^(b) = (# models containing j) / (# models)
    M <- length(models)
    if (M == 0) {
      prop_mat[b, ] <- 0
    } else {
      counts <- setNames(integer(p), scope)
      for (S in models) {
        for (j in S) counts[[j]] <- counts[[j]] + 1L
      }
      prop_mat[b, ] <- as.numeric(counts) / M
    }
  }

  pi <- colMeans(prop_mat, na.rm = TRUE)
  list(
    pi = pi,
    resample_stats = prop_mat,
    meta = list(B = B, resample = resample, m = m, p = p, scope = scope, build_args = build_args)
  )
}
