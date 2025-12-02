#' Internal: select plausible, stable models
#'
#' @param mp_full Output from \code{build_paths()}
#' @param stability_pi Named numeric stability vector
#' @param Delta AIC window
#' @param tau Stability threshold
#'
#' @return Data frame of plausible models
#' @keywords internal
select_plausible_models <- function(
    mp_full,
    stability_pi,
    Delta = 2,
    tau = 0.6
) {
  var_names   <- mp_full$var_names
  step_models <- mp_full$step_models
  step_AICs   <- mp_full$step_AICs

  model_key <- function(idx) {
    if (length(idx) == 0) return("")
    paste(sort(idx), collapse = ",")
  }

  all_models <- unlist(step_models, recursive = FALSE)
  all_AICs   <- unlist(step_AICs)

  if (!length(all_models)) {
    return(data.frame())
  }

  keys <- vapply(all_models, model_key, character(1))
  df <- data.frame(
    key = keys,
    AIC = all_AICs,
    stringsAsFactors = FALSE
  )

  agg <- aggregate(AIC ~ key, data = df, FUN = min)
  agg <- agg[agg$key != "", ]

  if (!nrow(agg)) {
    return(data.frame())
  }

  best_AIC      <- min(agg$AIC)
  agg$AIC_diff  <- agg$AIC - best_AIC
  plausible     <- agg[agg$AIC_diff <= Delta, ]

  if (!nrow(plausible)) {
    return(data.frame())
  }

  size      <- integer(nrow(plausible))
  avg_stab  <- numeric(nrow(plausible))
  vars_str  <- character(nrow(plausible))

  for (i in seq_len(nrow(plausible))) {
    key <- plausible$key[i]
    idx <- as.integer(strsplit(key, ",")[[1]])
    size[i] <- length(idx)

    these_vars <- var_names[idx]
    vars_str[i] <- paste(these_vars, collapse = " + ")
    avg_stab[i] <- mean(stability_pi[these_vars])
  }

  plausible$size          <- size
  plausible$vars          <- vars_str
  plausible$avg_stability <- avg_stab

  plausible_final <- plausible[plausible$avg_stability >= tau, ]
  if (!nrow(plausible_final)) {
    plausible_final <- plausible
  }

  plausible_final[order(plausible_final$AIC), ]
}

#' Plausible model selection using AIC and stability
#'
#' @param forest Object returned by \code{build_paths()}
#' @param pi Named numeric vector of stability scores
#' @param Delta AIC tolerance window
#' @param tau Minimum average stability
#'
#' @return Data frame of plausible models
#' @export
plausible_models <- function(
    forest,
    pi,
    Delta = 2,
    tau = 0.6
) {
  select_plausible_models(
    mp_full      = forest,
    stability_pi = pi,
    Delta        = Delta,
    tau          = tau
  )
}
