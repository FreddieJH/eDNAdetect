#' Predict Species Detection Probability
#'
#' Predicts the probability of species detection based on read counts using
#' a fitted GLM model. Returns point estimates with confidence intervals.
#'
#' @param reads Numeric vector of read counts. Must contain non-negative
#'   values representing sequencing read counts for different samples or taxa.
#' @param model_name Model to use - either DORY or DLOOP model
#'
#' @return A list containing:
#'   \describe{
#'     \item{fit}{Numeric vector of predicted detection probabilities}
#'     \item{lwr}{Numeric vector of lower 95% confidence interval bounds}
#'     \item{upr}{Numeric vector of upper 95% confidence interval bounds}
#'   }
#'
#' @details The function applies a centred log-ratio (CLR) transformation to
#'   the read counts and uses total sample reads as a predictor. Confidence
#'   intervals are calculated using ±1.96 standard errors on the link scale
#'   before back-transformation.
#'
#' @examples
#' \dontrun{
#' reads <- c(1500, 2000, 0, 0, 750, 3000, 0, 8000)
#' predictions <- predict_detection(reads)
#'
#' # Access results
#' predictions$fit  # Point estimates
#' predictions$lwr  # Lower 95% CI bounds
#' predictions$upr  # Upper 95% CI bounds
#' }
#'
#' @export
predict_detection <- function(
  reads,
  model_name = "DLOOP",
  ci = 0.8,
  pseudocount = 25
) {
  if (is.null(model_name)) {
    stop(
      "Please provide the model to be used; either 'model_name = DORY' or 'model_name = DLOOP'"
    )
  }

  if (class(reads) != "numeric") {
    stop(
      "The 'reads' argument must be a numeric value or vector."
    )
  }

  if (!is.numeric(ci) || ci < 0 || ci > 1) {
    stop("Credible Interval (ci) must be a numeric value between 0 and 1.")
  }

  if (model_name == "DORY") {
    mod <- dory_brms
  } else if (model_name == "DLOOP") {
    mod <- dloop_brms
  }

  ci_lwr <- (1 - ci) / 2
  ci_upr <- 1 - ci_lwr

  predication_table <-
    epred_draws(
      dory_brms,
      newdata = data.frame(reads_log = log(reads + pseudocount)),
      re_formula = NA
    ) |>
    summarise(
      fit = quantile(.epred, 0.5),
      lwr = quantile(.epred, ci_lwr),
      upr = quantile(.epred, ci_upr),
      .groups = "drop"
    )

  out <- list(
    fit = predication_table$fit,
    lwr = predication_table$lwr,
    upr = predication_table$upr
  )

  return(out)
}
