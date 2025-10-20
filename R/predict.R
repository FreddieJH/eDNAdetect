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
  pseudocount = 25
) {
  if (is.null(model_name)) {
    stop(
      "Please provide the model to be used; either 'model_name = DORY' or 'model_name = DLOOP'"
    )
  }

  if(model_name == "DORY"){
      mod <- dory_mod
  } else if (model_name == "DLOOP") {
      mod <- dloop_mod
  }

  preds_link <- suppressWarnings(stats::predict(
    mod,
    newdata = data.frame(
      reads_log = log(reads + pseudocount)
    ),
    se.fit = TRUE,
    re.form = NA
  ))

  link_f_spp <- stats::family(mod)$linkinv


  out <- list(fit = link_f_spp(preds_link$fit), 
  lwr = link_f_spp(preds_link$fit - 1.96 * preds_link$se.fit), 
  upr = link_f_spp(preds_link$fit + 1.96 * preds_link$se.fit))

  return(out)

}
