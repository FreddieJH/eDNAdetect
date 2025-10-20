#' Predict Species Detection Probability
#'
#' Predicts the probability of species detection based on read counts using
#' a fitted GLM model. Returns point estimates with confidence intervals.
#'
#' @param reads_vec Numeric vector of read counts. Must contain non-negative
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
#' predictions <- pred_detect(reads)
#'
#' # Access results
#' predictions$fit  # Point estimates
#' predictions$lwr  # Lower 95% CI bounds
#' predictions$upr  # Upper 95% CI bounds
#' }
#'
#' @importFrom compositions clr
#' @export
predict_detection <- function(
  reads_vec,
  model_name = c("DORY", "DLOOP"),
  fit_type = c("fit", "upr", "lwr"), 
  pseudocount = 25
) {
  if (is.null(mod)) {
    stop(
      "Please provide the fitted model object"
    )
  }

  if (length(fit_type) > 1) {
    stop(
      "Please specify your type of fit (Options: model estimate (fit), upper 95% confidence interval (lwr), lower 95% confidence interval (upr))"
    )
  }
  if(model_name == "DORY"){
      mod <- data("dory_mod", package = "eDNAdetect", envir = environment())
  } else if (model_name == "DLOOP") {
      mod <- data("dloop_mod", package = "eDNAdetect", envir = environment())
  }


  reads_log <- log(reads_vec + pseudocount)

  preds_link <- stats::predict(
    mod,
    newdata = list(
      reads_clr = reads_log
    ),
    se.fit = TRUE
  )

  link_f_spp <- stats::family(mod)$linkinv

  if (fit_type == "fit") {
    out <- link_f_spp(preds_link$fit)
  } else if (fit_type == "lwr") {
    out <- link_f_spp(preds_link$fit - 1.96 * preds_link$se.fit)
  } else if (fit_type == "upr") {
    out <- link_f_spp(preds_link$fit + 1.96 * preds_link$se.fit)
  }
  return(out)
}
