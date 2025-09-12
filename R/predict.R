#' Predict Species Detection Probability
#'
#' Predicts the probability of species detection based on read counts using
#' a fitted GLM model. Returns point estimates with confidence intervals.
#'
#' @param reads_vec Numeric vector of read counts. Must contain non-negative
#'   values representing sequencing read counts for different samples or taxa.
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
#' reads <- c(1500, 2000, 0, 0, 750, 3000, 0)
#' predictions <- pred_detect(reads)
#'
#' # Access results
#' predictions$fit  # Point estimates
#' predictions$lwr  # Lower 95% CI bounds
#' predictions$upr  # Upper 95% CI bounds
#' }
#'
#' @seealso \code{\link{get_spp_model}} for the underlying model
#'
#' @importFrom compositions clr
#' @export
pred_detect <- function(reads_vec) {
  n_sample_reads <- rep(sum(reads_vec), times = length(reads_vec))
  reads_clr <- as.numeric(compositions::clr(reads_vec))

  preds_link_spp <- predict(
    eDNAdetect::species_mod,
    newdata = list(
      reads_clr = reads_clr,
      sample_reads = n_sample_reads
    ),
    se.fit = TRUE
  )

  link_f_spp <- family(eDNAdetect::species_mod)$linkinv

  return(list(
    fit = link_f_spp(preds_link_spp$fit),
    lwr = link_f_spp(preds_link_spp$fit - 1.96 * preds_link_spp$se.fit),
    upr = link_f_spp(preds_link_spp$fit + 1.96 * preds_link_spp$se.fit)
  ))
}
