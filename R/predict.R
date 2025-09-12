#' Get Species Detection Model
#'
#' Loads the fitted GLM model for species detection from the package data.
#' The model is cached after first load for improved performance.
#'
#' @param reload Logical. If TRUE, forces reloading of the model from file.
#'   Default is FALSE.
#'
#' @return A fitted GLM object used for species detection predictions.
#'
#' @details The model file should be located in the package's extdata directory
#'   as 'model.rda'. The function will cache the model in memory after first
#'   load to improve performance on subsequent calls.
#'
#' @examples
#' \dontrun{
#' model <- get_spp_model()
#' }
#'
#' @export
get_spp_model <- function() {
  model_path <- system.file("extdata", "model.rda", package = "eDNAdetect")
  load(model_path)
  return(get(ls()[1]))
}

species_mod <- NULL
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
  species_mod <- get_spp_model()
  n_sample_reads <- rep(sum(reads_vec), times = length(reads_vec))
  reads_clr <- as.numeric(compositions::clr(reads_vec))

  preds_link_spp <- predict(
    species_mod,
    newdata = list(
      reads_clr = reads_clr,
      sample_reads = n_sample_reads
    ),
    se.fit = TRUE
  )

  link_f_spp <- family(species_mod)$linkinv

  return(list(
    fit = link_f_spp(preds_link_spp$fit),
    lwr = link_f_spp(preds_link_spp$fit - 1.96 * preds_link_spp$se.fit),
    upr = link_f_spp(preds_link_spp$fit + 1.96 * preds_link_spp$se.fit)
  ))
}
