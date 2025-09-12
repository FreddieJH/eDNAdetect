#' Fitted Generalised Linear Model
#'
#' A fitted model predicting the probability of detection based on the number of reads in given sample
#'
#' @format ## `species_mod`
#' A fiited GLM
#' \describe{
#'   \item{reads}{Number of reads of a given label}
#' \item{reads_clr}{Centered-log ratio transformation of reads}
#'   \item{sample_reads}{total number of sample reads}
#'   ...
#' }
"species_mod"
