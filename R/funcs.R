' Convert Numbers to Human-Readable Strings
#'
#' Converts numeric values to abbreviated string format using "K" (thousands)
#' and "M" (millions) suffixes for improved readability in plots and tables.
#'
#' @param number A numeric value to convert. Can handle `NA` values.
#'
#' @return A character string with the formatted number. Numbers < 1,000 are
#'   returned as-is. Numbers ≥ 1,000 and < 1,000,000 are formatted with "K"
#'   suffix (e.g., "1.5K"). Numbers ≥ 1,000,000 and < 1,000,000,000 are
#'   formatted with "M" suffix (e.g., "2.3M"). Larger numbers return as
#'   character without abbreviation. `NA` input returns `NA`.
#'
#' @examples
#' num2str(500)        # "500"
#' num2str(1500)       # "1.5K"
#' num2str(2500000)    # "2.5M"
#' num2str(NA)         # NA
#' num2str(-1500)      # "-1.5K"
#'
#' @export
num2str <- function(number) {
  if (is.na(number)) {
    return(NA)
  } else if (abs(number) < 1e3) {
    return(as.character(number))
  } else if (abs(number) >= 1e3 & abs(number) < 1e6) {
    return(paste0(round(number / 1e3, 1), "K"))
  } else if (abs(number) >= 1e6 & abs(number) < 1e9) {
    return(paste0(round(number / 1e6, 1), "M"))
  } else {
    # You can extend this for B (billions) or T (trillions) if needed
    return(as.character(number)) # Fallback for very large numbers
  }
}
