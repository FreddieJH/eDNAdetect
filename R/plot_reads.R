#' Create Stacked Bar Plots from Wide-Format Data
#'
#' Converts wide-format data to long format and creates horizontal stacked bar
#' plots showing either absolute values or proportions. Useful for visualising
#' compositional data such as species abundance, sequencing reads, or survey responses.
#'
#' @param data A data frame in wide format where the first column contains
#'   sample identifiers and subsequent columns contain numeric values to plot.
#' @param proportions Logical. If `TRUE` (default), creates a proportional
#'   stacked bar plot. If `FALSE`, creates an absolute value plot.
#' @param sample_col Character. Name of the sample column. If `NULL` (default),
#'   uses the first column name.
#'
#' @return A ggplot2 object showing horizontal stacked bars.
#'
#' @details The function automatically pivots wide data to long format,
#'   calculates sample totals and proportions, then creates a stacked bar plot.
#'   When `proportions = TRUE`, bars sum to 1.0 for each sample. When
#'   `proportions = FALSE`, bars show absolute values.
#'
#' @examples
#' # Create example data
#' example_data <- data.frame(
#'   sample = c("Site_A", "Site_B", "Site_C"),
#'   species_1 = c(150, 200, 100),
#'   species_2 = c(300, 150, 250),
#'   species_3 = c(50, 100, 75)
#' )
#'
#' # Proportional plot (default)
#' plot_stacked_bars(example_data)
#'
#' # Absolute values plot
#' plot_stacked_bars(example_data, proportions = FALSE)
#'
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @export
plot_reads <- function(data, proportions = TRUE, sample_col = NULL) {
  # Get sample column name
  if (is.null(sample_col)) {
    sample_col <- colnames(data)[1]
  }

  # Pivot to long format and calculate proportions
  plot_data <- data |>
    tidyr::pivot_longer(
      -!!sample_col,
      names_to = "label",
      values_to = "reads"
    ) |>
    dplyr::add_count(
      !!rlang::sym(sample_col),
      wt = reads,
      name = "sample_reads"
    ) |>
    dplyr::mutate(reads_prop = reads / sample_reads)

  # Create base plot
  p <- plot_data |>
    ggplot2::ggplot() +
    ggplot2::aes(y = !!rlang::sym(sample_col), fill = label) +
    ggplot2::geom_col()

  # Add x aesthetic based on proportions parameter
  if (proportions) {
    p <- p + ggplot2::aes(x = reads_prop)
  } else {
    p <- p + ggplot2::aes(x = reads)
  }

  return(p)
}
