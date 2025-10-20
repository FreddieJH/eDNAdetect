#' Plot model fit with confidence intervals and detection zones
#'
#' @param reads Numeric vector of read counts
#' @param labels Character string vector of the label names
#' @param model_name Model to use - either DORY or DLOOP model
#' @return ggplot2 object
#' @export
plot_model <- function(reads, labels = NULL, model_name = "DLOOP") {

  if (!is.null(labels)) {
    if(length(labels) != length(reads)){
      stop("The length of the 'reads' vector must be the same length as the 'labels' vector")
    }
  }

  model_preds <- predict_detection(reads, model_name)


  detection_conf_cols = c(
    "No detection" = "grey70",
    "Low confidence" = "#CA0020",
    "Medium confidence" = "#F4A582",
    "High confidence" = "#92C5DE"
  )

  confidence_thresholds <- c(0.25, 0.5, 0.75)

  p <-
    dplyr::tibble(
      id_num = 1:length(reads),
      id = if (!is.null(labels)) labels else paste0("Label ", seq_along(reads)),
      fit = model_preds$fit,
      lwr = model_preds$lwr,
      upr = model_preds$upr
    ) |>
    dplyr::mutate(
      id = forcats::fct_reorder(id, -id_num)
    ) |>
    ggplot2::ggplot() +
    ggplot2::geom_rect(
      xmin = 0,
      xmax = confidence_thresholds[1],
      ymin = -Inf,
      ymax = Inf,
      fill = detection_conf_cols[1],
      alpha = 0.5
    ) +
    ggplot2::geom_rect(
      xmin = confidence_thresholds[1],
      xmax = confidence_thresholds[2],
      ymin = -Inf,
      ymax = Inf,
      fill = detection_conf_cols[2],
      alpha = 0.5
    ) +
    ggplot2::geom_rect(
      xmin = confidence_thresholds[2],
      xmax = confidence_thresholds[3],
      ymin = -Inf,
      ymax = Inf,
      fill = detection_conf_cols[3],
      alpha = 0.5
    ) +
    ggplot2::geom_rect(
      xmin = confidence_thresholds[3],
      xmax = 1,
      ymin = -Inf,
      ymax = Inf,
      fill = detection_conf_cols[4],
      alpha = 0.5
    ) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(y = id, xmin = lwr, xmax = upr),
      width = 0.1
    ) +
    ggplot2::geom_point(ggplot2::aes(x = fit, y = id), size = 4) +
    ggplot2::scale_x_continuous(
      label = scales::label_percent(),
      limits = c(0, 1)
    ) +
    ggplot2::labs(x = "Pr(prescence)", y = "Label ID") +
    ggplot2::theme_classic(20)
  return(p)
}


# basic_reads <- c(1, 1000, 100, 10000)

# p_basic <- plot_model(
#   reads = basic_reads,
#   model_name = "DLOOP"
# ) +
#   ggplot2::ggtitle("Basic scenario")

# p_basic_m10 <- plot_model(
#   reads = basic_reads * 10,
#   model_name = "DLOOP"
# ) +
#   ggplot2::ggtitle("All samples x10")
# patchwork::wrap_plots(p_basic, p_basic_m10, ncol = 1)

# p_basic_m100 <- plot_model(
#   reads = basic_reads * 100,
#   model_name = "DLOOP"
# ) +
#   ggplot2::ggtitle("All samples x100")
# patchwork::wrap_plots(p_basic, p_basic_m100, ncol = 1)

# p_basic_m1000 <- plot_model(
#   reads = basic_reads * 1000,
#   model_name = "DLOOP"
# ) +
#   ggplot2::ggtitle("All samples x1000")
# patchwork::wrap_plots(p_basic, p_basic_m1000, ncol = 1)

# p_basic_addpos <- plot_model(
#   reads = c(basic_reads, 100, 1000, 100),
#   model_name = "DLOOP"
# ) +
#   ggplot2::ggtitle("Basic scenario + additional postive labels")
# patchwork::wrap_plots(p_basic, p_basic_addpos, ncol = 1)

# p_basic_addzeros <- plot_model(
#   reads = c(basic_reads, 0, 0, 0, 0),
#   model_name = "DLOOP"
# ) +
#   ggplot2::ggtitle("Basic scenario + additional zero labels")
# patchwork::wrap_plots(p_basic, p_basic_addzeros, ncol = 1)

# p_basic_addextreme <- plot_model(
#   reads = c(basic_reads, 10000000),
#   model_name = "DLOOP"
# ) +
#   ggplot2::ggtitle("Basic scenario + additional extreme read")
# patchwork::wrap_plots(p_basic, p_basic_addextreme, ncol = 1)

# p_basic_addmanyextreme <- plot_model(
#   reads = c(basic_reads, 10000000, 10000000, 10000000),
#   model_name = "DLOOP"
# ) +
#   ggplot2::ggtitle("Basic scenario + additional multiple extreme reads")
# patchwork::wrap_plots(p_basic, p_basic_addmanyextreme, ncol = 1)

# # DLOOP
# p_basic <- plot_model(
#   reads = basic_reads,
#   model_name = "DLOOP"
# ) +
#   ggplot2::ggtitle("Basic scenario")

# p_basic_m10 <- plot_model(
#   reads = basic_reads * 10,
#   mod = dloop
# ) +
#   ggplot2::ggtitle("All samples x10")

# p_basic_m100 <- plot_model(
#   reads = basic_reads * 100,
#   mod = dloop_mod
# ) +
#   ggplot2::ggtitle("All samples x100")

# p_basic_m1000 <- plot_model(
#   reads = basic_reads * 1000,
#   mod = dloop_mod
# ) +
#   ggplot2::ggtitle("All samples x1000")

# p_basic_addpos <- plot_model(
#   reads = c(basic_reads, 100, 1000, 100),
#   mod = dloop_mod
# ) +
#   ggplot2::ggtitle("Basic scenario + additional postive labels")
# patchwork::wrap_plots(p_basic, p_basic_addpos, ncol = 1)

# p_basic_addzeros <- plot_model(
#   reads = c(basic_reads, 0, 0, 0, 0),
#   mod = dloop_mod
# ) +
#   ggplot2::ggtitle("Basic scenario + additional zero labels")
# patchwork::wrap_plots(p_basic, p_basic_addzeros, ncol = 1)

# p_basic_addextreme <- plot_model(
#   reads = c(basic_reads, 10000000),
#   mod = dloop_mod
# ) +
#   ggplot2::ggtitle("Basic scenario + additional extreme read")
# patchwork::wrap_plots(p_basic, p_basic_addextreme, ncol = 1)

# p_basic_addmanyextreme <- plot_model(
#   reads = c(basic_reads, 10000000, 10000000, 10000000),
#   mod = dloop_mod
# ) +
#   ggplot2::ggtitle("Basic scenario + additional multiple extreme reads")
# patchwork::wrap_plots(p_basic, p_basic_addmanyextreme, ncol = 1)

# plot_model(
#   reads = c(0, 0, 0, 100000, 0, 0),
#   mod = dloop_mod
# ) +
#   ggplot2::ggtitle("extreme reads")
#   log(c(0, 0, 0, 100000, 0, 0) + 25)
