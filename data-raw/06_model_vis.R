log_reads_predication_table <- tidyr::expand_grid(
  reads_log = seq(3, 13, by = 0.01)
)

boot_function <- function(model) {
  predict(
    model,
    newdata = log_reads_predication_table,
    re.form = NA
  )
}

dory_boot_results <- lme4::bootMer(
  dory_mod,
  FUN = boot_function,
  nsim = 100,
  parallel = "multicore",
  ncpus = parallel::detectCores() - 1
)

dloop_boot_results <- lme4::bootMer(
  dloop_mod,
  FUN = boot_function,
  nsim = 100,
  parallel = "multicore",
  ncpus = parallel::detectCores() - 1
)

dory_prediction_table <- log_reads_predication_table |>
  dplyr::mutate(
    mod_fit = predict(
      dory_mod,
      newdata = log_reads_predication_table,
      re.form = NA
    ),
    mod_lwr = apply(boot_results$t, 2, quantile, probs = 0.025),
    mod_upr = apply(boot_results$t, 2, quantile, probs = 0.975)
  ) |>
  dplyr::mutate(
    fit = family(dory_mod)$linkinv(mod_fit),
    lwr = family(dory_mod)$linkinv(mod_lwr),
    upr = family(dory_mod)$linkinv(mod_upr)
  )

dloop_prediction_table <- log_reads_predication_table |>
  dplyr::mutate(
    mod_fit = predict(
      dloop_mod,
      newdata = log_reads_predication_table,
      re.form = NA
    ),
    mod_lwr = apply(boot_results$t, 2, quantile, probs = 0.025),
    mod_upr = apply(boot_results$t, 2, quantile, probs = 0.975)
  ) |>
  dplyr::mutate(
    fit = family(dloop_mod)$linkinv(mod_fit),
    lwr = family(dloop_mod)$linkinv(mod_lwr),
    upr = family(dloop_mod)$linkinv(mod_upr)
  )
# dory_prediction_table <-
#   tidyr::expand_grid(
#     reads_log = seq(
#       3,
#       13,
#       by = 0.01
#     )
#   ) |>
#   dplyr::mutate(
#     mod_fit = as.numeric(
#       predict(
#         dory_mod,
#         newdata = data.frame(
#           reads_log = reads_log
#         ),
#         se.fit = TRUE,
#         re.form = NA
#       )$fit
#     ),
#     mod_se = as.numeric(
#       predict(
#         dory_mod,
#         newdata = data.frame(
#           reads_log = reads_log
#         ),
#         se.fit = TRUE,
#         re.form = NA
#       )$se.fit
#     )
#   ) |>
#   dplyr::mutate(
#     fit = family(dory_mod)$linkinv(mod_fit),
#     lwr = family(dory_mod)$linkinv(mod_fit - 1.96 * mod_se),
#     upr = family(dory_mod)$linkinv(mod_fit + 1.96 * mod_se)
#   )

# dloop_prediction_table <-
#   tidyr::expand_grid(
#     reads_log = seq(
#       3,
#       13,
#       by = 0.01
#     )
#   ) |>
#   dplyr::mutate(
#     mod_fit = as.numeric(
#       stats::predict(
#         dloop_mod,
#         newdata = data.frame(
#           reads_log = reads_log
#         ),
#         se.fit = TRUE,
#         re.form = NA
#       )$fit
#     ),
#     mod_se = as.numeric(
#       predict(
#         dloop_mod,
#         newdata = data.frame(
#           reads_log = reads_log
#         ),
#         se.fit = TRUE,
#         re.form = NA
#       )$se.fit
#     )
#   ) |>
#   dplyr::mutate(
#     fit = family(dloop_mod)$linkinv(mod_fit),
#     lwr = family(dloop_mod)$linkinv(mod_fit - 1.96 * mod_se),
#     upr = family(dloop_mod)$linkinv(mod_fit + 1.96 * mod_se)
#   )

# dloop_prediction_table |>
#   # dory_prediction_table |>
#   dplyr::mutate(assay_reads_lab = purrr::map_chr(total_assay_reads, num2str)) |>
#   dplyr::mutate(
#     assay_reads_lab = forcats::fct_reorder(assay_reads_lab, total_assay_reads)
#   ) |>
#   ggplot2::ggplot() +
#   ggplot2::geom_point(
#     ggplot2::aes(
#       x = reads_prop,
#       y = possible_detection,
#       # size = log10(total_assay_reads),
#     ),
#     size = 4,
#     data = mod_dat_dloop,
#     pch = 21,
#     alpha = 0.3,
#   ) +
#   ggplot2::geom_ribbon(
#     ggplot2::aes(
#       x = reads_prop,
#       ymin = lwr,
#       ymax = upr,
#       fill = assay_reads_lab
#     ),
#     alpha = 0.2
#   ) +
#   ggplot2::geom_line(ggplot2::aes(
#     x = reads_prop,
#     y = fit,
#     col = assay_reads_lab
#   )) +
#   ggplot2::scale_y_continuous(label = scales::label_percent()) +
#   ggplot2::scale_x_continuous(label = scales::label_percent()) +
#   # ggplot2::scale_x_log10() +
#   ggplot2::labs(
#     x = "OTU proportion in sample",
#     y = "Pr(Presence)",
#     fill = "Sample reads",
#     col = "Sample reads"
#   ) +
#   ggplot2::theme_classic(20) +
#   ggplot2::theme(
#     legend.position = c(1, 0),
#     legend.justification = c(1, 0),
#     legend.background = ggplot2::element_rect(fill = "transparent")
#   )

# # CLR VERSION
# dloop_prediction_table <-
#   tidyr::expand_grid(
#     reads_clr = seq(
#       min(mod_dat_dloop$reads_clr),
#       max(mod_dat_dloop$reads_clr),
#       length.out = 10000
#     ),
#     sample_reads = c(1e3, 1e4, 1e5, 1e6)
#   ) |>
#   dplyr::mutate(
#     mod_fit = as.numeric(
#       stats::predict(
#         dloop_mod,
#         newdata = list(
#           reads_clr = reads_clr,
#           sample_reads = sample_reads
#         ),
#         se.fit = TRUE
#       )$fit
#     ),
#     mod_se = as.numeric(
#       predict(
#         dloop_mod,
#         newdata = list(
#           reads_clr = reads_clr,
#           sample_reads = sample_reads
#         ),
#         se.fit = TRUE
#       )$se.fit
#     )
#   ) |>
#   dplyr::mutate(
#     fit = family(dloop_mod)$linkinv(mod_fit),
#     lwr = family(dloop_mod)$linkinv(mod_fit - 1.96 * mod_se),
#     upr = family(dloop_mod)$linkinv(mod_fit + 1.96 * mod_se)
#   )

dory_mod_plot <-
  dory_prediction_table |>
  ggplot2::ggplot() +
  ggplot2::geom_point(
    ggplot2::aes(
      x = reads_log,
      y = possible_detection,
      # size = log10(total_assay_reads),
    ),
    size = 4,
    data = mod_dat_dloop |>
      filter((!str_detect(species, "Unidentified"))),
    pch = 21,
    alpha = 0.3,
  ) +
  ggplot2::geom_ribbon(
    ggplot2::aes(
      x = reads_log,
      ymin = lwr,
      ymax = upr
    ),
    alpha = 0.2
  ) +
  ggplot2::geom_line(ggplot2::aes(
    x = reads_log,
    y = fit
  )) +
  ggplot2::scale_y_continuous(label = scales::label_percent()) +
  # ggplot2::scale_x_log10() +
  ggplot2::labs(
    x = expression(log(OTU ~ reads + delta)),
    y = "Pr(Presence)"
  ) +
  ggplot2::theme_classic(20) +
  ggplot2::theme(
    legend.position = c(1, 0),
    legend.justification = c(1, 0),
    legend.background = ggplot2::element_rect(fill = "transparent")
  )

dloop_mod_plot <-
  dloop_prediction_table |>
  ggplot2::ggplot() +
  ggplot2::geom_point(
    ggplot2::aes(
      x = reads_log,
      y = possible_detection,
      # size = log10(total_assay_reads),
    ),
    size = 4,
    data = mod_dat_dory |>
      filter(
        (!str_detect(species, "Unidentified")) |
          species == "Unidentified Thunnus (genus)"
      ) |>
      filter(!str_detect(species, "Thunnus albacares")),
    pch = 21,
    alpha = 0.3,
  ) +
  ggplot2::geom_ribbon(
    ggplot2::aes(
      x = reads_log,
      ymin = lwr,
      ymax = upr
    ),
    alpha = 0.2
  ) +
  ggplot2::geom_line(ggplot2::aes(
    x = reads_log,
    y = fit
  )) +
  ggplot2::scale_y_continuous(label = scales::label_percent()) +
  # ggplot2::scale_x_log10() +
  ggplot2::labs(
    x = expression(log(OTU ~ reads + delta)),
    y = "Pr(Presence)"
  ) +
  ggplot2::theme_classic(20) +
  ggplot2::theme(
    legend.position = c(1, 0),
    legend.justification = c(1, 0),
    legend.background = ggplot2::element_rect(fill = "transparent")
  )

ggsave(
  filename = "data-raw/model-fit/dory.png",
  plot = dory_mod_plot,
  height = 10,
  width = 16
)

ggsave(
  filename = "data-raw/model-fit/dloop.png",
  plot = dloop_mod_plot,
  height = 10,
  width = 16
)
