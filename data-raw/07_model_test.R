source("data-raw/04_model_fit.R")

detection_conf_cols = c(
  "No detection" = "grey70",
  "Low confidence" = "#CA0020",
  "Medium confidence" = "#F4A582",
  "High confidence" = "#92C5DE"
)


mod_dat_dloop_withPreds <-
  modeldata_dloop |>
  dplyr::mutate(
    mod_fit = as.numeric(
      predict(
        dloop_mod,
        newdata = data.frame(
          reads_log = reads_log
        ),
        se.fit = TRUE,
        re.form = NA
      )$fit
    ),
    mod_se = as.numeric(
      predict(
        dloop_mod,
        newdata = data.frame(
          reads_log = reads_log
        ),
        se.fit = TRUE,
        re.form = NA
      )$se.fit
    )
  ) |>
  dplyr::mutate(
    fit = family(dloop_mod)$linkinv(mod_fit),
    lwr = family(dloop_mod)$linkinv(mod_fit - 1.96 * mod_se),
    upr = family(dloop_mod)$linkinv(mod_fit + 1.96 * mod_se)
  )


mod_dat_dory_withPreds <-
  modeldata_dory |>
  dplyr::mutate(
    mod_fit = as.numeric(
      predict(
        dory_mod,
        newdata = data.frame(
          reads_log = reads_log
        ),
        se.fit = TRUE,
        re.form = NA
      )$fit
    ),
    mod_se = as.numeric(
      predict(
        dory_mod,
        newdata = data.frame(
          reads_log = reads_log
        ),
        se.fit = TRUE,
        re.form = NA
      )$se.fit
    )
  ) |>
  dplyr::mutate(
    fit = family(dory_mod)$linkinv(mod_fit),
    lwr = family(dory_mod)$linkinv(mod_fit - 1.96 * mod_se),
    upr = family(dory_mod)$linkinv(mod_fit + 1.96 * mod_se)
  )


dloop_false_negatives <-
  mod_dat_dloop_withPreds |>
  filter(!str_detect(species, "Unidentified")) |>
  filter(possible_detection == 1, lwr < .25) |>
  pull(sample) |>
  unique()

dloop_false_positives <-
  mod_dat_dloop_withPreds |>
  filter(!str_detect(species, "Unidentified")) |>
  filter(possible_detection == 0, lwr > .25) |>
  pull(sample) |>
  unique()

dory_false_negatives <-
  mod_dat_dory_withPreds |>
  filter(!str_detect(species, "Unidentified")) |>
  filter(possible_detection == 1, lwr < .25) |>
  pull(sample) |>
  unique()

dory_false_positives <-
  mod_dat_dory_withPreds |>
  filter(!str_detect(species, "Unidentified")) |>
  filter(possible_detection == 0, lwr > .25) |>
  pull(sample) |>
  unique()


for (s in unique(mod_dat_dloop_withPreds$sample)) {
  p <-
    mod_dat_dloop_withPreds |>
    dplyr::filter(sample == s) |>
    dplyr::mutate(species = forcats::fct_reorder(species, reads)) |>
    dplyr::mutate(
      detection_col = ifelse(
        possible_detection == 1,
        "blue",
        "grey70"
      )
    ) |>
    dplyr::mutate(
      name = glue::glue("<p style='color:{detection_col}'>{species}</p>"),
      name = forcats::fct_reorder(name, reads_prop)
    ) |>
    dplyr::mutate(
      fit_lab = ifelse(
        fit > 0.001,
        paste0(
          signif(fit * 100, 2),
          "% (",
          signif(lwr * 100, 2),
          "-",
          signif(upr * 100, 2),
          "%)"
        ),
        NA
      )
    ) |>
    dplyr::mutate(
      conf = dplyr::case_when(
        lwr > 0.75 ~ detection_conf_cols[4],
        lwr > 0.50 ~ detection_conf_cols[3],
        lwr > 0.25 ~ detection_conf_cols[2],
        .default = detection_conf_cols[1]
      )
    ) |>
    ggplot2::ggplot() +
    ggplot2::aes(x = reads, y = name, fill = conf) +
    ggplot2::geom_col() +
    ggplot2::geom_text(ggplot2::aes(label = fit_lab), hjust = 0, size = 8) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.01, 0.5))
    ) +
    ggplot2::theme_classic(20) +
    ggplot2::theme(
      axis.text.y = ggtext::element_markdown(),
      axis.title.y = ggplot2::element_blank()
    )

  ggplot2::ggsave(
    p,
    filename = paste0("data-raw/model_ests/", s, ".png"),
    height = 10,
    width = 10
  )
}
plot_model_test <- function(
  sample_vector,
  type = c("dloop", "dory"),
  error_type = "false_negative"
) {
  get(paste0("mod_dat_", type, "_withPreds")) |>
    dplyr::filter(sample == i) |>
    dplyr::mutate(species = forcats::fct_reorder(species, reads)) |>
    dplyr::mutate(
      detection_col = ifelse(
        possible_detection == 1,
        "blue",
        "grey70"
      )
    ) |>
    dplyr::mutate(
      name = glue::glue("<p style='color:{detection_col}'>{species}</p>"),
      name = forcats::fct_reorder(name, reads_prop)
    ) |>
    dplyr::mutate(
      fit_lab = ifelse(
        fit > 0.001,
        paste0(
          signif(fit * 100, 2),
          "% (",
          signif(lwr * 100, 2),
          "-",
          signif(upr * 100, 2),
          "%)"
        ),
        NA
      )
    ) |>
    dplyr::mutate(
      conf = dplyr::case_when(
        lwr > 0.75 ~ detection_conf_cols[4],
        lwr > 0.50 ~ detection_conf_cols[3],
        lwr > 0.25 ~ detection_conf_cols[2],
        .default = detection_conf_cols[1]
      )
    ) |>
    ggplot2::ggplot() +
    ggplot2::aes(x = reads, y = name, fill = conf) +
    ggplot2::geom_col() +
    ggplot2::geom_text(ggplot2::aes(label = fit_lab), hjust = 0, size = 8) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.01, 0.5))
    ) +
    ggplot2::theme_classic(20) +
    ggplot2::theme(
      axis.text.y = ggtext::element_markdown(),
      axis.title.y = ggplot2::element_blank(),
    )
  ggsave(
    paste0("data-raw/model-mistakes/", type, "/", error_type, "/", i, ".png"),
    height = 10,
    width = 16
  )
}

for (i in dloop_false_negatives) {
  plot_model_test(i, "dloop", "false_negatives")
}

for (i in dloop_false_positives) {
  plot_model_test(i, "dloop", "false_positives")
}

for (i in dloop_false_negatives) {
  plot_model_test(i, "dloop", "false_negatives")
}

for (i in dloop_false_positives) {
  plot_model_test(i, "dloop", "false_positives")
}

for (i in dory_false_negatives) {
  plot_model_test(i, "dory", "false_negatives")
}

for (i in dory_false_positives) {
  plot_model_test(i, "dory", "false_positives")
}
