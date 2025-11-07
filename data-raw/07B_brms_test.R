library(dplyr)
library(modelr)
library(wesanderson)
library(tidybayes)
library(ggplot2)

if (!file.exists("data/dory_brms.rda")) {
  source("data-raw/04B_brms.R")
}

if (!file.exists("data/dloop_brms.rda")) {
  source("data-raw/04B_brms.R")
}

if (!file.exists("data-raw/files/cleaned/modeldata_dory.csv")) {
  source("data-raw/02_model_prep.R")
}

if (!file.exists("data-raw/files/cleaned/modeldata_dloop.csv")) {
  source("data-raw/02_model_prep.R")
}

load("data/dory_brms.rda")
load("data/dloop_brms.rda")
modeldata_dory <- read_csv("data-raw/files/cleaned/modeldata_dory.csv")
modeldata_dloop <- read_csv("data-raw/files/cleaned/modeldata_dloop.csv")
mypal <- wes_palette("Zissou1")
detection_conf_cols = c(
  "No detection" = "grey70",
  "Low confidence" = mypal[5],
  "Medium confidence" = mypal[3],
  "High confidence" = mypal[2]
)


dory_preds <- epred_draws(
  dory_brms,
  newdata = modeldata_dory |> select(reads_log),
  re_formula = NA
) |>
  group_by(reads_log) |>
  summarise(
    fit = quantile(.epred, 0.5),
    lwr = quantile(.epred, 0.1),
    upr = quantile(.epred, 0.9)
  )

dory_plot_data <-
  modeldata_dory %>%
  left_join(dory_preds, by = "reads_log") %>%
  mutate(
    model_status = case_when(
      possible_detection & lwr > 0.25 ~ "True positive",
      !possible_detection & lwr < 0.25 ~ "True negative",
      possible_detection & lwr < 0.25 ~ "False negative",
      !possible_detection & lwr > 0.25 ~ "False positive",
      TRUE ~ NA
    )
  ) |>
  mutate(
    actual = if_else(possible_detection == 1, "Detection", "No detection"),
    predicted = if_else(lwr > 0.25, "Predicted yes", "Predicted no")
  ) %>%
  count(actual, predicted, model_status) %>%
  mutate(
    actual = factor(actual, levels = c("No detection", "Detection")),
    predicted = factor(predicted, levels = c("Predicted no", "Predicted yes"))
  )


dory_pred_marginals <-
  dory_plot_data %>%
  summarise(pred_n = sum(n), .by = predicted) %>%
  arrange(predicted) %>%
  mutate(
    xmin = lag(cumsum(pred_n), default = 0) / sum(plot_data$n),
    xmax = cumsum(pred_n) / sum(plot_data$n)
  )

percentages_plot_dory <-
  dory_plot_data %>%
  left_join(dory_pred_marginals, by = "predicted") %>%
  arrange(predicted, actual) %>%
  mutate(
    prop_within_pred = n / pred_n,
    ymin = lag(cumsum(prop_within_pred), default = 0),
    ymax = cumsum(prop_within_pred),
    .by = predicted
  ) %>%
  mutate(
    xmin = xmin,
    xmax = xmax,
    label = paste0(
      model_status,
      " n = ",
      n,
      " (",
      percent(prop_within_pred, accuracy = 0.1),
      ")"
    ),
    cx = (xmin + xmax) / 2,
    cy = (ymin + ymax) / 2
  ) |>
  mutate(correct = str_detect(model_status, "True")) |>
  ggplot() +
  geom_rect(
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = correct),
    colour = "white",
    size = 0.6
  ) +
  scale_fill_discrete(palette = mypal[2:3]) +
  geom_text(aes(x = cx, y = cy, label = label), size = 5) +
  theme_void(20) +
  theme(legend.position = "none")


dloop_preds <- epred_draws(
  dloop_brms,
  newdata = modeldata_dloop |> select(reads_log),
  re_formula = NA
) |>
  group_by(reads_log) |>
  summarise(
    fit = quantile(.epred, 0.5),
    lwr = quantile(.epred, 0.1),
    upr = quantile(.epred, 0.9)
  )

dloop_plot_data <-
  modeldata_dloop %>%
  left_join(dloop_preds, by = "reads_log") %>%
  mutate(
    model_status = case_when(
      possible_detection & lwr > 0.25 ~ "True positive",
      !possible_detection & lwr < 0.25 ~ "True negative",
      possible_detection & lwr < 0.25 ~ "False negative",
      !possible_detection & lwr > 0.25 ~ "False positive",
      TRUE ~ NA
    )
  ) |>
  mutate(
    actual = if_else(possible_detection == 1, "Detection", "No detection"),
    predicted = if_else(lwr > 0.25, "Predicted yes", "Predicted no")
  ) %>%
  count(actual, predicted, model_status) %>%
  mutate(
    actual = factor(actual, levels = c("No detection", "Detection")),
    predicted = factor(predicted, levels = c("Predicted no", "Predicted yes"))
  )


dloop_pred_marginals <-
  dloop_plot_data %>%
  summarise(pred_n = sum(n), .by = predicted) %>%
  arrange(predicted) %>%
  mutate(
    xmin = lag(cumsum(pred_n), default = 0) / sum(plot_data$n),
    xmax = cumsum(pred_n) / sum(plot_data$n)
  )

percentages_plot_dloop <-
  dloop_plot_data %>%
  left_join(dloop_pred_marginals, by = "predicted") %>%
  arrange(predicted, actual) %>%
  mutate(
    prop_within_pred = n / pred_n,
    ymin = lag(cumsum(prop_within_pred), default = 0),
    ymax = cumsum(prop_within_pred),
    .by = predicted
  ) %>%
  mutate(
    xmin = xmin,
    xmax = xmax,
    label = paste0(
      model_status,
      " n = ",
      n,
      " (",
      percent(prop_within_pred, accuracy = 0.1),
      ")"
    ),
    cx = (xmin + xmax) / 2,
    cy = (ymin + ymax) / 2
  ) |>
  mutate(correct = str_detect(model_status, "True")) |>
  ggplot() +
  geom_rect(
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = correct),
    colour = "white",
    size = 0.6
  ) +
  scale_fill_discrete(palette = mypal[2:3]) +
  geom_text(aes(x = cx, y = cy, label = label), size = 5) +
  theme_void(20) +
  theme(legend.position = "none")

ggsave(
  percentages_plot_dory,
  filename = "data-raw/model-fit/dory_brms_percentages.png",
  height = 10,
  width = 15
)
ggsave(
  percentages_plot_dloop,
  filename = "data-raw/model-fit/dloop_brms_percentages.png",
  height = 10,
  width = 15
)

mod_dat_dloop_withPreds <-
  modeldata_dloop |>
  dplyr::mutate(
    mod_fit = as.numeric(
      predict(
        dloop_brms,
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

plot_model_test <- function(tbl) {
  tbl |>
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
}


dory_withpreds <-
  modeldata_dory |>
  left_join(dory_preds, by = "reads_log") %>%
  mutate(
    model_status = case_when(
      possible_detection & lwr > 0.25 ~ "True positive",
      !possible_detection & lwr < 0.25 ~ "True negative",
      possible_detection & lwr < 0.25 ~ "False negative",
      !possible_detection & lwr > 0.25 ~ "False positive",
      TRUE ~ NA
    )
  )

dloop_withpreds <-
  modeldata_dloop |>
  left_join(dloop_preds, by = "reads_log") %>%
  mutate(
    model_status = case_when(
      possible_detection & lwr > 0.25 ~ "True positive",
      !possible_detection & lwr < 0.25 ~ "True negative",
      possible_detection & lwr < 0.25 ~ "False negative",
      !possible_detection & lwr > 0.25 ~ "False positive",
      TRUE ~ NA
    )
  )


TP_dory <- dory_withpreds |>
  filter(model_status == "True positive") |>
  pull(sample) |>
  unique()
TN_dory <- dory_withpreds |>
  filter(model_status == "True negative") |>
  pull(sample) |>
  unique()
FP_dory <- dory_withpreds |>
  filter(model_status == "False positive") |>
  pull(sample) |>
  unique()
FN_dory <- dory_withpreds |>
  filter(model_status == "False negative") |>
  pull(sample) |>
  unique()


TP_dloop <- dloop_withpreds |>
  filter(model_status == "True positive") |>
  pull(sample) |>
  unique()
TN_dloop <- dloop_withpreds |>
  filter(model_status == "True negative") |>
  pull(sample) |>
  unique()
FP_dloop <- dloop_withpreds |>
  filter(model_status == "False positive") |>
  pull(sample) |>
  unique()
FN_dloop <-
  dloop_withpreds |>
  filter(model_status == "False negative") |>
  pull(sample) |>
  unique()


for (i in FP_dory) {
  p <- plot_model_test(dory_withpreds |> filter(sample == i))
  ggsave(
    plot = p,
    paste0("data-raw/model-mistakes/dory/falsePositive/", i, ".png"),
    height = 10,
    width = 16
  )
}

for (i in FN_dory) {
  p <- plot_model_test(dory_withpreds |> filter(sample == i))
  ggsave(
    plot = p,
    paste0("data-raw/model-mistakes/dory/falseNegative/", i, ".png"),
    height = 10,
    width = 16
  )
}

for (i in FP_dloop) {
  p <- plot_model_test(dloop_withpreds |> filter(sample == i))
  ggsave(
    plot = p,
    paste0("data-raw/model-mistakes/dloop/falsePositive/", i, ".png"),
    height = 10,
    width = 16
  )
}

for (i in FN_dloop) {
  p <- plot_model_test(dloop_withpreds |> filter(sample == i))
  ggsave(
    plot = p,
    paste0("data-raw/model-mistakes/dloop/falseNegative/", i, ".png"),
    height = 10,
    width = 16
  )
}

for (i in unique(dory_withpreds$sample)) {
  p <- plot_model_test(dory_withpreds |> filter(sample == i))
  ggsave(
    plot = p,
    paste0("data-raw/model-ests/dory/", i, ".png"),
    height = 10,
    width = 16
  )
}

for (i in unique(dloop_withpreds$sample)) {
  p <- plot_model_test(dloop_withpreds |> filter(sample == i))
  ggsave(
    plot = p,
    paste0("data-raw/model-ests/dloop/", i, ".png"),
    height = 10,
    width = 16
  )
}
