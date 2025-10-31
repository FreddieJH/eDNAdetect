library(wesanderson)
library(tidybayes)
library(dplyr)
library(ggplot2)
library(readr)
library(scales)

if (!file.exists("data/dory_brms.rda")) {
  source("data-raw/04B_brms.R")
}

if (!file.exists("data/dloop_brms.rda")) {
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

pred_detect <- function(reads_vec, model, ci) {
  ci_lwr <- (1 - ci) / 2
  ci_upr <- 1 - ci_lwr

  epred_draws(
    dory_brms,
    newdata = data.frame(reads_log = reads_vec),
    re_formula = NA
  ) |>
    summarise(
      fit = quantile(.epred, 0.5),
      lwr = quantile(.epred, ci_lwr),
      upr = quantile(.epred, ci_upr),
      .groups = "drop"
    )
}

pred_detect(seq(3, 13, 0.1), dory_brms, ci = 0.95) |>
  ggplot(aes(exp(reads_log), fit)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, fill = mypal[1]) +
  geom_line(aes(y = fit), col = mypal[1]) +
  geom_point(
    aes(y = possible_detection),
    pch = 21,
    data = modeldata_dory
  ) +
  scale_y_continuous(labels = label_percent()) +
  scale_x_log10(labels = label_log(), guide = "axis_logticks") +
  labs(
    x = expression(log(OTU ~ reads + delta)),
    y = "Pr(Presence)"
  ) +
  theme_classic(20)

pred_detect(seq(3, 13, 0.1), dloop_brms, ci = 0.95) |>
  ggplot(aes(exp(reads_log), fit)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, fill = mypal[3]) +
  geom_line(aes(y = fit), col = mypal[3]) +
  geom_point(
    aes(y = possible_detection),
    pch = 21,
    data = modeldata_dloop
  ) +
  scale_y_continuous(labels = label_percent()) +
  scale_x_log10(labels = label_log(), guide = "axis_logticks") +
  labs(
    x = expression(log(OTU ~ reads + delta)),
    y = "Pr(Presence)"
  ) +
  theme_classic(20)
