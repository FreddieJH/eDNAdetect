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

model_vis_dory <-
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

model_vis_dloop <-
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
    x = expression(OTU ~ reads + delta),
    y = "Pr(Presence)"
  ) +
  theme_classic(20)

ggsave(
  model_vis_dory,
  filename = "data-raw/model-fit/dory_brms.png",
  height = 10,
  width = 15
)
ggsave(
  model_vis_dloop,
  filename = "data-raw/model-fit/dloop_brms.png",
  height = 10,
  width = 15
)

num2str <- function(tx) {
  div <- findInterval(
    as.numeric(gsub("\\,", "", tx)),
    c(0, 1e3, 1e6, 1e9, 1e12)
  )
  paste0(
    round(as.numeric(gsub("\\,", "", tx)) / 10^(3 * (div - 1)), 2),
    c("", "K", "M", "B", "T")[div]
  )
}

posterior_vis_dory <-
  epred_draws(
    dloop_brms,
    newdata = data.frame(
      reads_log = log(c(1e1, 1e2, 1e3, 5e3, 1e4, 5e4, 1e5, 5e5))
    ),
    re_formula = NA
  ) |>
  mutate(reads = num2str(exp(reads_log))) |>
  mutate(reads = fct_reorder(reads, reads_log)) |>
  ggplot(aes(.epred, col = reads)) +
  geom_density(linewidth = 1.5) +
  xlim(c(0, 1)) +
  scale_x_continuous(label = label_percent()) +
  labs(
    x = "Pr(Presence)",
    y = 'Posterior probability density',
    col = "Number of reads"
  ) +
  theme_classic(20) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = c(0, 1),
    legend.justification = c(-0.1, 1.1),
    legend.background = element_rect(fill = "transparent")
  )

posterior_vis_dloop <-
  epred_draws(
    dloop_brms,
    newdata = data.frame(
      reads_log = log(c(1e1, 1e2, 1e3, 5e3, 1e4, 5e4, 1e5, 5e5))
    ),
    re_formula = NA
  ) |>
  mutate(reads = num2str(exp(reads_log))) |>
  mutate(reads = fct_reorder(reads, reads_log)) |>
  ggplot(aes(.epred, col = reads)) +
  geom_density(linewidth = 1.5) +
  xlim(c(0, 1)) +
  scale_x_continuous(label = label_percent()) +
  labs(
    x = "Pr(Presence)",
    y = 'Posterior probability density',
    col = "Number of reads"
  ) +
  theme_classic(20) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = c(0, 1),
    legend.justification = c(-0.1, 1.1),
    legend.background = element_rect(fill = "transparent")
  )


ggsave(
  posterior_vis_dory,
  filename = "data-raw/model-fit/dory_brms_posterior.png",
  height = 10,
  width = 15
)
ggsave(
  posterior_vis_dloop,
  filename = "data-raw/model-fit/dloop_brms_posterior.png",
  height = 10,
  width = 15
)
