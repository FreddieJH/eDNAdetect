library(tidybayes)
library(dplyr)

if (!file.exists("data/dory_brms.rda")) {
  source("data-raw/04B_brms.R")
}

if (!file.exists("data/dloop_brms.rda")) {
  source("data-raw/04B_brms.R")
}

load("data/dory_brms.rda")
load("data/dloop_brms.rda")

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
