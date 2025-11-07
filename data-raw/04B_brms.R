library(brms)
library(readr)
library(modelr)

if (!file.exists("data-raw/files/cleaned/modeldata_dory.csv")) {
  source("data-raw/02_model_prep.R")
}

if (!file.exists("data-raw/files/cleaned/modeldata_dloop.csv")) {
  source("data-raw/02_model_prep.R")
}

if (!file.exists("data-raw/files/cleaned/catch_lookup.csv")) {
  source("data-raw/01_import_clean.R")
}

catch_lookup <- read_csv("data-raw/files/cleaned/catch_lookup.csv")
modeldata_dory <- read_csv("data-raw/files/cleaned/modeldata_dory.csv")
modeldata_dloop <- read_csv("data-raw/files/cleaned/modeldata_dloop.csv")

dory_simple <- brm(
  possible_detection ~ 1 + reads_log,
  data = modeldata_dory,
  family = bernoulli(),
  chains = 4,
  cores = 4
)

dory_randomeffect <- brm(
  possible_detection ~ 1 + reads_log + (1 | vessel_trip),
  data = modeldata_dory,
  family = bernoulli(),
  chains = 4,
  cores = 4
)

dory_nestedrandomeffect <- brm(
  possible_detection ~ 1 + reads_log + (1 | vessel_trip:hold),
  data = modeldata_dory,
  family = bernoulli(),
  chains = 4,
  cores = 4
)

dory_nested2randomeffect <- brm(
  possible_detection ~ 1 + reads_log + (1 | vessel_trip:hold:replicate),
  data = modeldata_dory,
  family = bernoulli(),
  chains = 4,
  cores = 4
)

# Nested, but without replicate nest, is the best model
loo(
  dory_simple,
  dory_randomeffect,
  dory_nestedrandomeffect,
  dory_nested2randomeffect,
  compare = TRUE
)

dory_brms <- dory_nestedrandomeffect

dloop_simple <- brm(
  possible_detection ~ 1 + reads_log,
  data = modeldata_dloop,
  family = bernoulli(),
  chains = 4,
  cores = 4
)

dloop_randomeffect <- brm(
  possible_detection ~ 1 + reads_log + (1 | vessel_trip),
  data = modeldata_dloop,
  family = bernoulli(),
  chains = 4,
  cores = 4
)

dloop_nestedrandomeffect <- brm(
  possible_detection ~ 1 + reads_log + (1 | vessel_trip:hold),
  data = modeldata_dloop,
  family = bernoulli(),
  chains = 4,
  cores = 4
)

dloop_nested2randomeffect <- brm(
  possible_detection ~ 1 + reads_log + (1 | vessel_trip:hold:replicate),
  data = modeldata_dloop,
  family = bernoulli(),
  chains = 4,
  cores = 4
)

# Nested, but without replicate nest, is the best model
loo(
  dloop_simple,
  dloop_randomeffect,
  dloop_nestedrandomeffect,
  dloop_nested2randomeffect,
  compare = TRUE
)

dloop_brms <- dloop_nestedrandomeffect


usethis::use_data(dory_brms, overwrite = TRUE)
usethis::use_data(dloop_brms, overwrite = TRUE)


modeldata_dory %>%
  data_grid(reads_log = seq_range(reads_log, n = 101)) %>%
  add_epred_rvars(xx3, re_formula = NA) %>%
  ggplot(aes(x = reads_log)) +
  stat_lineribbon(aes(dist = .epred), .width = c(.95), alpha = 1 / 3) +
  geom_point(
    aes(y = possible_detection),
    pch = 21,
    data = modeldata_dory
  ) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Set2") +
  theme_classic(20)

tidybayes::add_epred_draws(xx3)
data.frame(reads_log = seq(3, 13, by = 0.1)) |>
  add_epred_rvars(xx3, re_formula = NA)
posterior::as_draws_df(xx3) |>
  as_tibble()


epred_draws(
  xx3,
  newdata = data.frame(reads_log = seq(3, 13, by = 0.1)),
  re_formula = NA
) |>
  summarise(
    fit = quantile(.epred, 0.5),
    lwr = quantile(.epred, 0.025),
    upr = quantile(.epred, 0.975)
  ) |>
  ggplot(aes(reads_log, fit)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.1)
