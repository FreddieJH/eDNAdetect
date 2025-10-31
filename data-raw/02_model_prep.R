library(compositions)
library(readr)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)

if (!file.exists("data-raw/files/cleaned/reads_dory_clean.csv")) {
  source("data-raw/01_import_clean.R")
}

if (!file.exists("data-raw/files/cleaned/reads_dloop_clean.csv")) {
  source("data-raw/01_import_clean.R")
}

if (!file.exists("data-raw/files/cleaned/catch_lookup.csv")) {
  source("data-raw/01_import_clean.R")
}

# centered log-ratio but with multiplicative pseudocount
clr_ps <- function(reads_vec) {
  min_nonzero <- min(reads_vec[reads_vec > 0])
  pseudocount <- 0.01 * min_nonzero

  as.numeric(clr(reads_vec + pseudocount))
}

log_plus_halfSmallest <- function(
  reads_vec,
  pseudocount = 25,
  log_base = exp(1)
) {
  as.numeric(log(reads_vec + pseudocount, base = log_base))
}

prep_data <- function(df) {
  df |>
    left_join(
      read_csv("data-raw/files/cleaned/catch_lookup.csv") |>
        nest(.by = c(vessel, trip, hold), .key = "catch"),
      by = c("vessel", "trip", "hold")
    ) |>
    dplyr::mutate(
      spp_in_catch = purrr::map2_dbl(
        .x = species,
        .y = catch,
        .f = function(x, y) {
          x %in% y$catch
        }
      ),
      genus_in_catch = purrr::map2_dbl(
        .x = species,
        .y = catch,
        .f = function(x, y) {
          x %in% y$genus_catch
        }
      ),
      family_in_catch = purrr::map2_dbl(
        .x = species,
        .y = catch,
        .f = function(x, y) {
          x %in% y$family_catch
        }
      ),
    ) |>
    dplyr::mutate(
      possible_detection = as.numeric(
        rowSums(dplyr::pick(spp_in_catch, genus_in_catch, family_in_catch)) > 0
      )
    ) |>
    dplyr::mutate(
      sample_reads = sum(reads),
      reads_clr = clr_ps(reads),
      reads_prop = reads / sample_reads,
      reads_log = log_plus_halfSmallest(reads, min(df$reads[df$reads > 0]) / 2),
      .by = sample
    ) |>
    select(-catch) |>
    mutate(vessel_trip = paste(vessel, trip)) |>
    mutate(
      vessel_trip = as.factor(vessel_trip),
      hold = as.factor(hold),
      replicate = as.factor(replicate)
    )
}

write_csv(
  x = read_csv("data-raw/files/cleaned/reads_dory_clean.csv") |>
    prep_data() |>
    filter(
      (!str_detect(species, "Unidentified")) |
        species == "Unidentified Thunnus (genus)"
    ) |>
    filter(!str_detect(species, "Thunnus albacares")),
  file = "data-raw/files/cleaned/modeldata_dory.csv"
)
write_csv(
  x = read_csv("data-raw/files/cleaned/reads_dloop_clean.csv") |>
    prep_data() |>
    filter(!str_detect(species, "Unidentified")),
  file = "data-raw/files/cleaned/modeldata_dloop.csv"
)
