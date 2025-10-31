library(ggplot2)
library(readr)
library(dplyr)
library(purrr)

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


plot_samples <- function(df, prop = TRUE) {
  df |>
    ggplot2::ggplot() +
    {
      if (prop) ggplot2::aes(x = reads / sample_reads)
    } +
    {
      if (!prop) ggplot2::aes(x = reads)
    } +
    ggplot2::aes(y = sample, fill = species) +
    ggplot2::geom_col() +
    ggplot2::labs(x = "Reads", y = NULL)
}

if (FALSE) {
  plot_samples(modeldata_dory |> dplyr::filter(vessel == "EM"), prop = TRUE)
  plot_samples(modeldata_dloop |> dplyr::filter(vessel == "EM"), prop = TRUE)
  plot_samples(modeldata_dory |> dplyr::filter(vessel == "EM"), prop = FALSE)
  plot_samples(modeldata_dloop |> dplyr::filter(vessel == "EM"), prop = FALSE)
}


plot_sample <- function(df, prop = FALSE) {
  axis_colours <-
    df %>%
    distinct(species, possible_detection) %>%
    arrange(species) %>%
    mutate(colour = ifelse(possible_detection == 1, "blue", "black")) %>%
    pull(colour)

  df |>
    # arrange(desc(reads)) |>
    ggplot2::ggplot() +
    {
      if (prop) ggplot2::aes(x = reads / sample_reads)
    } +
    {
      if (!prop) ggplot2::aes(x = reads)
    } +
    ggplot2::aes(y = species, fill = as.factor(possible_detection)) +
    ggplot2::geom_col() +
    ggplot2::labs(
      x = "Reads",
      y = NULL,
      title = unique(df$sample),
      subtitle = catch_lookup |>
        dplyr::filter(
          vessel %in% df$vessel,
          trip %in% df$trip,
          hold %in% df$hold
        ) |>
        dplyr::pull(catch) |>
        paste(collapse = ", ")
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      axis.text.y = element_text(colour = axis_colours)
    )
}
if (FALSE) {
  plot_sample(
    modeldata_dory |> dplyr::filter(sample == "EM_53.1a"),
    prop = FALSE
  )
  plot_sample(
    modeldata_dloop |> dplyr::filter(sample == "EM_53.1a"),
    prop = FALSE
  )
}

list(
  dory = modeldata_dory,
  dloop = modeldata_dloop
) |>
  imap(function(data, name) {
    outdir <- paste0("data-raw/data-vis/", name, "_bysample/")

    unique(data$sample) |>
      walk(function(s) {
        filename <- paste0(outdir, s, ".png")
        if (!file.exists(filename)) {
          p <- plot_sample(filter(data, sample == s), prop = FALSE)
          ggsave(plot = p, filename = filename, height = 6, width = 6)
        }
      })
  })
