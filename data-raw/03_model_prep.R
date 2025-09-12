source("data-raw/02_clean.R")

pkgs <- c("compositions")
new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) {
  install.packages(new_pkgs)
}

mod_dat_spp <-
  data_species |>
  dplyr::left_join(
    catch_lookup |> tidyr::nest(.by = sample, .key = "catch"),
    by = dplyr::join_by(sample)
  ) |>
  dplyr::mutate(
    possible_detection = purrr::map2_dbl(
      .x = species,
      .y = catch,
      .f = function(x, y) {
        x %in% y$species
      }
    )
  ) |>
  dplyr::filter(species != 'Other') |>
  dplyr::mutate(
    sample_reads = sum(species_reads),
    reads_clr = as.numeric(compositions::clr(
      species_reads + min(species_reads[species_reads > 0]) / 2
    )),
    .by = sample
  )

cat(
  "\nModelling data produced. Objects created:\n",
  paste0(c("mod_dat_spp"), sep = "\n")
)
