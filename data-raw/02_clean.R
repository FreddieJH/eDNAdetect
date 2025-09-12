source("data-raw/01_import.R")

# 95 samples
sample_lookup <-
  catch_raw |>
  dplyr::select(sample = Sample) |>
  dplyr::mutate(
    vessel = stringr::str_extract(sample, "^[a-zA-Z]+"),
    trip = stringr::str_extract(sample, "(?<=_)[0-9]"),
    hold = stringr::str_extract(sample, "(?<=_[0-9]{1})[0-9]"),
    replicate = stringr::str_extract(sample, "[1-9a-z]+$")
  )

# 20 labels
label_lookup <-
  label_raw |>
  dplyr::select(label = Label, family, genus, species) |>
  dplyr::mutate(
    genus = dplyr::case_when(
      stringr::str_detect(genus, "(family)") ~ "Other",
      .default = genus
    ),
    species = dplyr::case_when(
      stringr::str_detect(species, "Thunnus") ~ "Thunnus spp", # includes genus and spp level thuns
      stringr::str_detect(species, "(family)|(genus)") ~ "Other",
      .default = species
    )
  )


data_long <-
  data_raw |>
  dplyr::select(-Sequence) |>
  dplyr::rename(label = Label) |>
  tidyr::pivot_longer(
    cols = -label,
    names_to = "sample",
    values_to = "reads"
  ) |>
  dplyr::left_join(
    label_lookup |> dplyr::select(label, species),
    by = dplyr::join_by(label)
  ) |>
  dplyr::filter(!(species %in% c("Homo sapiens", "Canis lupus"))) |>
  dplyr::mutate(sample_reads = sum(reads), .by = sample) |>
  dplyr::mutate(reads_prop = reads / sample_reads) |>
  dplyr::mutate(species = as.factor(species), sample = as.factor(sample))

catch_lookup_specific <-
  catch_raw |>
  dplyr::select(sample = Sample, catch_species_specific = Catch) |>
  tidyr::separate_longer_delim(cols = catch_species_specific, delim = ",") |>
  dplyr::mutate(
    catch_species_specific = stringr::str_trim(catch_species_specific),
    catch_species = dplyr::case_when(
      stringr::str_detect(catch_species_specific, "Thunnus") ~ "Thunnus spp",
      .default = catch_species_specific
    )
  )

genus_family <-
  dplyr::tibble(
    genus = c(
      "Xiphias",
      "Thunnus",
      "Brama",
      "Kajikia",
      "Tetrapturus",
      "Coryphaena",
      "Carcharhinus"
    ),
    family = c(
      "Xiphiidae",
      "Scombridae",
      "Scombridae",
      "Istiophoridae",
      "Istiophoridae",
      "Coryphaenidae",
      "Carcharhinidae"
    )
  )


catch_lookup <-
  catch_raw |>
  dplyr::select(sample = Sample, catch_species_specific = Catch) |>
  tidyr::separate_longer_delim(cols = catch_species_specific, delim = ",") |>
  dplyr::mutate(
    catch_species_specific = stringr::str_trim(catch_species_specific)
  ) |>
  dplyr::mutate(
    species = dplyr::case_when(
      stringr::str_detect(catch_species_specific, "Thunnus") ~ "Thunnus spp",
      .default = catch_species_specific
    )
  ) |>
  dplyr::select(sample, species) |>
  dplyr::distinct() |>
  dplyr::mutate(genus = stringr::str_extract(species, "^\\w+")) |>
  dplyr::left_join(genus_family, by = dplyr::join_by(genus))


spp_lookup <-
  dplyr::tibble(species = unique(catch_lookup$species)) |>
  dplyr::mutate(genus = stringr::str_extract(species, "^\\w+")) |>
  dplyr::left_join(genus_family, by = dplyr::join_by(genus))

# This will give us all the labels that might be detected for each of the catch species
detection_lookup <-
  label_lookup |>
  dplyr::rename_with(.fn = ~ paste0("label_", .x), .cols = -label) |>
  tidyr::expand_grid(spp_lookup) |>
  dplyr::mutate(
    detect = dplyr::case_when(
      species == label_species ~ TRUE,
      genus == label_genus ~ TRUE,
      family == label_family ~ TRUE,
      .default = FALSE
    )
  ) |>
  dplyr::filter(detect == TRUE)


# Grouping data by family, genus or species ------------------------------------------

data_family <-
  data_long |>
  dplyr::left_join(
    label_lookup |> dplyr::rename(label_species = species),
    by = dplyr::join_by(label)
  ) |>
  dplyr::summarise(
    family_reads = sum(reads),
    .by = c(sample, sample_reads, family)
  )


data_genus <-
  data_long |>
  dplyr::left_join(
    label_lookup |> dplyr::rename(label_species = species),
    by = dplyr::join_by(label)
  ) |>
  dplyr::mutate(
    genus = dplyr::case_when(
      stringr::str_detect(genus, "(family)") ~ "Other",
      .default = genus
    )
  ) |>
  dplyr::summarise(
    genus_reads = sum(reads),
    .by = c(sample, sample_reads, genus)
  )
# unique(data_genus$genus)

data_species <-
  data_long |>
  dplyr::left_join(
    label_lookup |> dplyr::rename(label_species = species),
    by = dplyr::join_by(label)
  ) |>
  dplyr::mutate(
    species = dplyr::case_when(
      stringr::str_detect(label_species, "Thunnus") ~ "Thunnus spp", # includes genus and spp level thuns
      stringr::str_detect(label_species, "(family)|(genus)") ~ "Other",
      .default = label_species
    )
  ) |>
  dplyr::summarise(
    species_reads = sum(reads),
    .by = c(sample, sample_reads, species)
  )
# unique(data_species$species)

rm(data_raw, label_raw, catch_raw, data_long)

cat("\nData successfully cleaned. Objects created:\n", paste0(ls(), sep = "\n"))
