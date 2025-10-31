library(readr)
library(purrr)
library(tools)
library(dplyr)
library(tidyr)
library(stringr)

# OTU to species lookup
otu_lookup <-
  bind_rows(
    read_csv("data-raw/files/labels_dory.csv") |>
      select(label = Label, family, genus, species) |>
      mutate(label = paste0(label, "_DORY")),
    read_csv("data-raw/files/labels_dloop.csv") |>
      select(label, family, genus, species) |>
      mutate(label = paste0(label, "_DLOOP"))
  )

write_csv(otu_lookup, "data-raw/files/cleaned/otu_lookup.csv")

# OTUs that correspond to Unidentified Genera
unidentified_genus <-
  otu_lookup |>
  select(-label) |>
  filter(str_detect(species, "(genus)")) |>
  distinct()

# OTUs that correspond to Unidentified Family
unidentified_family <-
  otu_lookup |>
  select(-label) |>
  filter(str_detect(species, "(family)")) |>
  distinct()

# what was caught on each vessell/trip/hold?
catch_lookup <-
  read_csv("data-raw/files/catch.csv") |>
  separate_longer_delim(cols = catch, delim = ",") |>
  mutate(catch = str_trim(catch)) |>
  mutate(across(everything(), as.character)) |>
  left_join(
    otu_lookup |>
      select(family, genus, catch = species) |>
      distinct()
  ) |>
  left_join(
    unidentified_genus |> select(family, genus, genus_catch = species)
  ) |>
  left_join(unidentified_family |> select(family, family_catch = species)) |>
  select(-c(family, genus))

write_csv(catch_lookup, "data-raw/files/cleaned/catch_lookup.csv")

# have to rename the samples so they match the catch data
sample_replacements <- c(
  "Mock_a" = "MOCK_11_1a",
  "Mock_b" = "MOCK_11_2b",
  "Mock_c" = "MOCK_11_3c",
  "Mock_d" = "MOCK_11_4d",
  "Yellowfin_tuna" = "MOCK_21_1a",
  "Southern_bluefin_tuna" = "MOCK_31_1a",
  "Yellowfin_tuna_2R" = "MOCK_21_2b",
  "Southern_bluefin_tuna_2R" = "MOCK_31_2b",
  "Mock_SBT_YFT" = "MOCK_41_1a"
)

all_samples <-
  c(
    c(names(
      read_csv("data-raw/files/reads_dory.csv") |> select(-c(Label, Sequence))
    )),
    c(names(read_csv("data-raw/files/reads_dloop.csv") |> select(-c(label))))
  ) |>
  unique()

# sample name to vessell/trip/hold/replicate
sample_lookup <-
  tibble(
    sample = all_samples
  ) |>
  mutate(sample = recode(sample, !!!sample_replacements)) |>
  mutate(
    vessel = str_extract(sample, "^[a-zA-Z]+"),
    trip = str_extract(sample, "(?<=_)[0-9]"),
    hold = str_extract(sample, "(?<=_[0-9]{1})[0-9]"),
    replicate = str_extract(sample, "[1-9a-z]+$")
  )

write_csv(sample_lookup, "data-raw/files/cleaned/sample_lookup.csv")

# Cleanup reads data
# ... and convert to long format for analysis
clean_reads <- function(df, suffix) {
  df |>
    rename_with(~"label", .cols = matches("^[Ll]abel$")) |>
    mutate(label = paste0(label, suffix)) |>
    select(-any_of("Sequence")) |>
    pivot_longer(
      cols = -label,
      names_to = "sample",
      values_to = "reads"
    ) |>
    left_join(otu_lookup, by = "label") |>
    filter(!species %in% c("Homo sapiens", "Canis lupus")) |>
    mutate(sample = recode(sample, !!!sample_replacements)) |>
    left_join(sample_lookup, by = "sample")
}

write_csv(
  x = read_csv("data-raw/files/reads_dory.csv") |> clean_reads("_DORY"),
  file = "data-raw/files/cleaned/reads_dory_clean.csv"
)
write_csv(
  x = read_csv("data-raw/files/reads_dloop.csv") |> clean_reads("_DLOOP"),
  file = "data-raw/files/cleaned/reads_dloop_clean.csv"
)
