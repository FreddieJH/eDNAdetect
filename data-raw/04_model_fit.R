source("data-raw/03_model_prep.R")

species_mod <- glm(
  possible_detection ~ reads_clr * log(sample_reads),
  family = "binomial",
  data = mod_dat_spp
)

saveRDS(species_mod, "extdata/species_model.rds")

cat(
  "\nModel fit. Object(s) created:\n",
  paste0(c("species_mod"), sep = "\n")
)
