library(stringr)

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

modeldata_dory_cleaned <-
  modeldata_dory |>
  filter(
    (!str_detect(species, "Unidentified")) |
      species == "Unidentified Thunnus (genus)"
  ) |>
  filter(!str_detect(species, "Thunnus albacares")) |>
  mutate(vessel_trip = paste(vessel, trip)) |>
  mutate(
    vessel_trip = as.factor(vessel_trip),
    hold = as.factor(hold),
    replicate = as.factor(replicate)
  )

modeldata_dloop_cleaned <-
  modeldata_dloop |>
  filter(!str_detect(species, "Unidentified")) |>
  mutate(vessel_trip = paste(vessel, trip)) |>
  mutate(
    vessel_trip = as.factor(vessel_trip),
    hold = as.factor(hold),
    replicate = as.factor(replicate)
  )


dory_mod_log <- glm(
  possible_detection ~ reads_log,
  family = "binomial",
  data = modeldata_dory_cleaned
)

dory_mod_log_nestedEffects <- lme4::glmer(
  possible_detection ~ reads_log + (1 | vessel_trip),
  family = "binomial",
  data = modeldata_dory_cleaned
)

dory_mod_log_nestedEffects2 <- lme4::glmer(
  possible_detection ~ reads_log + (1 | vessel_trip:hold),
  family = binomial,
  data = modeldata_dory_cleaned
)

dory_mod_prop <- glm(
  possible_detection ~ reads_prop,
  family = "binomial",
  data = modeldata_dory_cleaned
)

dory_mod_randomEffect <- lme4::glmer(
  possible_detection ~ reads_log + (1 | sample_reads),
  family = "binomial",
  data = modeldata_dory_cleaned
)

dory_mod_randomEffect2 <- lme4::glmer(
  possible_detection ~ reads_prop + (1 | sample_reads),
  family = "binomial",
  data = modeldata_dory_cleaned
)

AIC(
  dory_mod_log,
  dory_mod_prop,
  dory_mod_randomEffect,
  dory_mod_randomEffect2,
  dory_mod_log_nestedEffects,
  dory_mod_log_nestedEffects2
) |>
  arrange(AIC)

dory_mod <- dory_mod_log_nestedEffects


dloop_mod_log <- glm(
  possible_detection ~ reads_log,
  family = "binomial",
  data = modeldata_dloop |> filter(!str_detect(species, "Unidentified"))
)

dloop_mod_log_nestedEffects <- lme4::glmer(
  possible_detection ~ reads_log + (1 | vessel_trip),
  family = "binomial",
  data = modeldata_dloop |>
    filter(!str_detect(species, "Unidentified")) |>
    mutate(vessel_trip = paste(vessel, trip))
)

dloop_mod_log_prop <- glm(
  possible_detection ~ reads_log + reads_prop,
  family = "binomial",
  data = modeldata_dloop |> filter(!str_detect(species, "Unidentified"))
)

dloop_mod_prop <- glm(
  possible_detection ~ reads_prop,
  family = "binomial",
  data = modeldata_dloop |> filter(!str_detect(species, "Unidentified"))
)
dloop_mod_randomEffect <- lme4::glmer(
  possible_detection ~ reads_log + (1 | sample_reads_factor),
  family = "binomial",
  data = modeldata_dloop |>
    filter(!str_detect(species, "Unidentified")) |>
    mutate(sample_reads_factor = as.factor(sample_reads))
)

dloop_mod_randomEffect2 <- lme4::glmer(
  possible_detection ~ reads_log + (reads_log | sample_reads_factor),
  family = "binomial",
  data = modeldata_dloop |>
    filter(!str_detect(species, "Unidentified")) |>
    mutate(sample_reads_factor = as.factor(sample_reads))
)


AIC(
  dloop_mod_log,
  dloop_mod_log_prop,
  dloop_mod_prop,
  dloop_mod_randomEffect,
  dloop_mod_log_nestedEffects
) |>
  arrange(AIC)

dloop_mod <- dloop_mod_log_nestedEffects

# THIS ONE FOR CLR
dloop_mod_additive <- glm(
  possible_detection ~ reads_clr + log(sample_reads),
  family = "binomial",
  data = modeldata_dloop
)
dloop_mod_interaction <- glm(
  possible_detection ~ reads_clr * log(sample_reads),
  family = "binomial",
  data = modeldata_dloop
)

anova(dloop_mod_additive, dloop_mod_interaction, test = "Chisq")
AIC(dloop_mod, dloop_mod_additive, dloop_mod_interaction)

dloop_mod <- dloop_mod_log_nestedEffects

usethis::use_data(dory_mod, overwrite = TRUE)
usethis::use_data(dloop_mod, overwrite = TRUE)

cat(
  "\nModel fit. Object(s) created:\n",
  paste0(c("dory_mod", "dloop_mod"), sep = "\n")
)

# modeldata_dloop |>
#   dplyr::arrange(desc(reads_clr)) |>
#   ggplot2::ggplot(ggplot2::aes(reads_clr, possible_detection)) +
#   ggplot2::geom_point(pch = 21, alpha = .5)

# compositions::clr(c(100, 100, 100))
# compositions::clr(c(1000, 100, 100))

# compositions::clr(c(5000, 5000, 100, 0))
# compositions::clr(c(5000, 5000, 100))
# compositions::clr(c(5001000, 5000000, 4999000, 0, 0, 1))
# compositions::clr(c(5001, 5000, 4999, 1, 100, 100, 10, 0, 0, 1000))

# rclr(c(5000, 5000, 100))
# rclr(c(5001, 5000, 4999) * 1000)
# rclr(c(5001, 5000, 4999, 1, 100, 100, 10))

# rclr_transform(c(5001, 5000, 4999))
# rclr_transform(c(5001, 5000, 4999, 1, 100, 100, 10))
# # Load the vegan package
# library(vegan)

# # Example data (replace with your compositional data)
# data(varespec)

# # Perform robust CLR transformation
# # impute = TRUE handles zeros by only using observed taxa for geometric mean calculations
# varespec.rclr <- decostand(varespec, "rclr", impute = TRUE)

# # View the transformed data
# head(varespec.rclr)

# modeldata_dloop |>
#   dplyr::filter(species == "Thunnus obesus") |>
#   ggplot2::ggplot() +
#   ggplot2::aes(log(reads), possible_detection) +
#   ggplot2::geom_point()

# modeldata_dloop |>
#   dplyr::filter(species == "Thunnus obesus") |>
#   ggplot2::ggplot() +
#   ggplot2::aes(reads_prop, possible_detection) +
#   ggplot2::geom_point()

# modeldata_dloop |>
#   dplyr::filter(species == "Thunnus obesus") |>
#   ggplot2::ggplot() +
#   ggplot2::aes(reads_clr, possible_detection) +
#   ggplot2::geom_point()

# modeldata_dloop |>
#   dplyr::filter(species == "Thunnus obesus") |>
#   ggplot2::ggplot() +
#   ggplot2::aes(, possible_detection) +
#   ggplot2::geom_point()

# devtools::install_github("antagomir/vegan")
# library(vegan)

# # Test data
# set.seed(252)
# testdata <- matrix(round(runif(1000, 0, 100)), nrow = 20)
# testdata <- testdata - 50
# testdata[testdata < 0] <- 0
# rownames(testdata) <- paste0("row", seq_len(nrow(testdata)))
# colnames(testdata) <- paste0("col", seq_len(ncol(testdata)))

# # Aitchison equals to CLR + Euclid (pseudocount is necessary with clr)
# a1 <- vegan::vegdist(testdata + 1, method = "aitchison")
# a2 <- vegan::vegdist(
#   vegan::decostand(testdata + 1, "clr"),
#   method = "euclidean"
# )
# max(abs(a1 - a2)) < 1e-6 # Tolerance

# # Robust aitchison equals to rCLR + Euclid
# # and works without pseudocount
# a1 <- vegan::vegdist(testdata, method = "robust.aitchison")
# a2 <- vegan::vegdist(vegan::decostand(testdata, "rclr"), method = "euclidean")
# max(abs(a1 - a2)) < 1e-6 # Tolerance

# # Robust aitchison and aitchison are equal when there are no zeroes
# a1 <- vegan::vegdist(testdata + 1, method = "robust.aitchison")
# a2 <- vegan::vegdist(testdata + 1, method = "aitchison")
# max(abs(a1 - a2)) < 1e-6 # Tolerance

# rclr <- function(x, na.rm = TRUE) {
#   # Error with negative values
#   if (any(x < 0, na.rm = na.rm)) {
#     stop("'rclr' cannot be used with negative data", call. = FALSE)
#   }
#   # Log transform
#   clog <- log(x)
#   # Convert zeros to NAs in rclr
#   clog[is.infinite(clog)] <- NA
#   # Calculate log of geometric mean for every sample, ignoring the NAs
#   mean_clog <- mean(clog, na.rm = na.rm)
#   # Divide all values by their sample-wide geometric means
#   # Log and transpose back to original shape
#   xx <- log(x) - mean_clog
#   # If there were zeros, there are infinite values after logarithmic transform.
#   # Convert those to zero.
#   xx[is.infinite(xx)] <- 0
#   attr(xx, "parameters") <- list("means" = mean_clog)
#   xx
# }
