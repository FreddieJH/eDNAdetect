source("data-raw/00_pkgs.R")

clean_load <- function(file) {
  x <- readr::read_csv(
    paste0("data-raw/files/", file),
    show_col_types = FALSE
  )
  cat(paste(file, "loaded successfully\n"))
  return(x)
}

label_raw <- clean_load("tax_table_full.csv")
data_raw <- clean_load("otu_table_full.csv")
catch_raw <- clean_load("Catch_data_full.csv")

rm(clean_load)
