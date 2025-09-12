# Packages common to all the scripts
pkgs <- c(
  "dplyr",
  "readr",
  "purrr",
  "tidyr",
  "stringr"
)

new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) {
  install.packages(new_pkgs)
}
