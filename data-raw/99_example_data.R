example_data <- dplyr::tibble(
  sample = paste0("sample_", sprintf("%02d", 1:10)),
  label_01 = c(0, 450000, 1300, 5, 4893, 0, 0, 28000, 0, 64045),
  label_02 = c(89000, 0, 0, 5400, 0, 0, 3645, 0, 0, 0),
  label_03 = c(0, 0, 67000, 0, 237, 4200, 0, 20, 93000, 0),
  label_04 = c(0, 270, 0, 81078, 0, 0, 10, 5800, 0, 0),
  label_05 = c(88207, 67, 0, 1000, 0, 75000, 0, 1000, 0, 84720)
)

usethis::use_data(example_data, overwrite = TRUE)
