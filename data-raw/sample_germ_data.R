## code to prepare `myGermData` dataset goes here

library(readr)
sample_germ_data <- read_csv("data-raw/sample-germ-data.csv")

usethis::use_data(sample_germ_data, overwrite = TRUE)
