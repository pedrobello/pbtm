## code to prepare `myPrimingData` dataset goes here

library(readr)
sample_priming_data <- read_csv("data-raw/sample-priming-data.csv")

usethis::use_data(sample_priming_data, overwrite = TRUE)
