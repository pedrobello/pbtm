## code to prepare `myGermData` dataset goes here

library(readr)
myGermData <- read_csv("data-raw/sample-germ-data.csv")

usethis::use_data(myGermData, overwrite = TRUE)
