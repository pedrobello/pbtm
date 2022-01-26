## code to prepare `myPrimingData` dataset goes here

library(readr)
myPrimingData <- read_csv("data-raw/sample-priming-data.csv")

usethis::use_data(myPrimingData, overwrite = TRUE)
