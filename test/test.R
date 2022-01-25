library(tidyverse)

sample_data <- read_csv("test/sample-germ-data.csv")
sample_priming_data <- read_csv("sample-priming-data.csv")

myGermData <- sample_data
myPrimingData <- sample_priming_data
save(myGermData, myPrimingData, file = "./data/sample_data.RData")

load("data/sample_data.RData")


## PBT MODELS ##

# thermaltime
calcTTSubOModel(sample_data)
calcTTSubOModel(sample_data, plot = F)
tt_results <- calcTTSubOModel(sample_data)
plotTTSubOModel(sample_data, tt_results)


# hydrotime
calcHTModel(sample_data)
calcHTModel(sample_data, plot = F)
ht_results <- calcHTModel(sample_data)
plotHTModel(sample_data, ht_results)


# hydro thermal time
calcHTTModel(sample_data)
calcHTTModel(sample_data, plot = F)
htt_results <- calcHTTModel(sample_data)
plotHTTModel(sample_data, htt_results)


## PRIMING MODELS ##

# hydro priming
sample_priming_data
calcSpeed(sample_priming_data, treatments = c("PrimingWP", "PrimingDuration"))
ht_speed_data <- calcSpeed(sample_priming_data, treatments = c("PrimingWP", "PrimingDuration"))
calcHPModel(ht_speed_data)
hp_results <- calcHPModel(ht_speed_data)
plotHPModel(speed_data, hp_results)



# hydro thermal priming
sample_priming_data
calcSpeed(sample_priming_data, treatments = c("PrimingWP", "PrimingTemp", "PrimingDuration"))
htp_speed_data <- calcSpeed(sample_priming_data, treatments = c("PrimingWP", "PrimingTemp", "PrimingDuration"))
calcHTPModel(htp_speed_data)
htp_results <- calcHTPModel(htp_speed_data)
plotHTPModel(htp_speed_data, htp_results)
