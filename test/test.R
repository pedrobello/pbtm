library(tidyverse)

sample_data <- read_csv("test/sample-germ-data.csv")
sample_priming_data <- read_csv("sample-priming-data.csv")

myGermData <- sample_data
myPrimingData <- sample_priming_data
save(myGermData, myPrimingData, file = "./data/sample_data.RData")

load("data/sample_data.RData")



## BASIC FUNCTIONS ##

# calcSpeed

calcSpeed(sample_data, treatments = c("GermWP", "GermTemp"), f = 0.5)
calcSpeed(sample_data, treatments = c("GermWP", "GermTemp"), f = c(0.25, 0.5))
calcSpeed(sample_data, treatments = c("GermWP", "GermTemp"), f = c(0.10, 0.16, 0.5))

f <- 0.5

f <- c(0.1, 0.5, 0.3, 0.3)
f <- sort(unique(f))


# when 1 fraction
sample_data %>%
  dplyr::group_by_at(c("GermWP", "GermTemp")) %>%
  dplyr::summarise(
    n = dplyr::n(),
    !!eval(paste0("T", f * 100)) := stats::approx(.data$CumFraction, .data$CumTime, xout = f, ties = "ordered", rule = 2)$y,
    !!eval(paste0("GR", f * 100)) := 1 / .data[[paste0("T", f * 100)]],
    .groups = "drop"
  )

# when multiple fractions
sample_data %>%
  dplyr::group_by_at(c("GermWP", "GermTemp")) %>%
  summarise(
    {
      df <- approx(CumFraction, CumTime, xout = f, ties = "ordered", rule = 2)
      names(df) <- c("Frac", "T")
      df <- as_tibble(df)
      drop_na(df)
    },
    .groups = "drop"
  ) %>%
  mutate(
    GR = 1 / T,
    Frac = Frac * 100) %>%
  pivot_longer(
    cols = c("T", "GR")
  ) %>%
  mutate(name = paste0(name, Frac)) %>%
  select(-"Frac") %>%
  pivot_wider()


for (f in fraction) {
  print(f)
}

sample_data %>%
  dplyr::group_by_at(c("GermWP", "GermTemp")) %>%
  dplyr::summarise(
    n = dplyr::n(),
    fracs = fraction
  )


build_approx <- function(x, y, f) {
  tibble(
    !!eval(paste0("T", f * 100)) := approx(x, y, f, ties = "ordered", rule = 2)$y
  )
}

build_approx(sample_data$CumFraction, sample_data$CumTime, 0.5)

build_approx <- function(x, y, f) {
  tibble(
    !!eval(paste0("T", f * 100)) := approx(x, y, f, ties = "ordered", rule = 2)$y
  )
}


# mergeTrts

sample_data %>%
  mutate(GermTemp = as.factor(GermTemp)) %>%
  ggplot(aes(x = CumTime, y = CumFraction, color = GermTemp)) +
  geom_smooth(method = "glm", method.args = list(family = "quasibinomial"), se = F) +
  geom_point(data = sample_data, size = 1, color = "grey") +
  geom_point()

sample_data %>%
  filter(GermWP == 0) %>%
  mutate(GermTemp = as.factor(GermTemp)) %>%
  ggplot(aes(x = CumTime, y = CumFraction, color = GermTemp)) +
  geom_smooth(method = "glm", method.args = list(family = "quasibinomial"), se = F) +
  geom_point(data = sample_data, size = 1, color = "grey") +
  geom_point()


# merge on temperature
sample_data %>%
  mergeTrts(keep.trts = "GermTemp") %>%
  mutate(GermTemp = as.factor(GermTemp)) %>%
  ggplot(aes(x = CumTime, y = CumFraction, color = GermTemp)) +
  geom_smooth(method = "loess", se = F) +
  geom_point(data = sample_data, size = 1, color = "grey") +
  geom_point()

# merge on water potential
sample_data %>%
  mergeTrts(keep.trts = "GermWP") %>%
  mutate(GermWP = as.factor(GermWP)) %>%
  ggplot(aes(x = CumTime, y = CumFraction, color = GermWP)) +
  geom_smooth(method = "loess", se = F) +
  geom_point(data = sample_data, size = 1, color = "grey") +
  geom_point()



# calcCumFrac
sample_seed_data <- read_csv("test/sample-seed-data.csv")

sample_seed_data %>%
  group_by(TrtID) %>%
  arrange(CumTime, nGerminated) %>%
  distinct(CumTime, .keep_all = T) %>%
  mutate(CumFraction = nGerminated / nTotal) %>%
  select(-c("nTotal", "nGerminated")) %>%
  ggplot(aes(x = CumTime, y = CumFraction, color = GermTemp, group = TrtID)) +
  geom_smooth(method = "glm", method.args = list(family = "quasibinomial"), se = F) +
  geom_point()

calcCumFrac(sample_seed_data)

calcCumFrac(sample_seed_data) %>%
  ggplot(aes(x = CumTime, y = CumFraction, color = as.factor(GermTemp), group = TrtID)) +
  geom_smooth(method = "glm", method.args = list(family = "quasibinomial"), se = F) +
  geom_point()




## BASIC PLOTS ##

speed_data <- calcSpeed(sample_data, treatments = "TrtDesc")
plotRateVsTrt(speed_data, x = "TrtDesc")



## PBT MODELS ##

# thermaltime
calcTTSubOModel(sample_data)

calcTTSubOModel(sample_data, plot = F)
tt_results <- calcTTSubOModel(sample_data)
plotTTSubOModel(sample_data, tt_results)

sample_data %>%
  filter(GermWP == 0) %>%
  calcTTSubOModel()

tt_results <- sample_data %>%
  filter(GermWP == 0) %>%
  calcTTSubOModel()
tt_results$Plot

sample_data %>%
  mergeTrts(keep.trts = "GermTemp") %>%
  calcTTSubOModel()



# hydrotime
calcHTModel(sample_data)

sample_data %>%
  filter(GermTemp == 25) %>%
  calcHTModel()

calcHTModel(sample_data, plot = F)
ht_results <- calcHTModel(sample_data)
plotHTModel(sample_data, ht_results)

sample_data %>%
  mergeTrts(keep.trts = "GermWP") %>%
  calcHTModel()


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

sample_priming_data %>%
  calcSpeed(treatments = c("PrimingWP", "PrimingDuration")) %>%
  calcHPModel()

sample_priming_data %>%
  mergeTrts(keep.trts = c("PrimingWP", "PrimingDuration")) %>%
  calcSpeed(treatments = c("PrimingWP", "PrimingDuration")) %>%
  calcHPModel()


# should there be a thermal priming model?


# hydro thermal priming
sample_priming_data
calcSpeed(sample_priming_data, treatments = c("PrimingWP", "PrimingTemp", "PrimingDuration"))
htp_speed_data <- calcSpeed(sample_priming_data, treatments = c("PrimingWP", "PrimingTemp", "PrimingDuration"))
calcHTPModel(htp_speed_data)
htp_results <- calcHTPModel(htp_speed_data)
plotHTPModel(htp_speed_data, htp_results)


sample_priming_data %>%
  calcSpeed(treatments = c("PrimingWP", "PrimingTemp", "PrimingDuration")) %>%
  calcHTPModel()
