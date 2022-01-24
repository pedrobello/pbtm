
#' A Function to calculate the ThermalTime model parameters. ------------------------
#'
#' This function calculates the temperature base (Tb), the thermal time to 50% of the population (ThetaT50) and the standard deviation (sigma).
#' @param data Time course and cumulative dataset to be used in the Thermaltime model. The original dataframe template should be used or column names should be modified similarly to the template. A column with time in hours (CumTime) + a column with cumulative fractions (CumFraction) and the experiment temperature (Germ.temp) are required. Filter the dataframe to only have treatments with temperature equal or under to  optimal temperature level.
#' @param germ.temp Column containing germination temperature values.
#' @param cum.time Column containing cumulative elapsed time.
#' @param cum.frac Column containing cumulative fraction germinated.
#' @param max.cum.frac Sets the ceiling cumulative fraction for the model when treatment at optimal condition displays a lower maximum cumulative fraction. Use it on your own discretion.
#' @param plot Should the model results be plotted?
#' @keywords Thermal time model parameters
#' @export
#' @examples calcTTSubOModel(MyData)
#' calcTTSubOModel(MyData)

calcTTSubOModel <- function(data, germ.temp = "GermTemp", cum.time = "CumTime", cum.frac = "CumFraction", max.cum.frac = 1, plot = TRUE) {

  # data and argument checks
  if (!is.data.frame(data)) stop("Data is not a valid data frame.")

  # check validity of columns defs
  if (!is.element(germ.temp, names(data))) stop("Germination temperature column '", germ.temp, "' not found in data frame.")
  if (!is.element(cum.time, names(data))) stop("Cumulative time column '", cum.time, "' not found in data frame.")
  if (!is.element(cum.frac, names(data))) stop("Cumulative fraction column '", cum.frac, "' not found in data frame.")

  # check validity of fraction argument
  if (!is.numeric(max.cum.frac)) stop("Non-numeric fraction value specified.")
  if (max.cum.frac < 0 || max.cum.frac > 1) stop("Fraction value ", max.cum.frac, " is out of range, must be between 0 and 1.")

  temp <- data[[germ.temp]]
  time <- data[[cum.time]]
  germ <- data[[cum.frac]]

  # Calculate Thermaltime Suboptimal Model Parameters- nls plus algorithm port used to add constraints on the parameters
  model <- stats::nls(
    formula = germ ~ max.cum.frac * stats::pnorm(
      log(time, base = 10),
      mean = thetaT50 - log(temp - Tb, base = 10),
      sd = sigma,
      log = FALSE),
    start = list(
      Tb = 6,
      thetaT50 = 3,
      sigma = 0.09),
    lower = list(
      Tb = 0,
      thetaT50 = 0.5,
      sigma = 0.0001),
    upper = list(
      Tb = 15,
      thetaT50 = 50,
      sigma = 0.5),
    algorithm = "port")

  # get some estimation of goodness of fit
  Corr <- stats::cor(germ, stats::predict(model)) ^ 2

  # passing fitted Hydrotime Model Parameters
  Tb <- summary(model)$coefficients[[1]]
  ThetaT50 <- summary(model)$coefficients[[2]]
  Sigma <- summary(model)$coefficients[[3]]

  message("ThermalTime Suboptimal Model Summary:")
  show(summary(model))

  results <- list(
    Type = "ThermalTime Suboptimal",
    Model = model,
    Plot = NULL,
    MaxCumFrac = max.cum.frac,
    Tb = Tb,
    ThetaT50 = ThetaT50,
    Sigma = Sigma,
    Correlation = Corr
  )

  if (plot == TRUE) {
    plt <- plotTTSubOModel(
      data,
      results,
      germ.temp = germ.temp,
      cum.time = cum.time,
      cum.frac = cum.frac
    )
    results$Plot <- plt
    show(plt)
  }

  return(results)
}



#' A Function to calculate the Hydrotime model parameters.
#'
#' This function calculates the hydrotime constant (HT), the median water potential base (Psib50) and the standard deviation (sigma).
#' @param data time course and cumulative dataset to be used in the hydrotime model. The original dataframe template should be used or column names should be modified similarly to the template. A column with time in hours CumTime) + a column with cumulative fractions (CumFraction) and the experiment water potential (Germ.wp) are required. Filter the dataframe to only have treatments with water potential at the same temperature level.
#' @param germ.wp Column containing germination water potential.
#' @param cum.time Column containing cumulative elapsed time.
#' @param cum.frac Column containing cumulative fraction germinated.
#' @param max.cum.frac Sets the ceiling cumulative fraction for the model when treatment at optimal condition displays a lower maximum cumulative fraction. Use it on your own discretion.
#' @param plot Should the model results be plotted?
#' @keywords hydrotime model parameters
#' @export
#' @examples calcHTModel(MyData)
#' calcHTModel(MyData)

calcHTModel <- function(data, germ.wp = "GermWP", cum.time = "CumTime", cum.frac = "CumFraction", max.cum.frac = 1, plot = TRUE) {

  # data and argument checks
  if (!is.data.frame(data)) stop("Data is not a valid data frame.")

  # check validity of columns defs
  if (!is.element(germ.wp, names(data))) stop("Germination water potential column '", germ.wp, "' not found in data frame.")
  if (!is.element(cum.time, names(data))) stop("Cumulative time column '", cum.time, "' not found in data frame.")
  if (!is.element(cum.frac, names(data))) stop("Cumulative fraction column '", cum.frac, "' not found in data frame.")

  # check validity of fraction argument
  if (!is.numeric(max.cum.frac)) stop("Non-numeric fraction value specified.")
  if (max.cum.frac < 0 || max.cum.frac > 1) stop("Fraction value ", max.cum.frac, " is out of range, must be between 0 and 1.")

  wp <- data[[germ.wp]]
  time <- data[[cum.time]]
  germ <- data[[cum.frac]]

  # Calculate Hydrotime Model Parameters- nls plus algorithm port used to add constraints on the parameters
  model <- stats::nls(
    formula = germ ~ max.cum.frac * stats::pnorm(
      wp - (HT / time),
      Psib50,
      Sigma,
      log = FALSE),
    start = list(
      HT = 60,
      Psib50 = -0.8,
      Sigma = 0.2),
    lower = list(
      HT = 1,
      Psib50 = -5,
      Sigma = 0.0001),
    upper = list(
      HT = 1000,
      Psib50 = -0.000000001,
      Sigma = 2),
    algorithm = "port")

  #get some estimation of goodness of fit
  corr <- cor(germ, stats::predict(model)) ^ 2

  # Passing fitted Hydrotime Model Parameters
  HT <- summary(model)$coefficients[[1]]
  Psib50 <- summary(model)$coefficients[[2]]
  Sigma <- summary(model)$coefficients[[3]]

  message("HydroTime Model Summary:")
  show(summary(model))

  results <- list(
    Type = "HydroTime",
    Model = model,
    Plot = NULL,
    MaxCumFrac = max.cum.frac,
    HT = HT,
    Psib50 = Psib50,
    Sigma = Sigma,
    Correlation = corr)

  if (plot == TRUE) {
    plt <- plotHTModel(
      data,
      results,
      germ.wp = germ.wp,
      cum.time = cum.time,
      cum.frac = cum.frac
    )
    results$Plot <- plt
    show(plt)
  }

  return(results)
}



#' A Function to calculate the Hydrothermal time model parameters.
#'
#' This function calculates the hydrohermal time constant (HT), the temperature base (tb), the median water potential base (Psib50) and the standard deviation (sigma).
#' @param data time course and cumulative dataset to be used in the hydrothermal time model. The original dataframe template should be used or column names should be modified similarly to the template. A column with time (CumTime) + a column with cumulative fractions (CumFraction), the experimental water potential (Germ.wp) and the temperature (Germ.temp) are required.
#' @param germ.wp Column containing germination water potential.
#' @param germ.temp Column containing germination temperature.
#' @param cum.time Column containing cumulative elapsed time.
#' @param cum.frac Column containin cumulative fraction germinated.
#' @param max.cum.frac Sets the ceiling cumulative fraction for the model when treatment at optimal condition displays a lower maximum cumulative fraction. Use it on your own discretion.
#' @param base.temp is the temperature base that can be provided and have a fixed value or left in blank to be calculated with other hydrothermal parameters.
#' @param plot Should the model results be plotted?
#' @keywords hydrothermal time model parameters
#' @export
#' @examples calcHTTModel(MyData)
#' calcHTTModel(MyData)

calcHTTModel <- function(data, germ.wp = "GermWP", germ.temp = "GermTemp", cum.time = "CumTime", cum.frac = "CumFraction", max.cum.frac = 1, base.temp = NULL, plot = TRUE) {

  # data and argument checks
  if (!is.data.frame(data)) stop("Data is not a valid data frame.")

  # check validity of columns defs
  if (!is.element(germ.wp, names(data))) stop("Germination water potential column '", germ.wp, "' not found in data frame.")
  if (!is.element(germ.temp, names(data))) stop("Germination temperature column '", germ.temp, "' not found in data frame.")
  if (!is.element(cum.time, names(data))) stop("Cumulative time column '", cum.time, "' not found in data frame.")
  if (!is.element(cum.frac, names(data))) stop("Cumulative fraction column '", cum.frac, "' not found in data frame.")

  # check validity of fraction argument
  if (!is.numeric(max.cum.frac)) stop("Non-numeric fraction value specified.")
  if (max.cum.frac < 0 || max.cum.frac > 1) stop("Fraction value ", max.cum.frac, " is out of range, must be between 0 and 1.")

  # check validity of base temperature
  if (!is.null(base.temp) && !is.numeric(base.temp) && !length(base.temp) == 1) stop("Invalid base temperature '", base.temp, "'.")

  wp <- data[[germ.wp]]
  temp <- data[[germ.temp]]
  time <- data[[cum.time]]
  germ <- data[[cum.frac]]

  # model conditions
  start <- list(
    HT = 800,
    psib50 = -1,
    sigma = 0.4)
  lower <- list(
    HT = 1,
    psib50 = -5,
    sigma = 0.0001)
  upper <- list(
    HT = 5000,
    psib50 = 0,
    sigma = 10)

  if (is.null(base.temp)) {
    start$Tb <- 1
    lower$Tb <- 0
    upper$Tb <- 15
  } else {
    Tb <- base.temp
  }

  # Calculate Hydrotime Model Parameters and Tb - nls plus algorithm port used to add constraints on the parameters
  model <- stats::nls(
    germ ~ max.cum.frac * stats::pnorm(
      wp - (HT / ((temp - Tb) * time)),
      psib50,
      sigma,
      log = FALSE),
    start = start,
    lower = lower,
    upper = upper,
    algorithm = "port")

  # get model coefficients
  if (is.null(base.temp)) {
    HT <- summary(model)$coefficients[[1]]
    Tb <- summary(model)$coefficients[[2]]
    Psib50 <- summary(model)$coefficients[[3]]
    Sigma <- summary(model)$coefficients[[4]]
  } else {
    HT <- summary(model)$coefficients[[1]]
    Tb <- base.temp
    Psib50 <- summary(model)$coefficients[[2]]
    Sigma <- summary(model)$coefficients[[3]]
  }

  # estimate goodness of fit
  corr <- stats::cor(germ, stats::predict(model)) ^ 2

  message("HydroThermalTime model summary:")
  show(summary(model))

  results <- list(
    Type = "HydroThermalTime",
    Model = model,
    Plot = NULL,
    MaxCumFrac = max.cum.frac,
    HT = HT,
    Tb = Tb,
    Psib50 = Psib50,
    Sigma = Sigma,
    Correlation = corr)

  if (plot == TRUE) {
    plt <- plotHTTModel(
      data,
      results,
      germ.wp = germ.wp,
      germ.temp = germ.temp,
      cum.time = cum.time,
      cum.frac = cum.frac
    )
    results$Plot <- plt
    show(plt)
  }

  return(results)
}



# HTTModel test stuff -----------------------------------------------------

# Hydrothermal time Model - Create table to plot treatments with predicted model lines
# TreatData$WPFactor <<- with(TreatData, (as.factor(TreatData$Germ.wp)))
# Factor1 <<- TreatData$WPFactor
# Factor1Title <<- "Water \n Potential"
# TreatmentsWP <<- distinct(TreatData, Germ.wp, .keep_all = FALSE)
# TreatData$TempFactor <<- with(TreatData, (as.factor(TreatData$Germ.temp)))
# Factor2 <<- TreatData$TempFactor
# Factor2Title <<- "Temperature"
# TreatmentsTemp <<- distinct(TreatData, Germ.temp, .keep_all = FALSE)
# TreatmentHTT <<- distinct(TreatData, Germ.temp, Germ.wp, .keep_all = FALSE)
#
#
# Passing fitted Hydrotime Model Parameters for plot legend
# ModPar1Label <<- "HT =="
# ModPar2Label <<- "T[b]=="
# ModPar3Label <<- "psi[b](50)=="
# ModPar4Label <<- "sigma == "
# ModPar5Label <<- "R^2 == "
#
# ModPar1 <<- round(summary(HTTModel)$coefficients[[1]],2)
# ModPar2 <<- round(Tb[1],2)
# ModPar3 <<- round(psib50[1],3)
# ModPar4 <<- round(sigma[1],3)
# ModPar5 <<- round(Correlation[1],2)
#
# Dt <<- data.frame(Treatments$Germ.temp,Treatments$Germ.wp)
# tmps <<- Dt$Treatments.Germ.temp
# wps <<- Dt$Treatments.Germ.wp
#
# Function to plot all predicted treatments by the HYDROTHERMAL time model
# modellines <<-
#   mapply(function(Temp1, WP1) {
#     stat_function(fun=function(x){pnorm((+WP1-(HT/((Temp1-Tb)*x))),Psib50,Sigma, log= FALSE)*MaxGerm}, aes_(colour = factor(Temp1), alpha = factor(WP1)))
#   }, tmps, wps)
#
# Create columns on TreatDataClean with predicted and Normalized time values from model
# TreatDataClean <<-TreatDataClean %>% as_tibble() %>% mutate(
#   ModelPredicted = pnorm((+Germ.wp-(HT/((Germ.temp-Tb)*Germ.time.hours))),Psib50,Sigma, log= FALSE)*MaxGerm,
#   NormalizedTime = +((1-Germ.wp/(+Germ.wp-(HT/(Germ.time.hours*(Germ.temp-Tb)))))*Germ.time.hours*(Germ.temp-Tb))
# )
#
# Function to plot Normalized time using the HYDROTHERMAL time model
# modelNormalized <<-
#   alply(as.matrix(TreatmentsTemp), 1, function(Temp) {
#     stat_function(fun=function(x){pnorm((+0-(HT/((1)*x))),Psib50,Sigma, log= FALSE)*MaxGerm}, color="blue3")
#   })
#  NormalizedAxisTitle <<- "Normalized thermal time (°h)"
