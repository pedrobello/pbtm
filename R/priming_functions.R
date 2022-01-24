#' A Function to calculate the Hydropriming model parameters.
#'
#' This function calculates the minimal water potential for priming effects (Psib50min).
#' @param data Data rame with the calculated rates with treatments to be used in the Hydropriming model. The output of the calcSpeed function can be directly used here with the desired treatments. Priming water potential and priming duration treatments are required.
#' @param priming.wp Column containing the priming water potential treatments.
#' @param priming.duration Column containing the priming duration treatments.
#' @param rate Column containing the growth rate data from the calcSpeed function.
#' @param plot Should the model results be plotted?
#' @keywords hydropriming model parameters
#' @export
#' @examples

calcHPModel <- function(data, priming.wp = "PrimingWP", priming.duration = "PrimingDuration", rate = "GR50", plot = TRUE) {

  # data and argument checks
  if (!is.data.frame(data)) stop("Data is not a valid data frame.")

  # check validity of columns defs
  if (!is.element(priming.wp, names(data))) stop("Priming water potential column '", priming.wp, "' not found in data frame.")
  if (!is.element(priming.duration, names(data))) stop("Priming duration column '", priming.duration, "' not found in data frame.")
  if (!is.element(rate, names(data))) stop("Growth rate column '", rate, "' not found in data frame.")

  primingWP <- data[[priming.wp]]
  primingDuration <- data[[priming.duration]]
  growthRate <- data[[rate]]

  # initialize values for model calculation
  psiMin50 <- -1

  # HydroPriming - Calculate PsiMin50, y intercept and slope using NLS
  nls_model <- stats::nls(
    formula = growthRate ~ intercept + slope * ((primingWP - psiMin50) * primingDuration),
    start = list(
      intercept = 0.001,
      psiMin50 = -1,
      slope = 0.1),
    lower = list(
      intercept = 0.0000001,
      psiMin50 = -10,
      slope = 0.000000001),
    upper = list(
      intercept = 0.1,
      psiMin50 = -1,
      slope = 1),
    algorithm = "port")

  Intercept <- summary(nls_model)$coefficients[[1]]
  PsiMin50 <- summary(nls_model)$coefficients[[2]]
  Slope <- summary(nls_model)$coefficients[[3]]

  message("HydroPriming nonlinear least-squares model summary:")
  methods::show(summary(nls_model))

  #Hydropriming model linear regression
  l_model <- stats::lm(
    formula = growthRate ~ (primingWP - psiMin50) * primingDuration,
    data = data.frame(
      growthRate = growthRate,
      primingWP = primingWP,
      psiMin50 = PsiMin50,
      primingDuration = primingDuration
    ))

  message("HydroPriming linear model summary:")
  methods::show(summary(l_model))

  #Future refine of the model
  #lmer(GR50 ~ ((Treat.priming.wp-psiMin50)*Treat.priming.duration), PrimingTreats, start = c(psiMin50 = -1) )
  #fthetaHP <- function(PM50,Pwp,Pdur){(Pwp-PM50)*Pdur}
  #HPlModel <- lm(GR50 ~ fthetaHP(PsiMin50,Treat.priming.wp,Treat.priming.duration), PrimingTreats)
  #summary(HPlModel)

  #get some estimation of goodness of fit
  corr <- stats::cor(growthRate, stats::predict(l_model)) ^ 2
  RSquared <- round(corr[1], 2)

  results <- list(
    Type = "HydroPriming",
    Model = nls_model,
    LinearModel = l_model,
    Plot = NULL,
    Rate = rate,
    PsiMin50 = PsiMin50,
    Intercept = Intercept,
    Slope = Slope,
    RSquared = RSquared)

  if (plot == TRUE) {
    plt <- plotHPModel(
      data,
      results,
      priming.wp = priming.wp,
      priming.duration = priming.duration
    )
    results$Plot <- plt
    methods::show(plt)
  }

  return(results)
}


#' A Function to calculate the Hydrothermal priming model parameters.
#'
#' This function calculates the minimal water potential for priming effects (Psibmin50) and minimal temperature (Tmin).
#' @param data object with the calculated rates with treatments to be used in the Hydrothermal priming model. The output of the CalcSpeed function can be directly used here with the desired treatments. The fields with Treat.priming.wp, Treat.priming.temp and Treat.priming.duration need to be informed in the data file.
#' @param priming.wp Column containing the priming water potential treatments.
#' @param priming.temp Column containing the priming temperature treatments.
#' @param priming.duration Column containing the priming duration treatments.
#' @param rate Column containing the growth rate calculations from the calcSpeed function.
#' @param plot Should the model results be plotted?
#' @keywords hydrothermal priming model parameters
#' @export
#' @examples

calcHTPModel <- function(data, priming.wp = "PrimingWP", priming.temp = "PrimingTemp", priming.duration = "PrimingDuration", rate = "GR50", plot = TRUE) {

  # data and argument checks
  if (!is.data.frame(data)) stop("Data is not a valid data frame.")

  # check validity of columns defs
  if (!is.element(priming.wp, names(data))) stop("Priming water potential column '", priming.wp, "' not found in data frame.")
  if (!is.element(priming.temp, names(data))) stop("Priming temperature column '", priming.temp, "' not found in data frame.")
  if (!is.element(priming.duration, names(data))) stop("Priming duration column '", priming.duration, "' not found in data frame.")
  if (!is.element(rate, names(data))) stop("Growth rate column '", rate, "' not found in data frame.")

  primingWP <- data[[priming.wp]]
  primingTemp <- data[[priming.temp]]
  primingDuration <- data[[priming.duration]]
  growthRate <- data[[rate]]

  # HydroPriming - Calculate PsiMin50, y intercept and slope using NLS
  nls_model <- stats::nls(
    formula = growthRate ~ intercept + slope * ((primingWP - psiMin50) * (primingTemp - tMin) * primingDuration),
    start = list(
      intercept = 0.001,
      psiMin50 = -1,
      tMin = 12,
      slope = 0.1),
    lower = list(
      intercept = 0.0000001,
      psiMin50 = -10,
      tMmin = 0.5,
      slope = 0.000000001),
    upper = list(
      intercept = 0.1,
      psiMin50 = -1,
      tMin = 20,
      slope = 1),
    algorithm = "port")

  message("HydroThermalPriming nonlinear least-squares model summary:")
  methods::show(summary(nls_model))

  Intercept <- summary(nls_model)$coefficients[[1]]
  PsiMin50 <- summary(nls_model)$coefficients[[2]]
  Tmin <- summary(nls_model)$coefficients[[3]]
  Slope <- summary(nls_model)$coefficients[[4]]

  # Hydropriming model linear regression
  l_model <- stats::lm(
    formula = growthRate ~ ((primingWP - psiMin50) * (primingTemp - tMin) * primingDuration),
    data = data.frame(
      growthRate = growthRate,
      primingWP = primingWP,
      psiMin50 = PsiMin50,
      primingTemp = primingTemp,
      tMin = Tmin,
      primingDuration = primingDuration
    ))

  message("HydroThermalPriming linear model summary:")
  methods::show(summary(l_model))

  # get some estimation of goodness of fit
  corr <- stats::cor(growthRate, stats::predict(l_model)) ^ 2
  RSquared <- round(corr[1], 2)

  results <- list(
    Type = "HydroThermalPriming",
    Model = nls_model,
    LinearModel = l_model,
    Plot = NULL,
    Rate = rate,
    PsiMin50 = PsiMin50,
    Tmin = Tmin,
    Intercept = Intercept,
    Slope = Slope,
    RSquared = RSquared)

  if (plot == TRUE) {
    plt <- plotHTPModel(
      data,
      results,
      priming.wp = priming.wp,
      priming.temp = priming.temp,
      priming.duration = priming.duration
    )
    results$Plot <- plt
    methods::show(plt)
  }

  return(results)
}
