#-----------------Tested Functions - MAc OS and Windows
#
#' A Function to plot the both priming models.
#'
#' This function plots the priming models and calculated parameters.
#' @param data object with the calculated rates with treatments to be used in the Hydrothermal priming model. The output of the CalcSpeed function can be directly used here with the desired treatments. The fields with Treat.priming.wp, Treat.priming.temp and Treat.priming.duration need to be informed in the data file.
#' @param model is data object resulting from the CalcHPModel() or CalcHTPModel() functions containing the model information and parameter results.
#' @param priming.wp Column name for the priming water potential.
#' @param priming.duration Column name for the duration of priming.
#' @keywords plot priming model hydropriming hydrothermal priming
#' @importFrom rlang .data
#' @importFrom dplyr %>%
#' @export
#' @examples
#' "foo"

plotHPModel <- function(data, model, priming.wp = "PrimingWP", priming.duration = "PrimingDuration") {

  modelName <- "HydroPriming"
  modelParams <- c("Type", "Rate", "PsiMin50", "Intercept", "Slope")

  # data and argument checks
  if (!is.data.frame(data)) stop("Data is not a valid data frame.")
  if (!is.list(model)) stop("Model results must be in list format as output from any of the PBT model functions.")

  # check validity of columns defs

  if (!is.element(priming.wp, names(data))) stop("Priming water potential column '", priming.wp, "' not found in data frame.")
  if (!is.element(priming.duration, names(data))) stop("Priming duration column '", priming.duration, "' not found in data frame.")

  # check for presence of model results
  lapply(modelParams, function(m) {
    if (!is.element(m, names(model))) stop("Required param '", m, "' missing from supplied model results.")
  })
  if (model$Type != modelName) stop("Model type must be '", modelName, "' for this plot function.")
  if (!is.element(model$Rate, names(data))) stop("Model rate '", model$Rate, "' does not match a rate column from the supplied speed data.")

  # set local vars
  rate <- model$Rate
  psimin50 <- model$PsiMin50
  intercept <- model$Intercept
  slope <- model$Slope
  r2 <- model$RSquared

  data <- data %>%
    dplyr::mutate(Theta = (.data[[priming.wp]] - psimin50) * .data[[priming.duration]])

  # model params
  par1 <- paste("Psi[min](50)==", psimin50)
  par2 <- paste("Intercept==", intercept)
  par3 <- paste("Slope==", slope)
  par4 <- paste("R^2==", r2)

  # generate plot
  plt <- data %>%
    ggplot2::ggplot(ggplot2::aes(
      x = .data$Theta,
      y = .data[[rate]],
      color = as.factor(.data[[priming.wp]]),
      shape = as.factor(.data[[priming.duration]]))) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_abline(intercept = intercept, slope = slope, color = "blue") +
    ggplot2::scale_x_continuous(expand = c(0,0)) +
    ggplot2::labs(
      title = paste(modelName, "Model"),
      x = "Hydro Priming Time",
      y = "Germination Rate",
      color = "Priming Water Potential",
      shape = "Priming Duration"
    ) +
    ggplot2::annotate("text", x = -Inf, y = Inf, label = "Model Parameters", color = "grey0", hjust = -0.1, vjust = 1.5) +
    ggplot2::annotate("text", x = -Inf, y = Inf, label = par1, color = "grey0", parse = TRUE, hjust = -0.09, vjust = 2.5) +
    ggplot2::annotate("text", x = -Inf, y = Inf, label = par2, color = "grey0", parse = TRUE, hjust = -0.12, vjust = 4.5) +
    ggplot2::annotate("text", x = -Inf, y = Inf, label = par3, color = "grey0", parse = TRUE, hjust = -0.11, vjust = 5.8) +
    ggplot2::annotate("text", x = -Inf, y = Inf, label = par4, color = "grey0", parse = TRUE, hjust = -0.2, vjust = 6.3) +
    theme_scatter_plot

  return(plt)
}



#' A Function to plot the both priming models.
#'
#' This function plots the priming models and calculated parameters.
#' @param data object with the calculated rates with treatments to be used in the Hydrothermal priming model. The output of the CalcSpeed function can be directly used here with the desired treatments.
#' @param model is data object resulting from the calcHPModel or calcHTPModel functions containing the model information and parameter results.
#' @param priming.wp Column name for the priming water potential.
#' @param priming.temp Column name for the priming temperature.
#' @param priming.duration Column name for the duration of priming.
#' @keywords plot priming model hydropriming hydrothermal priming
#' @importFrom rlang .data
#' @importFrom dplyr %>%
#' @export
#' @examples
#' "foo"

plotHTPModel <- function(data, model, priming.wp = "PrimingWP", priming.temp = "PrimingTemp", priming.duration = "PrimingDuration") {

  modelName <- "HydroThermalPriming"
  modelParams <- c("Type", "Rate", "PsiMin50", "Tmin", "Intercept", "Slope")

  # data and argument checks
  if (!is.data.frame(data)) stop("Data is not a valid data frame.")
  if (!is.list(model)) stop("Model results must be in list format as output from any of the PBT model functions.")

  # check validity of columns defs

  if (!is.element(priming.wp, names(data))) stop("Priming water potential column '", priming.wp, "' not found in data frame.")
  if (!is.element(priming.temp, names(data))) stop("Priming temperature column '", priming.temp, "' not found in data frame.")
  if (!is.element(priming.duration, names(data))) stop("Priming duration column '", priming.duration, "' not found in data frame.")

  # check for presence of model results
  lapply(modelParams, function(m) {
    if (!is.element(m, names(model))) stop("Required param '", m, "' missing from supplied model results.")
  })
  if (model$Type != modelName) stop("Model type must be '", modelName, "' for this plot function.")
  if (!is.element(model$Rate, names(data))) stop("Model rate '", model$Rate, "' does not match a rate column from the supplied speed data.")

  # set local vars
  rate <- model$Rate
  psimin50 <- model$PsiMin50
  tmin <- model$Tmin
  intercept <- model$Intercept
  slope <- model$Slope
  r2 <- model$RSquared

  data <- data %>%
    dplyr::mutate(Theta = (.data[[priming.wp]] - psimin50) * (.data[[priming.temp]] - tmin) * .data[[priming.duration]])

  # model params
  par1 <- paste("Psi[min](50)==", round(psimin50, 3))
  par2 <- paste("T[min]==", round(tmin, 2))
  par3 <- paste("Intercept==", round(intercept, 4))
  par4 <- paste("Slope==", round(slope, 6))
  par5 <- paste("R^2==", round(r2, 2))

  # generate plot
  plt <- data %>%
    ggplot2::ggplot(ggplot2::aes(
      x = .data$Theta,
      y = .data[[rate]],
      color = as.factor(.data[[priming.wp]]),
      shape = as.factor(.data[[priming.duration]]),
      alpha = as.factor(.data[[priming.temp]]))) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_abline(intercept = intercept, slope = slope, color = "blue") +
    ggplot2::scale_x_continuous(expand = c(0,0)) +
    ggplot2::labs(
      title = paste(modelName, "Model"),
      x = "Hydro Thermal Priming Time",
      y = "Germination Rate",
      color = "Priming Water Potential",
      shape = "Priming Duration",
      alpha = "Priming Temperature"
    ) +
    ggplot2::annotate("text", x = -Inf, y = Inf, label = "Model Parameters", color = "grey0", hjust = -0.1, vjust = 1.5) +
    ggplot2::annotate("text", x = -Inf, y = Inf, label = par1, color = "grey0", parse = TRUE, hjust = -0.09, vjust = 2.5) +
    ggplot2::annotate("text", x = -Inf, y = Inf, label = par2, color = "grey0", parse = TRUE, hjust = -0.12, vjust = 4.5) +
    ggplot2::annotate("text", x = -Inf, y = Inf, label = par3, color = "grey0", parse = TRUE, hjust = -0.11, vjust = 5.8) +
    ggplot2::annotate("text", x = -Inf, y = Inf, label = par4, color = "grey0", parse = TRUE, hjust = -0.1, vjust = 7.2) +
    ggplot2::annotate("text", x = -Inf, y = Inf, label = par5, color = "grey0", parse = TRUE, hjust = -0.2, vjust = 7.8) +
    theme_scatter_plot

  return(plt)
}
