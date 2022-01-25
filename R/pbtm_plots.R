#' A Function to plot the selected calculated model and parameters and predicitions.
#'
#' This function plots the selected model and calculated parameters.
#' @param data A data frame containing time course and cumulative germination fractions to be used in the ThermalTime model. A column with time in hours, a column with cumulative fractions, and the experiment temperature are required.
#' @param model is the results list returned from from running any of the PBT models.
#' @param germ.temp Name of the column for the experimental temperature.
#' @param cum.time Name of the column for cumulative time.
#' @param cum.frac Name of the the column for cumulative fraction germinated.
#' @keywords plot population-based threshold model hydrotime
#' @importFrom rlang .data
#' @importFrom dplyr %>%
#' @export
#' @examples
#' "foo"

plotTTSubOModel <- function(data, model, germ.temp = "GermTemp", cum.time = "CumTime", cum.frac = "CumFraction") {

  modelName <- "ThermalTime Suboptimal"
  modelParams <- c("Type", "MaxCumFrac", "Tb", "ThetaT50", "Sigma", "Correlation")

  # data and argument checks
  if (!is.data.frame(data)) stop("Data is not a valid data frame.")
  if (!is.list(model)) stop("Model results must be in list format as output from any of the PBT model functions.")

  # check validity of columns defs
  if (!is.element(cum.time, names(data))) stop("Cumulative time column '", cum.time, "' not found in data frame.")
  if (!is.element(cum.frac, names(data))) stop("Cumulative fraction column '", cum.frac, "' not found in data frame.")

  # check for presence of model results
  lapply(modelParams, function(m) {
    if (!is.element(m, names(model))) stop("Required param '", m, "' missing from supplied model results.")
  })
  if (model$Type != modelName) stop("Model type must be '", modelName, "' for this plot function.")

  # model params
  maxCumFrac <- model$MaxCumFrac
  tb <- model$Tb
  thetaT50 <- model$ThetaT50
  sigma <- model$Sigma
  corr <- model$Correlation

  par1 <- paste("T[b] ==", round(tb, 1))
  par2 <- paste("ThetaT(50) ==", round(thetaT50, 3))
  par3 <- paste("sigma ==", round(sigma, 3))
  par4 <- paste("R^2 ==", round(corr, 2))

  # Plot all predicted treatments by the thermal time model
  df <- data %>% dplyr::distinct(.data[[germ.temp]], .keep_all = FALSE)

  modelLines <- mapply(function(temp) {
    ggplot2::stat_function(
      fun = function(x) {
        stats::pnorm(log(x, base = 10), thetaT50 - log(temp - tb, base = 10),  sigma, log = FALSE)
      },
      ggplot2::aes(color = as.factor(temp))
    )
  },
    df[[germ.temp]]
  )

  # generate the plot
  plt <- data %>%
    ggplot2::ggplot(ggplot2::aes(x = .data[[cum.time]], y = .data[[cum.frac]], color = as.factor(.data[[germ.temp]]))) +
    ggplot2::geom_point(shape = 19, size = 2) +
    modelLines +
    ggplot2::scale_y_continuous(labels = scales::percent, expand = c(0, 0), limits = c(0, 1.02)) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::expand_limits(x = 0, y = 0) +
    ggplot2::labs(
      title = paste(modelName, "Model"),
      x = "Time",
      y = "Cumulative fraction germinated (%)",
      color = "Temperature") +
    ggplot2::guides(color = ggplot2::guide_legend(reverse = T, order = 1)) +
    theme_scatter_plot +
    ggplot2::annotate("text", x = -Inf, y = 0.95, label = paste("Model Parameters"), color = "grey0", hjust = -0.1) +
    ggplot2::annotate("text", x = -Inf, y = 0.9, label = par1, color = "grey0", hjust = -0.2, parse = TRUE) +
    ggplot2::annotate("text", x = -Inf, y = 0.85, label = par2, color = "grey0", hjust = -0.1, parse = TRUE) +
    ggplot2::annotate("text", x = -Inf, y = 0.8, label = par3, color = "grey0", hjust = -0.2, parse = TRUE) +
    ggplot2::annotate("text", x = -Inf, y = 0.75, label = par4, color = "grey0", hjust = -0.2, parse = TRUE)

  return(plt)
}


#' A Function to plot the selected calculated model and parameters and predicitions.
#'
#' This function plots the selected model and calculated parameters.
#' @param data A data frame containing time course and cumulative germination fractions to be used in the ThermalTime model. A column with time in hours, a column with cumulative fractions, and the experiment temperature are required.
#' @param model is the results list returned from from running any of the PBT models.
#' @param germ.wp Name of the column for the water potential.
#' @param cum.time Name of the column for cumulative time.
#' @param cum.frac Name of the the column for cumulative fraction germinated.
#' @keywords plot population-based threshold model hydrotime
#' @importFrom rlang .data
#' @importFrom dplyr %>%
#' @export
#' @examples
#' "foo"

plotHTModel <- function(data, model, germ.wp = "GermWP", cum.time = "CumTime", cum.frac = "CumFraction") {

  modelName <- "HydroTime"
  modelParams <- c("Type", "MaxCumFrac", "HT", "Psib50", "Sigma", "Correlation")

  # data and argument checks
  if (!is.data.frame(data)) stop("Data is not a valid data frame.")
  if (!is.list(model)) stop("Model results must be in list format as output from any of the PBT model functions.")

  # check validity of columns defs
  if (!is.element(cum.time, names(data))) stop("Cumulative time column '", cum.time, "' not found in data frame.")
  if (!is.element(cum.frac, names(data))) stop("Cumulative fraction column '", cum.frac, "' not found in data frame.")

  # check for presence of model results
  lapply(modelParams, function(m) {
    if (!is.element(m, names(model))) stop("Required param '", m, "' missing from supplied model results.")
  })
  if (model$Type != modelName) stop("Model type must be '", modelName, "' for this plot function.")

  # set local vars
  maxCumFrac <- model$MaxCumFrac
  ht <- model$HT
  psib50 <- model$Psib50
  sigma <- model$Sigma
  corr <- model$Correlation

  # model params to display
  par1 <- paste("HT ==", round(ht, 2))
  par2 <- paste("Psi[b](50)==", round(psib50, 3))
  par3 <- paste("sigma == ", round(sigma, 3))
  par4 <- paste("R^2 == ", round(corr, 2))

  # plot all predicted treatments by the hydro time model
  df <- dplyr::distinct(data, .data[[germ.wp]], .keep_all = FALSE)

  modelLines <- mapply(function(wp) {
    ggplot2::stat_function(
      fun = function(x) {
        maxCumFrac * stats::pnorm(wp - (ht / x), psib50, sigma, log = FALSE)
      },
      ggplot2::aes(color = as.factor(wp))
    )
  },
    df[[germ.wp]]
  )

  # generate the plot
  plt <- data %>%
    ggplot2::ggplot(ggplot2::aes(x = .data[[cum.time]], y = .data[[cum.frac]], color = as.factor(.data[[germ.wp]]))) +
    ggplot2::geom_point(shape = 19, size = 2) +
    modelLines +
    ggplot2::scale_y_continuous(labels = scales::percent, expand = c(0, 0), limits = c(0, 1.02)) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::expand_limits(x = 0, y = 0) +
    ggplot2::labs(
      title = paste(modelName, "Model"),
      x = "Time",
      y = "Cumulative fraction germinated (%)",
      color = "Water Potential") +
    ggplot2::guides(color = ggplot2::guide_legend(reverse = T, order = 1)) +
    ggplot2::annotate("text", x = -Inf, y = 0.95, label = paste("Model Parameters"), color = "grey0", hjust = -0.1) +
    ggplot2::annotate("text", x = -Inf, y = 0.9, label = par1, color = "grey0", hjust = -0.2, parse = TRUE) +
    ggplot2::annotate("text", x = -Inf, y = 0.85, label = par2, color = "grey0", hjust = -0.1, parse = TRUE) +
    ggplot2::annotate("text", x = -Inf, y = 0.8, label = par3, color = "grey0", hjust = -0.2, parse = TRUE) +
    ggplot2::annotate("text", x = -Inf, y = 0.75, label = par4, color = "grey0", hjust = -0.2, parse = TRUE) +
    theme_scatter_plot

  return(plt)
}


#' A Function to plot the selected calculated model and parameters and predicitions.
#'
#' This function plots the selected model and calculated parameters.
#' @param data A data frame containing time course and cumulative germination fractions to be used in the ThermalTime model. A column with time in hours, a column with cumulative fractions, and the experiment temperature are required.
#' @param model is the results list returned from from running any of the PBT models.
#' @param germ.wp Name of the column for the water potential.
#' @param germ.temp Name of the column for the experimental temperature.
#' @param cum.time Name of the column for cumulative time.
#' @param cum.frac Name of the the column for cumulative fraction germinated.
#' @keywords plot population-based threshold model hydrotime
#' @importFrom rlang .data
#' @importFrom dplyr %>%
#' @export
#' @examples
#' "foo"

plotHTTModel <- function(data, model, germ.wp = "GermWP", germ.temp = "GermTemp", cum.time = "CumTime", cum.frac = "CumFraction") {

  modelName <- "HydroThermalTime"
  modelParams <- c("Type", "MaxCumFrac", "HT", "Psib50", "Sigma", "Correlation")

  # data and argument checks
  if (!is.data.frame(data)) stop("Data is not a valid data frame.")
  if (!is.list(model)) stop("Model results must be in list format as output from any of the PBT model functions.")

  # check validity of columns defs
  if (!is.element(cum.time, names(data))) stop("Cumulative time column '", cum.time, "' not found in data frame.")
  if (!is.element(cum.frac, names(data))) stop("Cumulative fraction column '", cum.frac, "' not found in data frame.")

  # check for presence of model results
  lapply(modelParams, function(m) {
    if (!is.element(m, names(model))) stop("Required param '", m, "' missing from supplied model results.")
  })
  if (model$Type != modelName) stop("Model type must be '", modelName, "' for this plot function.")

  # set local vars
  maxCumFrac <- model$MaxCumFrac
  ht <- model$HT
  psib50 <- model$Psib50
  tb <- model$Tb
  sigma <- model$Sigma
  corr <- model$Correlation

  # model params
  par1 <- paste("HT ==", round(ht, 2))
  par2 <- paste("T[b]==", round(tb, 2))
  par3 <- paste("psi[b](50)==", round(psib50,3))
  par4 <- paste("sigma == ", round(sigma, 3))
  par5 <- paste("R^2 == ", round(corr, 2))

  # function to plot all predicted treatments by the hydro thermal time model
  df <- data %>%
    dplyr::distinct(.data[[germ.wp]], .data[[germ.temp]], .keep_all = F) %>%
    dplyr::arrange(.data[[germ.wp]], .data[[germ.temp]])

  modelLines <- mapply(function(wp, temp) {
    ggplot2::stat_function(
      fun = function(x) {
        maxCumFrac * stats::pnorm(
          wp - (ht / ((temp - tb) * x)),
          psib50,
          sigma,
          log = FALSE
        )
      },
      ggplot2::aes(color = as.factor(wp), alpha = as.factor(temp))
    )
  },
    df[[germ.wp]],
    df[[germ.temp]]
  )

  # generate the plot
  plt <- data %>%
    ggplot2::ggplot(ggplot2::aes(
      x = .data[[cum.time]],
      y = .data[[cum.frac]],
      color = as.factor(.data[[germ.wp]]),
      alpha = as.factor(.data[[germ.temp]]))) +
    ggplot2::geom_point(size = 2) +
    modelLines +
    ggplot2::scale_y_continuous(labels = scales::percent, expand = c(0, 0), limits = c(0, 1.02)) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::expand_limits(x = 0, y = 0) +
    ggplot2::labs(
      title = paste(modelName, "Model"),
      x = "Time",
      y = "Cumulative fraction germinated (%)",
      color = "Water Potential",
      alpha = "Temperature") +
    ggplot2::guides(color = ggplot2::guide_legend(reverse = T, order = 1)) +
    ggplot2::annotate("text", x = -Inf, y = 0.95, label = paste("Model Parameters"), color = "grey0", hjust = -0.1) +
    ggplot2::annotate("text", x = -Inf, y = 0.9, label = par1, color = "grey0", hjust = -0.2, parse = TRUE) +
    ggplot2::annotate("text", x = -Inf, y = 0.85, label = par2, color = "grey0", hjust = -0.1, parse = TRUE) +
    ggplot2::annotate("text", x = -Inf, y = 0.8, label = par3, color = "grey0", hjust = -0.2, parse = TRUE) +
    ggplot2::annotate("text", x = -Inf, y = 0.75, label = par4, color = "grey0", hjust = -0.2, parse = TRUE) +
    ggplot2::annotate("text", x = -Inf, y = 0.7, label = par5, color = "grey0", hjust = -0.2, parse = TRUE) +
    theme_scatter_plot

  return(plt)
}
