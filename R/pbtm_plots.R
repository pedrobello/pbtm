#' A Function to plot the selected calculated model and parameters and predicitions.
#'
#' This function plots the selected model and calculated parameters.
#' @param data A data frame containing time course and cumulative germination fractions to be used in the ThermalTime model. A column with time in hours (CumTime) + a column with cumulative fractions (CumFraction) and the experiment temperature (GermTemp) are required.
#' @param model is the results list returned from from running any of the PBT models.
#' @keywords plot population-based threshold model
#' @importFrom rlang .data
#' @importFrom dplyr %>%
#' @export
#' @examples plotPBTModel(MyData, MyModelResults)
#' plotPBTModel(MyData, MyModelResults)

plotTTSubOModel <- function(data, model, germ.temp = "GermTemp", cum.time = "CumTime", cum.frac = "CumFraction") {

  # data and argument checks
  if (!is.data.frame(data)) stop("Data is not a valid data frame.")
  if (!is.list(model)) stop("Model results must be in list format as output from any of the PBT model functions.")

  # check validity of columns defs
  if (!is.element(cum.time, names(data))) stop("Cumulative time column '", cum.time, "' not found in data frame.")
  if (!is.element(cum.frac, names(data))) stop("Cumulative fraction column '", cum.frac, "' not found in data frame.")

  # check for presence of model results
  modelParams <- c("Type", "MaxCumFrac", "BaseTemp", "ThetaT50", "Sigma", "Correlation")
  lapply(modelParams, function(m) {
    if (!is.element(m, names(model))) stop("Required param '", m, "' missing from supplied model results.")
  })
  if (model$Type != "ThermalTime Suboptimal") stop("Model type must be 'ThermalTime Suboptimal' for this plot function.")

  # model params
  maxCumFrac <- model$MaxCumFrac
  Tb <- model$BaseTemp
  thetaT50 <- model$ThetaT50
  sigma <- model$Sigma
  corr <- model$Correlation

  par1 <- paste("T[b] ==", round(Tb, 1))
  par2 <- paste("ThetaT(50) ==", round(thetaT50, 3))
  par3 <- paste("sigma ==", round(sigma, 3))
  par4 <- paste("R^2 ==", round(corr, 2))

  # Function to plot all predicted treatments by the Thermaltime model
  df <- data %>% dplyr::distinct(.data[[germ.temp]], .keep_all = FALSE)
  modelLines <- plyr::alply(as.matrix(df), 1, function(temp) {
    ggplot2::stat_function(
      fun = function(x) {
        stats::pnorm(
          log(x, base = 10),
          thetaT50 - log(temp - Tb, base = 10),
          sigma,
          log = FALSE
        )
      },
      aes(color = as.factor(temp))
    )
  })

  plt <- data %>%
    ggplot2::ggplot(ggplot2::aes(x = .data[[cum.time]], y = .data[[cum.frac]], color = as.factor(.data[[germ.temp]]))) +
    geom_point(shape = 19, size = 2) +
    modelLines +
    scale_y_continuous(labels = scales::percent, expand = c(0, 0), limits = c(0, 1.02)) +
    scale_x_continuous(expand = c(0, 0)) +
    expand_limits(x = 0, y = 0) +
    labs(
      title = "HydroThermal Time Model",
      x = "Time",
      y = "Cumulative fraction germinated (%)",
      color = "Temperature") +
    guides(color = guide_legend(reverse = T, order = 1)) +
    theme_scatter_plot +
    annotate("text", x = -Inf, y = 0.95, label = paste("Model Parameters"), color = "grey0", hjust = -0.1) +
    annotate("text", x = -Inf, y = 0.9, label = par1, color = "grey0", hjust = -0.2, parse = TRUE) +
    annotate("text", x = -Inf, y = 0.85, label = par2, color = "grey0", hjust = -0.1, parse = TRUE) +
    annotate("text", x = -Inf, y = 0.8, label = par3, color = "grey0", hjust = -0.2, parse = TRUE) +
    annotate("text", x = -Inf, y = 0.75, label = par4, color = "grey0", hjust = -0.2, parse = TRUE)
  plt
}


plotHTModel <- function(data, model, ) {

#  else if (Model == "HT") { #Hydrotime suboptimal model selected  ------------------

    #Create columns on TreatDataClean with predicted values and Normalized time from model
    #TreatDataClean <<-TreatDataClean %>% as_tibble() %>% mutate(
    #  ModelPredicted = pnorm(Germ.wp-(HT/(Germ.time.hours)),Psib50,Sigma,log=FALSE)*MaxGerm,
    #  NormalizedTime = (1-(Germ.wp/(Germ.wp-(HT/Germ.time.hours))))*Germ.time.hours
    #)

    #Function to plot Normalized time using the HYDROTIME model
    # modelNormalized <<-
    #  alply(as.matrix(TreatmentsWP), 1, function(WP) {
    #   stat_function(fun=function(x){pnorm(0-(HT/(x)),Psib50,Sigma,log=FALSE)*MaxGerm}, color="blue3")
    # })
    #NormalizedAxisTitle <<- "Normalized time (h)"

    WP <- TreatData$Germ.wp
    HT <- ModelResults$HT
    psib50 <- ModelResults$psib50
    sigma <- ModelResults$sigma

    TreatFactor1 <- (as.factor(TreatData$Germ.wp))
    TreatFactor2 <- NA
    TreatFactor3 <- NA

    #Label for legends
    LegendTitleFactor1 <- "Water Potential"
    LegendTitleFactor2 <- NA
    LegendTitleFactor3 <- NA

    #Passing fitted Hydrotime Model Parameters for plot legend
    ModPar1Label <- "HT =="
    ModPar2Label <- "Psi[b](50)=="
    ModPar3Label <- "sigma == "
    ModPar4Label <- "R^2 == "
    ModPar5Label <- ""

    ModPar1 <- round(HT[1],2)
    ModPar2 <- round(psib50[1],3)
    ModPar3 <- round(sigma[1],3)
    ModPar4 <- round(Correlation[1],2)
    ModPar5 <- ""

    TreatmentsWP <- distinct(TreatData, Germ.wp, .keep_all = FALSE)

    #Function to plot all predicted treatments by the HYDROTIME model
    modellines <-
      plyr::alply(as.matrix(TreatmentsWP), 1, function(WP) {
        stat_function(fun=function(x){pnorm(WP-(HT/(x)),psib50,sigma,log=FALSE)*MaxCumFraction}, aes_(colour = factor(WP)))
      })

  } else if (Model == "HTT") { #Hydrothermal time suboptimal model selected  ------------------

    #Hydrothermal time Model - Create table to plot treatments with predicted model lines
    #TreatData$WPFactor <<- with(TreatData, (as.factor(TreatData$Germ.wp)))
    #Factor1 <<- TreatData$WPFactor
    #Factor1Title <<- "Water \n Potential"
    #TreatmentsWP <<- distinct(TreatData, Germ.wp, .keep_all = FALSE)
    #TreatData$TempFactor <<- with(TreatData, (as.factor(TreatData$Germ.temp)))
    #Factor2 <<- TreatData$TempFactor
    #Factor2Title <<- "Temperature"
    #TreatmentsTemp <<- distinct(TreatData, Germ.temp, .keep_all = FALSE)
    #TreatmentHTT <<- distinct(TreatData, Germ.temp, Germ.wp, .keep_all = FALSE)

    #Dt <<- data.frame(Treatments$Germ.temp,Treatments$Germ.wp)
    #tmps <<- Dt$Treatments.Germ.temp
    #wps <<- Dt$Treatments.Germ.wp

    WP <- TreatData$Germ.wp
    Temp <- TreatData$Germ.temp
    HT <- ModelResults$HT
    psib50 <- ModelResults$psib50
    sigma <- ModelResults$sigma
    Tb <- ModelResults$Tb

    TreatFactor1 <- (as.factor(TreatData$Germ.wp))
    TreatFactor2 <- (as.factor(TreatData$Germ.temp))
    TreatFactor3 <- NA

    #Label for legends
    LegendTitleFactor1 <- "Water Potential"
    LegendTitleFactor2 <- "Temperature"
    LegendTitleFactor3 <- NA


    #Passing fitted Hydrothermal time Model Parameters for plot legend
    ModPar1Label <- "HT =="
    ModPar2Label <- "T[b]=="
    ModPar3Label <- "psi[b](50)=="
    ModPar4Label <- "sigma == "
    ModPar5Label <- "R^2 == "

    ModPar1 <- round(HT[1],2)
    ModPar2 <- round(Tb[1],2)
    ModPar3 <- round(psib50[1],3)
    ModPar4 <- round(sigma[1],3)
    ModPar5 <- round(Correlation[1],2)

    TreatmentsWP <<- distinct(TreatData, Germ.wp, .keep_all = FALSE)
    TreatmentsTemp <<- distinct(TreatData, Germ.temp, .keep_all = FALSE)

    #Sort data with higher temperature and higher wp
    #TreatData <<- TreatData[order(TreatData$Treat.ID,-TreatData$Germ.temp, -TreatData$Germ.wp, TreatData$Germ.time.hours),]
    Treatments <- TreatData[order(TreatData$Treat.ID,-TreatData$Germ.temp, -TreatData$Germ.wp),]

    #still need to fix this. check mapply approach

    Dt <- data.frame(Treatments$Germ.temp,TreatData$Germ.wp)
    tmps <- Dt$Treatments.Germ.temp
    wps <- Dt$Treatments.Germ.wp

    #Dt <- CalcSpeed(TreatData, "Germ.wp","Germ.temp")
    #tmps <- Dt$TreatData.Germ.temp
    #wps <- Dt$TreatData.Germ.wp

    #Dt <- data.frame(TreatmentsTemp,TreatmentsWP)
    #tmps <- Dt$TreatData.Germ.temp
    #wps <- Dt$TreatData.Germ.wp

    #Function to plot all predicted treatments by the HYDROTHERMAL time model
    modellines <<-
      mapply(function(Temp1, WP1) {
        stat_function(fun=function(x){pnorm((+WP1-(HT/((Temp1-Tb)*x))),psib50,sigma, log= FALSE)*MaxCumFraction}, aes_(colour = factor(Temp1), alpha = factor(WP1)))
      }, tmps, wps)

  }




#-------
#Plot provided data with predicted lines indicated above.
  p <- ggplot(data=TreatData, aes(x=Time, y=Germ,color=TreatFactor1, alpha = TreatFactor2)) + geom_point(shape=19, size=2) + xlab("Time") + ylab("Cumulative (%)") +
    modellines + scale_alpha_discrete(range = c(0.5, 1.0)) +
    scale_y_continuous(labels = scales::percent, expand = c(0,0), limits = c(0,1.02)) +
    scale_x_continuous(expand = c(0,0)) + expand_limits(x = 0, y = 0) +
    guides(color=guide_legend(reverse=T, title=LegendTitleFactor1, order = 1),
           alpha=guide_legend(reverse=T, title=LegendTitleFactor2, order = 2)) + theme_scatter_plot +
    annotate("text", x = -Inf, y = 0.95, label = paste("Model Parameters"), color = "grey0", hjust = -0.1) +
    annotate("text", x = -Inf, y = 0.9, label = paste(ModPar1Label, ModPar1), color = "grey0", hjust = -0.15, parse = TRUE) +
    annotate("text", x = -Inf, y = 0.85, label = paste(ModPar2Label, ModPar2), color = "grey0", hjust = -0.12, parse = TRUE) +
    annotate("text", x = -Inf, y = 0.8, label = paste(ModPar3Label, ModPar3), color = "grey0", hjust = -0.15, parse = TRUE) +
    annotate("text", x = -Inf, y = 0.75, label = paste(ModPar4Label, ModPar4), color = "grey0", hjust = -0.15, parse = TRUE) +
    annotate("text", x = -Inf, y = 0.7, label = paste(ModPar5Label, ModPar5), color = "grey0", hjust = -0.2, parse = TRUE)
  p

}

