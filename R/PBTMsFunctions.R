#-----------------Tested Functions - MAc OS and Windows
#
#' A Function to calculate the Thermaltime model parameters.
#'
#' This function calculates the temperature base (Tb), the thermal time to 50% of the population (ThetaT50) and the standard deviation (sigma).
#' @param Data time course and cumulative dataset to be used in the Thermaltime model. The original dataframe template should be used or column names should be modified similarly to the template. A column with time in hours (CumTime) + a column with cumulative fractions (CumFract) and the experiment temperature (Germ.temp) are required. Filter the dataframe to only have treatments with temperature equal or under to  optimal temperature level.
#' @param MaxCumFract sets the ceiling cumulative fraction for the model when treatment at optimal condition displays a lower maximum cumulative fraction. Use it on your own discretion.
#' @keywords Thermal time model parameters
#' @export
#' @examples CalcTTSubOModel(myData)
#' CalcTTSubOModel(myData)
#'
CalcTTSubOModel <- function(Data, MaxCumFract)
{
  TreatData <- Data
  Germ <- TreatData$CumFract
  Time <- TreatData$CumTime
  Temp <- TreatData$Germ.temp

  if (missing(MaxCumFract)) { #MaxCumFract not informed
    MaxCumFract <- 1
  } else {
    MaxCumFract <- MaxCumFract
  }

  #Inform intial and limit values for the Hydrotime Model parameters
  # Initials
  iTb <- 6
  ithetaT50 <- 3
  iSigma <- 0.09
  #lower limits
  lTb <- 0
  lthetaT50 <- 0.5
  lSigma <- 0.0001
  #upper limits
  uTb <- 15
  uthetaT50 <- 50
  uSigma <- 0.5

  #Calculate Thermaltime Suboptimal Model Parameters- nls plus algorithm port used to add constraints on the parameters
  TTSubOModel <- nls(Germ ~ pnorm(log(Time, base = 10),mean = thetaT50-log(Temp-Tb, base = 10), sd = sigma, log= FALSE)*MaxCumFract, start=list(Tb=iTb,thetaT50=ithetaT50,sigma=iSigma),lower=list(Tb=lTb,thetaT50=lthetaT50,sigma=lSigma),upper=list(Tb=uTb,thetaT50=uthetaT50,sigma=uSigma), algorithm ="port")
  summary(TTSubOModel)

  #get some estimation of goodness of fit
  Correlation <- cor(Germ,predict(TTSubOModel))^2

  #Passing fitted Hydrotime Model Parameters
  Tb <- summary(TTSubOModel)$coefficients[[1]]
  thetaT50 <- summary(TTSubOModel)$coefficients[[2]]
  sigma <- summary(TTSubOModel)$coefficients[[3]]

  Model <- "TTsuboptimal"
  HTPModelResults <- data.frame(Model,Tb,thetaT50,sigma,MaxCumFract,Correlation)
  return(HTPModelResults)
}

#----------------------New Development - Under Testing




#' A Function to plot the selected calculated model and parameters and predicitions.
#'
#' This function plots the selected model and calculated parameters.
#' @param Data time course and cumulative dataset to be used in the Thermaltime model. The original dataframe template should be used or column names should be modified similarly to the template. A column with time in hours (CumTime) + a column with cumulative fractions (CumFract) and the experiment temperature (Germ.temp) are required. Filter the dataframe to only have treatments with temperature equal or under to  optimal temperature level.
#' @param ModelResults is data object resulting from the Calc PBTM functions containing the model information and parameter results.
#' @importFrom plyr alply
#' @keywords plot population-based threshold model
#' @export
#' @examples PlotPBTMModel(mydata, myModelResults)
#' PlotPBTMModel(mydata, myModelResults)
PlotPBTMModel <- function (Data, ModelResults)
{
  TreatData <- Data
  Germ <- TreatData$CumFract
  Time <- TreatData$CumTime
  MaxCumFract <- ModelResults$MaxCumFract
  Correlation <- ModelResults$Correlation
  Model <- ModelResults$Model

  if (Model == "TTsuboptimal")
  {

    Temp <- TreatData$Germ.temp
    Tb <- ModelResults$Tb
    thetaT50 <- ModelResults$thetaT50
    sigma <- ModelResults$sigma

    TreatFactor1 <- (as.factor(TreatData$Germ.temp))
    TreatFactor2 <- NA
    TreatFactor3 <- NA

    #Label for legends
    LegendTitleFactor1 <- "Temperature"
    LegendTitleFactor2 <- NA
    LegendTitleFactor3 <- NA

    #Passing fitted Hydrotime Model Parameters for plot legend
    ModPar1Label <- "T[b] =="
    ModPar2Label <- "thetaT(50)=="
    ModPar3Label <- "sigma == "
    ModPar4Label <- "R^2 == "
    ModPar5Label <- ""

    ModPar1 <- round(Tb[1],1)
    ModPar2 <- round(thetaT50[1],3)
    ModPar3 <- round(sigma[1],3)
    ModPar4 <- round(Correlation[1],2)
    ModPar5 <- ""

    TreatmentsTemp <<- distinct(TreatData, Germ.temp, .keep_all = FALSE)

    #Function to plot all predicted treatments by the Thermaltime model
    modellines <<-
      alply(as.matrix(TreatmentsTemp), 1, function(Temp) {
        stat_function(fun=function(x){pnorm(log(x, base = 10),thetaT50-log(Temp-Tb, base = 10),sigma,log=FALSE)*MaxCumFract}, aes_(colour = factor(Temp)))
      })
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

