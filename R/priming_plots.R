#-----------------Tested Functions - MAc OS and Windows
#
#' A Function to plot the both priming models.
#'
#' This function plots the priming models and calculated parameters.
#' @param myData object with the calculated rates with treatments to be used in the Hydrothermal priming model. The output of the CalcSpeed function can be directly used here with the desired treatments. The fields with Treat.priming.wp, Treat.priming.temp and Treat.priming.duration need to be informed in the data file.
#' @param ModelResults is data object resulting from the CalcHPModel() or CalcHTPModel() functions containing the model information and parameter results.
#' @param GR is the column name for the rate to be used in the model and needs to be informed in case the data object contains a GR different than GR50.
#' @keywords plot priming model hydropriming hydrothermal priming
#' @importFrom tibble as_tibble
#' @export
#' @examples PlotPrimingModel(myData, HPModelResults)
#' PlotPrimingModel(myData, HPModelResults)
PlotPrimingModel <- function(Data, ModelResults, GR)
{
  Treatments <- Data
  if (missing(GR)) { #GR not informed
    grColName <- "GR50"
  } else {
    grColName <- GR
  }
  grUsed <- eval(parse(text=paste("Treatments$",grColName, sep = "")))

  if (ModelResults$Model == "HP") { #HP Model identified and update Theta Hydropriming values
    Treatments <-Treatments %>% as_tibble() %>% dplyr::mutate(
      Theta = (Treatments$Treat.priming.wp-ModelResults$PsiMin50)*Treatments$Treat.priming.duration)
    #Treatment factor for plot
    TreatFactor1 <- (as.factor(Treatments$Treat.priming.wp))
    TreatFactor2 <- (as.factor(Treatments$Treat.priming.duration))
    TreatFactor3 <- NA

    factor2lab <- "Priming \n duration"
    factor3lab <- NA

    #Pass parameters for plot
    ModPar1Label <- "Psi[min](50)=="
    ModPar2Label <- "Intercept=="
    ModPar3Label <- "Slope=="
    xAxisTitlePriming <- "Hydropriming Time"

    ModPar1 <- ModelResults$PsiMin50 # PsiMin50 Value
    ModPar2 <- ModelResults$Intercept # Intercept
    ModPar3 <- ModelResults$Slope # Slope


  } else { #HTP identified and update Theta Hydrothermal priming values
    Treatments <-Treatments %>% as_tibble() %>% dplyr::mutate(
      Theta = (Treatments$Treat.priming.wp-ModelResults$PsiMin50)*(Treatments$Treat.priming.temp-ModelResults$Tmin)*Treatments$Treat.priming.duration)
    #Treatment factor for plot
    TreatFactor1 <- (as.factor(Treatments$Treat.priming.wp))
    TreatFactor2 <- (as.factor(Treatments$Treat.priming.temp))
    TreatFactor3 <- (as.factor(Treatments$Treat.priming.duration))

    factor2lab <- "Priming \n temperature"
    factor3lab <- "Priming \n duration"

    #Pass parameters for plot
    ModPar1Label <- "Psi[min](50)=="
    ModPar2Label <- "T[min]=="
    ModPar3Label <- "Intercept=="
    ModPar4Label <- "Slope=="
    xAxisTitlePriming <- "Hydrothermal priming Time"

    ModPar1 <- ModelResults$PsiMin50 # PsiMin50 Value
    ModPar2 <- ModelResults$Tmin # Tmin Value
    ModPar3 <- ModelResults$Intercept # Intercept
    ModPar4 <- ModelResults$Slope # Slope

  }

  pPM <- ggplot(data=Treatments, aes(x=Theta, y=grUsed, color=TreatFactor1, shape = TreatFactor2, alpha = TreatFactor3)) + geom_point(size=2) + xlab(xAxisTitlePriming) + ylab("Germination Rate") +
    scale_x_continuous(expand = c(0,0)) + geom_abline(intercept = ModelResults$Intercept, slope = ModelResults$Slope, color = "blue") +
    annotate("text", x = -Inf, y = Inf, label = paste("Model Parameters"), color = "grey0", hjust = -0.1, vjust = 1.5) +
    annotate("text", x = -Inf, y = Inf, label = paste(ModPar1Label, ModPar1), color = "grey0", parse = TRUE, hjust = -0.09, vjust = 2.5) +
    annotate("text", x = -Inf, y = Inf, label = paste(ModPar2Label, ModPar2), color = "grey0", parse = TRUE, hjust = -0.12, vjust = 4.5) +
    annotate("text", x = -Inf, y = Inf, label = paste(ModPar3Label, ModPar3), color = "grey0", parse = TRUE, hjust = -0.11, vjust = 5.8) +
    annotate("text", x = -Inf, y = Inf, label = paste("R^2 == ", ModelResults$RSquared), color = "grey0", parse = TRUE, hjust = -0.2, vjust = 6.3) +
    theme_scatter_plot + labs(color = "Priming WP", shape = factor2lab, alpha = factor3lab )

  pPM
}

#----------------------New Development - Under Testing




