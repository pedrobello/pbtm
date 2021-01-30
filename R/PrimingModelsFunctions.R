#-----------------Tested Functions - MAc OS and Windows
#
#' A Function to calculate the Hydropriming model parameters.
#'
#' This function calculates the minimal water potential for priming effects (Psib50min).
#' @param myData object with the calculated rates with treatments to be used in the Hydropriming model. The output of the CalcSpeed function can be directly used here with the desired treatments. The fields with Treat.priming.wp and Treat.priming.duration need to be informed in the data file.
#' @param GR is the column name for the rate to be used in the model and needs to be informed in case the data object contains a GR different than GR50.
#' @keywords hydropriming model parameters
#' @export
#' @examples CalcHPModel(myData, "GR20")
#' CalcHPModel(myData)
CalcHPModel <- function(Data, GR)
{
  Treatments <- Data
  if (missing(GR)) { #GR not informed
    grColName <- "GR50"
  } else {
    grColName <- GR
  }

  #Initilaze values for model calculation
  psiMin50 <- -1
  grUsed <- eval(parse(text=paste("Treatments$",grColName, sep = ""))) #grUsed <- PrimingTreats$GR50
  Pwp <- Treatments$Treat.priming.wp #Pwp <- PrimingTreats$Treat.priming.wp
  Pdur <- Treatments$Treat.priming.duration #Pdur <- PrimingTreats$Treat.priming.duration

  #Function to calculate Theta Hydro Priming
  fthetaHP <- function(PM50){(Pwp-PM50)*Pdur}

  # HydroPriming - Calculate PsiMin50, y intercept and slope using NLS
  GetPsiMin50 <- nls(grUsed ~ inTer + fthetaHP(psiMin50) * sLope, algorithm="port",
                     start=c(inTer=0.001,psiMin50=-1,sLope=0.1),lower=c(inTer=0.0000001,psiMin50=-10, sLope=0.000000001),upper=c(inTer=0.1,psiMin50=-1, sLope=1))

  PsiMin50 <- round(summary(GetPsiMin50)$coefficients[[2]],3) # PsiMin50 Value
  Intercept <- round(summary(GetPsiMin50)$coefficients[[1]],4) # Intercept
  Slope <- round(summary(GetPsiMin50)$coefficients[[3]],6) # Slope

  #Future refine of the model
  #lmer(GR50 ~ ((Treat.priming.wp-psiMin50)*Treat.priming.duration), PrimingTreats, start = c(psiMin50 = -1) )
  #fthetaHP <- function(PM50,Pwp,Pdur){(Pwp-PM50)*Pdur}
  #HPlModel <- lm(GR50 ~ fthetaHP(PsiMin50,Treat.priming.wp,Treat.priming.duration), PrimingTreats)
  #summary(HPlModel)

  #Hydropriming model linear regression
  HPlModel <- lm(grUsed ~ fthetaHP(PsiMin50))

  #get some estimation of goodness of fit
  Correlation <- cor(grUsed,predict(HPlModel))^2
  RSquared <- round(Correlation[1],2)

  Model <- "HP"
  HPModelResults <- data.frame(Model,PsiMin50,Intercept,Slope,RSquared)
  return(HPModelResults)
}


#' A Function to calculate the Hydrothermal priming model parameters.
#'
#' This function calculates the minimal water potential for priming effects (Psibmin50) and minimal temperature (Tmin).
#' @param myData object with the calculated rates with treatments to be used in the Hydrothermal priming model. The output of the CalcSpeed function can be directly used here with the desired treatments. The fields with Treat.priming.wp, Treat.priming.temp and Treat.priming.duration need to be informed in the data file.
#' @param GR is the column name for the rate to be used in the model and needs to be informed in case the data object contains a GR different than GR50.
#' @keywords hydrothermal priming model parameters
#' @export
#' @examples CalcHTPModel(myData, "GR90")
#' CalcHTPModel(myData)
CalcHTPModel <- function(Data, GR)
{
  Treatments <- Data
  if (missing(GR)) { #GR not informed
    grColName <- "GR50"
  } else {
    grColName <- GR
  }
  #Initilaze values for model calculation
  psiMin50i <- -1
  Tmini <- 12
  grUsed <- eval(parse(text=paste("Treatments$",grColName, sep = ""))) #grUsed <- PrimingTreats$GR50
  Pwp <- Treatments$Treat.priming.wp
  Ptemp <- Treatments$Treat.priming.temp
  Pdur <- Treatments$Treat.priming.duration

  #Function to calculate Theta Hydro Priming
  fthetaHTP <- function(PM50,TMIN){(Pwp-PM50)*(Ptemp-TMIN)*Pdur}

  # HydroPriming - Calculate PsiMin50, y intercept and slope using NLS
  GetPsiMin50Tmin <- nls(grUsed ~ inTer + fthetaHTP(psiMin50,Tmin) * sLope, algorithm="port",
                         start=c(inTer=0.001,psiMin50=psiMin50i,Tmin=Tmini,sLope=0.1),lower=c(inTer=0.0000001,psiMin50=-10,Tmin=0.5,sLope=0.000000001),upper=c(inTer=0.1,psiMin50=-1,Tmin=20,sLope=1))

  PsiMin50 <- round(summary(GetPsiMin50Tmin)$coefficients[[2]],3) # PsiMin50 Value
  Tmin <- round(summary(GetPsiMin50Tmin)$coefficients[[3]],3) # Tmin Value
  Intercept <- round(summary(GetPsiMin50Tmin)$coefficients[[1]],4) # Intercept
  Slope <- round(summary(GetPsiMin50Tmin)$coefficients[[4]],6) # Slope

  #Hydropriming model linear regression
  HTPlModel <- lm(grUsed ~ fthetaHTP(PsiMin50,Tmin))

  #get some estimation of goodness of fit
  Correlation <- cor(grUsed,predict(HTPlModel))^2
  RSquared <- round(Correlation[1],2)

  Model <- "HTP"
  HTPModelResults <- data.frame(Model,PsiMin50,Tmin,Intercept,Slope,RSquared)
  return(HTPModelResults)
}

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




