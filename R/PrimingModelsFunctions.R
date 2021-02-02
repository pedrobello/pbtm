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



#----------------------New Development - Under Testing




