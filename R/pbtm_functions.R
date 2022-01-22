#-----------------Tested Functions - MAc OS and Windows
#
#' A Function to calculate the Thermaltime model parameters. ------------------------
#'
#' This function calculates the temperature base (Tb), the thermal time to 50% of the population (ThetaT50) and the standard deviation (sigma).
#' @param Data time course and cumulative dataset to be used in the Thermaltime model. The original dataframe template should be used or column names should be modified similarly to the template. A column with time in hours (CumTime) + a column with cumulative fractions (CumFraction) and the experiment temperature (Germ.temp) are required. Filter the dataframe to only have treatments with temperature equal or under to  optimal temperature level.
#' @param MaxCumFraction sets the ceiling cumulative fraction for the model when treatment at optimal condition displays a lower maximum cumulative fraction. Use it on your own discretion.
#' @keywords Thermal time model parameters
#' @export
#' @examples CalcTTSubOModel(myData)
#' CalcTTSubOModel(myData)
#'
CalcTTSubOModel <- function(Data, MaxCumFraction)
{
  TreatData <- Data
  Germ <- TreatData$CumFraction
  Time <- TreatData$CumTime
  Temp <- TreatData$Germ.temp

  if (missing(MaxCumFraction)) { #MaxCumFraction not informed
    MaxCumFraction <- 1
  } else {
    MaxCumFraction <- MaxCumFraction
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
  TTSubOModel <- nls(Germ ~ pnorm(log(Time, base = 10),mean = thetaT50-log(Temp-Tb, base = 10), sd = sigma, log= FALSE)*MaxCumFraction, start=list(Tb=iTb,thetaT50=ithetaT50,sigma=iSigma),lower=list(Tb=lTb,thetaT50=lthetaT50,sigma=lSigma),upper=list(Tb=uTb,thetaT50=uthetaT50,sigma=uSigma), algorithm ="port")
  summary(TTSubOModel)

  #get some estimation of goodness of fit
  Correlation <- cor(Germ,predict(TTSubOModel))^2

  #Passing fitted Hydrotime Model Parameters
  Tb <- summary(TTSubOModel)$coefficients[[1]]
  thetaT50 <- summary(TTSubOModel)$coefficients[[2]]
  sigma <- summary(TTSubOModel)$coefficients[[3]]

  Model <- "TTsuboptimal"
  HTPModelResults <- data.frame(Model,Tb,thetaT50,sigma,MaxCumFraction,Correlation)
  return(HTPModelResults)
}

#' A Function to calculate the Hydrotime model parameters.  ------------------------
#'
#' This function calculates the hydrotime constant (HT), the median water potential base (Psib50) and the standard deviation (sigma).
#' @param Data time course and cumulative dataset to be used in the hydrotime model. The original dataframe template should be used or column names should be modified similarly to the template. A column with time in hours CumTime) + a column with cumulative fractions (CumFraction) and the experiment water potential (Germ.wp) are required. Filter the dataframe to only have treatments with water potential at the same temperature level.
#' @param MaxCumFraction sets the ceiling cumulative fraction for the model when treatment at optimal condition displays a lower maximum cumulative fraction. Use it on your own discretion.
#' @keywords hydrotime model parameters
#' @export
#' @examples CalcHTModel(Data)
#' CalcHTModel(Data, MaxCumFraction)
CalcHTModel <- function(Data, MaxCumFraction)
{
  TreatData <- Data
  Germ <- TreatData$CumFraction
  Time <- TreatData$CumTime
  WP <- TreatData$Germ.wp

  if (missing(MaxCumFraction)) { #MaxCumFraction not informed
    MaxCumFraction <- 1
  } else {
    MaxCumFraction <- MaxCumFraction
  }
  #Inform intial and limit values for the Hydrotime Model parameters
  # Initials
  iHT <- 60
  iPsib50 <- -0.8
  iSigma <- 0.2
  #lower limits
  lHT <- 1
  lPsib50 <- -5
  lSigma <- 0.0001
  #upper limits
  uHT <- 1000
  uPsib50 <- -0.000000001
  uSigma <- 2

  #Calculate Hydrotime Model Parameters- nls plus algorithm port used to add constraints on the parameters
  HTModel <<- nls(Germ ~ pnorm(+WP-(HT/Time), Psib50, Sigma, log= FALSE)*MaxCumFraction, start=list(HT=iHT,Psib50=iPsib50,Sigma=iSigma),lower=list(HT=lHT,Psib50=lPsib50,Sigma=lSigma),upper=list(HT=uHT,Psib50=uPsib50,Sigma=uSigma), algorithm ="port")
  summary(HTModel)

  #get some estimation of goodness of fit
  Correlation <<- cor(Germ,predict(HTModel))^2

  #Passing fitted Hydrotime Model Parameters
  HT <- summary(HTModel)$coefficients[[1]]
  psib50 <- summary(HTModel)$coefficients[[2]]
  sigma <- summary(HTModel)$coefficients[[3]]

  Model <- "HT"
  HTModelResults <- data.frame(Model,HT,psib50,sigma,MaxCumFraction,Correlation)
  return(HTModelResults)
}

#----------------------New Development - Under Testing



#' A Function to calculate the Hydrothermal time model parameters.
#'
#' This function calculates the hydrohermal time constant (HT), the temperature base (tb), the median water potential base (Psib50) and the standard deviation (sigma).
#' @param Data time course and cumulative dataset to be used in the hydrothermal time model. The original dataframe template should be used or column names should be modified similarly to the template. A column with time (CumTime) + a column with cumulative fractions (CumFraction), the experimental water potential (Germ.wp) and the temperature (Germ.temp) are required.
#' @param MaxCumFraction sets the ceiling cumulative fraction for the model when treatment at optimal condition displays a lower maximum cumulative fraction. Use it on your own discretion.
#' @param Tb is the temperature base that can be provided and have a fixed value or left in blank to be calculated with other hydrothermal parameters.
#' @keywords hydrothermal time model parameters
#' @export
#' @examples CalcHTTModel(mydata)
#' CalcHTTModel(mydata)
CalcHTTModel <- function(Data, MaxCumFraction, Tb)
{
  TreatData <- Data
  Germ <- TreatData$CumFraction
  Time <- TreatData$CumTime
  WP <- TreatData$Germ.wp
  Temp <- TreatData$Germ.temp

  if (missing(MaxCumFraction)) { #MaxCumFraction not informed
    MaxCumFraction <- 1
  } else {
    MaxCumFraction <- MaxCumFraction
  }
  #Inform intial and limit values for the Hydrothermal time model parameters
  # Initials
  iHT <- 800
  iTb <- 1
  iPsib50 <- -1
  iSigma <- 0.4
  #lower limits
  lHT <- 1
  lTb <- 0
  lPsib50 <- -5
  lSigma <- 0.0001
  #upper limits
  uHT <- 5000
  uTb <- 15
  uPsib50 <- 0
  uSigma <- 10

  if (missing(Tb)){
    #Calculate Hydrotime Model Parameters and Tb - nls plus algorithm port used to add constraints on the parameters
    HTTModel <<- nls(Germ ~ pnorm(WP-(HT/((Temp-Tb)*Time)),psib50,sigma, log= FALSE)*MaxCumFraction, start=list(HT=iHT,Tb=iTb,psib50=iPsib50,sigma=iSigma),lower=list(HT=lHT,Tb=lTb,psib50=lPsib50,sigma=lSigma),upper=list(HT=uHT,Tb=uTb,psib50=uPsib50,sigma=uSigma), algorithm ="port")
    summary(HTTModel)
    #Passing fitted Hydrotime Model Parameters
    HT <- summary(HTTModel)$coefficients[[1]]
    Tb <- summary(HTTModel)$coefficients[[2]]
    psib50 <- summary(HTTModel)$coefficients[[3]]
    sigma <- summary(HTTModel)$coefficients[[4]]
  } else {
    #Calculate Hydrotime Model Parameters- with fixed Tb
    HTTModel <<- nls(Germ ~ pnorm((+WP-(HT/((Temp-MyTb)*Time))),psib50,sigma, log= FALSE)*MaxCumFraction, start=list(HT=iHT,psib50=iPsib50,sigma=iSigma),lower=list(HT=lHT,psib50=lPsib50,sigma=lSigma),upper=list(HT=uHT,psib50=uPsib50,sigma=uSigma), algorithm ="port")
    summary(HTTModel)
    HT <<- summary(HTTModel)$coefficients[[1]]
    Tb <- Tb
    psib50 <<- summary(HTTModel)$coefficients[[2]]
    sigma <<- summary(HTTModel)$coefficients[[3]]

  }

  #get some estimation of goodness of fit
  Correlation <<- cor(Germ,predict(HTTModel))^2

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


  #Passing fitted Hydrotime Model Parameters for plot legend
  #ModPar1Label <<- "HT =="
  #ModPar2Label <<- "T[b]=="
  #ModPar3Label <<- "psi[b](50)=="
  #ModPar4Label <<- "sigma == "
  #ModPar5Label <<- "R^2 == "

  #ModPar1 <<- round(summary(HTTModel)$coefficients[[1]],2)
  #ModPar2 <<- round(Tb[1],2)
  #ModPar3 <<- round(psib50[1],3)
  #ModPar4 <<- round(sigma[1],3)
  #ModPar5 <<- round(Correlation[1],2)

  #Dt <<- data.frame(Treatments$Germ.temp,Treatments$Germ.wp)
  #tmps <<- Dt$Treatments.Germ.temp
  #wps <<- Dt$Treatments.Germ.wp

  #Function to plot all predicted treatments by the HYDROTHERMAL time model
  #modellines <<-
  #  mapply(function(Temp1, WP1) {
  #    stat_function(fun=function(x){pnorm((+WP1-(HT/((Temp1-Tb)*x))),Psib50,Sigma, log= FALSE)*MaxGerm}, aes_(colour = factor(Temp1), alpha = factor(WP1)))
  #  }, tmps, wps)

  #Create columns on TreatDataClean with predicted and Normalized time values from model
  #TreatDataClean <<-TreatDataClean %>% as_tibble() %>% mutate(
  #  ModelPredicted = pnorm((+Germ.wp-(HT/((Germ.temp-Tb)*Germ.time.hours))),Psib50,Sigma, log= FALSE)*MaxGerm,
  #  NormalizedTime = +((1-Germ.wp/(+Germ.wp-(HT/(Germ.time.hours*(Germ.temp-Tb)))))*Germ.time.hours*(Germ.temp-Tb))
  #)

  #Function to plot Normalized time using the HYDROTHERMAL time model
  #modelNormalized <<-
  #  alply(as.matrix(TreatmentsTemp), 1, function(Temp) {
  #    stat_function(fun=function(x){pnorm((+0-(HT/((1)*x))),Psib50,Sigma, log= FALSE)*MaxGerm}, color="blue3")
  #  })
  # NormalizedAxisTitle <<- "Normalized thermal time (°h)"


  Model <- "HTT"
  HTTModelResults <- data.frame(Model,HT,psib50,sigma,Tb,MaxCumFraction,Correlation)
  return(HTTModelResults)

}
