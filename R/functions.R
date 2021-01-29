#-----------------Tested Functions - MAc OS and Windows
#
#' A function that calculates the required time and rate to a desired fraction in cumulative curves
#'
#' This function allows you to calculate the time to a desired cumulative fraction and respective germination rate. Use this function on raw data to avoid loss of points closer to the desired cumulative fraction.
#' @param Data time course and cumulative dataset. Several treatments can be used at once as long as it respects the template and column names provided. A column with time in hours (CumTime) + a column with cumulative fractions (CumFract) are required with at least one additional column for relevant treatment (e.g., germination temperature or water potential)
#' @param Fraction from 0 to 1 used to calculate the time required for that level to be obtained in the cumulative time course. Standard value is 0.5 (50 percent), to calculate the time to 50 percent germination (T50) and respective germination rate (GR50). Fraction level can be entered and be used for calculation and change column name.
#' @param Treat1,Treat2,Treat3,Treat4,Treat5 are the names of the treatment columns to separate the dataset. The time course cumulative curves will be grouped for each distinct treatment that should be informed here. These column names do not need to be informed in case the provided template file is used to organize the data.
#' @keywords Tx, GRx, germination speed, germination rate
#' @importFrom dplyr group_by_at
#' @importFrom dplyr tally
#' @importFrom magrittr %>%
#' @export
#' @examples
#' CalcSpeed(MyData)
CalcSpeed <- function(Data, Fraction, Treat1, Treat2, Treat3, Treat4, Treat5)
{
  #Define the informed treatments (columns) to be grouped when calling the function or set standard to all columns besides CumFract and CumTime.
  if (missing(Treat5)) { #Treat5 not informed
    if (missing(Treat4)) { #Treat4 not informed
      if (missing(Treat3)) { #Treat3 not informed
        if (missing(Treat2)) { #Treat2 not informed
          if (missing(Treat1)) { #Treat1 not informed
            TreatColNames <- c("Treat.ID","Treat.desc","Treat.aging.time","Treat.priming.wp","Treat.priming.temp","Treat.priming.duration","Germ.wp","Germ.temp", "Germ.promoter.dosage", "Germ.inhibitor.dosage")
          } else {TreatColNames <- c(Treat1)}
        } else {TreatColNames <- c(Treat1,Treat2)}
      } else {TreatColNames <- c(Treat1,Treat2,Treat3)}
    } else {TreatColNames <- c(Treat1,Treat2,Treat3,Treat4)}
  } else {TreatColNames <- c(Treat1,Treat2,Treat3,Treat4,Treat5)}

  #Define the informed fraction for time and rate calculation. Standard is 0.5 -> 50 percent. Also populates new column name base in the fraction used.
  if (missing(Fraction)) { #Fraction not informed
    Frac <- 0.5
    FracSpeedLbl <- "T50"
    FracRateLbl <- "GR50"
  } else {
    Frac <- Fraction
    FracSpeedLbl <- paste("T",(Frac*100), sep = "")
    FracRateLbl <- paste("GR",(Frac*100), sep = "")
  }

  # Calculate the time to the fraction selected using linear interpolation (approx function) after it calculates the inverse of time to generate the rate.
  Treatments <- Data %>% group_by_at(TreatColNames) %>%
    dplyr::mutate(Tx = approx(CumFract,CumTime, xout=Frac, ties="ordered")$y,
                  GRx = 1/approx(CumFract,CumTime, xout=Frac, ties="ordered")$y)

  TreatColNames <- c(TreatColNames, "Tx", "GRx")

  # Separate all treatments without germination time courses
  Treatments <- Treatments %>% group_by_at(TreatColNames) %>% tally()

  #Replaces the column name with the used fraction instead.
  names(Treatments)[names(Treatments) == "Tx"] <- FracSpeedLbl
  names(Treatments)[names(Treatments) == "GRx"] <- FracRateLbl

  return(Treatments)
}


#' A Function to plot rate vs desired treatment.
#'
#' This function plots rates against the desired treatment.
#' @param Data should inform the table that resulted from the CalcSpeed function. It should have the summarized treatments and respective GR values.
#' @param x should indicate the treatment column name for the x axis (e.g., "Germ.temp", "Germ.wp" or others).
#' @param y should indicate the rate column name for the y axis if different than GR50 (e.g., "GR90", "GR10", etc).
#' @keywords plot rates Temperature
#' @import ggplot2
#' @export
#' @examples PlotRateVsTreat(MyCalcSpeedData,"Germ.temp")
#' PlotRateVsTreat(MyCalcSpeedData,"Germ.temp")
PlotRateVsTreat <- function (Data, x, y)
{
  Treatments <- Data
  if (missing(x)) { #x/treatment not informed
    print("Informed treatment for x axis.")
  } else {
    Treat <- x
  }
  if (missing(y)) { #y/rate not informed
    rate <- "GR50"
  } else {
    rate <- y
  }
  pGR <- ggplot(data=Treatments, aes_string(x=Treat, y=rate, color=Treat)) + geom_point(shape=19, size=2) +
    expand_limits(x = 0, y = 0) + theme_scatter_plot
  pGR
}


#' A Function to plot the cumulative raw data selected.
#'
#' This function plots the raw data previously selected. It will plot raw cumulative curves over time.
#' @param Data with time course data for one or more treatments.
#' @param Treat1,Treat2 are the desired treatment factors to be informed as column names. The first informed treatment will be separated as color and the second as shape.
#' @keywords plot raw data
#' @import ggplot2
#' @export
#' @examples PlotRawDt(MyData,"Germ.temp")
#' PlotRawDt(MyRawData,"Germ.temp")
PlotRawDt <- function(Data, Treat1, Treat2)
{
  TreatData <- Data
  #MaxTime <- TreatData$CumTime[which.max(TreatData$CumTime)] #Gets the longest time measurement in the dataset provided.
  #PlotTime <- ceiling(MaxTime/5)*5 #make maxTime a multiple of 5 - ceiling() rounds up
  #Increment <- round(PlotTime/5, digits = 0) #define tick mark separation

  gp <- "geom_point(shape=19, size=2)" #Add standard fixed shape in case shape is not used as the second factor

  if (missing(Treat1)) { #treatment 1 not informed
    print("Informed treatment for factor.")
  } else {
    eval(parse(text=paste("TreatData$",Treat1, " <- (factor(TreatData$",Treat1,"))", sep = "")))
    T1 <- Treat1
  }
  if (missing(Treat2)) { #treatment 2 not informed
    T2 <- NA
  } else {
    gp <- "geom_point(size=2)" #Add fixed shape in case shape is not used as the second factor
    eval(parse(text=paste("TreatData$",Treat2, " <- (factor(TreatData$",Treat2,"))", sep = "")))
    T2 <- Treat2
  }
  #if (missing(Treat3)) { #treatment 3 not informed - possible in the near future
  #  T3 <- NA
  #  alph <- ""
  #} else {
  # eval(parse(text=paste("(as.factor(TreatData$",Treat3,"))", sep = "")))
  #  T3 <- Treat3
  #  alph <- "color=T1, shape=T2, alpha=T3"
  #}

  #Plot All Treatments with fitted equation (Whole data plot here, including repetitive percentages)
  pRaw <<- ggplot(data=TreatData, aes_string(x="CumTime", y="CumFract", color=T1, shape=T2 )) +
    eval(parse(text=gp)) + geom_line() + xlab("Time") + ylab("Cumulative (%)") +
    scale_y_continuous(labels = scales::percent, expand = c(0,0), limits = c(0,1.02)) +
    scale_x_continuous(expand = c(0,0)) +
    theme_scatter_plot
  pRaw
}

#' A Function to clean cumulative curves on dataset.
#'
#' This function removes repetitive cumulative fractions/percentages, keeping only the initial presence of the value
#' @param Data object with the raw cumulative data that needs to be removed.
#' @param Treat1,Treat2,Treat3,Treat4,Treat5 are the names of the treatment columns to separate the dataset. The time course cumulative curves will be grouped for each distinct treatment that should be informed here.
#' @keywords Clean cumulative fraction repetitive percentage
#' @importFrom dplyr distinct
#' @export
#' @examples CleanData(mydata,"Treat.desc")
#' CleanData(mydata,"Treat.desc")
CleanData <- function(Data, Treat1, Treat2, Treat3, Treat4, Treat5) {
  #Clean Repetitive Percentages (keeps only initial presence of value)

  if (missing(Data)) { #data object not informed
    print("Informe the data object.")
  } else {
    TreatData <- Data
    if (missing(Treat1)) { #treatment 1 not informed
      print("Informed treatment for factor.")
    } else {
      T1 <-  Treat1
      if (missing(Treat2)) { #treatment 2 not informed
        T2 <- ""
      } else {
        T2 <- paste(",",Treat2, sep = "")
      }
      if (missing(Treat3)) { #treatment 3 not informed
        T3 <- ""
      } else {
        T3 <- paste(",",Treat3, sep = "")
      }
      if (missing(Treat4)) { #treatment 2 not informed
        T4 <- ""
      } else {
        T4 <- paste(",",Treat4, sep = "")
      }
      if (missing(Treat5)) { #treatment 2 not informed
        T5 <- ""
      } else {
        T5 <- paste(",",Treat5, sep = "")
      }

      TreatDataClean <- eval(parse(text=paste("distinct(TreatData,",T1,T2,T3,T4,T5,",CumFract, .keep_all = TRUE)", sep="")))

      return(TreatDataClean)
    }
  }
}

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


#ggplot package theme ----------------------------------------------------------------------------------
theme_scatter_plot <- theme(
  legend.background = element_blank(),
  legend.key = element_blank(),
  legend.title = element_text(size=12, color ="black"),
  legend.text = element_text(size=12, color ="black"),
  panel.background = element_blank(),
  panel.grid.minor.y=element_blank(),
  panel.grid.major.x=element_blank(),
  panel.grid.major.y=element_blank(),
  panel.grid.minor.x= element_blank(),
  panel.border = element_rect(colour = "grey50", fill=NA, size=0.5),
  strip.background = element_rect(colour="black", fill="white"),
  axis.ticks = element_line(color = "black", size =0.5),
  axis.text = element_text(size=12, color ="black"),
  axis.title = element_text(size=14, color ="black",face = "bold"),
  axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
  plot.title = element_blank())


#----------------------New Development - Under Testing


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

  #Thermaltime Suboptimal Model - Create table to plot treatments with predicted model lines
  #TreatData$TempFactor <<- with(TreatData, (as.factor(TreatData$Germ.temp)))
  #Factor1 <<- TreatData$TempFactor
  #Factor1Title <<- "Temperature"
  #TreatmentsTemp <<- distinct(TreatData, Germ.temp, .keep_all = FALSE)

  #Passing fitted Hydrotime Model Parameters
  Tb <- summary(TTSubOModel)$coefficients[[1]]
  thetaT50 <- summary(TTSubOModel)$coefficients[[2]]
  sigma <- summary(TTSubOModel)$coefficients[[3]]

  #Passing fitted Hydrotime Model Parameters for plot legend
  #ModPar1Label <<- "T[b] =="
  #ModPar2Label <<- "θT(50)=="
  #ModPar3Label <<- "σ == "
  #ModPar4Label <<- "R^2 == "
  #ModPar5Label <<- ""

  #ModPar1 <<- round(Tb[1],1)
  #ModPar2 <<- round(thetaT50[1],3)
  #ModPar3 <<- round(Sigma[1],3)
  #ModPar4 <<- round(Correlation[1],2)
  #ModPar5 <<- ""


  #Function to plot all predicted treatments by the Thermaltime model
  #modellines <<-
  #  alply(as.matrix(TreatmentsTemp), 1, function(Temp) {
  #    stat_function(fun=function(x){pnorm(log(x, base = 10),thetaT50-log(Temp-Tb, base = 10),Sigma,log=FALSE)*MaxGerm}, aes_(colour = factor(Temp)))
  #  })

  #Create column on TreatDataClean with predicted values from model
  #TreatDataClean <<-TreatDataClean %>% as_tibble() %>% mutate(
  #  ModelPredicted = pnorm(log(CumTime, base = 10),thetaT50-log(Germ.temp-Tb, base = 10),Sigma,log=FALSE)*MaxGerm,
  #  NormalizedTime = (Germ.temp-Tb)*CumTime
  #)

  #Function to plot Normalized time using the HYDROTIME model
  #modelNormalized <<-
  #  alply(as.matrix(TreatmentsTemp), 1, function(Temp) {
  #    stat_function(fun=function(x){pnorm(log(x, base = 10),thetaT50-log(1, base = 10),Sigma,log=FALSE)*MaxGerm}, color="blue3")
  #  })
  #NormalizedAxisTitle <<- "Normalized thermal time (°h)"

  #Plot raw data and predicted model with parameters
  #PlotPBTMModel()

  Model <- "TTsuboptimal"
  HTPModelResults <- data.frame(Model,Tb,thetaT50,sigma,MaxCumFract,Correlation)
  return(HTPModelResults)
}


#' A Function to plot the selected calculated model and parameters and predicitions.
#'
#' This function plots the selected model and calculated parameters.
#' @param Data time course and cumulative dataset to be used in the Thermaltime model. The original dataframe template should be used or column names should be modified similarly to the template. A column with time in hours (CumTime) + a column with cumulative fractions (CumFract) and the experiment temperature (Germ.temp) are required. Filter the dataframe to only have treatments with temperature equal or under to  optimal temperature level.
#' @param ModelResults is data object resulting from the Calc PBTM functions containing the model information and parameter results.
#' @importFrom plyr alply
#' @keywords plot population-based threshold model
#' @export
#' @examples prePlotPBTMModel()
#' prePlotPBTMModel()
prePlotPBTMModel <- function (Data, ModelResults)
{
  TreatData <- Data
  Germ <- TreatData$CumFract
  Time <- TreatData$CumTime
  Temp <- TreatData$Germ.temp
  Tb <- ModelResults$Tb
  thetaT50 <- ModelResults$thetaT50
  sigma <- ModelResults$sigma
  MaxCumFract <- ModelResults$MaxCumFract
  Correlation <- ModelResults$Correlation

  TreatFactor1 <- (as.factor(TreatData$Germ.temp))
  TreatFactor2 <- NA
  TreatFactor3 <- NA

  #Label for legends
  LegendTitleFactor1 <- "Temperature"
  LegendTitleFactor2 <- NA
  LegendTitleFactor3 <- NA

  #Passing fitted Hydrotime Model Parameters for plot legend
  ModPar1Label <- "T[b] =="
  ModPar2Label <- "θT(50)=="
  ModPar3Label <- "σ == "
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

  p <- ggplot(data=TreatData, aes(x=Time, y=Germ,color=TreatFactor1, alpha = TreatFactor2)) + geom_point(shape=19, size=2) + xlab("Time (hours)") + ylab("Germination (%)") +
    modellines + scale_alpha_discrete(range = c(0.5, 1.0)) +
    scale_y_continuous(labels = scales::percent, expand = c(0,0), limits = c(0,1.02)) +
    scale_x_continuous(expand = c(0,0)) +
    guides(color=guide_legend(reverse=T, title=LegendTitleFactor1, order = 1),
           alpha=guide_legend(reverse=T, title=LegendTitleFactor2, order = 2)) + theme_scatter_plot +
    annotate("text", x = -Inf, y = Inf, label = paste("Model \n Parameters"), color = "grey0") +
    annotate("text", x = -Inf, y = Inf, label = paste(ModPar1Label, ModPar1), color = "grey0", parse = TRUE) +
    annotate("text", x = -Inf, y = Inf, label = paste(ModPar2Label, ModPar2), color = "grey0", parse = TRUE) +
    annotate("text", x = -Inf, y = Inf, label = paste(ModPar3Label, ModPar3), color = "grey0", parse = TRUE) +
    annotate("text", x = -Inf, y = Inf, label = paste(ModPar4Label, ModPar4), color = "grey0", parse = TRUE) +
    annotate("text", x = -Inf, y = Inf, label = paste(ModPar5Label, ModPar5), color = "grey0", parse = TRUE)
  p

}



#' A Function to plot the selected calculated model and parameters and predicitions.
#'
#' This function plots the selected model and calculated parameters.
#' @param
#' @keywords plot population-based threshold model
#' @export
#' @examples PlotModl()
#' PlotModl()
PlotModl <- function ()
{
  #Plot All Treatments with fitted equation (Whole data plot here, including repetitive percentages)
  p <- ggplot(data=TreatData, aes(x=Time, y=Germ,color=TreatFactor1, alpha = TreatFactor2)) + geom_point(shape=19, size=2) + xlab("Time (hours)") + ylab("Germination (%)") +
    modellines + scale_alpha_discrete(range = c(0.5, 1.0)) +
    scale_y_continuous(labels = scales::percent, expand = c(0,0), limits = c(0,1.02)) +
    scale_x_continuous(expand = c(0,0)) +
    guides(color=guide_legend(reverse=T, title=LegendTitleFactor1, order = 1),
           alpha=guide_legend(reverse=T, title=LegendTitleFactor2, order = 2)) + theme_scatter_plot +
    annotate("text", x = -Inf, y = Inf, label = paste("Model \n Parameters"), color = "grey0") +
    annotate("text", x = -Inf, y = Inf, label = paste(ModPar1Label, ModPar1), color = "grey0", parse = TRUE) +
    annotate("text", x = -Inf, y = Inf, label = paste(ModPar2Label, ModPar2), color = "grey0", parse = TRUE) +
    annotate("text", x = -Inf, y = Inf, label = paste(ModPar3Label, ModPar3), color = "grey0", parse = TRUE) +
    annotate("text", x = -Inf, y = Inf, label = paste(ModPar4Label, ModPar4), color = "grey0", parse = TRUE) +
    annotate("text", x = -Inf, y = Inf, label = paste(ModPar5Label, ModPar5), color = "grey0", parse = TRUE)
  p
}
