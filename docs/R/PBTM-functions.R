#' A T50 and GR50 Function
#'
#' This function allows you to calculate the time to 50% germination (T50) and respective germination rate (GR50).
#' @param CalcT50nGR50() to calculate T50 and GR50 on the TreatData dataset.
#' @keywords T50, GR50, germination speed, germination rate
#' @export
#' @examples
#' CalcT50nGR50()
CalcT50nGR50 <- function()
{
  # Calculate Time to 50% Germination (T50) (calculate on raw data to avoid loss of points closer to 50% germination) + GR50
  Treatments <<- TreatData %>% group_by(Treat.ID, Treat.desc, Treat.aging.time, Treat.priming.wp, Treat.priming.temp,Treat.priming.duration,Germ.wp,Germ.temp, Germ.promoter.dosage, Germ.inhibitor.dosage) %>%
    dplyr::mutate(T50 = approx(Germ.fraction,Germ.time.hours, xout=0.5, ties="ordered")$y,
                  GR50 = 1/approx(Germ.fraction,Germ.time.hours, xout=0.5, ties="ordered")$y)

  # Separate all treatments without germination time courses
  Treatments <<- Treatments %>% group_by(Treat.ID, Treat.desc, Treat.aging.time, Treat.priming.wp, Treat.priming.temp,Treat.priming.duration,Germ.wp,Germ.temp, Germ.promoter.dosage, Germ.inhibitor.dosage, T50, GR50) %>% tally()
}


#' A Function to define the population-based threshold model(PBTM) you want to work with.
#'
#' This function allows you to choose the PBT model to be calculated. Choose the model that you want to work on (6.1)Hydropriming model; (2.1)Suboptimal Hydrotime; (2.2)Supra-optimal Hydrotime; (1.1)Thermaltime; (3.1)Suboptimal Hydrothermal Time; (3.2)Supra-optimal Hydrothermal Time.
#' @param DefineModel(1.1) selects the Thermal time Sub-optimal Model
#' @param DefineModel(2.1) selects the Suboptimal Hydrotime Model
#' @keywords PBTM
#' @return Sets up the package and print the selected model.
#' @export
#' @examples
#' DefineModel(2.1)
DefineModel <- function(ModelCode)
{
  CalcT50nGR50()

  if (ModelCode == 1.1) { #Thermal time suboptimal Model Selected
  ActiveModel <<- 1.1

  #Sort data with higher temperature
  TreatData <<- TreatData[order(TreatData$Treat.ID,-TreatData$Germ.temp, TreatData$Germ.time.hours),]
  Treatments <<- Treatments[order(Treatments$Treat.ID,-Treatments$Germ.temp),]

  #Treatment factor for plot
  TreatFactor1 <<- (as.factor(TreatData$Germ.temp))
  TreatFactor2 <<- NA
  TreatFactor3 <<- NA

  #Label for legends
  LegendTitleFactor1 <<- "Temperature"
  LegendTitleFactor2 <<- NA
  LegendTitleFactor3 <<- NA

  #Calculates how many temp and wp treatments for plotting purposes
  TempTreatNumber <<- length(unique(Treatments$Germ.temp))
  if (TempTreatNumber < 5)
  {
    my_colors <<- scale_colour_manual(values=rev(c("red", "darkorange","forestgreen","blue3")))
  }
  else
  {
    my_colors <<- scale_colour_brewer(palette = "Spectral")
  }

  print("Thermal time suboptimal model selected")


  } else if (ModelCode == 2.1) { #Hydrotime suboptimal model selected
    ActiveModel <<- 2.1

    #Sort data with higher wp
    TreatData <<- TreatData[order(TreatData$Treat.ID,-TreatData$Germ.wp, TreatData$Germ.time.hours),]
    Treatments <<- Treatments[order(Treatments$Treat.ID,-Treatments$Germ.wp),]

    #Treatment factor for plot
    TreatFactor1 <<- (as.factor(TreatData$Germ.wp))
    TreatFactor2 <<- NA
    TreatFactor3 <<- NA
    TreatFactorSum3 <<- NA

    #Label for legends
    LegendTitleFactor1 <<- "Water\nPotential"
    LegendTitleFactor2 <<- NA
    LegendTitleFactor3 <<- NA

    #Calculates how many temp and wp treatments for plotting purposes
    WPTreatNumber <<- length(unique(Treatments$Germ.wp))
    if (WPTreatNumber < 5)
      {
      my_colors <<- scale_colour_manual(values=rev(c("blue3", "forestgreen", "darkorange", "red")))
      }
    else
      {
      my_colors <<- scale_colour_brewer(palette = "Spectral")
      }

    print("Hydrotime suboptimal model selected")

  } else if (ModelCode == 2.2) { #Hydrotime supra-optimal model selected
    ActiveModel <<- 2.2

    DefinedHTo <<- FALSE

    #Sort data with higher wp
    TreatData <<- TreatData[order(TreatData$Treat.ID,-TreatData$Germ.wp, TreatData$Germ.time.hours),]
    Treatments <<- Treatments[order(Treatments$Treat.ID,-Treatments$Germ.wp),]

    #Treatment factor for plot
    TreatFactor1 <<- (as.factor(TreatData$Germ.wp))
    TreatFactor2 <<- NA
    TreatFactor3 <<- NA

    #Label for legends
    LegendTitleFactor1 <<- "Water\nPotential"
    LegendTitleFactor2 <<- NA
    LegendTitleFactor3 <<- NA


    #Calculates how many temp and wp treatments for plotting purposes
    WPTreatNumber <<- length(unique(Treatments$Germ.wp))
    if (WPTreatNumber < 5)
    {
      my_colors <<- scale_colour_manual(values=rev(c("blue3", "forestgreen", "darkorange", "red")))
    }
    else
    {
      my_colors <<- scale_colour_brewer(palette = "Spectral")
    }

    print("Hydrotime supra-optimal model selected")

  } else if (ModelCode == 3.1) { #Hydrothermal time suboptimal model selected
    ActiveModel <<- 3.1

    #Sort data with higher temperature and higher wp
    TreatData <<- TreatData[order(TreatData$Treat.ID,-TreatData$Germ.temp, -TreatData$Germ.wp, TreatData$Germ.time.hours),]
    Treatments <<- Treatments[order(Treatments$Treat.ID,-Treatments$Germ.temp, -Treatments$Germ.wp),]

    #Treatment factor for plot
    TreatFactor1 <<- (as.factor(TreatData$Germ.temp))
    TreatFactor2 <<- (as.factor(TreatData$Germ.wp))
    TreatFactor3 <<- NA

    #Label for legends
    LegendTitleFactor1 <<- "Temperature"
    LegendTitleFactor2 <<- "Water\nPotential"
    LegendTitleFactor3 <<- NA

    #Calculates how many temp and wp treatments for plotting purposes
    TempTreatNumber <<- length(unique(Treatments$Germ.temp))
    WPTreatNumber <<- length(unique(Treatments$Germ.wp))

    FixTb <<- FALSE

    #Calculates how many temp and wp treatments for plotting purposes
    TempTreatNumber <<- length(unique(Treatments$Germ.temp))
    if (TempTreatNumber < 5)
    {
      my_colors <<- scale_colour_manualmy_colors <<- scale_colour_manual(values=rev(c("red","darkorange","forestgreen","blue3")))
    }
    else
    {
      my_colors <<- scale_colour_brewer(palette = "Spectral")
    }

    print("Hydrothermal time suboptimal model selected")

  } else if (ModelCode == 3.2) { #Hydrothermal time supra-optimal model selected
    ActiveModel <<- 3.2

    #Sort data with higher temperature and higher wp
    TreatData <<- TreatData[order(TreatData$Treat.ID,-TreatData$Germ.temp, -TreatData$Germ.wp, TreatData$Germ.time.hours),]
    Treatments <<- Treatments[order(Treatments$Treat.ID,-Treatments$Germ.temp, -Treatments$Germ.wp),]

    #Treatment factor for plot
    TreatFactor1 <<- (as.factor(TreatData$Germ.temp))
    TreatFactor2 <<- (as.factor(TreatData$Germ.wp))
    TreatFactor3 <<- NA


    #Label for legends
    LegendTitleFactor1 <<- "Temperature"
    LegendTitleFactor2 <<- "Water\nPotential"
    LegendTitleFactor3 <<- NA

    #Calculates how many temp and wp treatments for plotting purposes
    TempTreatNumber <<- length(unique(Treatments$Germ.temp))
    WPTreatNumber <<- length(unique(Treatments$Germ.wp))

    FixTo <<- FALSE

    #Calculates how many temp and wp treatments for plotting purposes
    TempTreatNumber <<- length(unique(Treatments$Germ.temp))
    if (TempTreatNumber < 5)
    {
      my_colors <<- scale_colour_manualmy_colors <<- scale_colour_manual(values=rev(c("red","darkorange","forestgreen","blue3")))
    }
    else
    {
      my_colors <<- scale_colour_brewer(palette = "Spectral")
    }

    print("Hydrothermal time supra-optimal model selected")

  } else if (ModelCode == 6.1) { #Hydropriming model selected
      ActiveModel <<- 6.1

      #Label for legends
      LegendTitleFactor1 <<- "Priming\nWater\nPotential"
      LegendTitleFactor2 <<- "Priming\nDuration"
      LegendTitleFactor3 <<- NA

      #Treatment factor for plot
      TreatFactor1 <<- (as.factor(TreatData$Treat.priming.wp))
      TreatFactorSum1 <<- (as.factor(Treatments$Treat.priming.wp))
      TreatFactor2 <<- (as.factor(TreatData$Treat.priming.duration))
      TreatFactorSum2 <<- (as.factor(Treatments$Treat.priming.duration))
      TreatFactor3 <<- NA
      TreatFactorSum3 <<- NA

      #Calculates how many temp and wp treatments for plotting purposes
      WPTreatNumber <<- length(unique(Treatments$Treat.priming.wp))
      if (WPTreatNumber < 5)
      {
        my_colors <<- scale_colour_manual(values=rev(c("blue3", "forestgreen", "darkorange", "red")))
      }
      else
      {
        my_colors <<- scale_colour_brewer(palette = "Spectral")
      }

      print("Hydropriming model selected")

    } else if (ModelCode == 6.2) { #HydrotherMal priming model selected
      ActiveModel <<- 6.2

      #Label for legends
      LegendTitleFactor1 <<- "Priming\nWater\nPotential"
      LegendTitleFactor2 <<- "Priming\nTemperature"
      LegendTitleFactor3 <<- "Priming\nDuration"

      #Treatment factor for plot
      TreatFactor1 <<- (as.factor(TreatData$Treat.priming.wp))
      TreatFactorSum1 <<- (as.factor(Treatments$Treat.priming.wp))
      TreatFactor2 <<- (as.factor(TreatData$Treat.priming.temp))
      TreatFactorSum2 <<- (as.factor(Treatments$Treat.priming.temp))
      TreatFactor3 <<- (as.factor(TreatData$Treat.priming.duration))
      TreatFactorSum3 <<- (as.factor(Treatments$Treat.priming.duration))

      #Calculates how many temp and wp treatments for plotting purposes
      WPTreatNumber <<- length(unique(Treatments$Treat.priming.wp))
      if (WPTreatNumber < 5)
      {
        my_colors <<- scale_colour_manual(values=rev(c("blue3", "forestgreen", "darkorange", "red")))
      }
      else
      {
        my_colors <<- scale_colour_brewer(palette = "Spectral")
      }

      print("Hydrothermal priming model selected")

  } else {
    ActiveModel <<- 0
    print("No model selected") }
}

#' A Function to calculate and plot the PBTM selected.
#'
#' This function triggers the model calculation and outputs calculated parameters and graph prediction.
#' @param CALCnPLOTModel() calculates the selected model
#' @keywords PBTM
#' @export
#' @examples
CALCnPLOTModel <- function()
{
  CleanGermData()
  if (ActiveModel == 1.1) {
    #Calculate and Plot Thermaltime Model
    CalcTTSubOModel()

  } else if (ActiveModel == 2.1) {
    #Calculate and Plot Hydrotime Model
    CalcHTModel()

  } else if (ActiveModel == 2.2) {
    #Calculate and Plot Supra-optimal Hydrotime Model
    CalcSupraHTModel()

  } else if (ActiveModel == 3.1) {
    #Calculate and Plot Hydrothermal time Model
    CalcHTTModel()

  } else if (ActiveModel == 3.2) {
    #Calculate and Plot Supra-optimal Hydrothermal time Model
    CalcHTTModelSupra()

  } else if (ActiveModel == 6.1) {
    #Calculate and Plot Hydropriming Model
    CalcHPModel()

  } else if (ActiveModel == 6.2) {
    #Calculate and Plot Hydrothermal priming Model
    CalcHTPModel()

  } else {
    print("No model selected") }
}


#' A Function to get maximum germination percentage on dataset.
#'
#' This function normalizes all data to the percentage informed here. This is useful when a known fraction of the population is not viable.
#' @param DefineMaxGerm(0.85) will define germination to a max. limit of 85%
#' @keywords Maximum germination
#' @export
#' @examples DefineMaxGerm()
#' DefineMaxGerm()
DefineMaxGerm <- function(TreatID, GermWP, GermTEMP){
  DtSet <<- subset(mydata, Treat.ID == TreatID & Germ.wp == GermWP & Germ.temp == GermTEMP)
  MaxGerm <<- DtSet$Germ.fraction[which.max(DtSet$Germ.fraction)]
  print(paste0("Maximum germination to be used is ", MaxGerm))
}


#' A Function to clean cumulative curves on dataset.
#'
#' This function removes repetitive percentages, keeping only the intial presence of the value
#' @param CleanGermData() will remove repetitive points
#' @keywords Clean germination repetitive percentage
#' @export
#' @examples CleanGermData()
#' CleanGermData()
CleanGermData <- function() {
  #Clean Repetitive Percentages (keeps only initial presence of value)
  TreatDataClean <<- distinct(TreatData, Treat.ID, Germ.temp, Germ.wp, Germ.fraction, .keep_all = TRUE)

  #Pass clean data for model parameters
  Temp <<- TreatDataClean$Germ.temp
  WP <<- TreatDataClean$Germ.wp
  Time.hours <<- TreatDataClean$Germ.time.hours
  Germ <<- TreatDataClean$Germ.fraction

      if (ActiveModel == 6.1) { #Hydropriming Model Selected
        #Treatment factor for plot
        CleanTreatFactor1 <<- (as.factor(TreatDataClean$Treat.priming.wp))
        CleanTreatFactor2 <<- (as.factor(TreatDataClean$Treat.priming.duration))

      } else if (ActiveModel == 6.2) { #Hydrothermal priming Model Selected
        #Treatment factor for plot
        CleanTreatFactor1 <<- (as.factor(TreatDataClean$Treat.priming.wp))
        CleanTreatFactor2 <<- (as.factor(TreatDataClean$Treat.priming.temp))
        CleanTreatFactor3 <<- (as.factor(TreatDataClean$Treat.priming.duration))

      } else if (ActiveModel == 2.1) { #Hydrotime Model Selected
        #Treatment factor for plot
        CleanTreatFactor1 <<- (as.factor(TreatDataClean$Germ.wp))
        CleanTreatFactor2 <<- NA

      } else if (ActiveModel == 2.2) { #Supra-optimal Hydrotime Model Selected
        #Treatment factor for plot
        CleanTreatFactor1 <<- (as.factor(TreatDataClean$Germ.wp))
        CleanTreatFactor2 <<- NA

      } else if (ActiveModel == 1.1) { #Thermaltime Model Selected
        #Treatment factor for plot
        CleanTreatFactor1 <<- (as.factor(TreatDataClean$Germ.temp))
        CleanTreatFactor2 <<- NA

      } else if (ActiveModel == 3.1) { #Hydrothermal time Model Selected
        #Treatment factor for plot
        CleanTreatFactor1 <<- (as.factor(TreatDataClean$Germ.temp))
        CleanTreatFactor2 <<- (as.factor(TreatDataClean$Germ.wp))

      } else if (ActiveModel == 3.2) { #Supra-optimal Hydrothermal time Model Selected
        #Treatment factor for plot
        CleanTreatFactor1 <<- (as.factor(TreatDataClean$Germ.temp))
        CleanTreatFactor2 <<- (as.factor(TreatDataClean$Germ.wp))

      } else {
        ActiveModel <<- 0
        print("No model selected") }

  PlotCleanData()
}

#' A Function to manually define the temperature base (Tb) for model calculation.
#'
#' This function defines the Tb value and avoid it to be calculated together with other values.
#' @param FixedTb(11) will fix the Tb to 11C for further calculation of the model.
#' @keywords Known temperature base fixed Tb
#' @export
#' @examples FixedTb(11)
#' FixedTb()
FixedTb <- function(Tb1)
{
  FixTb <<- TRUE
  MyTb <<-Tb1
  print(paste0("Temperature base manually defined at ", MyTb))
}

#' A Function to manually define the optimal temperature (To) for model calculation.
#'
#' This function defines the To value and avoid it to be calculated together with other values.
#' @param FixedTo(22.5) will fix optimal temperature to 22.5C for further model calculation.
#' @keywords Known otpimal temperature fixed To
#' @export
#' @examples FixedTo(22.5)
#' FixedTo()
FixedTo <- function(To1)
{
  FixTo <<- TRUE
  MyTo <<-To1
  print(paste0("Optimal temperature manually defined at ", MyTo))
}

#' A Function to manually define the Hydrotime constant at optimal temperature (HTo) for model calculation.
#'
#' This function defines the HTo value and avoid it to be calculated together with other values.
#' @param DefineHTo(50) will define hydrotime constant to 50 for further supra-optimal hydrotime model calculation.
#' @keywords Known Hydrotime constant supra-optimal hydrotime temperature fixed HTo
#' @export
#' @examples DefineHTo(50)
#' DefineHTo()
DefineHTo <- function(HTo1)
{
  DefinedHTo <<- TRUE
  HTo <<-HTo1
  print(paste0("Hydrotime constant was manually defined at ", HTo))
}


#Plotting Settings----------------------------------------------------------------------------------
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

#' A Function to plot the raw data selected.
#'
#' This function plots the raw data previously selected.
#' @param PlotRawData() will plot raw cumulative curves over time.
#' @keywords plot raw data
#' @export
#' @examples PlotRawData()
#' PlotRawData()
PlotRawData <- function()
{
  MaxTime <<- TreatData$Germ.time.hours[which.max(TreatData$Germ.time.hours)]
  PlotTime <<- ceiling(MaxTime/5)*5
  Increment <<- round(PlotTime/5, digits = 0)

  #Plot All Treatments with fitted equation (Whole data plot here, including repetitive percentages)
  pRaw <<- ggplot(data=TreatData, aes(x=Germ.time.hours, y=Germ.fraction, color=TreatFactor1, alpha = TreatFactor2, shape = TreatFactor3)) +
    geom_point(shape=19, size=2) + geom_line() + xlab("Time (hours)") + ylab("Germination (%)") + my_colors + scale_alpha_discrete(range = c(0.3, 1.0)) +
    scale_y_continuous(labels = scales::percent, expand = c(0,0), limits = c(0,1.02)) +
    scale_x_continuous(expand = c(0,0), limits = c(0,PlotTime+(Increment/5))) +
    guides(color=guide_legend(reverse=T, title=LegendTitleFactor1, order = 1),
           alpha=guide_legend(reverse=T, title=LegendTitleFactor2, order = 2),
           alpha=guide_legend(reverse=T, title=LegendTitleFactor3, order = 3)) +
    theme_scatter_plot


  pRaw
}


#' A Function to plot the cleaned data selected.
#'
#' This function plots the cleaned data previously selected after removing repetitive percentages. Only intial presence of the percentage is displayed.
#' @param PlotCleanData() will plot cleaned cumulative curves over time.
#' @keywords plot cleaned data
#' @export
#' @examples PlotCleanData()
#' PlotCleanData()
PlotCleanData <- function()
{
  CleanMaxTime <<- TreatDataClean$Germ.time.hours[which.max(TreatDataClean$Germ.time.hours)]
  CleanPlotTime <<- ceiling(CleanMaxTime/5)*5
  CleanIncrement <<- round(CleanPlotTime/5, digits = 0)

  #Plot All Treatments with fitted equation (Whole data plot here, including repetitive percentages)
  pClean <<- ggplot(data=TreatDataClean, aes(x=Germ.time.hours, y=Germ.fraction, color=CleanTreatFactor1, alpha = CleanTreatFactor2)) +
    geom_point(shape=19, size=2) + geom_line() + xlab("Time (hours)") + ylab("Germination (%)") + my_colors + scale_alpha_discrete(range = c(0.3, 1.0)) +
    scale_y_continuous(labels = scales::percent, expand = c(0,0), limits = c(0,1.02)) +
    scale_x_continuous(expand = c(0,0), limits = c(0,CleanPlotTime+(CleanIncrement/5))) +
    guides(color=guide_legend(reverse=T, title=LegendTitleFactor1, order = 1),
           alpha=guide_legend(reverse=T, title=LegendTitleFactor2, order = 2)) +
    theme_scatter_plot
  pClean
}


#' A Function to plot the both priming models.
#'
#' This function plots the priming models and calculated parameters.
#' @param
#' @keywords plot priming model hydropriming hydrothermal priming
#' @export
#' @examples PlotPrimingModel()
#' PlotPrimingModel()
PlotPrimingModel <- function ()
{
  pPM <- ggplot(data=Treatments, aes(x=Theta, y=GR50, color=TreatFactorSum1, alpha = TreatFactorSum2)) + geom_point(shape=19, size=2) + xlab(xAxisTitlePriming) + ylab("Germination Rate (GR50)") +
  scale_x_continuous(expand = c(0,0), limits = c(0,PlotTheta+(IncrementTheta/5))) + my_colors + scale_alpha_discrete(range = c(0.3, 1.0)) +
  geom_abline(intercept = Inter, slope = Slope, color = "blue") +
  annotate("text", x = -Inf, y = Inf, label = paste("Model Parameters"), color = "grey0", hjust = -0.1, vjust = 1.5) +
  annotate("text", x = -Inf, y = Inf, label = paste(ModPar1Label, ModPar1), color = "grey0", parse = TRUE, hjust = -0.09, vjust = 2.5) +
  annotate("text", x = -Inf, y = Inf, label = paste(ModPar2Label, ModPar2), color = "grey0", parse = TRUE, hjust = -0.12, vjust = 4.5) +
  annotate("text", x = -Inf, y = Inf, label = paste(ModPar3Label, ModPar3), color = "grey0", parse = TRUE, hjust = -0.11, vjust = 5.8) +
  annotate("text", x = -Inf, y = Inf, label = paste("R^2 == ", RSquaredPlot), color = "grey0", parse = TRUE, hjust = -0.2, vjust = 6.3) +
  theme_scatter_plot + theme(legend.position="none")

  #Plot Hydropriming Model with two columns
  grid.arrange(pRaw,pPM, ncol=2)
}

#' A Function to plot the selected calculated model and parameters and predicitions.
#'
#' This function plots the selected model and calculated parameters.
#' @param
#' @keywords plot population-based threshold model
#' @export
#' @examples PlotPBTMModel()
#' PlotPBTMModel()
PlotPBTMModel <- function ()
{
  #Plot All Treatments with fitted equation (Whole data plot here, including repetitive percentages)
  p <- ggplot(data=TreatData, aes(x=Germ.time.hours, y=Germ.fraction,color=TreatFactor1, alpha = TreatFactor2)) + geom_point(shape=19, size=2) + xlab("Time (hours)") + ylab("Germination (%)") +
    modellines + my_colors + scale_alpha_discrete(range = c(0.5, 1.0)) +
    scale_y_continuous(labels = scales::percent, expand = c(0,0), limits = c(0,1.02)) +
    scale_x_continuous(expand = c(0,0), limits = c(0,PlotTime+(Increment/5))) +
    guides(color=guide_legend(reverse=T, title=LegendTitleFactor1, order = 1),
           alpha=guide_legend(reverse=T, title=LegendTitleFactor2, order = 2)) + theme_scatter_plot +
    annotate("text", x =PlotTime*.12, y = .96, label = paste("Model \n Parameters"), color = "grey0") +
    annotate("text", x =PlotTime*.12, y = .88, label = paste(ModPar1Label, ModPar1), color = "grey0", parse = TRUE) +
    annotate("text", x =PlotTime*.12, y = .83, label = paste(ModPar2Label, ModPar2), color = "grey0", parse = TRUE) +
    annotate("text", x =PlotTime*.12, y = .78, label = paste(ModPar3Label, ModPar3), color = "grey0", parse = TRUE) +
    annotate("text", x =PlotTime*.12, y = .73, label = paste(ModPar4Label, ModPar4), color = "grey0", parse = TRUE) +
    annotate("text", x =PlotTime*.12, y = .68, label = paste(ModPar5Label, ModPar5), color = "grey0", parse = TRUE)
  p
}

#' A Function to plot the selected calculated model and parameters and predicitions with the cleaned data.
#'
#' This function plots the selected model and calculated parameters with the cleaned data.
#' @param
#' @keywords plot population-based threshold model
#' @export
#' @examples PlotCleanModel()
#' PlotCleanModel()
PlotCleanModel <- function ()
{
  #Plot All Treatments with fitted equation (Whole data plot here, including repetitive percentages)
  p <- ggplot(data=TreatDataClean, aes(x=Germ.time.hours, y=Germ.fraction,color=CleanTreatFactor1, alpha = CleanTreatFactor2)) + geom_point(shape=19, size=2) + xlab("Time (hours)") + ylab("Germination (%)") +
    modellines + my_colors + scale_alpha_discrete(range = c(0.5, 1.0)) +
    scale_y_continuous(labels = scales::percent, expand = c(0,0), limits = c(0,1.02)) +
    scale_x_continuous(expand = c(0,0), limits = c(0,CleanPlotTime+(CleanIncrement/5))) +
    guides(color=guide_legend(reverse=T, title=LegendTitleFactor1, order = 1),
           alpha=guide_legend(reverse=T, title=LegendTitleFactor2, order = 2)) + theme_scatter_plot +
    annotate("text", x =CleanPlotTime*.12, y = .96, label = paste("Model \n Parameters"), color = "grey0") +
    annotate("text", x =CleanPlotTime*.12, y = .88, label = paste(ModPar1Label, ModPar1), color = "grey0", parse = TRUE) +
    annotate("text", x =CleanPlotTime*.12, y = .83, label = paste(ModPar2Label, ModPar2), color = "grey0", parse = TRUE) +
    annotate("text", x =CleanPlotTime*.12, y = .78, label = paste(ModPar3Label, ModPar3), color = "grey0", parse = TRUE) +
    annotate("text", x =CleanPlotTime*.12, y = .73, label = paste(ModPar4Label, ModPar4), color = "grey0", parse = TRUE) +
    annotate("text", x =CleanPlotTime*.12, y = .68, label = paste(ModPar5Label, ModPar5), color = "grey0", parse = TRUE)
  p
}

#' A Function to plot GR50 vs temperature.
#'
#' This function plots germination rates against temperature.
#' @param
#' @keywords plot GR50 Temperature
#' @export
#' @examples PlotGR50vsTemp()
#' PlotGR50vsTemp()
PlotGR50vsTemp <- function ()
{

  pGR50 <- ggplot(data=Treatments, aes(x=Germ.temp, y=GR50)) + geom_point(shape=19, size=2) + xlab("Temperature (°C)") +
    ylab(bquote('Germination rate ('*h^-1*')')) +
    expand_limits(x = 0, y = 0) + theme_scatter_plot
  pGR50
}

#' A Function to plot Psib50 vs temperature.
#'
#' This function plots median water potential base against temperature .
#' @param
#' @keywords plot Psib50 Temperature
#' @export
#' @examples PlotPsib50vsTemp()
#' PlotPsib50vsTemp()
PlotPsib50vsTemp <- function ()
{
  #Still needs work.

  pGR50 <- ggplot(data=Treatments, aes(x=Germ.temp, y=GR50)) + geom_point(shape=19, size=2) + xlab("Temperature (°C)") +
    ylab(bquote('Germination rate ('*h^-1*')')) +
    expand_limits(x = 0, y = 0) + theme_scatter_plot
  pGR50
}

#' A Function to plot actual and model predicted germination.
#'
#' This function plots germination time points predicted by the model with actual germination percentage with treatments.
#' @param
#' @keywords plot actual predicted model data
#' @export
#' @examples PlotActualvsModelPredictedTreatments()
#' PlotActualvsModelPredictedTreatments()
PlotActualvsModelPredictedTreatments <- function ()
{
  FitActualPredicted <- lm(Germ.fraction ~ ModelPredicted, data = TreatDataClean)
  summary(FitActualPredicted)

  #Plot All Treatments with fitted equation (Whole data plot here, including repetitive percentages)
  pM <- ggplot(data=TreatDataClean, aes(x=Germ.fraction, y=ModelPredicted,color=CleanTreatFactor1, alpha = CleanTreatFactor2)) + geom_point(shape=19, size=2) + xlab("Actual Germination (%)") + ylab("Model Prediction (%)") +
    my_colors + scale_alpha_discrete(range = c(0.5, 1.0)) + theme_scatter_plot + abline(FitActualPredicted)
  pM
}

#' A Function to plot actual and model predicted germination.
#'
#' This function plots germination time points predicted by the model with actual germination percentage.
#' @param
#' @keywords plot actual predicted model data
#' @export
#' @examples PlotActualvsModelPredicted()
#' PlotActualvsModelPredicted()
PlotActualvsModelPredicted <- function ()
{

  #Plot All Treatments with fitted equation (Whole data plot here, including repetitive percentages)
  pM <- ggplot(data=TreatDataClean, aes(x=Germ.fraction, y=ModelPredicted)) + geom_point(shape=19, size=2) + xlab("Actual Germination (%)") + ylab("Model Prediction (%)") +
    my_colors + scale_alpha_discrete(range = c(0.5, 1.0)) + theme_scatter_plot +stat_summary(fun.data=mean_cl_normal) +
    geom_smooth(method='lm',formula=y~x) +
    annotate(x=0.2, y=0.9,
           label=paste("R^2 == ", round(cor(TreatDataClean$Germ.fraction, TreatDataClean$ModelPredicted)^2,2)),
           geom="text", size=4, parse = TRUE)
  pM
}

#' A Function to plot model normalized time with actual germination data.
#'
#' This function plots normalized time predicted by the model with actual germination percentage.
#' @param
#' @keywords plot normalized time actual predicted model data
#' @export
#' @examples PlotNormalizedTime()
#' PlotNormalizedTime()
PlotNormalizedTime <- function ()
{

  #Plot All Treatments with fitted equation (Whole data plot here, including repetitive percentages)
  pM <- ggplot(data=TreatDataClean, aes(x=NormalizedTime, y=Germ.fraction)) + geom_point(shape=19, size=2) + xlab(NormalizedAxisTitle) + ylab("Germination (%)") +
    my_colors + scale_alpha_discrete(range = c(0.5, 1.0)) + theme_scatter_plot + xlim(0, NA) +
    modelNormalized
  pM
}

#' A Function to calculate the Hydrotime model parameters.
#'
#' This function calculates the hydrotime constant (HT), the median water potential base (Psib50) and the standard deviation (sigma).
#' @param
#' @keywords hydrotime model parameters
#' @export
#' @examples CalcHTModel()
#' CalcHTModel()
CalcHTModel <- function()
{
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
  HTModel <<- nls(Germ ~ pnorm(+WP-(HT/Time.hours), Psib50, Sigma, log= FALSE)*MaxGerm, start=list(HT=iHT,Psib50=iPsib50,Sigma=iSigma),lower=list(HT=lHT,Psib50=lPsib50,Sigma=lSigma),upper=list(HT=uHT,Psib50=uPsib50,Sigma=uSigma), algorithm ="port")
  summary(HTModel)

  #get some estimation of goodness of fit
  Correlation <<- cor(Germ,predict(HTModel))^2

  #Hydrotime Model - Create table to plot treatments with predicted model lines
  TreatData$WPFactor <<- with(TreatData, (as.factor(TreatData$Germ.wp)))
  Factor1 <<- TreatData$WPFactor
  Factor1Title <<- "Water \n Potential"
  TreatmentsWP <<- distinct(TreatData, Germ.wp, .keep_all = FALSE)

  #Passing fitted Hydrotime Model Parameters
  HT <<- summary(HTModel)$coefficients[[1]]
  Psib50 <<- summary(HTModel)$coefficients[[2]]
  Sigma <<- summary(HTModel)$coefficients[[3]]

  #Passing fitted Hydrotime Model Parameters for plot legend
  ModPar1Label <<- "HT =="
  ModPar2Label <<- "Psi[b](50)=="
  ModPar3Label <<- "σ == "
  ModPar4Label <<- "R^2 == "
  ModPar5Label <<- ""

  ModPar1 <<- round(summary(HTModel)$coefficients[[1]],2)
  ModPar2 <<- round(Psib50[1],3)
  ModPar3 <<- round(Sigma[1],3)
  ModPar4 <<- round(Correlation[1],2)
  ModPar5 <<- ""

  RSquaredPlot <<-  paste("R^2 == ", round(Correlation[1],2 ))

  #Function to plot all predicted treatments by the HYDROTIME model
  modellines <<-
    alply(as.matrix(TreatmentsWP), 1, function(WP) {
      stat_function(fun=function(x){pnorm(WP-(HT/(x)),Psib50,Sigma,log=FALSE)*MaxGerm}, aes_(colour = factor(WP)))
    })

  #Create columns on TreatDataClean with predicted values and Normalized time from model
  TreatDataClean <<-TreatDataClean %>% as_tibble() %>% mutate(
    ModelPredicted = pnorm(Germ.wp-(HT/(Germ.time.hours)),Psib50,Sigma,log=FALSE)*MaxGerm,
    NormalizedTime = (1-(Germ.wp/(Germ.wp-(HT/Germ.time.hours))))*Germ.time.hours
  )

  #Function to plot Normalized time using the HYDROTIME model
  modelNormalized <<-
    alply(as.matrix(TreatmentsWP), 1, function(WP) {
      stat_function(fun=function(x){pnorm(0-(HT/(x)),Psib50,Sigma,log=FALSE)*MaxGerm}, color="blue3")
    })
  NormalizedAxisTitle <<- "Normalized time (h)"

  #Plot raw data and predicted model with parameters
  PlotPBTMModel()
}

#' A Function to calculate the supra-optimal Hydrotime model parameters.
#'
#' This function calculates the supra-optimal hydrotime constant (HT), the median water potential base (Psib50) and the variation (sigma).
#' @param
#' @keywords supra-optimal hydrotime model parameters
#' @export
#' @examples CalcSupraHTModel()
#' CalcSupraHTModel()
CalcSupraHTModel <- function()
{
  if (DefinedHTo == TRUE) { #Check if Hydrotime constant was informed before calculation model

  #Inform intial and limit values for the Hydrotime Model parameters
  # Initials
  iPsib50 <- -0.6
  iSigma <- 0.15
  #lower limits
  lPsib50 <- -4
  lSigma <- 0.01
  #upper limits
  uPsib50 <- -0.001
  uSigma <- 2

  #Calculate Hydrotime Model Parameters- nls plus algorithm port used to add constraints on the parameters
  SupraHTModel <<- nls(Germ ~ pnorm(WP-(HTo/Time.hours), Psib50, Sigma, log= FALSE)*MaxGerm, start=list(Psib50=iPsib50,Sigma=iSigma),lower=list(Psib50=lPsib50,Sigma=lSigma),upper=list(Psib50=uPsib50,Sigma=uSigma), algorithm ="port")
  summary(SupraHTModel)

  #get some estimation of goodness of fit
  Correlation <<- cor(Germ,predict(SupraHTModel))^2

  #Hydrotime Model - Create table to plot treatments with predicted model lines
  TreatData$WPFactor <<- with(TreatData, (as.factor(TreatData$Germ.wp)))
  Factor1 <<- TreatData$WPFactor
  Factor1Title <<- "Water \n Potential"
  TreatmentsWP <<- distinct(TreatData, Germ.wp, .keep_all = FALSE)

  #Passing fitted Hydrotime Model Parameters
  HT <<- HTo
  Psib50 <<- summary(SupraHTModel)$coefficients[[1]]
  Sigma <<- summary(SupraHTModel)$coefficients[[2]]

  #Passing fitted Hydrotime Model Parameters for plot legend
  ModPar1Label <<- "HT =="
  ModPar2Label <<- "Psi[b](50)=="
  ModPar3Label <<- "σ == "
  ModPar4Label <<- "R^2 == "
  ModPar5Label <<- ""

  ModPar1 <<- round(HTo[1],3)
  ModPar2 <<- round(Psib50[1],3)
  ModPar3 <<- round(Sigma[1],3)
  ModPar4 <<- round(Correlation[1],2)
  ModPar5 <<- ""

  RSquaredPlot <<-  paste("R^2 == ", round(Correlation[1],2 ))

  #Function to plot all predicted treatments by the HYDROTIME model
  modellines <<-
    alply(as.matrix(TreatmentsWP), 1, function(WP) {
      stat_function(fun=function(x){pnorm(WP-(HTo/(x)),Psib50,Sigma,log=FALSE)*MaxGerm}, aes_(colour = factor(WP)))
    })

  #Create columns on TreatDataClean with predicted values and Normalized time from model
  TreatDataClean <<-TreatDataClean %>% as_tibble() %>% mutate(
    ModelPredicted = pnorm(Germ.wp-(HTo/(Germ.time.hours)),Psib50,Sigma,log=FALSE)*MaxGerm,
    NormalizedTime = (1-(Germ.wp/(Germ.wp-(HTo/Germ.time.hours))))*Germ.time.hours
  )

  #Function to plot Normalized time using the HYDROTIME model
  modelNormalized <<-
    alply(as.matrix(TreatmentsWP), 1, function(WP) {
      stat_function(fun=function(x){pnorm(0-(HTo/(x)),Psib50,Sigma,log=FALSE)*MaxGerm}, color="blue3")
    })
  NormalizedAxisTitle <<- "Normalized time (h)"

  #Plot raw data and predicted model with parameters
  PlotPBTMModel()

  } else
    print("Hydrotime constant must be defined first.")
}


#' A Function to calculate the Hydropriming model parameters.
#'
#' This function calculates the minimal water potential for priming effects (Psib50min).
#' @param
#' @keywords hydropriming model parameters
#' @export
#' @examples CalcHPModel()
#' CalcHPModel()
CalcHPModel <- function()
{
  #Initilaze values for model calculation
  ΨMin50 <- -1
  gr50 <- Treatments$GR50
  Pwp <- Treatments$Treat.priming.wp
  Pdur <- Treatments$Treat.priming.duration

  #Function to calculate Theta Hydro Priming
  fθHP <- function(PM50){(Pwp-PM50)*Pdur}

  # HydroPriming - Calculate PsiMin50, y intercept and slope using NLS
  GetPsiMin50 <<- nls(gr50 ~ inTer + fθHP(ΨMin50) * sLope, algorithm="port",
                      start=c(inTer=0.001,ΨMin50=-1,sLope=0.1),lower=c(inTer=0.0000001,ΨMin50=-10, sLope=0.000000001),upper=c(inTer=0.1,ΨMin50=-1, sLope=1))

  #Pass PsiMin50 to Treatments table
  Treatments <<-Treatments %>% as_tibble() %>% mutate(
    PsiMin50 = summary(GetPsiMin50)$coefficients[[2]])

  #Treatments$PsiMin50 <<- summary(GetPsiMin50)$coefficients[[2]]

  #Hydropriming model linear regression
  HPlModel <<- lm(gr50 ~ fθHP(Treatments$PsiMin50))

  #Pass parameters for plot
  ModPar1Label <<- "Psi[min](50)=="
  ModPar2Label <<- "Intercept=="
  ModPar3Label <<- "Slope=="
  xAxisTitlePriming <<- "Hydropriming Time"

  ModPar1 <<- round(summary(GetPsiMin50)$coefficients[[2]],3) # PsiMin50 Value
  ModPar2 <<- round(summary(GetPsiMin50)$coefficients[[1]],4) # Intercept
  ModPar3 <<- round(summary(GetPsiMin50)$coefficients[[3]],6) # Slope

  Inter <<- summary(GetPsiMin50)$coefficients[[1]]
  Slope <<- summary(GetPsiMin50)$coefficients[[3]]

  #get some estimation of goodness of fit
  Correlation <<- cor(gr50,predict(HPlModel))^2
  RSquaredPlot <<- round(Correlation[1],2)

  #Update Theta Hydropriming values
  Treatments <<-Treatments %>% as_tibble() %>% mutate(
    Theta = (Treatments$Treat.priming.wp-Treatments$PsiMin50)*Treatments$Treat.priming.duration)

  #Treatments$Theta <<- (Treatments$Treat.priming.wp-Treatments$PsiMin50)*Treatments$Treat.priming.duration

  #Set legend position at bottom right corner
  LegPosV <<- 0.01
  LegPosH <<- 0.80

  #Get higher Hydropriming time (theta HP) calculated on the dataset
  MaxTheta <<- Treatments$Theta[which.max(Treatments$Theta)]
  PlotTheta <<- (round(MaxTheta/10, digits = 0)+1)*10
  IncrementTheta <<- round(PlotTheta/50, digits = 0)*10

  PlotPrimingModel()
}

#' A Function to calculate the Hydrothermal priming model parameters.
#'
#' This function calculates the minimal water potential for priming effects (Psibmin50) and minimal temperature (Tmin).
#' @param
#' @keywords hydrothermal priming model parameters
#' @export
#' @examples CalcHTPModel()
#' CalcHTPModel()
CalcHTPModel <- function()
{
  #Initilaze values for model calculation
  ΨMin50i <- -1
  Tmini <- 12
  gr50 <- Treatments$GR50
  Pwp <- Treatments$Treat.priming.wp
  Ptemp <- Treatments$Treat.priming.temp
  Pdur <- Treatments$Treat.priming.duration

  #Function to calculate Theta Hydro Priming
  fθHTP <- function(PM50,TMIN){(Pwp-PM50)*(Ptemp-TMIN)*Pdur}

  # HydroPriming - Calculate PsiMin50, y intercept and slope using NLS
  GetPsiMin50Tmin <<- nls(gr50 ~ inTer + fθHTP(ΨMin50,Tmin) * sLope, algorithm="port",
                      start=c(inTer=0.001,ΨMin50=ΨMin50i,Tmin=Tmini,sLope=0.1),lower=c(inTer=0.0000001,ΨMin50=-10,Tmin=0.5,sLope=0.000000001),upper=c(inTer=0.1,ΨMin50=-1,Tmin=20,sLope=1))

  #Pass PsiMin50 to Treatments table
  Treatments <<-Treatments %>% as_tibble() %>% mutate(
    PsiMin50 = summary(GetPsiMin50Tmin)$coefficients[[2]],
    Tmin = summary(GetPsiMin50Tmin)$coefficients[[3]])

  #Treatments$PsiMin50 <<- summary(GetPsiMin50Tmin)$coefficients[[2]]
  #Treatments$Tmin <<- summary(GetPsiMin50Tmin)$coefficients[[2]]

  #Hydropriming model linear regression
  HTPlModel <<- lm(gr50 ~ fθHTP(Treatments$PsiMin50,Treatments$Tmin))

  #Pass parameters for plot
  ModPar1Label <<- "Psi[min](50)=="
  ModPar2Label <<- "T[min]=="
  ModPar3Label <<- "Intercept=="
  ModPar4Label <<- "Slope=="
  xAxisTitlePriming <<- "Hydrothermal priming Time"

  ModPar1 <<- round(summary(GetPsiMin50Tmin)$coefficients[[2]],3) # PsiMin50 Value
  ModPar2 <<- round(summary(GetPsiMin50Tmin)$coefficients[[3]],3) # Tmin Value
  ModPar3 <<- round(summary(GetPsiMin50Tmin)$coefficients[[1]],4) # Intercept
  ModPar4 <<- round(summary(GetPsiMin50Tmin)$coefficients[[4]],4) # Slope

  Inter <<- summary(GetPsiMin50Tmin)$coefficients[[1]]
  Slope <<- summary(GetPsiMin50Tmin)$coefficients[[4]]

  #get some estimation of goodness of fit
  Correlation <<- cor(gr50,predict(HTPlModel))^2
  RSquaredPlot <<- round(Correlation[1],2)

  #Update Theta Hydropriming values
  Treatments <<-Treatments %>% as_tibble() %>% mutate(
    Theta = (Treatments$Treat.priming.wp-Treatments$PsiMin50)*(Treatments$Treat.priming.temp-Tmin)*Treatments$Treat.priming.duration,
    Tmin = Tmin)

  #Treatments$Theta <<- (Treatments$Treat.priming.wp-Treatments$PsiMin50)*Treatments$Treat.priming.duration
  #Treatments$Tmin <<- (Treatments$Treat.priming.wp-Treatments$PsiMin50)*Treatments$Treat.priming.duration

  #Set legend position at bottom right corner
  LegPosV <<- 0.01
  LegPosH <<- 0.80

  #Get higher Hydropriming time (theta HP) calculated on the dataset
  MaxTheta <<- Treatments$Theta[which.max(Treatments$Theta)]
  PlotTheta <<- (round(MaxTheta/10, digits = 0)+1)*10
  IncrementTheta <<- round(PlotTheta/50, digits = 0)*10

  PlotPrimingModel()
}

#' A Function to calculate the Thermaltime model parameters.
#'
#' This function calculates the temperature base (Tb), the thermal time to 50% of the population (ThetaT50) and the standard deviation (sigma).
#' @param
#' @keywords Thermal time model parameters
#' @export
#' @examples CalcTTSubOModel()
#' CalcTTSubOModel()
#'
CalcTTSubOModel <- function()
{
  #Inform intial and limit values for the Hydrotime Model parameters
  # Initials
  iTb <<- 6
  iθT50 <<- 3
  iSigma <<- 0.09
  MaxGerm <<- 1
  #lower limits
  lTb <<- 0
  lθT50 <<- 0.5
  lSigma <<- 0.0001
  #upper limits
  uTb <<- 15
  uθT50 <<- 50
  uSigma <<- 0.5

  #Calculate Thermaltime Suboptimal Model Parameters- nls plus algorithm port used to add constraints on the parameters
  TTSubOModel <<- nls(Germ ~ pnorm(log(Time.hours, base = 10),mean = θT50-log(Temp-Tb, base = 10), sd = Sigma, log= FALSE)*MaxGerm, start=list(Tb=iTb,θT50=iθT50,Sigma=iSigma),lower=list(Tb=lTb,θT50=lθT50,Sigma=lSigma),upper=list(Tb=uTb,θT50=uθT50,Sigma=uSigma), algorithm ="port")
  summary(TTSubOModel)

  #get some estimation of goodness of fit
  Correlation <<- cor(Germ,predict(TTSubOModel))^2

  #Thermaltime Suboptimal Model - Create table to plot treatments with predicted model lines
  TreatData$TempFactor <<- with(TreatData, (as.factor(TreatData$Germ.temp)))
  Factor1 <<- TreatData$TempFactor
  Factor1Title <<- "Temperature"
  TreatmentsTemp <<- distinct(TreatData, Germ.temp, .keep_all = FALSE)

  #Passing fitted Hydrotime Model Parameters
  Tb <<- summary(TTSubOModel)$coefficients[[1]]
  θT50 <<- summary(TTSubOModel)$coefficients[[2]]
  Sigma <<- summary(TTSubOModel)$coefficients[[3]]

  #Passing fitted Hydrotime Model Parameters for plot legend
  ModPar1Label <<- "T[b] =="
  ModPar2Label <<- "θT(50)=="
  ModPar3Label <<- "σ == "
  ModPar4Label <<- "R^2 == "
  ModPar5Label <<- ""

  ModPar1 <<- round(Tb[1],1)
  ModPar2 <<- round(θT50[1],3)
  ModPar3 <<- round(Sigma[1],3)
  ModPar4 <<- round(Correlation[1],2)
  ModPar5 <<- ""


  #Function to plot all predicted treatments by the Thermaltime model
  modellines <<-
    alply(as.matrix(TreatmentsTemp), 1, function(Temp) {
      stat_function(fun=function(x){pnorm(log(x, base = 10),θT50-log(Temp-Tb, base = 10),Sigma,log=FALSE)*MaxGerm}, aes_(colour = factor(Temp)))
    })

  #Create column on TreatDataClean with predicted values from model
  TreatDataClean <<-TreatDataClean %>% as_tibble() %>% mutate(
    ModelPredicted = pnorm(log(Germ.time.hours, base = 10),θT50-log(Germ.temp-Tb, base = 10),Sigma,log=FALSE)*MaxGerm,
    NormalizedTime = (Germ.temp-Tb)*Germ.time.hours
  )

  #Function to plot Normalized time using the HYDROTIME model
  modelNormalized <<-
    alply(as.matrix(TreatmentsTemp), 1, function(Temp) {
      stat_function(fun=function(x){pnorm(log(x, base = 10),θT50-log(1, base = 10),Sigma,log=FALSE)*MaxGerm}, color="blue3")
    })
  NormalizedAxisTitle <<- "Normalized thermal time (°h)"

  #Plot raw data and predicted model with parameters
  PlotPBTMModel()
}

#' A Function to calculate the Hydrothermal time model parameters.
#'
#' This function calculates the hydrohermal time constant (HT), the temperature base (tb), the median water potential base (Psib50) and the standard deviation (sigma).
#' @param
#' @keywords hydrothermal time model parameters
#' @export
#' @examples CalcHTTModel()
#' CalcHTTModel()
CalcHTTModel <- function()
{
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
  uSigma <<- 10

  if (FixTb == FALSE) {
  #Calculate Hydrotime Model Parameters- nls plus algorithm port used to add constraints on the parameters
  HTTModel <<- nls(Germ ~ pnorm(WP-(HT/((Temp-Tb)*Time.hours)),Psib50,Sigma, log= FALSE)*MaxGerm, start=list(HT=iHT,Tb=iTb,Psib50=iPsib50,Sigma=iSigma),lower=list(HT=lHT,Tb=lTb,Psib50=lPsib50,Sigma=lSigma),upper=list(HT=uHT,Tb=uTb,Psib50=uPsib50,Sigma=uSigma), algorithm ="port")
  summary(HTTModel)
  #Passing fitted Hydrotime Model Parameters
  HT <<- summary(HTTModel)$coefficients[[1]]
  Tb <<- summary(HTTModel)$coefficients[[2]]
  Psib50 <<- summary(HTTModel)$coefficients[[3]]
  Sigma <<- summary(HTTModel)$coefficients[[4]]
  } else {
  #Calculate Hydrotime Model Parameters- with fixed Tb
  HTTModel <<- nls(Germ ~ pnorm((+WP-(HT/((Temp-MyTb)*Time.hours))),Psib50,Sigma, log= FALSE)*MaxGerm, start=list(HT=iHT,Psib50=iPsib50,Sigma=iSigma),lower=list(HT=lHT,Psib50=lPsib50,Sigma=lSigma),upper=list(HT=uHT,Psib50=uPsib50,Sigma=uSigma), algorithm ="port")
  summary(HTTModel)
  HT <<- summary(HTTModel)$coefficients[[1]]
  Tb <<- MyTb
  Psib50 <<- summary(HTTModel)$coefficients[[2]]
  Sigma <<- summary(HTTModel)$coefficients[[3]]

  }

  #get some estimation of goodness of fit
  Correlation <<- cor(Germ,predict(HTTModel))^2

  #Hydrothermal time Model - Create table to plot treatments with predicted model lines
  TreatData$WPFactor <<- with(TreatData, (as.factor(TreatData$Germ.wp)))
  Factor1 <<- TreatData$WPFactor
  Factor1Title <<- "Water \n Potential"
  TreatmentsWP <<- distinct(TreatData, Germ.wp, .keep_all = FALSE)
  TreatData$TempFactor <<- with(TreatData, (as.factor(TreatData$Germ.temp)))
  Factor2 <<- TreatData$TempFactor
  Factor2Title <<- "Temperature"
  TreatmentsTemp <<- distinct(TreatData, Germ.temp, .keep_all = FALSE)
  TreatmentHTT <<- distinct(TreatData, Germ.temp, Germ.wp, .keep_all = FALSE)


  #Passing fitted Hydrotime Model Parameters for plot legend
  ModPar1Label <<- "HT =="
  ModPar2Label <<- "T[b]=="
  ModPar3Label <<- "Psi[b](50)=="
  ModPar4Label <<- "σ == "
  ModPar5Label <<- "R^2 == "

  ModPar1 <<- round(summary(HTTModel)$coefficients[[1]],2)
  ModPar2 <<- round(Tb[1],2)
  ModPar3 <<- round(Psib50[1],3)
  ModPar4 <<- round(Sigma[1],3)
  ModPar5 <<- round(Correlation[1],2)

  Dt <<- data.frame(Treatments$Germ.temp,Treatments$Germ.wp)
  tmps <<- Dt$Treatments.Germ.temp
  wps <<- Dt$Treatments.Germ.wp

  #Function to plot all predicted treatments by the HYDROTHERMAL time model
  modellines <<-
    mapply(function(Temp1, WP1) {
      stat_function(fun=function(x){pnorm((+WP1-(HT/((Temp1-Tb)*x))),Psib50,Sigma, log= FALSE)*MaxGerm}, aes_(colour = factor(Temp1), alpha = factor(WP1)))
    }, tmps, wps)

  #Create columns on TreatDataClean with predicted and Normalized time values from model
  TreatDataClean <<-TreatDataClean %>% as_tibble() %>% mutate(
    ModelPredicted = pnorm((+Germ.wp-(HT/((Germ.temp-Tb)*Germ.time.hours))),Psib50,Sigma, log= FALSE)*MaxGerm,
    NormalizedTime = +((1-Germ.wp/(+Germ.wp-(HT/(Germ.time.hours*(Germ.temp-Tb)))))*Germ.time.hours*(Germ.temp-Tb))
  )

  #Function to plot Normalized time using the HYDROTHERMAL time model
  modelNormalized <<-
    alply(as.matrix(TreatmentsTemp), 1, function(Temp) {
      stat_function(fun=function(x){pnorm((+0-(HT/((1)*x))),Psib50,Sigma, log= FALSE)*MaxGerm}, color="blue3")
    })
  NormalizedAxisTitle <<- "Normalized thermal time (°h)"

  #Plot raw data and predicted model with parameters
  PlotPBTMModel()
}

#' A Function to calculate the Supra-optimal hydrothermal time model parameters.
#'
#' This function calculates the supra-optimal hydrothermal time constant (HT), the optimal temperature (To), the median water potential base (Psib50), the standard deviation (sigma) and the slope of the relationship between Psib50 and temperature (kt).
#' @param
#' @keywords supra-optimal hydrothermal time model parameters
#' @export
#' @examples CalcHTTModelSupra()
#' CalcHTTModelSupra()
#'
CalcHTTModelSupra <- function()
{
  #Inform intial and limit values for the Hydrothermal time model parameters
  # Initials
  iHT <- 10
  iTo <- 25
  iPsib50 <- -0.8
  iSigma <- 0.2
  iKt <- 0.5
  #lower limits
  lHT <- 1
  lTo <- 24
  lPsib50 <- -3
  lSigma <- 0.001
  lKt <- 0.001
  #upper limits
  uHT <- 500
  uTo <- 29
  uPsib50 <- -0.01
  uSigma <<- 0.99
  uKt <- 10

  if (FixTo == FALSE) {
    #Calculate Hydrothermaltime Supra-optimal Model Parameters- nls plus algorithm port used to add constraints on the parameters
    HTTModelSupra <<- nls(Germ ~ pnorm(((+WP-Kt*(Temp-To))-(HT/((Temp-MyTb)*Time.hours))),Psib50,Sigma, log= FALSE)*MaxGerm, start=list(HT=iHT,To=iTo,Psib50=iPsib50,Sigma=iSigma,Kt=iKt),lower=list(HT=lHT,To=lTo,Psib50=lPsib50,Sigma=lSigma,Kt=lKt),upper=list(HT=uHT,To=uTo,Psib50=uPsib50,Sigma=uSigma,Kt=uKt), algorithm ="port")
    summary(HTTModelSupra)
    #Passing fitted Hydrothermaltime Supra-optimal Model Parameters
    HT <<- summary(HTTModelSupra)$coefficients[[1]]
    To <<- summary(HTTModelSupra)$coefficients[[2]]
    Psib50 <<- summary(HTTModelSupra)$coefficients[[3]]
    Sigma <<- summary(HTTModelSupra)$coefficients[[4]]
    Kt <<- summary(HTTModelSupra)$coefficients[[5]]
  } else {
    #Calculate Hydrothermaltime Supra-optimal Model Parameters- with fixed Tb
    HTTModelSupra <<- nls(Germ ~ pnorm(((+WP-Kt*(Temp-MyTo))-(HT/((Temp-MyTb)*Time.hours))),Psib50,Sigma, log= FALSE)*MaxGerm, start=list(HT=iHT,Psib50=iPsib50,Sigma=iSigma,Kt=iKt),lower=list(HT=lHT,Psib50=lPsib50,Sigma=lSigma,Kt=lKt),upper=list(HT=uHT,Psib50=uPsib50,Sigma=uSigma,Kt=uKt), algorithm ="port")
    summary(HTTModelSupra)
    HT <<- summary(HTTModelSupra)$coefficients[[1]]
    To <<- MyTo
    Psib50 <<- summary(HTTModelSupra)$coefficients[[2]]
    Sigma <<- summary(HTTModelSupra)$coefficients[[3]]
    Kt <<- summary(HTTModelSupra)$coefficients[[4]]

  }

  #get some estimation of goodness of fit
  Correlation <<- cor(Germ,predict(HTTModelSupra))^2

  #Hydrothermaltime Supra-optimal Model - Create table to plot treatments with predicted model lines
  TreatData$WPFactor <<- with(TreatData, (as.factor(TreatData$Germ.wp)))
  Factor1 <<- TreatData$WPFactor
  Factor1Title <<- "Water \n Potential"
  TreatmentsWP <<- distinct(TreatData, Germ.wp, .keep_all = FALSE)
  TreatData$TempFactor <<- with(TreatData, (as.factor(TreatData$Germ.temp)))
  Factor2 <<- TreatData$TempFactor
  Factor2Title <<- "Temperature"
  TreatmentsTemp <<- distinct(TreatData, Germ.temp, .keep_all = FALSE)
  TreatmentHTT <<- distinct(TreatData, Germ.temp, Germ.wp, .keep_all = FALSE)

  #Passing fitted Hydrothermaltime Supra-optimal Model Parameters for plot legend
  ModPar1Label <<- "HT =="
  ModPar2Label <<- "T[o]=="
  ModPar3Label <<- "Psi[b](50)=="
  ModPar4Label <<- "σ == "
  ModPar5Label <<- "R^2 == "

  ModPar1 <<- round(summary(HTTModelSupra)$coefficients[[1]],2)
  ModPar2 <<- round(To[1],2)
  ModPar3 <<- round(Psib50[1],3)
  ModPar4 <<- round(Sigma[1],3)
  ModPar5 <<- round(Correlation[1],2)

  Dt <<- data.frame(Treatments$Germ.temp,Treatments$Germ.wp)
  tmps <<- Dt$Treatments.Germ.temp
  wps <<- Dt$Treatments.Germ.wp

  #Function to plot all predicted treatments by the Supra-optimal Hydrothermaltime model
  modellines <<-
    mapply(function(Temp1, WP1) {
      stat_function(fun=function(x){pnorm(((+WP1-Kt*(Temp1-To))-(HT/((Temp1-MyTb)*x))),Psib50,Sigma,log=FALSE)*MaxGerm}, aes_(colour = factor(Temp1), alpha = factor(WP1)))
    }, tmps, wps)

  #Create column on TreatDataClean with predicted and normalized values from model
  TreatDataClean <<-TreatDataClean %>% as_tibble() %>% mutate(
    ModelPredicted = pnorm(((+Germ.wp-Kt*(Germ.temp-To))-(HT/((Germ.temp-MyTb)*Germ.time.hours))),Psib50,Sigma, log= FALSE)*MaxGerm,
    NormalizedTime = +((1-Germ.wp/(+Germ.wp-(HT/(Germ.time.hours*(Germ.temp-Tb)))))*Germ.time.hours*(Germ.temp-Tb))
    )


  #Function to plot Normalized time using the HYDROTHERMAL time model
  modelNormalized <<-
    alply(as.matrix(TreatmentsTemp), 1, function(Temp) {
      stat_function(fun=function(x){pnorm((+0-(HT/((1)*x))),Psib50,Sigma, log= FALSE)*MaxGerm}, color="blue3")
    })
  NormalizedAxisTitle <<- "Normalized thermal time (°h)"

  #Plot raw data and predicted model with parameters
  PlotPBTMModel()
}

#' A Function to plots the hydropriming and hydrothermal priming models in a 3D surface graph.
#'
#' This function plots the germination rate (GR50) for each temperature (when more than one temperature was used) in relation to priming duration and water potential.
#' @param
#' @keywords priming models hydropriming hydrothermal priming 3D surface
#' @export
#' @examples PlotPrimingMatrix()
#' PlotPrimingMatrix()
#'
PlotPrimingMatrix <- function()
{
  #Load treatments table and calculate GR50 in case that was not performed earlier.
  CalcT50nGR50()

  #Checking amount of priming temperatures used.


  #Order priming treatments in a proper form for priming matrix
  TreatsPriming <<- Treatments %>% group_by(Treat.priming.wp, Treat.priming.temp,Treat.priming.duration, GR50) %>% tally()
  TreatsPriming <<- TreatsPriming[order(TreatsPriming$Treat.priming.temp,TreatsPriming$Treat.priming.wp,TreatsPriming$Treat.priming.duration),]
  TreatsPriming <- subset(TreatsPriming, Treat.priming.wp < 0)

  #Separate all treatments by priming temperature, water potential and duration and calculate amount of each treatment.
  PrTemps <- unique(TreatsPriming$Treat.priming.temp)
  PrTempAmt <- length(unique(TreatsPriming$Treat.priming.temp))
  PrWPs <- unique(TreatsPriming$Treat.priming.wp)
  PrWPAmt <- length(unique(TreatsPriming$Treat.priming.wp))
  PrDurs <- unique(TreatsPriming$Treat.priming.duration)
  PrDurAmt <- length(unique(TreatsPriming$Treat.priming.duration))

  #List and load all the priming temperature used. Up to 10 at this time. (Build loop in the future)
  TreatsPrTemp1 <- subset(TreatsPriming, Treat.priming.temp == unique(TreatsPriming$Treat.priming.temp)[1])
  TreatsPrTemp2 <- subset(TreatsPriming, Treat.priming.temp == unique(TreatsPriming$Treat.priming.temp)[2])
  TreatsPrTemp3 <- subset(TreatsPriming, Treat.priming.temp == unique(TreatsPriming$Treat.priming.temp)[3])
  TreatsPrTemp4 <- subset(TreatsPriming, Treat.priming.temp == unique(TreatsPriming$Treat.priming.temp)[4])
  TreatsPrTemp5 <- subset(TreatsPriming, Treat.priming.temp == unique(TreatsPriming$Treat.priming.temp)[5])
  TreatsPrTemp6 <- subset(TreatsPriming, Treat.priming.temp == unique(TreatsPriming$Treat.priming.temp)[6])
  TreatsPrTemp7 <- subset(TreatsPriming, Treat.priming.temp == unique(TreatsPriming$Treat.priming.temp)[7])
  TreatsPrTemp8 <- subset(TreatsPriming, Treat.priming.temp == unique(TreatsPriming$Treat.priming.temp)[8])
  TreatsPrTemp9 <- subset(TreatsPriming, Treat.priming.temp == unique(TreatsPriming$Treat.priming.temp)[9])
  TreatsPrTemp10 <- subset(TreatsPriming, Treat.priming.temp == unique(TreatsPriming$Treat.priming.temp)[10])

  #Organize and build priming matrix for all temperatures used for priming.
  Temp1 <- matrix(TreatsPrTemp1$GR50, nrow = PrDurAmt, ncol = PrWPAmt)
  Temp2 <- matrix(TreatsPrTemp2$GR50, nrow = PrDurAmt, ncol = PrWPAmt)
  Temp3 <- matrix(TreatsPrTemp3$GR50, nrow = PrDurAmt, ncol = PrWPAmt)
  Temp4 <- matrix(TreatsPrTemp4$GR50, nrow = PrDurAmt, ncol = PrWPAmt)
  Temp5 <- matrix(TreatsPrTemp5$GR50, nrow = PrDurAmt, ncol = PrWPAmt)
  Temp6 <- matrix(TreatsPrTemp6$GR50, nrow = PrDurAmt, ncol = PrWPAmt)
  Temp7 <- matrix(TreatsPrTemp7$GR50, nrow = PrDurAmt, ncol = PrWPAmt)
  Temp8 <- matrix(TreatsPrTemp8$GR50, nrow = PrDurAmt, ncol = PrWPAmt)
  Temp9 <- matrix(TreatsPrTemp9$GR50, nrow = PrDurAmt, ncol = PrWPAmt)
  Temp10 <- matrix(TreatsPrTemp10$GR50, nrow = PrDurAmt, ncol = PrWPAmt)

  #Plot all priming temperature matrixes. Plotly package used here.
  SurfPlotPrimingMatrix <- plot_ly(showscale = FALSE ) %>%
    add_surface(x = PrWPs , y = PrDurs, z = Temp1, opacity = 1, colorscale = list(c(0,1),c("rgb(19,50,117)","rgb(0,183,255)")), name = paste0(unique(TreatsPriming$Treat.priming.temp)[1], "C")) %>%
    add_surface(x = PrWPs , y = PrDurs, z = Temp2, opacity = 1, colorscale = list(c(0,1),c("rgb(23,46,17)","rgb(45,209,0)")), name = paste0(unique(TreatsPriming$Treat.priming.temp)[2], "C")) %>%
    add_surface(x = PrWPs , y = PrDurs, z = Temp3, opacity = 1, colorscale = list(c(0,1),c("rgb(232,116,49)","rgb(191,204,43)")), name = paste0(unique(TreatsPriming$Treat.priming.temp)[3], "C")) %>%
    add_surface(x = PrWPs , y = PrDurs, z = Temp4, opacity = 1, colorscale = list(c(0,1),c("rgb(209,14,4","rgb(230,97,90)")), name = paste0(unique(TreatsPriming$Treat.priming.temp)[4], "C")) %>%
    add_surface(x = PrWPs , y = PrDurs, z = Temp5, opacity = 1, colorscale = list(c(0,1),c("rgb(77,4,0)","rgb(161,129,127)")), name = paste0(unique(TreatsPriming$Treat.priming.temp)[5], "C")) %>%
    add_surface(x = PrWPs , y = PrDurs, z = Temp6, opacity = 1, colorscale = list(c(0,1),c("rgb(19,50,117)","rgb(0,183,255)")), name = paste0(unique(TreatsPriming$Treat.priming.temp)[6], "C")) %>%
    add_surface(x = PrWPs , y = PrDurs, z = Temp7, opacity = 1, colorscale = list(c(0,1),c("rgb(23,46,17)","rgb(45,209,0)")), name = paste0(unique(TreatsPriming$Treat.priming.temp)[7], "C")) %>%
    add_surface(x = PrWPs , y = PrDurs, z = Temp8, opacity = 1, colorscale = list(c(0,1),c("rgb(232,116,49)","rgb(191,204,43)")), name = paste0(unique(TreatsPriming$Treat.priming.temp)[8], "C")) %>%
    add_surface(x = PrWPs , y = PrDurs, z = Temp9, opacity = 1, colorscale = list(c(0,1),c("rgb(23,46,17)","rgb(45,209,0)")), name = paste0(unique(TreatsPriming$Treat.priming.temp)[9], "C")) %>%
    add_surface(x = PrWPs , y = PrDurs, z = Temp10, opacity = 1, colorscale = list(c(0,1),c("rgb(232,116,49)","rgb(191,204,43)")), name = paste0(unique(TreatsPriming$Treat.priming.temp)[10], "C")) %>%

    layout(
      ## title = "Layout options in a 3d scatter plot",
      scene = list(
        surfacecolor = "rgb(244, 244, 248)",
        xaxis = list(title = "WP"),
        yaxis = list(title = "Duration"),
        zaxis = list(title = "GR50"),
        camera = list(eye = list(x = -1.25, y = 1.25, z = 0.75), center = list(x = 0, y = 0, z = -0.1), up = list(x = 0, y = 0, z = 1)),
        type = "perspective"
      ))


  #Plot priming matrix
  SurfPlotPrimingMatrix
}

