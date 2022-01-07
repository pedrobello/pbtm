#-----------------Tested Functions - Mac OS and Windows
#
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
  pRaw <<- ggplot(data=TreatData, aes_string(x="CumTime", y="CumFraction", color=T1, shape=T2 )) +
    eval(parse(text=gp)) + geom_line() + xlab("Time") + ylab("Cumulative (%)") +
    scale_y_continuous(labels = scales::percent, expand = c(0,0), limits = c(0,1.02)) +
    scale_x_continuous(expand = c(0,0)) + expand_limits(x = 0, y = 0) +
    theme_scatter_plot
  pRaw
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
