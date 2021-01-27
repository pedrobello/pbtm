#-----------------Tested Functions - MAc OS and Windows
#
#' A function that calculates the required time and rate to a desired fraction in cumulative curves
#'
#' This function allows you to calculate the time to a desired cumulative fraction and respective germination rate. Use this function on raw data to avoid loss of points closer to the desired cumulative fraction.
#' @param Data time course and cumulative dataset. Several treatments can be used at once as long as it respects the template and column names provided. A column with time in hours (Germ.time.hours) + a column with cumulative fractions (Germ.fraction) are required with at least one additional column for relevant treatment (e.g., germination temperature or water potential)
#' @param Fraction from 0 to 1 used to calculate the time required for that level to be obtained in the cumulative time course. Standard value is 0.5 (50 percent), to calculate the time to 50 percent germination (T50) and respective germination rate (GR50). Fraction level can be entered and be used for calculation and change column name.
#' @param T1ColName,T2ColName,T3ColName,T4ColName,T5ColName are the names of the treatment columns to separate the dataset. The time course cumulative curves will be grouped for each distinct treatment that should be informed here. These column names do not need to be informed in case the provided template file is used to organize the data.
#' @keywords Tx, GRx, germination speed, germination rate
#' @importFrom dplyr group_by_at
#' @importFrom dplyr tally
#' @importFrom magrittr %>%
#' @export
#' @examples
#' CalcSpeed(MyData)
CalcSpeed <- function(Data, Fraction, T1ColName, T2ColName, T3ColName, T4ColName, T5ColName)
{
  #Define the informed treatments (columns) to be grouped when calling the function or set standard to all columns besides CumFract and CumTime.
  if (missing(T5ColName)) { #T5ColName not informed
    if (missing(T4ColName)) { #T4ColName not informed
      if (missing(T3ColName)) { #T3ColName not informed
        if (missing(T2ColName)) { #T2ColName not informed
          if (missing(T1ColName)) { #T1ColName not informed
            TreatColNames <- c("Treat.ID","Treat.desc","Treat.aging.time","Treat.priming.wp","Treat.priming.temp","Treat.priming.duration","Germ.wp","Germ.temp", "Germ.promoter.dosage", "Germ.inhibitor.dosage")
          } else {TreatColNames <- c(T1ColName)}
        } else {TreatColNames <- c(T1ColName,T2ColName)}
      } else {TreatColNames <- c(T1ColName,T2ColName,T3ColName)}
    } else {TreatColNames <- c(T1ColName,T2ColName,T3ColName,T4ColName)}
  } else {TreatColNames <- c(T1ColName,T2ColName,T3ColName,T4ColName,T5ColName)}

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

#----------------------New Development - Under Testing

#' A Function to plot rate vs treatment temperature.
#'
#' This function plots rates against the desired treatment.
#' @param Data should inform the table with time course data to be plotted
#' @param x should indicate the treatment column name for the x axis.
#' @param y should indicate the rate column name for the y axis if different than GR50.
#' @keywords plot GR50 Temperature
#' @importFrom ggplot2 ggplot
#' @export
#' @examples PlotRateVsTreat()
#' PlotRateVsTreat()
PlotRateVsTreat <- function (Data, x, y)
{
  Treatments <- Data
  Treat <- x
  if (missing(y)) { #y/rate not informed
    rate <- "GR50"
  } else {
    rate <- y
  }
  pGR <- ggplot(data=Treatments, aes(x=Treat, y=rate)) + geom_point(shape=19, size=2) + xlab("Temperature (°C)") +
    ylab(bquote('Rate ('*h^-1*')')) +
    expand_limits(x = 0, y = 0) + theme_scatter_plot
  pGR
}
