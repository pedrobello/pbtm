#' A Hello, world! Function
#'
#' This is an example function named 'hello'
#' @param hello()
#' @export
#' @examples
#' hello()
hello <- function() {
  print("Hello, world this is a test function!")
}

#' A T50 and GR50 Function
#'
#' This function allows you to calculate the time to 50% germination (T50) and respective germination rate (GR50).
#' @param Data time course and cumulative dataset. A column with time in hours (Germ.time.hours) + a column with cumulative fractions (Germ.fraction) are required with at least one additional column for revelant treatment (e.g., germination temperature or water potential)
#' @param T1ColName,T2ColName,T3ColName,T4ColName,T5ColName are the names of the treatment columns to separate the dataset. The time course cumulative curves will be grouped for each distinct treatment that should be informed here. These column names do not need to be informed in case the provided template file is used to organize the data.
#' @keywords T50, GR50, germination speed, germination rate
#' @importFrom dplyr group_by_at
#' @importFrom dplyr tally
#' @importFrom magrittr %>%
#' @export
#' @examples
#' CalcT50nGR50(MyData)
CalcT50nGR50 <- function(Data, T1ColName, T2ColName, T3ColName, T4ColName, T5ColName)
{
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

  # Calculate Time to 50% Germination (T50) (calculate on raw data to avoid loss of points closer to 50% germination) + GR50
  Treatments <- Data %>% group_by_at(TreatColNames) %>%
    dplyr::mutate(T50 = approx(CumFract,CumTime, xout=0.5, ties="ordered")$y,
                  GR50 = 1/approx(CumFract,CumTime, xout=0.5, ties="ordered")$y)

  TreatColNames <- c(TreatColNames, "T50", "GR50")

  # Separate all treatments without germination time courses
  Treatments <- Treatments %>% group_by_at(TreatColNames) %>% tally()
  return(Treatments)
}

#' A Tx and GRx Function
#'
#' This function allows you to calculate the time to a desired germination fraction (Tx) and respective germination rate (GRx).
#' @param Data time course and cumulative dataset. A column with time in hours (Germ.time.hours) + a column with cumulative fractions (Germ.fraction) are required with at least one additional column for revelant treatment (e.g., germination temperature or water potential)
#' @param T1ColName,T2ColName,T3ColName,T4ColName,T5ColName are the names of the treatment columns to separate the dataset. The time course cumulative curves will be grouped for each distinct treatment that should be informed here. These column names do not need to be informed in case the provided template file is used to organize the data.
#' @keywords Tx, GRx, germination speed, germination rate
#' @importFrom dplyr group_by_at
#' @importFrom dplyr tally
#' @importFrom magrittr %>%
#' @export
#' @examples
#' CalcTxnGRx(MyData)
CalcTxnGRx <- function(Data, Fraction, T1ColName, T2ColName, T3ColName, T4ColName, T5ColName)
{
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

  if (missing(Fraction)) { #Fraction not informed
    Frac <- 0.5
    FracSpeedLbl <- "T50"
    FracRateLbl <- "GR50"
  } else {
    Frac <- Fraction
    FracSpeedLbl <- paste("T",(Frac*100))
    FracRateLbl <- paste("GR",(Frac*100))
  }

  # Calculate Time to 50% Germination (T50) (calculate on raw data to avoid loss of points closer to 50% germination) + GR50
  Treatments <- Data %>% group_by_at(TreatColNames) %>%
    dplyr::mutate(Tx = approx(CumFract,CumTime, xout=Frac, ties="ordered")$y,
                  GRx = 1/approx(CumFract,CumTime, xout=Frac, ties="ordered")$y)


  TreatColNames <- c(TreatColNames, "Tx", "GRx")

  # Separate all treatments without germination time courses
  Treatments <- Treatments %>% group_by_at(TreatColNames) %>% tally()

  names(Treatments)[names(Treatments) == "Tx"] <- FracSpeedLbl
  names(Treatments)[names(Treatments) == "GRx"] <- FracRateLbl

  return(Treatments)
}
