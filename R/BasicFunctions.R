#-----------------Tested Functions - Mac OS and Windows
#
#' A function that calculates the required time to a desired fraction in cumulative curves and respective rate.
#'
#' This function allows you to calculate the time to a desired cumulative fraction and respective germination rate. Use this function on raw data to avoid loss of points closer to the desired cumulative fraction.
#' @param Data time course and cumulative dataset. Several treatments can be used at once as long as it respects the template and column names provided. A column with time in hours (CumTime) + a column with cumulative fractions (CumFraction) are required with at least one additional column for relevant treatment (e.g., germination temperature or water potential)
#' @param Fraction from 0 to 1 used to calculate the time required for that level to be obtained in the cumulative time course. Standard value is 0.5 (50 percent), to calculate the time to 50 percent germination (T50) and respective germination rate (GR50). Fraction level can be entered and be used for calculation and change column name.
#' @param Treat1,Treat2,Treat3,Treat4,Treat5 are the names of the treatment columns to separate the dataset. The time course cumulative curves will be grouped for each distinct treatment that should be informed here. These column names do not need to be informed in case the provided template file is used to organize the data.
#' @keywords Tx, GRx, germination speed, germination rate
#' @importFrom dplyr group_by_at
#' @importFrom dplyr mutate
#' @importFrom dplyr tally
#' @importFrom magrittr %>%
#' @export
#' @examples
#' CalcSpeed(MyData)
CalcSpeed <- function(Data, Fraction, Treat1, Treat2, Treat3, Treat4, Treat5)
{
  #Define the informed treatments (columns) to be grouped when calling the function or set standard to all columns besides CumFraction and CumTime.
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
  Treatments <- Data %>% dplyr::group_by_at(TreatColNames) %>%
    dplyr::mutate(Tx = approx(CumFraction,CumTime, xout=Frac, ties="ordered")$y,
                  GRx = 1/approx(CumFraction,CumTime, xout=Frac, ties="ordered")$y)

  TreatColNames <- c(TreatColNames, "Tx", "GRx")

  # Separate all treatments without germination time courses
  Treatments <- Treatments %>% dplyr::group_by_at(TreatColNames) %>% dplyr::tally()

  #Replaces the column name with the used fraction instead.
  names(Treatments)[names(Treatments) == "Tx"] <- FracSpeedLbl
  names(Treatments)[names(Treatments) == "GRx"] <- FracRateLbl

  return(Treatments)
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

      TreatDataClean <- eval(parse(text=paste("dplyr::distinct(TreatData,",T1,T2,T3,T4,T5,",CumFraction, .keep_all = TRUE)", sep="")))

      return(TreatDataClean)
    }
  }
}

#----------------------New Development - Under Testing
