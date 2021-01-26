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
#' @param TimeColName,CumFractColName are the names of the data columns that will be used to calculate values. In the case the column names for time (Germ.time.hours) and cumulative fraction (Germ.fraction) used are equal to the template provided, these column names to not need to be provided.
#' @keywords T50, GR50, germination speed, germination rate
#' @importFrom dplyr group_by_at
#' @importFrom dplyr tally
#' @importFrom magrittr %>%
#' @export
#' @examples
#' CalcT50nGR50(MyData)
CalcT50nGR50 <- function(Data, T1ColName, T2ColName, T3ColName, T4ColName, T5ColName, TimeColName, CumFractColName)
{
  if (missing(T5ColName)) { #T5ColName not informed
    if (missing(T4ColName)) { #T4ColName not informed
      if (missing(T3ColName)) { #T3ColName not informed
        if (missing(T2ColName)) { #T2ColName not informed
          if (missing(T1ColName)) { #T1ColName not informed
            TreatColNames <- c("Treat.ID","Treat.desc","Treat.aging.time","Treat.priming.wp","Treat.priming.temp","Treat.priming.duration","Germ.wp","Germ.temp", "Germ.promoter.dosage", "Germ.inhibitor.dosage")
            } else {TreatColNames <- T1ColName}
        } else {TreatColNames <- paste(T1ColName,",",T2ColName)}
      } else {TreatColNames <- paste(T1ColName,",",T2ColName,",",T3ColName)}
    } else {TreatColNames <- paste(T1ColName,",",T2ColName,",",T3ColName,",",T4ColName)}
  } else {TreatColNames <- paste(T1ColName,",",T2ColName,",",T3ColName,",",T4ColName,",",T5ColName)}

  print(TreatColNames)

  if (missing(TimeColName)) { #TimeColName not informed
    Germ.time.hours <- "Germ.time.hours"
  } else {
    Germ.time.hours <- TimeColName
  }

  if (missing(CumFractColName)) { #CumFractColName not informed
    Germ.fraction <- "Germ.fraction"
  } else {
    Germ.fraction <- CumFractColName
  }

  # Calculate Time to 50% Germination (T50) (calculate on raw data to avoid loss of points closer to 50% germination) + GR50
  Treatments <- Data %>% group_by_at(TreatColNames) %>%
    dplyr::mutate(T50 = approx(Germ.fraction,Germ.time.hours, xout=0.5, ties="ordered")$y,
                  GR50 = 1/approx(Germ.fraction,Germ.time.hours, xout=0.5, ties="ordered")$y)

  # Separate all treatments without germination time courses
  Treatments <- Treatments %>% group_by_at(TreatColNames, "T50", "GR50") %>% tally()
  return(Treatments)
}

