#' A function that calculates the required time to a desired fraction in cumulative curves and respective rate.
#'
#' This function allows you to calculate the time to a desired cumulative fraction and respective germination rate. Use this function on raw data to avoid loss of points closer to the desired cumulative fraction.
#'
#' @param data time course and cumulative dataset. Several treatments can be used at once as long as it respects the template and column names provided. A column with time in hours (CumTime) + a column with cumulative fractions (CumFraction) are required with at least one additional column for relevant treatment (e.g., germination temperature or water potential)
#' @param fraction from 0 to 1 used to calculate the time required for that level to be obtained in the cumulative time course. Default value is 0.5 (50 percent), to calculate the time to 50 percent germination (T50) and respective germination rate (GR50). Fraction level can be entered and be used for calculation and change column name.
#' @param treatments are the names of the treatment columns to separate the dataset. The time course cumulative curves will be grouped for each distinct treatment that should be informed here. These column names do not need to be informed in case the provided template file is used to organize the data.
#' @param cum.time is the column name in the dataset for cumulative time elapsed during the germination trial.
#' @param cum.frac is the column name in the data for cumulative fraction germinated.
#'
#' @return A dataframe with germination speed and germination rate calculated from specified treatment columns.
#' @keywords Tx, GRx, germination speed, germination rate
#' @importFrom rlang .data
#' @importFrom dplyr %>%
#' @export
#' @examples
#' CalcSpeed(MyData)

CalcSpeed <- function(data, fraction = 0.5, treatments = c("TrtID"), cum.time = "CumTime", cum.frac = "CumFraction") {

  # Argument checks ----

  # check if data is valid
  if (!is.data.frame(data)) stop("Data is not a valid data frame.")

  # check validity of fraction argument
  if (!is.numeric(fraction)) stop("Non-numeric fraction value specified.")
  if (fraction < 0 || fraction > 1) stop("Fraction value ", fraction, " is out of range, must be between 0 and 1.")

  # check validity of treatments list
  if (!is.character(treatments)) stop("Invalid or no treatments list specified.")
  trts <- intersect(names(data), treatments)
  if (length(trts) == 0) stop("Specified treatments do not occur in the data frame.")
  if (!setequal(treatments, trts)) warning("Treatments '", paste(setdiff(treatments, names(data)), collapse = ", "), "' do not appear in the data frame.")

  # check validity of observation columns
  if (!is.element(cum.time, names(data))) stop("Cumulative time column '", cum.time, "' not found in data frame.")
  names(data)[names(data) == cum.time] <- "CumTime"

  if (!is.element(cum.frac, names(data))) stop("Cumulative fraction column '", cum.frac, "' not found in data frame.")
  names(data)[names(data) == cum.frac] <- "CumFraction"


  # Compute germination time and growth rate ----

  # regenerate cumulative fractions depending on grouping trts
  df <- data %>%
    dplyr::group_by_at(trts) %>%
    dplyr::arrange(.data$CumTime, .data$CumFraction) %>%
    dplyr::mutate(FracDiff = .data$CumFraction - dplyr::lag(.data$CumFraction, default = 0))

  # merge values that occur at the same timepoint
  df <- df %>%
    dplyr::group_by(.data$CumTime, .add = T) %>%
    dplyr::summarise(
      FracDiff = sum(.data$FracDiff),
      n = dplyr::n(),
      .groups = "drop_last") %>%
    dplyr::mutate(CumFraction = cumsum(.data$FracDiff) / sum(.data$FracDiff))

  # Calculate the time to the fraction selected using linear interpolation (approx function).
  df <- df %>%
    dplyr::group_by_at(trts) %>%
    dplyr::summarise(
      n = dplyr::n(),
      Tx = stats::approx(.data$CumFraction, .data$CumTime, xout = fraction, ties = "ordered", rule = 2)$y,
      GRx = 1 / .data$Tx,
      .groups = "drop"
    )

  # Rename time and rate columns based on specified fraction.
  names(df)[names(df) == "Tx"] <- paste0("Tx", fraction * 100)
  names(df)[names(df) == "GRx"] <- paste0("GRx", fraction * 100)

  return(df)
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
    print("Missing the data object.")
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

      TreatDataClean <- eval(
        parse(
          text = paste("dplyr::distinct(TreatData,", T1, T2, T3, T4, T5, ",CumFraction, .keep_all = TRUE)", sep = "")
        )
      )

      return(TreatDataClean)
    }
  }
}
