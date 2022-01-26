#' A function that calculates the required time to a desired fraction in cumulative curves and respective rate.
#'
#' This function allows you to calculate the time to a desired cumulative fraction and respective germination rate. Use this function on raw data to avoid loss of points closer to the desired cumulative fraction.
#'
#' @param data time course and cumulative dataset. Several treatments can be used at once as long as it respects the template and column names provided. A column with time in hours (CumTime) + a column with cumulative fractions (CumFraction) are required with at least one additional column for relevant treatment (e.g., germination temperature or water potential)
#' @param treatments are the names of the treatment columns to separate the dataset. The time course cumulative curves will be grouped for each distinct treatment that should be informed here. These column names do not need to be informed in case the provided template file is used to organize the data.
#' @param f A single value or vector of values from 0 to 1 used to calculate the time required for that level to be obtained in the cumulative time course. Default value is 0.5 (50 percent), to calculate the time to 50 percent germination (T50) and respective germination rate (GR50).
#' @param cum.time is the column name in the dataset for cumulative time elapsed during the germination trial.
#' @param cum.frac is the column name in the data for cumulative fraction germinated.
#' @return A dataframe with germination speed and germination rate calculated from specified treatment columns.
#' @keywords Tx, GRx, germination speed, germination rate
#' @importFrom rlang .data
#' @importFrom rlang :=
#' @importFrom dplyr %>%
#' @export
#' @examples
#' "foo"

calcSpeed <- function(data, treatments = c("TrtID"), f = 0.5, cum.time = "CumTime", cum.frac = "CumFraction") {

  # check if data is valid
  if (!is.data.frame(data)) stop("Data is not a valid data frame.")

  # check validity of fraction argument
  if (!is.numeric(f)) stop("Non-numeric fraction value specified.")
  f <- sort(unique(f))
  lapply(f, function(i) {
    if (i < 0 || i > 1) stop("Fraction value ", i, " is out of range, must be between 0 and 1.")
  })

  # check validity of treatments list
  if (!is.character(treatments) || length(treatments) == 0) stop("Invalid or no treatments specified.")
  trts <- intersect(names(data), treatments)
  if (length(trts) == 0) stop("Specified treatments do not occur in the data frame.")
  if (!setequal(treatments, trts)) warning("Treatments '", paste(setdiff(treatments, names(data)), collapse = ", "), "' do not appear in the data frame.")

  # check validity of observation columns
  if (!is.element(cum.time, names(data))) stop("Cumulative time column '", cum.time, "' not found in data frame.")
  if (!is.element(cum.frac, names(data))) stop("Cumulative fraction column '", cum.frac, "' not found in data frame.")

  # Calculate the time to the fraction selected using linear interpolation (approx function).
  if (length(f) == 1) {
    df <- data %>%
      dplyr::group_by_at(trts) %>%
      dplyr::summarise(
        n = dplyr::n(),
        !!eval(paste0("T", f * 100)) := stats::approx(.data[[cum.frac]], .data[[cum.time]], xout = f, ties = "ordered", rule = 2)$y,
        !!eval(paste0("GR", f * 100)) := 1 / .data[[paste0("T", f * 100)]],
        .groups = "drop"
      )
  } else {
    df <- data %>%
      dplyr::group_by_at(trts) %>%
      dplyr::summarise(
        {
          df <- stats::approx(.data[[cum.frac]], .data[[cum.time]], xout = f, ties = "ordered", rule = 2)
          names(df) <- c("Frac", "T")
          df <- dplyr::as_tibble(df)
          tidyr::drop_na(df)
        },
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        GR = 1 / .data$T,
        Frac = .data$Frac * 100) %>%
      tidyr::pivot_longer(
        cols = c("T", "GR")
      ) %>%
      dplyr::mutate(name = paste0(.data$name, .data$Frac)) %>%
      dplyr::select(-"Frac") %>%
      tidyr::pivot_wider()
  }

  return(df)
}


# TODO: Seems like we could eliminate this function if all it does is run dplyr::distinct. The user can do that themselves.
#' A Function to clean cumulative curves on dataset.
#'
#' This function removes repetitive cumulative fractions/percentages, keeping only the initial presence of the value
#' @param data A data frame object with the raw cumulative data that needs to be cleaned.
#' @param treatments are the names of the treatment columns to separate the dataset. The time course cumulative curves will be grouped for each distinct treatment that should be informed here.
#' @param cum.frac Column name for cumulative fraction germinated.
#' @return A data frame with duplicate observations removed
#' @keywords clean cumulative fraction repetitive percentage
#' @importFrom rlang .data
#' @importFrom dplyr %>%
#' @export
#' @examples
#' "foo"

cleanData <- function(data, treatments = c("TrtID"), cum.frac = "CumFraction") {

  # check if data is valid
  if (!is.data.frame(data)) stop("Data is not a valid data frame.")

  # check validity of treatments list
  if (!is.character(treatments) || length(treatments) == 0) stop("Invalid or no treatments specified.")
  trts <- intersect(names(data), treatments)
  if (length(trts) == 0) stop("Specified treatments do not occur in the data frame.")
  if (!setequal(treatments, trts)) warning("Treatments '", paste(setdiff(treatments, names(data)), collapse = ", "), "' do not appear in the data frame.")

  # check the cumulative fraction column
  if (!is.element(cum.frac, names(data))) stop("Cumulative fraction column '", cum.frac, "' not found in data frame.")

  # clean the data
  df <- data %>%
    dplyr::distinct(dplyr::across(c(treatments, cum.frac)), .keep_all = T)

  message("Removed ", nrow(data) - nrow(df), " rows with duplicated cumulative fraction values.")

  return(df)
}


#' A function to generate cumulative germination fractions from raw seed counts.
#'
#' @param data A dataset containing seed germination observations.
#' @param trt.id Column name containing ids which reflect the individual treatments.
#' @param cum.time Column name containing cumulative elapsted time.
#' @param n.total Total number of seeds for the given treatment.
#' @param n.germinated Column name containing total number of seeds germinated through the current time.
#' @return a data frame with cumulative fractions generated from raw seed counts.
#' @keywords seed germination cumulative fraction
#' @importFrom rlang .data
#' @importFrom dplyr %>%
#' @export
#' @examples
#' "foo"

calcCumFrac <- function(data, trt.id = "TrtID", cum.time = "CumTime", n.total = "nTotal", n.germinated = "nGerminated") {

  # check if data is valid
  if (!is.data.frame(data)) stop("Data is not a valid data frame.")

  # check columns
  if (!is.element(trt.id, names(data))) stop("Treatment ID column '", trt.id, "' not found in the data frame.")
  if (!is.element(cum.time, names(data))) stop("Cumulative time column '", cum.time, "' not found in the data frame.")
  if (!is.element(n.total, names(data))) stop("Total seed count column '", n.total, "' not found in the data frame.")
  if (!is.element(n.germinated, names(data))) stop("Germinated seed count column '", n.germinated, "' not found in the data frame.")

  df <- data %>%
    dplyr::group_by(.data[[trt.id]]) %>%
    dplyr::arrange(.data[[cum.time]], .data[[n.germinated]]) %>%
    dplyr::distinct(.data[[cum.time]], .keep_all = T) %>%
    dplyr::mutate(CumFraction = .data[[n.germinated]] / .data[[n.total]]) %>%
    dplyr::select(-c(n.germinated, n.total))

  return(df)
}


#' A function to regenerate cumulative fractions by merging treatments.
#'
#' @param data A data frame object with the raw cumulative data that needs to be cleaned.
#' @param keep.trts One or more treatments to keep from the original set of treatments.
#' @param trt.id The column specifying the original set of treatments used to generate the cumulative fractions.
#' @param cum.time Column name for the cumulative time data.
#' @param cum.frac Column name for the cumulative fraction germinated.
#' @return A data frame with new cumulative fraction values.
#' @keywords cumulative fraction merge treatments
#' @importFrom rlang .data
#' @importFrom rlang :=
#' @importFrom dplyr %>%
#' @export
#' @examples
#' "foo"

mergeTrts <- function(data, keep.trts, trt.id = "TrtID", cum.time = "CumTime", cum.frac = "CumFraction") {

  # check if data is valid
  if (!is.data.frame(data)) stop("Data is not a valid data frame.")

  # check validity of treatments list
  if (!is.character(keep.trts) || length(keep.trts) == 0) stop("Invalid or no treatments specified.")
  trts <- intersect(names(data), keep.trts)
  if (length(trts) == 0) stop("Specified treatments do not occur in the data frame.")
  if (!setequal(keep.trts, trts)) warning("Treatments '", paste(setdiff(keep.trts, names(data)), collapse = ", "), "' do not appear in the data frame.")

  # check validity of observation columns
  if (!is.element(trt.id, names(data))) stop("Treatment ID column '", trt.id, "' not found in data frame.")
  if (!is.element(cum.time, names(data))) stop("Cumulative time column '", cum.time, "' not found in data frame.")
  if (!is.element(cum.frac, names(data))) stop("Cumulative fraction column '", cum.frac, "' not found in data frame.")


  # find the difference in cumulative fraction between data points
  df <- data %>%
    dplyr::group_by(.data[[trt.id]]) %>%
    dplyr::arrange(.data[[trt.id]], .data[[cum.time]], .data[[cum.frac]]) %>%
    dplyr::mutate(FracDiff = .data[[cum.frac]] - dplyr::lag(.data[[cum.frac]], default = 0))

  # generate new cumulative fractions based on kept treatments
  df <- df %>%
    dplyr::group_by(dplyr::across(c(keep.trts, cum.time))) %>%
    dplyr::summarise(FracDiff = sum(.data$FracDiff), .groups = "drop_last") %>%
    dplyr::mutate(!!cum.frac := cumsum(.data$FracDiff) / sum(.data$FracDiff)) %>%
    dplyr::select(-"FracDiff")

  return(df)
}
