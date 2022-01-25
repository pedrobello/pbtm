#' A function that calculates the required time to a desired fraction in cumulative curves and respective rate.
#'
#' This function allows you to calculate the time to a desired cumulative fraction and respective germination rate. Use this function on raw data to avoid loss of points closer to the desired cumulative fraction.
#'
#' @param data time course and cumulative dataset. Several treatments can be used at once as long as it respects the template and column names provided. A column with time in hours (CumTime) + a column with cumulative fractions (CumFraction) are required with at least one additional column for relevant treatment (e.g., germination temperature or water potential)
#' @param fraction from 0 to 1 used to calculate the time required for that level to be obtained in the cumulative time course. Default value is 0.5 (50 percent), to calculate the time to 50 percent germination (T50) and respective germination rate (GR50). Fraction level can be entered and be used for calculation and change column name.
#' @param treatments are the names of the treatment columns to separate the dataset. The time course cumulative curves will be grouped for each distinct treatment that should be informed here. These column names do not need to be informed in case the provided template file is used to organize the data.
#' @param cum.time is the column name in the dataset for cumulative time elapsed during the germination trial.
#' @param cum.frac is the column name in the data for cumulative fraction germinated.
#' @return A dataframe with germination speed and germination rate calculated from specified treatment columns.
#' @keywords Tx, GRx, germination speed, germination rate
#' @importFrom rlang .data
#' @importFrom dplyr %>%
#' @export
#' @examples
#' load("data/sample_data.RData")
#' mySpeedData <- calcSpeed(myGermData, treatments = c("GermWP", "GermTemp"))
#' print(mySpeedData)

calcSpeed <- function(data, fraction = 0.5, treatments = c("TrtID"), cum.time = "CumTime", cum.frac = "CumFraction") {

  # check if data is valid
  if (!is.data.frame(data)) stop("Data is not a valid data frame.")

  # check validity of fraction argument
  if (!is.numeric(fraction)) stop("Non-numeric fraction value specified.")
  if (fraction < 0 || fraction > 1) stop("Fraction value ", fraction, " is out of range, must be between 0 and 1.")

  # check validity of treatments list
  if (!is.character(treatments) || length(treatments) == 0) stop("Invalid or no treatments specified.")
  trts <- intersect(names(data), treatments)
  if (length(trts) == 0) stop("Specified treatments do not occur in the data frame.")
  if (!setequal(treatments, trts)) warning("Treatments '", paste(setdiff(treatments, names(data)), collapse = ", "), "' do not appear in the data frame.")

  # check validity of observation columns
  if (!is.element(cum.time, names(data))) stop("Cumulative time column '", cum.time, "' not found in data frame.")
  if (!is.element(cum.frac, names(data))) stop("Cumulative fraction column '", cum.frac, "' not found in data frame.")

  # TODO: Maybe should extract this section into its own function.
  # regenerate cumulative fractions depending on grouping trts
  df <- data %>%
    dplyr::group_by_at(trts) %>%
    dplyr::arrange(.data[[cum.time]], .data[[cum.frac]]) %>%
    dplyr::mutate(FracDiff = .data[[cum.frac]] - dplyr::lag(.data[[cum.frac]], default = 0))

  # merge values that occur at the same timepoint
  df <- df %>%
    dplyr::group_by(.data[[cum.time]], .add = T) %>%
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
      Tx = stats::approx(.data$CumFraction, .data[[cum.time]], xout = fraction, ties = "ordered", rule = 2)$y,
      GRx = 1 / .data$Tx,
      .groups = "drop"
    )

  # Rename time and rate columns based on specified fraction.
  names(df)[names(df) == "Tx"] <- paste0("T", fraction * 100)
  names(df)[names(df) == "GRx"] <- paste0("GR", fraction * 100)

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
  names(data)[names(data) == cum.frac] <- "CumFraction"

  # clean the data
  df <- data %>%
    dplyr::distinct(dplyr::across(c(treatments, cum.frac)), .keep_all = T)

  message("Removed ", nrow(data) - nrow(df), " rows with duplicated cumulative fraction values.")

  return(df)
}
