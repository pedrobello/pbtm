#' A T50 and GR50 Function
#'
#' This function allows you to calculate the time to 50% germination (T50) and respective germination rate (GR50).
#' @param CalcT50nGR50() to calculate T50 and GR50 on the TreatData dataset.
#' @keywords T50, GR50, germination speed, germination rate
#' @importFrom magrittr %>%
#' @import dplyr
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

