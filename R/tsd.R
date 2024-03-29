#' Create a S3 `aedseo_tsd` (time-series data) object from observed data
#' and corresponding dates.
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' This function takes observed data and the corresponding date vector and
#' converts it into a `aedseo_tsd` object, which is a time series
#' data structure that can be used for time series analysis.
#'
#' @param observed A numeric vector containing the observations.
#' @param time A date vector containing the corresponding dates.
#' @param time_interval A character vector specifying the time interval.
#' Choose between "day," "week," or "month."
#'
#' @return A `aedseo_tsd` object containing:
#'   - 'time': The time point for for where the observations is observed.
#'   - 'observed': The observed value at the time point.
#'
#' @export
#'
#' @examples
#' # Create a `aedseo_tsd` object from daily data
#' daily_tsd <- tsd(
#'   observed = c(10, 15, 20, 18),
#'   time = as.Date(
#'     c("2023-01-01", "2023-01-02", "2023-01-03", "2023-01-04")
#'   ),
#'   time_interval = "day"
#' )
#'
#' # Create a `aedseo_tsd` object from weekly data
#' weekly_tsd <- tsd(
#'   observed = c(100, 120, 130),
#'   time = as.Date(
#'     c("2023-01-01", "2023-01-08", "2023-01-15")
#'   ),
#'   time_interval = "week"
#' )
#'
#' # Create a `aedseo_tsd` object from monthly data
#' monthly_tsd <- tsd(
#'   observed = c(500, 520, 540),
#'   time = as.Date(
#'     c("2023-01-01", "2023-02-01", "2023-03-01")
#'   ),
#'   time_interval = "month"
#' )
#'
tsd <- function(observed, time, time_interval = c("day", "week", "month")) {
  # Throw an error if any of the inputs are not supported
  time_interval <- rlang::arg_match(time_interval)

  # Collect the input in a tibble
  tbl <- tibble::tibble(
    time = time,
    observed = observed
  )

  # Create the time series data object
  tsd <- tibble::new_tibble(
    x = tbl,
    class = "aedseo_tsd",
    time_interval = time_interval
  )

  return(tsd)
}
