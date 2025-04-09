## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, warning=FALSE, message=FALSE--------------------------------------
library(aedseo)

## ----include = FALSE----------------------------------------------------------
withr::local_seed(222)
# Construct an 'tsd' object with time series data
tsd_data <- generate_seasonal_data(
  years = 5,
  start_date = as.Date("2020-10-18"),
  trend_rate = 1.002,
  noise_overdispersion = 5,
  relative_epidemic_concentration = 3,
  time_interval = "week"
)

tsd_data <- tsd_data |>
  dplyr::filter(time <= "2025-02-01")


## ----dpi=300------------------------------------------------------------------
plot(tsd_data)

## ----dpi=300------------------------------------------------------------------
tsd_onset <- seasonal_onset(
  tsd = tsd_data,
  k = 5,
  family = "quasipoisson",
  na_fraction_allowed = 0.4,
  season_start = 21, # Season starts in week 21
  season_end = 20, # Season ends in week 20 the following year
  only_current_season = FALSE
)

consecutive_gr_warn <- consecutive_growth_warnings(
  onset_output = tsd_onset
)

autoplot(
  consecutive_gr_warn,
  k = 5,
  skip_current_season = TRUE
) +
  ggplot2::geom_vline(
    ggplot2::aes(xintercept = 25, linetype = "Threshold"),
    color = "black", linewidth = 0.6
  ) +
  ggplot2::scale_linetype_manual(
    name   = "",
    values = c("Threshold" = "dashed")
  )

## -----------------------------------------------------------------------------
consecutive_gr_warn |>
  dplyr::filter(!is.na(significant_counter)) |>
  dplyr::filter(season != max(consecutive_gr_warn$season)) |>
  dplyr::group_by(season) |>
  dplyr::filter(significant_counter == max(significant_counter)) |>
  dplyr::mutate(disease_threshold = sum_of_cases / 5,
                week = ISOweek::ISOweek(reference_time)) |>
  dplyr::select(season, week, disease_threshold)

## -----------------------------------------------------------------------------
seasonal_output <- combined_seasonal_output(
  tsd = tsd_data,
  disease_threshold = 25,
  method = "intensity_levels",
  family = "quasipoisson"
)

## -----------------------------------------------------------------------------
summary(seasonal_output$onset_output)

## -----------------------------------------------------------------------------
summary(seasonal_output$burden_output)

## ----dpi=300------------------------------------------------------------------
# Adjust y_lower_bound dynamically to remove noisy small values
disease_threshold <- 25
y_lower_bound <- ifelse(disease_threshold < 10, 1, 5)

plot(
  x = seasonal_output,
  y_lower_bound = y_lower_bound,
  time_interval = "3 weeks"
)

## ----dpi=300------------------------------------------------------------------
# Get `tsd_onset` object
tsd_onset <- seasonal_onset(
  tsd = tsd_data,
  disease_threshold = 25,
  family = "quasipoisson",
  season_start = 21,
  season_end = 20,
  only_current_season = FALSE
)

historical_summary(tsd_onset)

