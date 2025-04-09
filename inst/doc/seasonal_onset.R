## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, warning=FALSE, message=FALSE--------------------------------------
library(aedseo)

## -----------------------------------------------------------------------------
# Construct an 'tsd' object with time series data
set.seed(222)
tsd_data <- generate_seasonal_data(
  years = 3,
  start_date = as.Date("2020-10-18"),
  trend_rate = 1.002,
  noise_overdispersion = 3,
  relative_epidemic_concentration = 2,
  time_interval = "week"
)

## -----------------------------------------------------------------------------
seasonal_onset_results <- seasonal_onset(
  tsd = tsd_data,
  k = 5,
  level = 0.95,
  disease_threshold = 20,
  family = "quasipoisson",
  season_start = 21,
  season_end = 20,
  only_current_season = FALSE
)

## ----dpi=300------------------------------------------------------------------
plot(seasonal_onset_results)

## -----------------------------------------------------------------------------
prediction <- predict(seasonal_onset_results, n_step = 5)

## ----echo = FALSE-------------------------------------------------------------
prediction <- prediction |>
  dplyr::select(-t) |>
  dplyr::rename(observation = estimate)

# Extract two seasons for better visualisation:
seasonal_onset_short <- seasonal_onset_results |>
  dplyr::filter(season %in% c("2023/2024", "2024/2025"))

# Plot observations and predictions
autoplot(seasonal_onset_short)$observed +
  # Add the prediction ribbon
  ggplot2::geom_ribbon(
    data = prediction,
    ggplot2::aes(
      x = reference_time,
      ymin = lower,
      ymax = upper,
      fill = "Observations with \n Confidence intervals"
    ),
    alpha = 0.2,
    inherit.aes = FALSE
  ) +
  ggplot2::geom_line(
    data = prediction,
    ggplot2::aes(
      x = reference_time,
      y = observation,
      color = "Observations with \n Confidence intervals"
    ),
    linewidth = 1,
    inherit.aes = FALSE
  ) +
  ggplot2::scale_color_manual(
    name = "Predictions",
    values = c("Observations with \n Confidence intervals" = "red")
  ) +
  ggplot2::scale_fill_manual(
    name = "Predictions",
    values = c("Observations with \n Confidence intervals" = "red")
  ) +
  ggplot2::coord_cartesian(
    xlim = c(
      min(seasonal_onset_short$reference_time),
      max(prediction$reference_time)
    )
  )

## ----summary------------------------------------------------------------------
summary(seasonal_onset_results)

