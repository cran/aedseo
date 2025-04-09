## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, warning=FALSE, message=FALSE--------------------------------------
library(aedseo)

## -----------------------------------------------------------------------------
seasonal_wave_sim_weekly <- generate_seasonal_data(
  years = 3,
  start_date = as.Date("2021-05-26"),
  amplitude = 100,
  mean = 100,
  trend_rate = 1.003,
  time_interval = "week"
)

## ----dpi=300------------------------------------------------------------------
plot(seasonal_wave_sim_weekly, time_interval = "5 weeks")

## ----dpi=300------------------------------------------------------------------
seasonal_wave_sim_monthly <- generate_seasonal_data(
  years = 4,
  start_date = as.Date("2021-05-26"),
  amplitude = 50,
  mean = 50,
  trend_rate = 0.99,
  time_interval = "month"
)
plot(
  seasonal_wave_sim_monthly,
  time_interval = "3 months",
  y_label = "Monthly observations"
)

## ----dpi=300------------------------------------------------------------------
seasonal_wave_sim_daily <- generate_seasonal_data(
  years = 3,
  start_date = as.Date("2021-05-26"),
  amplitude = 50,
  mean = 50,
  time_interval = "day"
)
plot(
  seasonal_wave_sim_daily,
  time_interval = "50 days",
  y_label = "Daily observations"
)

## ----dpi=300------------------------------------------------------------------
seasonal_wave_sim_daily_phase_shift <- generate_seasonal_data(
  years = 3,
  start_date = as.Date("2021-05-26"),
  amplitude = 50,
  mean = 50,
  phase = 1,
  time_interval = "day"
)
plot(
  seasonal_wave_sim_daily_phase_shift,
  time_interval = "50 days",
  y_label = "Daily observations"
)

## ----dpi=300------------------------------------------------------------------
sim_no_noise <- generate_seasonal_data(
  years = 3,
  start_date = as.Date("2021-05-26"),
  amplitude = 100,
  mean = 100,
  noise_overdispersion = 0,
  time_interval = "week"
)
plot(
  sim_no_noise,
  time_interval = "5 weeks"
)

## ----dpi=300------------------------------------------------------------------
sim_poisson_noise <- generate_seasonal_data(
  years = 3,
  start_date = as.Date("2021-05-26"),
  amplitude = 100,
  mean = 100,
  noise_overdispersion = 1,
  time_interval = "week"
)
plot(
  sim_poisson_noise,
  time_interval = "5 weeks"
)

## ----dpi=300------------------------------------------------------------------
sim_nb_noise <- generate_seasonal_data(
  years = 3,
  start_date = as.Date("2021-05-26"),
  amplitude = 100,
  mean = 100,
  noise_overdispersion = 5,
  time_interval = "week"
)
plot(
  sim_nb_noise,
  time_interval = "5 weeks"
)

## ----dpi=300------------------------------------------------------------------
sim_sinus <- generate_seasonal_data(
  years = 2,
  start_date = as.Date("2021-05-26"),
  amplitude = 100,
  mean = 100,
  relative_epidemic_concentration = 1,
  time_interval = "week"
)
plot(
  sim_sinus,
  time_interval = "5 weeks"
)

## ----dpi=300------------------------------------------------------------------
sim_conc <- generate_seasonal_data(
  years = 2,
  start_date = as.Date("2021-05-26"),
  amplitude = 100,
  mean = 100,
  relative_epidemic_concentration = 4,
  time_interval = "week"
)
plot(
  sim_conc,
  time_interval = "5 weeks"
)

