## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, warning=FALSE, message=FALSE--------------------------------------
library(aedseo)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)

## -----------------------------------------------------------------------------
mu_t <- function(
    t,
    theta = 1,
    exp_beta = 1.001,
    gamma_sin = 1,
    gamma_cos = 1,
    trend = 1,
    j = 1,
    ...) {
  # Start construction of linear predictor
  linear_predictor <- theta
  # ... add a trend if the scenario request it
  if (trend == 1) {
    linear_predictor <- linear_predictor + log(exp_beta) * t
  }
  # ... add a seasonal component
  linear_predictor <- linear_predictor +
    gamma_sin * sin(2 * pi * t * j / 52) + gamma_cos * cos(2 * pi * t * j / 52)

  return(exp(linear_predictor))
}

simulate_from_nbinom <- function(...) {
  # Set some default values for the simulation
  default_pars <- list(
    "t" = 1,
    "theta" = 1,
    "exp_beta" = 1.001,
    "gamma_sin" = 1,
    "gamma_cos" = 1,
    "trend" = 1,
    "j" = 1,
    "phi" = 1,
    "seed" = 42
  )

  # Match call
  mc <- as.list(match.call())[-1]
  # ... and change parameters relative to the call
  index_mc <- !names(default_pars) %in% names(mc)
  mc <- append(mc, default_pars[index_mc])[names(default_pars)]

  # Set the seed
  set.seed(mc$seed)

  # Calculate the number of time points
  n <- length(eval(mc$t))
  # Calculate mu_t
  mu_t_scenario <- do.call(what = "mu_t", args = mc)
  # ... and compute the variance of the negative binomial distribution
  variance <- mu_t_scenario * mc$phi
  # ... and infer the size in the negative binomial distribution
  size <- (mu_t_scenario + mu_t_scenario^2) / variance
  # Plugin and simulate the data
  simulation <- rnbinom(n = n, mu = mu_t_scenario, size = size)

  return(list("mu_t" = mu_t_scenario, "simulation" = simulation, "pars" = mc))
}

# Define the number of years and the number of months within a year
years <- 3
weeks <- 52
# ... calculate the total number of observations
n <- years * weeks
# ... and a vector containing the overall time passed
time_overall <- 1:n
# Create arbitrary dates
dates <- seq.Date(from = as.Date("2010-01-01"), by = "week", length.out = n)


# Simulate the data
simulation <- simulate_from_nbinom(t = time_overall, theta = log(1000), phi = 5)

# Collect the data in a 'tibble'
sim_data <- tibble(
  Date = dates,
  mu_t = simulation$mu_t,
  y = simulation$simulation
)

## -----------------------------------------------------------------------------
# Construct an 'aedseo_tsd' object with the time series data
tsd_data <- tsd(
  observed = sim_data$y,
  time = sim_data$Date,
  time_interval = "week"
)

## -----------------------------------------------------------------------------
# Have a glance at the time varying mean and the simulated data
plot(tsd_data)

## -----------------------------------------------------------------------------
# Apply the 'aedseo' algorithm
aedseo_results <- aedseo(
  tsd = tsd_data,
  k = 5,
  level = 0.95,
  disease_threshold = 2000,
  family = "quasipoisson"
)

## -----------------------------------------------------------------------------
# Employing the S3 `plot()` method
plot(aedseo_results)

## -----------------------------------------------------------------------------
# Example: Predict growth rates for the next 5 time steps
(prediction <- predict(aedseo_results, n_step = 5))

## ----summary------------------------------------------------------------------
summary(aedseo_results)

