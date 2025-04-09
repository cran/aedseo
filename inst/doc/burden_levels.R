## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, warning=FALSE, message=FALSE--------------------------------------
library(aedseo)

## ----echo = FALSE, dpi=300----------------------------------------------------
withr::local_seed(222)
# Construct an 'tsd' object with time series data
tsd_data_clean <- generate_seasonal_data(
  years = 3,
  start_date = as.Date("2021-10-18"),
  noise_overdispersion = 1,
  relative_epidemic_concentration = 3,
  time_interval = "week"
)

# Run models
intensity_levels <- seasonal_burden_levels(
  tsd = tsd_data_clean,
  disease_threshold = 10,
  method = "intensity_levels",
  conf_levels = 0.95
)
peak_levels <- seasonal_burden_levels(
  tsd = tsd_data_clean,
  disease_threshold = 10,
  method = "peak_levels",
  conf_levels = c(0.4, 0.9, 0.975),
  n_peak = 8
)

# Create data frame
burden_levels_df <- data.frame(
  Level = names(
    c(intensity_levels$values,
      peak_levels$values
    )
  ),
  Threshold = c(
    intensity_levels$values,
    peak_levels$values
  ),
  Method = c(
    rep("Intensity Levels", 4),
    rep("Peak Levels", 4)
  )
)
burden_levels_df$Level <- factor(
  burden_levels_df$Level,
  levels = c("high", "medium", "low", "very low")
)

# Calculate y_tics
y_tics_log10 <- pretty(c(log10(burden_levels_df$Threshold)))
y_tics_levels <- 10^(y_tics_log10)

# For each tic, find the closest magnitude to round correctly
round_to_nearest <- function(x) {
  magnitude <- 10^floor(log10(x))
  plyr::round_any(x, accuracy = magnitude)
}
y_tics <- sapply(y_tics_levels, round_to_nearest)

# Round max value
rounded_max_threshold <- (round(max(burden_levels_df$Threshold) / 100) * 100) + 100
y_tics <- c(y_tics[-5], rounded_max_threshold)

# Create the plot
burden_levels_df |>
  ggplot2::ggplot(ggplot2::aes(x = 0, y = Threshold, linetype = Level, color = Level)) +
  ggplot2::geom_hline(
    ggplot2::aes(yintercept = Threshold, linetype = Level, color = Level),
    linewidth = 1
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Observations",
    linetype = "Aedseo levels",
    color = "Aedseo levels"
  ) +
  ggplot2::scale_linetype_manual(
    values = c(
      "very low" = "dotted",
      "low" = "dashed",
      "medium" = "longdash",
      "high" = "solid"
    )
  ) +
  ggplot2::scale_color_manual(
    values = c(
      "very low" = "#44ce1b",
      "low" = "#bbdb44",
      "medium" = "#f2a134",
      "high" = "#e51f1f"
    )
  ) +
  ggplot2::facet_wrap(~ Method, ncol = 2) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    legend.position = "right",
    legend.key.width = grid::unit(2, "cm"),
    legend.text = ggplot2::element_text(size = 11, color = "black"),
    strip.text = ggplot2::element_text(size = 11, color = "black"),
    axis.text = ggplot2::element_text(size = 9, color = "black"),
    axis.title.y = ggplot2::element_text(size = 11, color = "black")
  ) +
  ggplot2::scale_y_log10(
    breaks = y_tics,
    limits = range(y_tics),
    labels = scales::label_comma(big.mark = ".", decimal.mark = ",")
  ) +
  ggplot2::guides(linetype = "none")

## ----echo = FALSE, fig.width=10, fig.height=4, dpi=300------------------------
withr::local_seed(222)
# Construct 'tsd' objects with time series data
tsd_data_noise <- generate_seasonal_data(
  years = 5,
  start_date = as.Date("2020-10-18"),
  amplitude = 600,
  mean = 600,
  phase = 0,
  noise_overdispersion = 10,
  relative_epidemic_concentration = 3,
  time_interval = "week"
)

tsd_data_noise_and_pos_trend <- generate_seasonal_data(
  years = 5,
  start_date = as.Date("2020-10-18"),
  amplitude = 600,
  mean = 600,
  phase = 0,
  noise_overdispersion = 10,
  trend_rate = 1.002,
  relative_epidemic_concentration = 3,
  time_interval = "week"
)

tsd_data_noise_and_neg_trend <- generate_seasonal_data(
  years = 5,
  start_date = as.Date("2020-10-18"),
  amplitude = 600,
  mean = 600,
  phase = 0,
  noise_overdispersion = 10,
  trend_rate = 0.99,
  relative_epidemic_concentration = 3,
  time_interval = "week"
)

# Remove days after week 20 in last season to get 5 seasons data
tsd_data_all <- rbind(
  tsd_data_noise |>
    dplyr::mutate(Data = "Noise"),
  tsd_data_noise_and_pos_trend |>
    dplyr::mutate(Data = "Noise and positive trend"),
  tsd_data_noise_and_neg_trend |>
    dplyr::mutate(Data = "Noise and negative trend")
) |>
  dplyr::filter(time <= "2025-05-12") |>
  dplyr::mutate(
    Data = factor(Data, levels = c("Noise", "Noise and positive trend", "Noise and negative trend"))
  )

start_date <- min(tsd_data_all$time)
end_date <- max(tsd_data_all$time)

tsd_data_all |>
  ggplot2::ggplot(ggplot2::aes(
    x = time,
    y = observation,
    color = Data,
    group = Data
  )) +
  ggplot2::geom_line(linewidth = 0.7) +
  ggplot2::geom_point() +
  ggplot2::scale_color_manual(
    name   = "Seasonal data",
    values = c(
      "Noise" = "blue",
      "Noise and positive trend" = "green",
      "Noise and negative trend" = "red"
    )
  ) +
  aedseo:::time_interval_x_axis(
    start_date = start_date,
    end_date   = end_date,
    time_interval_step = "6 weeks"
  ) +
  ggplot2::labs(y = "Weekly observations") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text    = ggplot2::element_text(size = 9, color = "black", family = "sans"),
    axis.title.x = ggplot2::element_text(size = 11, color = "black", family = "sans"),
    axis.title.y = ggplot2::element_text(size = 11, color = "black", family = "sans")
  )

## -----------------------------------------------------------------------------
intensity_levels_n <- seasonal_burden_levels(
  tsd = tsd_data_noise,
  disease_threshold = 10,
  method = "intensity_levels",
  conf_levels = 0.975
)
summary(intensity_levels_n)

## ----include = FALSE----------------------------------------------------------
intensity_levels_n_pos_t <- seasonal_burden_levels(
  tsd = tsd_data_noise_and_pos_trend,
  disease_threshold = 20,
  method = "intensity_levels",
  conf_levels = 0.975
)
intensity_levels_n_neg_t <- seasonal_burden_levels(
  tsd = tsd_data_noise_and_neg_trend,
  disease_threshold = 5,
  method = "intensity_levels",
  conf_levels = 0.975
)

## -----------------------------------------------------------------------------
peak_levels_n <- seasonal_burden_levels(
  tsd = tsd_data_noise,
  disease_threshold = 10,
  method = "peak_levels",
  conf_levels = c(0.4, 0.9, 0.975),
  n_peak = 8
)
summary(peak_levels_n)

## ----include = FALSE----------------------------------------------------------
peak_levels_n_pos_t <- seasonal_burden_levels(
  tsd = tsd_data_noise_and_pos_trend,
  disease_threshold = 20,
  method = "peak_levels",
  conf_levels = c(0.4, 0.9, 0.975),
  n_peak = 8
)

peak_levels_n_neg_t <- seasonal_burden_levels(
  tsd = tsd_data_noise_and_neg_trend,
  disease_threshold = 5,
  method = "peak_levels",
  conf_levels = c(0.4, 0.9, 0.975),
  n_peak = 8
)

## -----------------------------------------------------------------------------
# Remove current season such as previous seasons predict for newest season
previous_seasons <- tsd_data_all |>
  dplyr::mutate(season = epi_calendar(time)) |>
  dplyr::filter(season != "2024/2025") |>
  dplyr::select(-season)

# Run mem algorithm
mem_thresholds <- previous_seasons |>
  dplyr::group_by(Data) |>
  dplyr::group_modify(~ {
    mem_data <- .x |>
      dplyr::mutate(season = aedseo::epi_calendar(time),
                    week = lubridate::isoweek(time)) |>
      dplyr::select(-time) |>
      tidyr::pivot_wider(names_from = season, values_from = observation) |>
      dplyr::select(-week)
    # Run mem
    mem_result <- mem::memmodel(mem_data)
    # Extract thresholds
    mem_thresholds <- tibble::tibble(
      `epidemic threshold \n (mem)` = mem_result$epidemic.thresholds[1],
      `medium` = mem_result$intensity.thresholds[1],
      `high` = mem_result$intensity.thresholds[2],
      `very high` = mem_result$intensity.thresholds[3]
    )
  })

## ----echo = FALSE, fig.width=10, fig.height=8, dpi=300------------------------
#### Create data frame
burden_levels_df <- tibble::tibble(
  Level = names(
    c(intensity_levels_n$values,
      intensity_levels_n_pos_t$values,
      intensity_levels_n_neg_t$values,
      peak_levels_n$values,
      peak_levels_n_pos_t$values,
      peak_levels_n_neg_t$values
    )
  ),
  Threshold = c(
    intensity_levels_n$values,
    intensity_levels_n_pos_t$values,
    intensity_levels_n_neg_t$values,
    peak_levels_n$values,
    peak_levels_n_pos_t$values,
    peak_levels_n_neg_t$values
  ),
  Method = c(
    rep("Intensity levels", 4),
    rep("Intensity levels", 4),
    rep("Intensity levels", 4),
    rep("Peak levels", 4),
    rep("Peak levels", 4),
    rep("Peak levels", 4)
  ),
  Data = c(
    rep("Noise", 4),
    rep("Noise and positive trend", 4),
    rep("Noise and negative trend", 4),
    rep("Noise", 4),
    rep("Noise and positive trend", 4),
    rep("Noise and negative trend", 4)
  )
)

mem_levels_df <- mem_thresholds |>
  tidyr::pivot_longer(cols = `epidemic threshold \n (mem)`:`very high`,
                      names_to = "Level",
                      values_to = "Threshold") |>
  dplyr::mutate(Method = "mem levels")

# Combine the threshold data frames from the two methods
levels_all <- dplyr::bind_rows(burden_levels_df, mem_levels_df) |>
  dplyr::mutate(
    Level = factor(Level, levels = c("very high", "high", "medium", "low", "very low",
                                     "epidemic threshold \n (mem)")),
    Method = factor(Method, levels = c("Intensity levels", "Peak levels", "mem levels")),
    Data = factor(Data, levels = c("Noise", "Noise and positive trend", "Noise and negative trend"))
  )

# Merge observations
all_observations <- tsd_data_all |>
  dplyr::filter(time <= "2025-04-01") |>
  dplyr::filter(time >= "2024-05-20") |>
  dplyr::filter(!is.na(observation), observation > 3)

# Calculate y_tics
y_tics_log10 <- pretty(c(log10(levels_all$Threshold)))
y_tics_levels <- 10^(y_tics_log10)

# For each tic, find the closest magnitude to round correctly
round_to_nearest <- function(x) {
  magnitude <- 10^floor(log10(x))
  plyr::round_any(x, accuracy = magnitude)
}
y_tics <- sapply(y_tics_levels, round_to_nearest)

# Create a combined plot
all_observations |>
  ggplot2::ggplot(ggplot2::aes(x = time, y = observation)) +
  ggplot2::geom_point(ggplot2::aes(color = "observations")) +
  ggplot2::geom_hline(
    data = levels_all,
    ggplot2::aes(yintercept = Threshold, linetype = Level, color = Level),
    linewidth = 1
  ) +
  ggplot2::scale_linetype_manual(
    values = c(
      "very low" = "dotted",
      "low" = "dashed",
      "medium" = "longdash",
      "high" = "solid",
      "very high" = "dotdash",
      "epidemic threshold \n (mem)" = "dotted"
    )
  ) +
  ggplot2::scale_color_manual(
    name = "",
    values = c(
      "observations" = "black",
      "very low" = "#44ce1b",
      "low" = "#bbdb44",
      "medium" = "#f2a134",
      "high" = "#e51f1f",
      "very high" = "#891212",
      "epidemic threshold \n (mem)" = "#44ce1b"
    )
  ) +
  ggplot2::facet_grid(Data ~ Method) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    legend.position = "right",
    legend.key.width = grid::unit(2, "cm"),
    legend.text  = ggplot2::element_text(size = 11, color = "black"),
    strip.text   = ggplot2::element_text(size = 11, color = "black"),
    axis.text    = ggplot2::element_text(size = 9, color = "black"),
    axis.title.x = ggplot2::element_text(size = 11, color = "black"),
    axis.title.y = ggplot2::element_text(size = 11, color = "black")
  ) +
  ggplot2::scale_y_log10(
    breaks   = y_tics,
    limits   = range(y_tics),
    labels   = scales::label_comma(big.mark = ".", decimal.mark = ",")
  ) +
  aedseo:::time_interval_x_axis(
    start_date         = min(all_observations$time),
    end_date           = as.Date("2025-05-12"),
    time_interval_step = "5 weeks"
  ) +
  ggplot2::guides(linetype = "none")

