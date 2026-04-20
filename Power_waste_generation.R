library(tidyverse)
library(triangle)
library(data.table)
library(truncnorm)
library(purrr)

# 1. Set directories and read input datasets ------------------------------
project_dir <- "project_dir"
power_input_dir <- file.path(project_dir, "data", "Power_sector_data")
output_dir <- file.path(project_dir, "results", "power_waste")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

coal_data <- read_csv("coal_power_generation.csv")
pv_data <- read_csv("PV_waste_input.csv")
wind_data <- read_csv("wind_waste_input.csv")
waste_coeff <- read_csv("coal_waste_coefficients.csv")


# 2. Set simulation parameters --------------------------------------------
n_sim <- 10000
output_years <- c(2025, 2030, 2035, 2040, 2045, 2050)

set.seed(42)

# 3. Define functions -----------------------------------------------------
weibull_retire_cdf <- function(age, shape, scale) {
  pweibull(age, shape = shape, scale = scale)
}

ci_half_width <- function(x, lower_prob = 0.025, upper_prob = 0.975) {
  lower_bound <- unname(quantile(x, lower_prob, na.rm = TRUE))
  upper_bound <- unname(quantile(x, upper_prob, na.rm = TRUE))
  (abs(upper_bound) - abs(lower_bound)) / 2
}

# 4. Monte Carlo simulation -----------------------------------------------
results_list <- vector("list", n_sim)
pb <- txtProgressBar(min = 0, max = n_sim, style = 3)

for (sim_id in seq_len(n_sim)) {
  # 4.1 Simulate coal-power waste generation
  coal_sample <- coal_data %>%
    mutate(
      coal_power_sim = pmap_dbl(
        list(
          a = coal_power_TWh * 0.8,
          b = coal_power_TWh * 1.2,
          c = coal_power_TWh
        ),
        ~ rtriangle(n = 1, a = ..1, b = ..2, c = ..3)
      )
    )

  waste_sample <- waste_coeff %>%
    mutate(
      coefficient_sim = rtruncnorm(
        n = n(),
        a = 0,
        b = Inf,
        mean = coefficient_wt_per_TWh,
        sd = sd
      )
    )

  coal_waste <- coal_sample %>%
    crossing(waste_sample) %>%
    mutate(
      waste_mt = coal_power_sim * coefficient_sim / 100,
      simulation = sim_id
    ) %>%
    select(
      simulation,
      year,
      region,
      scenario,
      waste_type,
      waste_mt
    )

  # 4.2 Simulate PV waste generation
  pv_sample <- pv_data %>%
    mutate(
      cumulative_cap_sim = pmap_dbl(
        list(
          a = cumulative_capacity_GW * 0.8,
          b = cumulative_capacity_GW * 1.2,
          c = cumulative_capacity_GW
        ),
        ~ rtriangle(n = 1, a = ..1, b = ..2, c = ..3)
      )
    )

  eta_early_pv <- rnorm(1, mean = 20.84, sd = 20.84 * 0.10)
  beta_early_pv <- rnorm(1, mean = 2.83, sd = 2.83 * 0.10)
  eta_regular_pv <- rnorm(1, mean = 30.00, sd = 30.00 * 0.10)
  beta_regular_pv <- rnorm(1, mean = 5.38, sd = 5.38 * 0.10)

  csi_share <- rtriangle(1, a = 0.76, b = 0.95, c = 0.95)
  cdte_share <- 1 - csi_share

  area_csi <- rtriangle(1, a = 1.46 * 0.8, b = 1.46 * 1.2, c = 1.46)
  mass_csi <- rtriangle(1, a = 23.0 * 0.8, b = 23.0 * 1.2, c = 23.0)
  power_csi <- rtriangle(
    1,
    a = (224 / 1.46) * 0.8,
    b = (224 / 1.46) * 1.2,
    c = (224 / 1.46)
  )

  area_cdte <- rtriangle(1, a = 0.72 * 0.8, b = 0.72 * 1.2, c = 0.72)
  mass_cdte <- rtriangle(1, a = 12.0 * 0.8, b = 12.0 * 1.2, c = 12.0)
  power_cdte <- rtriangle(
    1,
    a = (65 / 0.72) * 0.8,
    b = (65 / 0.72) * 1.2,
    c = (65 / 0.72)
  )

  panel_info <- tibble(
    tech = c("PV_CSi", "PV_CdTe"),
    market_share = c(csi_share, cdte_share),
    area_m2 = c(area_csi, area_cdte),
    mass_kg_m2 = c(mass_csi, mass_cdte),
    power_w_m2 = c(power_csi, power_cdte)
  ) %>%
    mutate(
      unit_mass_t_per_gw = market_share * (mass_kg_m2 / area_m2) / power_w_m2 * 1e6
    )

  pv_retire_detail <- pv_sample %>%
    arrange(region, scenario, year) %>%
    group_by(region, scenario) %>%
    mutate(
      previous_cap = lag(cumulative_cap_sim, default = 0),
      added_installation = cumulative_cap_sim - previous_cap,
      installation_year = year
    ) %>%
    ungroup() %>%
    select(region, scenario, installation_year, added_installation) %>%
    crossing(retirement_year = 2000:2050) %>%
    filter(retirement_year > installation_year) %>%
    mutate(
      age = retirement_year - installation_year,
      scale_use = if_else(installation_year <= 2013, eta_early_pv, eta_regular_pv),
      shape_use = if_else(installation_year <= 2013, beta_early_pv, beta_regular_pv),
      cdf_retirement = weibull_retire_cdf(age, shape = shape_use, scale = scale_use),
      cdf_previous = if_else(
        age - 1 >= 0,
        weibull_retire_cdf(age - 1, shape = shape_use, scale = scale_use),
        0
      ),
      retire_fraction = pmax(0, pmin(1, cdf_retirement - cdf_previous)),
      retired_gw = added_installation * retire_fraction
    ) %>%
    group_by(region, scenario, year = retirement_year) %>%
    summarise(
      retired_gw_total = sum(retired_gw, na.rm = TRUE),
      .groups = "drop"
    )

  pv_waste <- pv_retire_detail %>%
    crossing(panel_info) %>%
    mutate(
      waste_mt = retired_gw_total * unit_mass_t_per_gw / 1e6,
      simulation = sim_id
    ) %>%
    select(
      simulation,
      year,
      region,
      scenario,
      waste_type = tech,
      waste_mt
    ) %>%
    filter(year %in% output_years)

  # 4.3 Simulate wind-power waste generation
  wind_sample <- wind_data %>%
    mutate(
      cumulative_cap_sim = pmap_dbl(
        list(
          a = cumulative_capacity_GW * 0.8,
          b = cumulative_capacity_GW * 1.2,
          c = cumulative_capacity_GW
        ),
        ~ rtriangle(n = 1, a = ..1, b = ..2, c = ..3)
      )
    )

  scale_early_wind <- rnorm(1, mean = 15.2, sd = 15.2 * 0.10)
  shape_early_wind <- rnorm(1, mean = 2.2, sd = 2.2 * 0.10)
  scale_regular_wind <- rnorm(1, mean = 22.8, sd = 22.8 * 0.10)
  shape_regular_wind <- rnorm(1, mean = 4.1, sd = 4.1 * 0.10)

  gfrp_share <- rtriangle(1, a = 0.80 * 0.8, b = 0.80 * 1.2, c = 0.80)
  cfrp_share <- 1 - gfrp_share
  mass_gfrp <- rtriangle(1, a = 12350 * 0.8, b = 12350 * 1.2, c = 12350)
  mass_cfrp <- rtriangle(1, a = 8645 * 0.8, b = 8645 * 1.2, c = 8645)

  wind_info <- tibble(
    tech = c("Wind_GFRP", "Wind_CFRP"),
    market_share = c(gfrp_share, cfrp_share),
    mass_t_per_gw = c(mass_gfrp, mass_cfrp)
  ) %>%
    mutate(
      unit_mass_t_per_gw = market_share * mass_t_per_gw
    )

  wind_retire_detail <- wind_sample %>%
    arrange(region, scenario, year) %>%
    group_by(region, scenario) %>%
    mutate(
      previous_cap = lag(cumulative_cap_sim, default = 0),
      added_installation = cumulative_cap_sim - previous_cap,
      installation_year = year
    ) %>%
    ungroup() %>%
    select(region, scenario, installation_year, added_installation) %>%
    crossing(retirement_year = 2000:2050) %>%
    filter(retirement_year > installation_year) %>%
    mutate(
      age = retirement_year - installation_year,
      scale_use = if_else(
        installation_year <= 2013,
        scale_early_wind,
        scale_regular_wind
      ),
      shape_use = if_else(
        installation_year <= 2013,
        shape_early_wind,
        shape_regular_wind
      ),
      cdf_retirement = weibull_retire_cdf(age, shape = shape_use, scale = scale_use),
      cdf_previous = if_else(
        age - 1 >= 0,
        weibull_retire_cdf(age - 1, shape = shape_use, scale = scale_use),
        0
      ),
      retire_fraction = pmax(0, pmin(1, cdf_retirement - cdf_previous)),
      retired_gw = added_installation * retire_fraction
    ) %>%
    group_by(region, scenario, year = retirement_year) %>%
    summarise(
      retired_gw_total = sum(retired_gw, na.rm = TRUE),
      .groups = "drop"
    )

  wind_waste <- wind_retire_detail %>%
    crossing(wind_info) %>%
    mutate(
      waste_mt = retired_gw_total * unit_mass_t_per_gw / 1e6,
      simulation = sim_id
    ) %>%
    select(
      simulation,
      year,
      region,
      scenario,
      waste_type = tech,
      waste_mt
    ) %>%
    filter(year %in% output_years)

  # 4.4 Combine all waste types
  all_waste_sim <- bind_rows(
    coal_waste,
    pv_waste,
    wind_waste
  ) %>%
    group_by(simulation, year, region, scenario, waste_type) %>%
    summarise(
      total_waste_mt = sum(waste_mt, na.rm = TRUE),
      .groups = "drop"
    )

  results_list[[sim_id]] <- all_waste_sim
  setTxtProgressBar(pb, sim_id)
}

close(pb)

# 5. Combine all simulation outputs ---------------------------------------
all_results <- rbindlist(results_list)

save(
  all_results,
  file = file.path(output_dir, "power_waste_generation.Rdata")
)