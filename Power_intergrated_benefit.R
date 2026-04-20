library(tidyverse)
library(data.table)
library(purrr)
library(triangle)
library(truncnorm)

# 1. Set directories ------------------------------------------------------
project_dir <- "projec_dir"
input_dir <- file.path(project_dir, "data", "Power_sector_data")
output_dir <- file.path(project_dir, "output", "Power_integrated_benefit_results")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Input datasets -------------------------------------------------------
# Power generation data
coal_data <- read_csv("coal_power_generation.csv")
pv_data <- read_csv("PV_waste_input.csv")
wind_data <- read_csv("wind_waste_input.csv")

# Waste-generation and benefit parameter data
waste_coeff <- read_csv("coal_waste_coefficients.csv")
power_tech_percent <- read_csv("Power_tech_percent.csv")
tech_benefit <- read_csv("tech_benefit.csv")

carbon_price_df <- read_csv("carbon_price.csv") %>%
  rename(
    development_scenario = scenario,
    CO2_price = carbon_price_2020
  ) %>%
  select(year, development_scenario, CO2_price)

# 3. Set simulation parameters --------------------------------------------
n_sim <- 10000
disc_rate <- 0.03
base_year <- 2020
ci_lower_prob <- ci_lower_prob
ci_upper_prob <- ci_upper_prob

set.seed(42)

# 4. Define functions -----------------------------------------------------
weibull_retire_cdf <- function(age, shape, scale) {
  pweibull(age, shape = shape, scale = scale)
}

generate_triangle <- function(x, variation = 0.10) {
  if (is.na(x)) {
    return(NA_real_)
  }

  min_val <- min(x * (1 - variation), x * (1 + variation))
  max_val <- max(x * (1 - variation), x * (1 + variation))
  mode_val <- x

  rtriangle(1, a = min_val, b = max_val, c = mode_val)
}

ci_half_width <- function(lower, upper) {
  (abs(upper) - abs(lower)) / 2
}

# 5. Monte Carlo simulation -----------------------------------------------
results_list <- vector("list", n_sim)
lca_result <- vector("list", n_sim)

pb <- txtProgressBar(min = 0, max = n_sim, style = 3)

for (sim_id in seq_len(n_sim)) {

  # 5.1 Simulate coal-power waste generation
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
      waste_Mt = coal_power_sim * coefficient_sim / 100,
      simulation = sim_id
    ) %>%
    select(
      year, region, scenario,
      waste_type,
      waste_Mt,
      simulation
    )

  # 5.2 Simulate PV waste generation
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
  eta_reg_pv <- rnorm(1, mean = 30.00, sd = 30.00 * 0.10)
  beta_reg_pv <- rnorm(1, mean = 5.38, sd = 5.38 * 0.10)

  csi_share <- rtriangle(1, a = 0.76, b = 0.95, c = 0.95)
  cdte_share <- 1 - csi_share

  area_csi <- rtriangle(1, a = 1.46 * 0.8, b = 1.46 * 1.2, c = 1.46)
  mass_csi <- rtriangle(1, a = 23.0 * 0.8, b = 23.0 * 1.2, c = 23.0)
  pow_csi <- rtriangle(
    1,
    a = (224 / 1.46) * 0.8,
    b = (224 / 1.46) * 1.2,
    c = (224 / 1.46)
  )

  area_cdte <- rtriangle(1, a = 0.72 * 0.8, b = 0.72 * 1.2, c = 0.72)
  mass_cdte <- rtriangle(1, a = 12.0 * 0.8, b = 12.0 * 1.2, c = 12.0)
  pow_cdte <- rtriangle(
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
    power_W_m2 = c(pow_csi, pow_cdte)
  ) %>%
    mutate(
      unit_mass_t_per_GW = market_share * (mass_kg_m2 / area_m2) / power_W_m2 * 1e6
    )

  pv_retire_detail <- pv_sample %>%
    arrange(region, scenario, year) %>%
    group_by(region, scenario) %>%
    mutate(
      cap_prev = lag(cumulative_cap_sim, default = 0),
      add_install = cumulative_cap_sim - cap_prev,
      year_inst = year
    ) %>%
    ungroup() %>%
    select(region, scenario, year_inst, add_install) %>%
    crossing(year_ret = 2000:2050) %>%
    filter(year_ret > year_inst) %>%
    mutate(
      age = year_ret - year_inst,
      alpha_use = if_else(year_inst <= 2013, eta_early_pv, eta_reg_pv),
      beta_use = if_else(year_inst <= 2013, beta_early_pv, beta_reg_pv),
      cdf_ret = weibull_retire_cdf(age, shape = beta_use, scale = alpha_use),
      cdf_prev = if_else(
        age - 1 >= 0,
        weibull_retire_cdf(age - 1, shape = beta_use, scale = alpha_use),
        0
      ),
      retire_frac = pmax(0, pmin(1, cdf_ret - cdf_prev)),
      retired_GW = add_install * retire_frac
    ) %>%
    group_by(region, scenario, year = year_ret) %>%
    summarise(
      retired_GW_total = sum(retired_GW, na.rm = TRUE),
      .groups = "drop"
    )

  pv_waste <- pv_retire_detail %>%
    crossing(panel_info) %>%
    mutate(
      waste_Mt = retired_GW_total * unit_mass_t_per_GW / 1e6,
      simulation = sim_id
    ) %>%
    select(
      simulation,
      year, region, scenario,
      waste_type = tech,
      waste_Mt
    ) %>%
    filter(year %in% c(2025, 2030, 2035, 2040, 2045, 2050))

  # 5.3 Simulate wind-power waste generation
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

  alpha_early_w <- rnorm(1, mean = 15.2, sd = 15.2 * 0.10)
  beta_early_w <- rnorm(1, mean = 2.2, sd = 2.2 * 0.10)
  alpha_reg_w <- rnorm(1, mean = 22.8, sd = 22.8 * 0.10)
  beta_reg_w <- rnorm(1, mean = 4.1, sd = 4.1 * 0.10)

  gfrp_share <- rtriangle(1, a = 0.80 * 0.8, b = 0.80 * 1.2, c = 0.80)
  cfrp_share <- 1 - gfrp_share
  mass_gfrp <- rtriangle(1, a = 12350 * 0.8, b = 12350 * 1.2, c = 12350)
  mass_cfrp <- rtriangle(1, a = 8645 * 0.8, b = 8645 * 1.2, c = 8645)

  wind_info <- tibble(
    tech = c("Wind_GFRP", "Wind_CFRP"),
    market_share = c(gfrp_share, cfrp_share),
    mass_t_per_GW = c(mass_gfrp, mass_cfrp)
  ) %>%
    mutate(
      unit_mass_t_per_GW = market_share * mass_t_per_GW
    )

  wind_retire_detail <- wind_sample %>%
    arrange(region, scenario, year) %>%
    group_by(region, scenario) %>%
    mutate(
      prev_cap = lag(cumulative_cap_sim, default = 0),
      add_install = cumulative_cap_sim - prev_cap,
      year_inst = year
    ) %>%
    ungroup() %>%
    select(region, scenario, year_inst, add_install) %>%
    crossing(year_ret = 2000:2050) %>%
    filter(year_ret > year_inst) %>%
    mutate(
      age = year_ret - year_inst,
      alpha_use = if_else(year_inst <= 2013, alpha_early_w, alpha_reg_w),
      beta_use = if_else(year_inst <= 2013, beta_early_w, beta_reg_w),
      cdf_rt = weibull_retire_cdf(age, shape = beta_use, scale = alpha_use),
      cdf_prev = if_else(
        age - 1 >= 0,
        weibull_retire_cdf(age - 1, shape = beta_use, scale = alpha_use),
        0
      ),
      retire_frac = pmax(0, pmin(1, cdf_rt - cdf_prev)),
      retired_GW = add_install * retire_frac
    ) %>%
    group_by(region, scenario, year = year_ret) %>%
    summarise(
      retired_GW_total = sum(retired_GW, na.rm = TRUE),
      .groups = "drop"
    )

  wind_waste <- wind_retire_detail %>%
    crossing(wind_info) %>%
    mutate(
      waste_Mt = retired_GW_total * unit_mass_t_per_GW / 1e6,
      simulation = sim_id
    ) %>%
    select(
      simulation,
      year, region, scenario,
      waste_type = tech,
      waste_Mt
    ) %>%
    filter(year %in% c(2025, 2030, 2035, 2040, 2045, 2050))

  # 5.4 Combine all waste streams
  all_waste_i <- bind_rows(
    coal_waste,
    pv_waste,
    wind_waste
  ) %>%
    group_by(simulation, year, region, scenario, waste_type) %>%
    summarise(
      total_waste_Mt = sum(waste_Mt, na.rm = TRUE),
      .groups = "drop"
    )

  # 5.5 Allocate waste to treatment technologies and calculate benefits
  waste_with_tech <- all_waste_i %>%
    left_join(
      power_tech_percent,
      by = c("year", "region", "waste_type"),
      relationship = "many-to-many"
    ) %>%
    mutate(
      development_scenario = scenario,
      combined_scenario = paste0(development_scenario, "-", tech_scenario),
      waste_to_tech_Mt = total_waste_Mt * percent
    ) %>%
    select(
      simulation, year, region,
      development_scenario, tech_scenario, combined_scenario,
      waste_type, tech, waste_to_tech_Mt
    )

  tech_param_sim <- tech_benefit %>%
    mutate(
      EB_1_sim = map_dbl(EB_1, generate_triangle),
      EB_2_sim = map_dbl(EB_2_cost, generate_triangle),
      EB_3_sim = map_dbl(EB_3_cost, generate_triangle),
      EB_4_sim = map_dbl(EB_4_cost, generate_triangle),
      EB_5_sim = map_dbl(EB_5_cost, generate_triangle),
      EB_6_sim = map_dbl(EB_6_cost, generate_triangle),
      EB_7_sim = map_dbl(EB_7_cost, generate_triangle),
      EB_8_sim = map_dbl(EB_8_cost, generate_triangle),
      EB_9_sim = map_dbl(EB_9_cost, generate_triangle),
      EB_10_sim = map_dbl(EB_10_cost, generate_triangle),
      EB_11_sim = map_dbl(EB_11_cost, generate_triangle),
      EB_12_sim = map_dbl(EB_12_cost, generate_triangle),
      EB_13_sim = map_dbl(EB_13_cost, generate_triangle),
      EB_14_sim = map_dbl(EB_14_cost, generate_triangle),
      EB_15_sim = map_dbl(EB_15_cost, generate_triangle),
      EB_16_sim = map_dbl(EB_16_cost, generate_triangle),
      EB_17_sim = map_dbl(EB_17_cost, generate_triangle),
      EB_18_sim = map_dbl(EB_18_cost, generate_triangle),
      cost_sim = pmap_dbl(
        list(a = Cost * 0.9, b = Cost * 1.1, c = Cost),
        ~ rtriangle(1, ..1, ..2, ..3)
      ),
      PB_sim = pmap_dbl(
        list(a = PB * 0.9, b = PB * 1.1, c = PB),
        ~ rtriangle(1, ..1, ..2, ..3)
      )
    ) %>%
    select(
      tech,
      EB_1_sim, EB_2_sim, EB_3_sim, EB_4_sim, EB_5_sim,
      EB_6_sim, EB_7_sim, EB_8_sim, EB_9_sim, EB_10_sim,
      EB_11_sim, EB_12_sim, EB_13_sim, EB_14_sim, EB_15_sim,
      EB_16_sim, EB_17_sim, EB_18_sim,
      cost_sim, PB_sim
    )

  carbon_price_df_sim <- carbon_price_df %>%
    mutate(
      carbon_price_sim = pmap_dbl(
        list(a = CO2_price * 0.9, b = CO2_price * 1.1, c = CO2_price),
        ~ rtriangle(1, ..1, ..2, ..3)
      )
    )

  revenue_calc <- waste_with_tech %>%
    left_join(
      tech_param_sim,
      by = "tech"
    ) %>%
    filter(year >= 2025) %>%
    left_join(
      carbon_price_df_sim,
      by = c("year", "development_scenario")
    ) %>%
    mutate(
      discount_factor = 1 / ((1 + disc_rate) ^ (year - base_year)),
      cost_total = waste_to_tech_Mt * cost_sim / 1e3,
      cost_total_disc = cost_total * discount_factor,
      PB_total = waste_to_tech_Mt * PB_sim / 1e3,
      PB_total_disc = PB_total * discount_factor,
      CO2_decreasing_Gt = -EB_1_sim * waste_to_tech_Mt / 1e6,
      CO2_benefit_raw = -EB_1_sim * waste_to_tech_Mt * carbon_price_sim / 1e6,
      CO2_benefit_disc = CO2_benefit_raw * discount_factor,
      EB_2_benefit_disc = -EB_2_sim * waste_to_tech_Mt * discount_factor / 1e3,
      EB_3_benefit_disc = -EB_3_sim * waste_to_tech_Mt * discount_factor / 1e3,
      EB_4_benefit_disc = -EB_4_sim * waste_to_tech_Mt * discount_factor / 1e3,
      EB_5_benefit_disc = -EB_5_sim * waste_to_tech_Mt * discount_factor / 1e3,
      EB_6_benefit_disc = -EB_6_sim * waste_to_tech_Mt * discount_factor / 1e3,
      EB_7_benefit_disc = -EB_7_sim * waste_to_tech_Mt * discount_factor / 1e3,
      EB_8_benefit_disc = -EB_8_sim * waste_to_tech_Mt * discount_factor / 1e3,
      EB_9_benefit_disc = -EB_9_sim * waste_to_tech_Mt * discount_factor / 1e3,
      EB_10_benefit_disc = -EB_10_sim * waste_to_tech_Mt * discount_factor / 1e3,
      EB_11_benefit_disc = -EB_11_sim * waste_to_tech_Mt * discount_factor / 1e3,
      EB_12_benefit_disc = -EB_12_sim * waste_to_tech_Mt * discount_factor / 1e3,
      EB_13_benefit_disc = -EB_13_sim * waste_to_tech_Mt * discount_factor / 1e3,
      EB_14_benefit_disc = -EB_14_sim * waste_to_tech_Mt * discount_factor / 1e3,
      EB_15_benefit_disc = -EB_15_sim * waste_to_tech_Mt * discount_factor / 1e3,
      EB_16_benefit_disc = -EB_16_sim * waste_to_tech_Mt * discount_factor / 1e3,
      EB_17_benefit_disc = -EB_17_sim * waste_to_tech_Mt * discount_factor / 1e3,
      EB_18_benefit_disc = -EB_18_sim * waste_to_tech_Mt * discount_factor / 1e3
    )

  revenue_summary <- revenue_calc %>%
    group_by(
      simulation, year, region,
      development_scenario, tech_scenario, combined_scenario
    ) %>%
    summarise(
      total_waste_mt = sum(waste_to_tech_Mt, na.rm = TRUE),
      total_cost_all = sum(cost_total_disc, na.rm = TRUE),
      total_PB_all = sum(PB_total_disc, na.rm = TRUE),
      total_CO2_decreasing_Gt = sum(CO2_decreasing_Gt, na.rm = TRUE),
      GWP_benefit = sum(CO2_benefit_disc, na.rm = TRUE),
      total_EB_2_benefit = sum(EB_2_benefit_disc, na.rm = TRUE),
      total_EB_3_benefit = sum(EB_3_benefit_disc, na.rm = TRUE),
      total_EB_4_benefit = sum(EB_4_benefit_disc, na.rm = TRUE),
      total_EB_5_benefit = sum(EB_5_benefit_disc, na.rm = TRUE),
      total_EB_6_benefit = sum(EB_6_benefit_disc, na.rm = TRUE),
      total_EB_7_benefit = sum(EB_7_benefit_disc, na.rm = TRUE),
      total_EB_8_benefit = sum(EB_8_benefit_disc, na.rm = TRUE),
      total_EB_9_benefit = sum(EB_9_benefit_disc, na.rm = TRUE),
      total_EB_10_benefit = sum(EB_10_benefit_disc, na.rm = TRUE),
      total_EB_11_benefit = sum(EB_11_benefit_disc, na.rm = TRUE),
      total_EB_12_benefit = sum(EB_12_benefit_disc, na.rm = TRUE),
      total_EB_13_benefit = sum(EB_13_benefit_disc, na.rm = TRUE),
      total_EB_14_benefit = sum(EB_14_benefit_disc, na.rm = TRUE),
      total_EB_15_benefit = sum(EB_15_benefit_disc, na.rm = TRUE),
      total_EB_16_benefit = sum(EB_16_benefit_disc, na.rm = TRUE),
      total_EB_17_benefit = sum(EB_17_benefit_disc, na.rm = TRUE),
      total_EB_18_benefit = sum(EB_18_benefit_disc, na.rm = TRUE),
      Ecosystem_benefit =
        total_EB_2_benefit + total_EB_3_benefit + total_EB_4_benefit +
        total_EB_5_benefit + total_EB_6_benefit + total_EB_7_benefit +
        total_EB_8_benefit + total_EB_9_benefit + total_EB_15_benefit,
      Toxicity_benefit =
        total_EB_10_benefit + total_EB_12_benefit + total_EB_13_benefit +
        total_EB_14_benefit,
      Resources_benefit =
        total_EB_16_benefit + total_EB_17_benefit + total_EB_18_benefit,
      .groups = "drop"
    ) %>%
    filter(region != "world")

  # 5.6 Calculate LCA indicator results
  tech_lca_sim <- tech_benefit %>%
    mutate(
      lca_1_sim = map_dbl(EB_1, generate_triangle),
      lca_2_sim = map_dbl(EB_2, generate_triangle),
      lca_3_sim = map_dbl(EB_3, generate_triangle),
      lca_4_sim = map_dbl(EB_4, generate_triangle),
      lca_5_sim = map_dbl(EB_5, generate_triangle),
      lca_6_sim = map_dbl(EB_6, generate_triangle),
      lca_7_sim = map_dbl(EB_7, generate_triangle),
      lca_8_sim = map_dbl(EB_8, generate_triangle),
      lca_9_sim = map_dbl(EB_9, generate_triangle),
      lca_10_sim = map_dbl(EB_10, generate_triangle),
      lca_11_sim = map_dbl(EB_11, generate_triangle),
      lca_12_sim = map_dbl(EB_12, generate_triangle),
      lca_13_sim = map_dbl(EB_13, generate_triangle),
      lca_14_sim = map_dbl(EB_14, generate_triangle),
      lca_15_sim = map_dbl(EB_15, generate_triangle),
      lca_16_sim = map_dbl(EB_16, generate_triangle),
      lca_17_sim = map_dbl(EB_17, generate_triangle),
      lca_18_sim = map_dbl(EB_18, generate_triangle)
    )

  lca_calc <- waste_with_tech %>%
    left_join(
      tech_lca_sim,
      by = "tech"
    ) %>%
    filter(year >= 2025) %>%
    mutate(
      EB_1_benefit_lca = -lca_1_sim * waste_to_tech_Mt / 1e6,
      EB_2_benefit_lca = -lca_2_sim * waste_to_tech_Mt,
      EB_3_benefit_lca = -lca_3_sim * waste_to_tech_Mt,
      EB_4_benefit_lca = -lca_4_sim * waste_to_tech_Mt,
      EB_5_benefit_lca = -lca_5_sim * waste_to_tech_Mt,
      EB_6_benefit_lca = -lca_6_sim * waste_to_tech_Mt,
      EB_7_benefit_lca = -lca_7_sim * waste_to_tech_Mt,
      EB_8_benefit_lca = -lca_8_sim * waste_to_tech_Mt,
      EB_9_benefit_lca = -lca_9_sim * waste_to_tech_Mt,
      EB_10_benefit_lca = -lca_10_sim * waste_to_tech_Mt,
      EB_11_benefit_lca = -lca_11_sim * waste_to_tech_Mt,
      EB_12_benefit_lca = -lca_12_sim * waste_to_tech_Mt,
      EB_13_benefit_lca = -lca_13_sim * waste_to_tech_Mt,
      EB_14_benefit_lca = -lca_14_sim * waste_to_tech_Mt,
      EB_15_benefit_lca = -lca_15_sim * waste_to_tech_Mt,
      EB_16_benefit_lca = -lca_16_sim * waste_to_tech_Mt,
      EB_17_benefit_lca = -lca_17_sim * waste_to_tech_Mt,
      EB_18_benefit_lca = -lca_18_sim * waste_to_tech_Mt
    )

  lca_summary <- lca_calc %>%
    group_by(
      simulation, year, region,
      development_scenario, tech_scenario, combined_scenario
    ) %>%
    summarise(
      total_waste_mt = sum(waste_to_tech_Mt, na.rm = TRUE),
      total_EB_1_lca = sum(EB_1_benefit_lca, na.rm = TRUE),
      total_EB_2_lca = sum(EB_2_benefit_lca, na.rm = TRUE),
      total_EB_3_lca = sum(EB_3_benefit_lca, na.rm = TRUE),
      total_EB_4_lca = sum(EB_4_benefit_lca, na.rm = TRUE),
      total_EB_5_lca = sum(EB_5_benefit_lca, na.rm = TRUE),
      total_EB_6_lca = sum(EB_6_benefit_lca, na.rm = TRUE),
      total_EB_7_lca = sum(EB_7_benefit_lca, na.rm = TRUE),
      total_EB_8_lca = sum(EB_8_benefit_lca, na.rm = TRUE),
      total_EB_9_lca = sum(EB_9_benefit_lca, na.rm = TRUE),
      total_EB_10_lca = sum(EB_10_benefit_lca, na.rm = TRUE),
      total_EB_11_lca = sum(EB_11_benefit_lca, na.rm = TRUE),
      total_EB_12_lca = sum(EB_12_benefit_lca, na.rm = TRUE),
      total_EB_13_lca = sum(EB_13_benefit_lca, na.rm = TRUE),
      total_EB_14_lca = sum(EB_14_benefit_lca, na.rm = TRUE),
      total_EB_15_lca = sum(EB_15_benefit_lca, na.rm = TRUE),
      total_EB_16_lca = sum(EB_16_benefit_lca, na.rm = TRUE),
      total_EB_17_lca = sum(EB_17_benefit_lca, na.rm = TRUE),
      total_EB_18_lca = sum(EB_18_benefit_lca, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(region != "world")

  results_list[[sim_id]] <- revenue_summary
  lca_result[[sim_id]] <- lca_summary

  setTxtProgressBar(pb, sim_id)
}

close(pb)

# 6. Combine all simulation outputs and summarize uncertainty -------------
all_results <- rbindlist(results_list) %>%
  filter(region != "world", year >= 2025)

save(
  all_results,
  file = file.path(output_dir, "Power_all_benefits.RData")
)

save(
  lca_result,
  file = file.path(output_dir, "Power_lca_results.RData")
)

# 7. Summarize Monte Carlo uncertainty using 95% CI -----------------------
uncertainty_summary <- all_results %>%
  group_by(year, region, development_scenario, tech_scenario, combined_scenario) %>%
  summarise(
    mean_waste_Mt = mean(total_waste_mt),
    lower_95_waste_Mt = quantile(total_waste_mt, ci_lower_prob),
    upper_95_waste_Mt = quantile(total_waste_mt, ci_upper_prob),
    waste_Mt_half_width = ci_half_width(lower_95_waste_Mt, upper_95_waste_Mt),

    mean_CO2_decreasing_Gt = mean(total_CO2_decreasing_Gt),
    lower_95_CO2_decreasing_Gt = quantile(total_CO2_decreasing_Gt, ci_lower_prob),
    upper_95_CO2_decreasing_Gt = quantile(total_CO2_decreasing_Gt, ci_upper_prob),
    CO2_decreasing_Gt_half_width = ci_half_width(
      lower_95_CO2_decreasing_Gt,
      upper_95_CO2_decreasing_Gt
    ),

    mean_cost_billion_usd = mean(total_cost_all),
    lower_95_cost_billion_usd = quantile(total_cost_all, ci_lower_prob),
    upper_95_cost_billion_usd = quantile(total_cost_all, ci_upper_prob),
    cost_billion_usd_half_width = ci_half_width(
      lower_95_cost_billion_usd,
      upper_95_cost_billion_usd
    ),

    mean_product_benefit_billion_usd = mean(total_PB_all),
    lower_95_product_benefit_billion_usd = quantile(total_PB_all, ci_lower_prob),
    upper_95_product_benefit_billion_usd = quantile(total_PB_all, ci_upper_prob),
    product_benefit_billion_usd_half_width = ci_half_width(
      lower_95_product_benefit_billion_usd,
      upper_95_product_benefit_billion_usd
    ),

    mean_gwp_benefit_billion_usd = mean(GWP_benefit),
    lower_95_gwp_benefit_billion_usd = quantile(GWP_benefit, ci_lower_prob),
    upper_95_gwp_benefit_billion_usd = quantile(GWP_benefit, ci_upper_prob),
    gwp_benefit_billion_usd_half_width = ci_half_width(
      lower_95_gwp_benefit_billion_usd,
      upper_95_gwp_benefit_billion_usd
    ),

    mean_ecosystem_benefit_billion_usd = mean(Ecosystem_benefit),
    lower_95_ecosystem_benefit_billion_usd = quantile(Ecosystem_benefit, ci_lower_prob),
    upper_95_ecosystem_benefit_billion_usd = quantile(Ecosystem_benefit, ci_upper_prob),
    ecosystem_benefit_billion_usd_half_width = ci_half_width(
      lower_95_ecosystem_benefit_billion_usd,
      upper_95_ecosystem_benefit_billion_usd
    ),

    mean_toxicity_benefit_billion_usd = mean(Toxicity_benefit),
    lower_95_toxicity_benefit_billion_usd = quantile(Toxicity_benefit, ci_lower_prob),
    upper_95_toxicity_benefit_billion_usd = quantile(Toxicity_benefit, ci_upper_prob),
    toxicity_benefit_billion_usd_half_width = ci_half_width(
      lower_95_toxicity_benefit_billion_usd,
      upper_95_toxicity_benefit_billion_usd
    ),

    mean_resources_benefit_billion_usd = mean(Resources_benefit),
    lower_95_resources_benefit_billion_usd = quantile(Resources_benefit, ci_lower_prob),
    upper_95_resources_benefit_billion_usd = quantile(Resources_benefit, ci_upper_prob),
    resources_benefit_billion_usd_half_width = ci_half_width(
      lower_95_resources_benefit_billion_usd,
      upper_95_resources_benefit_billion_usd
    ),
    .groups = "drop"
  )

write.csv(
  uncertainty_summary,
  file.path(output_dir, "Power_summary_cost_benefit_all_benefits.csv"),
  row.names = FALSE
)
