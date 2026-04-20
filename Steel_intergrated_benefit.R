library(tidyverse)
library(triangle)
library(data.table)
library(progress)
library(truncnorm)


# 1. Set directories and global parameters --------------------------------
project_dir <- "projec_dir"
input_dir <- file.path(project_dir, "data", "Steel_sector_data")
output_dir <- file.path(project_dir, "output", "Steel_integrated_benefit_results")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setwd(project_dir)

disc_rate <- 0.03
base_year <- 2020
n_sim <- 10000
ci_lower_prob <- 0.025
ci_upper_prob <- 0.975

set.seed(42)

# 2. Input datasets -------------------------------------------------------
steel_production <- read_csv("steel_production.csv")
steel_route_share <- read_csv("steel_route_share.csv")
waste_factors <- read_csv("steel_waste_factors.csv")

steel_tech_percent <- read_csv("steel_tech_percent.csv")
tech_benefit <- read_csv("tech_benefit.csv")

carbon_price_df <- read_csv("carbon_price.csv") %>%
  rename(
    development_scenario = scenario,
    CO2_price = carbon_price_2020
  ) %>%
  select(year, development_scenario, CO2_price)


# 3. Define functions -----------------------------------------------------
simulate_triangle <- function(x, variation = 0.10) {
  if (is.na(x)) {
    return(NA_real_)
  }

  min_val <- min(x * (1 - variation), x * (1 + variation))
  max_val <- max(x * (1 - variation), x * (1 + variation))
  mode_val <- x

  rtriangle(1, a = min_val, b = max_val, c = mode_val)
}

# Use a truncated normal distribution to avoid negative waste factors.
simulate_truncnorm_positive <- function(mu, sigma, lower_bound = 1e-9) {
  if (is.na(sigma) || sigma <= 0) {
    return(max(mu, lower_bound))
  }

  rtruncnorm(1, a = lower_bound, b = Inf, mean = mu, sd = sigma)
}

calc_ci_half_width <- function(lower, upper) {
  (abs(upper) - abs(lower)) / 2
}

# 4. Initialize storage objects -------------------------------------------
results_list <- vector("list", n_sim)
lca_results_list <- vector("list", n_sim)

progress_bar <- txtProgressBar(min = 0, max = n_sim, style = 3)  # style = 3 shows percentage and progress bar

# 5. Monte Carlo simulation -----------------------------------------------
for (i in seq_len(n_sim)) {

  # 5.1 Simulate steel production
  steel_sample <- steel_production %>%
    mutate(
      steel_production_sim = pmap_dbl(
        list(
          steel_production * 0.8,
          steel_production * 1.2,
          steel_production
        ),
        ~ rtriangle(1, ..1, ..2, ..3)
      )
    )

  # 5.2 Simulate BOF/EAF shares
  route_sample <- steel_route_share %>%
    mutate(
      BOF_share_sim = pmap_dbl(
        list(
          a = BOF_share * 0.8,
          b = BOF_share * 1.2,
          c = BOF_share
        ),
        ~ pmin(pmax(rtriangle(1, ..1, ..2, ..3), 0), 1)
      ),
      EAF_share_sim = 1 - BOF_share_sim
    )

  # 5.3 Merge production and route shares
  steel_data_sim <- steel_sample %>%
    select(year, region, scenario, steel_production_sim) %>%
    left_join(
      route_sample %>%
        select(year, region, scenario, BOF_share_sim, EAF_share_sim),
      by = c("year", "region", "scenario")
    ) %>%
    mutate(
      BOF_production = steel_production_sim * BOF_share_sim,
      EAF_production = steel_production_sim * EAF_share_sim
    ) %>%
    select(year, region, scenario, BOF_production, EAF_production)

  # 5.4 Convert wide production data to long format
  steel_long_sim <- steel_data_sim %>%
    pivot_longer(
      cols = c(BOF_production, EAF_production),
      names_to = "route",
      values_to = "production"
    ) %>%
    mutate(
      route = if_else(route == "BOF_production", "BOF", "EAF")
    )

  # 5.5 Simulate waste factors
  waste_factor_sample <- waste_factors %>%
    mutate(
      grp = if_else(
        waste_type %in% c("BF_Slag", "Fly_Ash_Steel"),
        waste_type,
        paste(waste_type, route, sep = "|")
      )
    ) %>%
    group_by(grp) %>%
    mutate(
      waste_factor_sim = simulate_truncnorm_positive(
        mu = first(waste_factor),
        sigma = first(sd)
      )
    ) %>%
    ungroup() %>%
    select(waste_type, route, waste_factor_sim)

  # 5.6 Calculate waste output
  waste_output <- steel_long_sim %>%
    left_join(
      waste_factor_sample,
      by = "route",
      relationship = "many-to-many"
    ) %>%
    mutate(
      waste_amount_Mt = production * waste_factor_sim  # Mt
    ) %>%
    select(year, region, scenario, route, waste_type, waste_amount_Mt)

  # 5.7 Aggregate waste across BOF and EAF and tag the simulation ID
  waste_summary <- waste_output %>%
    group_by(year, region, scenario, waste_type) %>%
    summarise(
      total_waste_mt = sum(waste_amount_Mt, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(simulation = i)

  # 5.8 Expand waste flows to technology-specific pathways
  waste_with_tech <- waste_summary %>%
    left_join(
      steel_tech_percent,
      by = c("year", "region", "waste_type"),
      relationship = "many-to-many"
    ) %>%
    mutate(
      combined_scenario = paste0(scenario, "-", tech_scenario),
      waste_to_tech = total_waste_mt * percent
    ) %>%
    select(
      simulation,
      year,
      region,
      development_scenario = scenario,
      tech_scenario,
      combined_scenario,
      waste_type,
      tech,
      waste_to_tech
    )

  # 5.9 Simulate cost-benefit parameters
  tech_param_sim <- tech_benefit %>%
    mutate(
      EB_1_sim = map_dbl(EB_1, simulate_triangle),
      EB_2_sim = map_dbl(EB_2_cost, simulate_triangle),
      EB_3_sim = map_dbl(EB_3_cost, simulate_triangle),
      EB_4_sim = map_dbl(EB_4_cost, simulate_triangle),
      EB_5_sim = map_dbl(EB_5_cost, simulate_triangle),
      EB_6_sim = map_dbl(EB_6_cost, simulate_triangle),
      EB_7_sim = map_dbl(EB_7_cost, simulate_triangle),
      EB_8_sim = map_dbl(EB_8_cost, simulate_triangle),
      EB_9_sim = map_dbl(EB_9_cost, simulate_triangle),
      EB_10_sim = map_dbl(EB_10_cost, simulate_triangle),
      EB_11_sim = map_dbl(EB_11_cost, simulate_triangle),
      EB_12_sim = map_dbl(EB_12_cost, simulate_triangle),
      EB_13_sim = map_dbl(EB_13_cost, simulate_triangle),
      EB_14_sim = map_dbl(EB_14_cost, simulate_triangle),
      EB_15_sim = map_dbl(EB_15_cost, simulate_triangle),
      EB_16_sim = map_dbl(EB_16_cost, simulate_triangle),
      EB_17_sim = map_dbl(EB_17_cost, simulate_triangle),
      EB_18_sim = map_dbl(EB_18_cost, simulate_triangle),
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

  # 5.10 Simulate carbon prices
  carbon_price_df_sim <- carbon_price_df %>%
    mutate(
      carbon_price_sim = pmap_dbl(
        list(a = CO2_price * 0.9, b = CO2_price * 1.1, c = CO2_price),
        ~ rtriangle(1, ..1, ..2, ..3)
      )
    )

  # 5.11 Calculate discounted economic and environmental benefits
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

      # Cost and product benefits:
      # waste_to_tech (Mt) * parameter (USD/t) / 1e3 = billion USD
      cost_total = waste_to_tech * cost_sim / 1e3,
      cost_total_disc = cost_total * discount_factor,
      PB_total = waste_to_tech * PB_sim / 1e3,
      PB_total_disc = PB_total * discount_factor,

      # CO2 reduction and carbon benefits:
      # EB_1 unit = kg per tonne; waste_to_tech unit = Mt
      CO2_decreasing_Gt = -EB_1_sim * waste_to_tech / 1e6,
      CO2_benefit_raw = -EB_1_sim * waste_to_tech * carbon_price_sim / 1e6,
      CO2_benefit_disc = CO2_benefit_raw * discount_factor,

      # Discounted environmental benefits (billion USD)
      EB_2_benefit_disc = -EB_2_sim * waste_to_tech * discount_factor / 1e3,
      EB_3_benefit_disc = -EB_3_sim * waste_to_tech * discount_factor / 1e3,
      EB_4_benefit_disc = -EB_4_sim * waste_to_tech * discount_factor / 1e3,
      EB_5_benefit_disc = -EB_5_sim * waste_to_tech * discount_factor / 1e3,
      EB_6_benefit_disc = -EB_6_sim * waste_to_tech * discount_factor / 1e3,
      EB_7_benefit_disc = -EB_7_sim * waste_to_tech * discount_factor / 1e3,
      EB_8_benefit_disc = -EB_8_sim * waste_to_tech * discount_factor / 1e3,
      EB_9_benefit_disc = -EB_9_sim * waste_to_tech * discount_factor / 1e3,
      EB_10_benefit_disc = -EB_10_sim * waste_to_tech * discount_factor / 1e3,
      EB_11_benefit_disc = -EB_11_sim * waste_to_tech * discount_factor / 1e3,
      EB_12_benefit_disc = -EB_12_sim * waste_to_tech * discount_factor / 1e3,
      EB_13_benefit_disc = -EB_13_sim * waste_to_tech * discount_factor / 1e3,
      EB_14_benefit_disc = -EB_14_sim * waste_to_tech * discount_factor / 1e3,
      EB_15_benefit_disc = -EB_15_sim * waste_to_tech * discount_factor / 1e3,
      EB_16_benefit_disc = -EB_16_sim * waste_to_tech * discount_factor / 1e3,
      EB_17_benefit_disc = -EB_17_sim * waste_to_tech * discount_factor / 1e3,
      EB_18_benefit_disc = -EB_18_sim * waste_to_tech * discount_factor / 1e3
    )

  # 5.12 Summarize discounted benefits by year, region, scenario, and technology scenario
  revenue_summary <- revenue_calc %>%
    group_by(simulation, year, region, development_scenario, tech_scenario, combined_scenario) %>%
    summarise(
      total_waste_mt = sum(waste_to_tech, na.rm = TRUE),
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

  # 5.13 Calculate LCA-based environmental benefits
  tech_lca_sim <- tech_benefit %>%
    mutate(
      lca_1_sim = map_dbl(EB_1, simulate_triangle),
      lca_2_sim = map_dbl(EB_2, simulate_triangle),
      lca_3_sim = map_dbl(EB_3, simulate_triangle),
      lca_4_sim = map_dbl(EB_4, simulate_triangle),
      lca_5_sim = map_dbl(EB_5, simulate_triangle),
      lca_6_sim = map_dbl(EB_6, simulate_triangle),
      lca_7_sim = map_dbl(EB_7, simulate_triangle),
      lca_8_sim = map_dbl(EB_8, simulate_triangle),
      lca_9_sim = map_dbl(EB_9, simulate_triangle),
      lca_10_sim = map_dbl(EB_10, simulate_triangle),
      lca_11_sim = map_dbl(EB_11, simulate_triangle),
      lca_12_sim = map_dbl(EB_12, simulate_triangle),
      lca_13_sim = map_dbl(EB_13, simulate_triangle),
      lca_14_sim = map_dbl(EB_14, simulate_triangle),
      lca_15_sim = map_dbl(EB_15, simulate_triangle),
      lca_16_sim = map_dbl(EB_16, simulate_triangle),
      lca_17_sim = map_dbl(EB_17, simulate_triangle),
      lca_18_sim = map_dbl(EB_18, simulate_triangle)
    )

  lca_calc <- waste_with_tech %>%
    left_join(
      tech_lca_sim,
      by = "tech"
    ) %>%
    filter(year >= 2025) %>%
    mutate(
      EB_1_benefit_lca = -lca_1_sim * waste_to_tech / 1e6,
      EB_2_benefit_lca = -lca_2_sim * waste_to_tech,
      EB_3_benefit_lca = -lca_3_sim * waste_to_tech,
      EB_4_benefit_lca = -lca_4_sim * waste_to_tech,
      EB_5_benefit_lca = -lca_5_sim * waste_to_tech,
      EB_6_benefit_lca = -lca_6_sim * waste_to_tech,
      EB_7_benefit_lca = -lca_7_sim * waste_to_tech,
      EB_8_benefit_lca = -lca_8_sim * waste_to_tech,
      EB_9_benefit_lca = -lca_9_sim * waste_to_tech,
      EB_10_benefit_lca = -lca_10_sim * waste_to_tech,
      EB_11_benefit_lca = -lca_11_sim * waste_to_tech,
      EB_12_benefit_lca = -lca_12_sim * waste_to_tech,
      EB_13_benefit_lca = -lca_13_sim * waste_to_tech,
      EB_14_benefit_lca = -lca_14_sim * waste_to_tech,
      EB_15_benefit_lca = -lca_15_sim * waste_to_tech,
      EB_16_benefit_lca = -lca_16_sim * waste_to_tech,
      EB_17_benefit_lca = -lca_17_sim * waste_to_tech,
      EB_18_benefit_lca = -lca_18_sim * waste_to_tech
    )

  lca_summary <- lca_calc %>%
    group_by(simulation, year, region, development_scenario, tech_scenario, combined_scenario) %>%
    summarise(
      total_waste_mt = sum(waste_to_tech, na.rm = TRUE),
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

  results_list[[i]] <- revenue_summary
  lca_results_list[[i]] <- lca_summary
  setTxtProgressBar(progress_bar, i)
}

close(progress_bar)

# 6. Combine all simulation outputs ---------------------------------------
all_results <- rbindlist(results_list) %>%
  filter(region != "world")

save(all_results, file = file.path(output_dir, "Steel_all_benefits.Rdata"))
save(lca_results_list, file = file.path(output_dir, "Steel_lca_results.Rdata"))

# 7. Summarize Monte Carlo uncertainty using 95% CI -----------------------
uncertainty_summary <- all_results %>%
  group_by(year, region, development_scenario, tech_scenario, combined_scenario) %>%
  summarise(
    mean_waste = mean(total_waste_mt),
    lower_95_waste = quantile(total_waste_mt, ci_lower_prob),
    upper_95_waste = quantile(total_waste_mt, ci_upper_prob),
    waste_half_width = calc_ci_half_width(lower_95_waste, upper_95_waste),

    mean_CO2_decreasing_Gt = mean(total_CO2_decreasing_Gt),
    lower_95_CO2_decreasing_Gt = quantile(total_CO2_decreasing_Gt, ci_lower_prob),
    upper_95_CO2_decreasing_Gt = quantile(total_CO2_decreasing_Gt, ci_upper_prob),
    CO2_decreasing_Gt_half_width = calc_ci_half_width(
      lower_95_CO2_decreasing_Gt,
      upper_95_CO2_decreasing_Gt
    ),

    mean_cost = mean(total_cost_all),
    lower_95_cost = quantile(total_cost_all, ci_lower_prob),
    upper_95_cost = quantile(total_cost_all, ci_upper_prob),
    cost_half_width = calc_ci_half_width(lower_95_cost, upper_95_cost),

    mean_PB = mean(total_PB_all),
    lower_95_PB = quantile(total_PB_all, ci_lower_prob),
    upper_95_PB = quantile(total_PB_all, ci_upper_prob),
    PB_half_width = calc_ci_half_width(lower_95_PB, upper_95_PB),

    mean_GWP_benefit = mean(GWP_benefit),
    lower_95_GWP_benefit = quantile(GWP_benefit, ci_lower_prob),
    upper_95_GWP_benefit = quantile(GWP_benefit, ci_upper_prob),
    GWP_benefit_half_width = calc_ci_half_width(
      lower_95_GWP_benefit,
      upper_95_GWP_benefit
    ),

    mean_Ecosystem_benefit = mean(Ecosystem_benefit),
    lower_95_Ecosystem_benefit = quantile(Ecosystem_benefit, ci_lower_prob),
    upper_95_Ecosystem_benefit = quantile(Ecosystem_benefit, ci_upper_prob),
    Ecosystem_benefit_half_width = calc_ci_half_width(
      lower_95_Ecosystem_benefit,
      upper_95_Ecosystem_benefit
    ),

    mean_Toxicity_benefit = mean(Toxicity_benefit),
    lower_95_Toxicity_benefit = quantile(Toxicity_benefit, ci_lower_prob),
    upper_95_Toxicity_benefit = quantile(Toxicity_benefit, ci_upper_prob),
    Toxicity_benefit_half_width = calc_ci_half_width(
      lower_95_Toxicity_benefit,
      upper_95_Toxicity_benefit
    ),

    mean_Resources_benefit = mean(Resources_benefit),
    lower_95_Resources_benefit = quantile(Resources_benefit, ci_lower_prob),
    upper_95_Resources_benefit = quantile(Resources_benefit, ci_upper_prob),
    Resources_benefit_half_width = calc_ci_half_width(
      lower_95_Resources_benefit,
      upper_95_Resources_benefit
    ),
    .groups = "drop"
  )

write.csv(
  uncertainty_summary,
  file.path(output_dir, "steel_summary_cost_benefit_all_benefits.csv"),
  row.names = FALSE
)
