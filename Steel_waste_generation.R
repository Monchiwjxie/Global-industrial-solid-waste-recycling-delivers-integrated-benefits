library(tidyverse)
library(triangle)
library(data.table)
library(purrr)
library(truncnorm)

# 1. Configure directories and read input datasets ------------------------
project_dir <- "project_dir"
steel_input_dir <- file.path(project_dir, "data", "Steel_sector_data")
output_dir <- file.path(project_dir, "result", "steel_waste")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

steel_production <- read_csv("steel_production.csv")
steel_route_share <- read_csv("steel_route_share.csv")
waste_factors <- read_csv("steel_waste_factors.csv")

# 2. Set simulation parameters --------------------------------------------
n_sim <- 10000
set.seed(42)

# 3. Define helper functions ----------------------------------------------
sample_triangle <- function(x, lower_mult = 0.8, upper_mult = 1.2) {
  rtriangle(
    n = 1,
    a = x * lower_mult,
    b = x * upper_mult,
    c = x
  )
}

sample_truncnorm_positive <- function(mu, sigma, eps = 1e-9) {
  if (is.na(sigma) || sigma <= 0) {
    max(mu, eps)
  } else {
    rtruncnorm(
      n = 1,
      a = eps,
      b = Inf,
      mean = mu,
      sd = sigma
    )
  }
}

ci_half_width <- function(lower, upper) {
  (abs(upper) - abs(lower)) / 2
}

# 4. Monte Carlo simulation -----------------------------------------------
results_list <- vector("list", n_sim)

pb <- txtProgressBar(min = 0, max = n_sim, style = 3)

for (sim_id in seq_len(n_sim)) {

  # 4.1 Simulate steel production with triangular uncertainty
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

  # 4.2 Simulate BOF and EAF shares with triangular uncertainty
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

  # 4.3 Merge production and route shares
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

  # 4.4 Convert wide route data to long format
  steel_long_sim <- steel_data_sim %>%
    pivot_longer(
      cols = c(BOF_production, EAF_production),
      names_to = "route",
      values_to = "production"
    ) %>%
    mutate(
      route = if_else(route == "BOF_production", "BOF", "EAF")
    )

  # 4.5 Simulate waste factors with truncated normal uncertainty
  waste_factor_sample <- waste_factors %>%
    mutate(
      grouping_key = if_else(
        waste_type %in% c("BF_Slag", "Fly_Ash_Steel"),
        waste_type,
        paste(waste_type, route, sep = "|")
      )
    ) %>%
    group_by(grouping_key) %>%
    mutate(
      waste_factor_sim = sample_truncnorm_positive(
        mu = first(waste_factor),
        sigma = first(sd)
      )
    ) %>%
    ungroup() %>%
    select(waste_type, route, waste_factor_sim)

  # 4.6 Calculate waste generation by route and waste type
  waste_output <- steel_long_sim %>%
    left_join(
      waste_factor_sample,
      by = "route",
      relationship = "many-to-many"
    ) %>%
    mutate(
      waste_amount_mt = production * waste_factor_sim
    ) %>%
    select(year, region, scenario, route, waste_type, waste_amount_mt)

  # 4.7 Aggregate BOF and EAF results for each simulation
  waste_summary <- waste_output %>%
    group_by(year, region, scenario, waste_type) %>%
    summarise(
      total_waste_mt = sum(waste_amount_mt, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(simulation = sim_id)

  results_list[[sim_id]] <- waste_summary
  setTxtProgressBar(pb, sim_id)
}

close(pb)

# 5. Merge simulation outputs ---------------------------------------------
all_results <- rbindlist(results_list)

save(
  all_results,
  file = file.path(output_dir, "steel_waste_generation.Rdata")
)