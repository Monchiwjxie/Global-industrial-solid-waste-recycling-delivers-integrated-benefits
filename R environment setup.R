# ============================================================
# Reproducible R environment setup
# Expected R version: 4.4.3
# ============================================================

expected_r_version <- "4.4.3"

current_r_version <- paste(R.version$major, R.version$minor, sep = ".")
message("Current R version: ", current_r_version)
message("Expected R version: ", expected_r_version)

# ------------------------------------------------------------
# Required packages and target versions
# ------------------------------------------------------------
required_packages <- c(
  "tidyverse"  = "2.0.0",
  "triangle"   = "1.0",
  "data.table" = "1.17.2",
  "progress"   = "1.2.3",
  "truncnorm"  = "1.0-9",
  "purrr"      = "1.0.4"
)

# ------------------------------------------------------------
# Install/check all required packages
# ------------------------------------------------------------
for (pkg in names(required_packages)) {
  install_if_needed(pkg, required_packages[[pkg]])
}

# ------------------------------------------------------------
# Load packages
# ------------------------------------------------------------
library(tidyverse)
library(triangle)
library(data.table)
library(progress)
library(truncnorm)
library(purrr)
