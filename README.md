# Global industrial solid-waste recycling delivers integrated climate, environmental and economic benefits

This repository provides the original data and code for the manuscript entitled **Global industrial solid-waste recycling delivers integrated climate, environmental and economic benefits**.

## Overview

All scripts were developed and tested in **R 4.4.3**.  
The required R packages and environment configuration are described in **R environment setup.R**.

This repository includes code for estimating industrial solid-waste generation and integrated recycling benefits in the **steel** and **power** sectors across six regions:

- Brazil
- China
- India
- EU
- Middle East
- USA

## Repository contents

### 1. Waste generation scripts

- **Steel_waste_generation.R**  
  Calculates industrial solid-waste generation in the **steel sector** for the six study regions and performs **10,000 Monte Carlo simulations**.

- **Power_waste_generation.R**  
  Calculates industrial solid-waste generation in the **power sector** for the six study regions and performs **10,000 Monte Carlo simulations**.

### 2. Integrated benefit scripts

- **Steel_intergrated_benefit.R**  
  Calculates the **integrated benefits** of industrial solid-waste recycling in the steel sector, including:
  - climate benefits
  - environmental benefits
  - direct economic benefits

- **Power_intergrated_benefit.R**  
  Calculates the **integrated benefits** of industrial solid-waste recycling in the power sector, including:
  - climate benefits
  - environmental benefits
  - direct economic benefits

### 3. Data

- **Data/**  
  Contains the original input data used in this study.

### 4. Reference tables

- **Environmental indicator comparison table.docx**  
  Provides the correspondence between the **18 midpoint environmental impact indicators** and their abbreviated names used in the scripts.

- **Recycling technologies comparison table.docx**  
  Provides the correspondence between **industrial solid-waste recycling technologies** and their abbreviated names used in the scripts.

  
## Notes

1. Please make sure that all required packages are installed before running the scripts.
2. The scripts are designed to reproduce the analysis presented in the manuscript.
3. File names are kept consistent with the original project structure.


## Contact

For questions regarding the data or code, please contact the corresponding author(s) of the manuscript.
