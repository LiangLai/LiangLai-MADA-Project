# Overview
This is Liang Lai's class project repository for analyzing the relationship between BMIV and cognitive decline in older adults using the Health and Retirement Study (HRS) data.

## Project Status
###  Completed
- Data loading and merging (RAND HRS Longitudinal File 2020 + Langa-Weir Cognitive Classification Dataset)
- Data cleaning and wrangling
- Variable creation and recoding
- Final processed dataset saved as `hrs_cognition_final.rds`

###  In Progress
- Exploratory Data Analysis (EDA)

###  Upcoming
- Statistical modeling and regression analysis
- Manuscript writing


## Reproduction Instructions
1. For **full reproduction** (including data processing): Download raw data files as described in `code/processing-code/README.md`, then run `processingfile-hrs.qmd`
2. For **EDA and analysis only**: The processed dataset `hrs_cognition_final.rds` is available in `data/processed-data/`. Run the EDA scripts directly.