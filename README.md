# Association between BMI Variability and Incident Dementia in Older Adults
**Author:** Liang Lai

**Data:** Health and Retirement Study (HRS), 1996-2020

## Project Overview
This project examines whether long-term BMI variability, measured as the coefficient of variation across repeated BMI measurements, is associated with incident dementia in older adults. The main analysis uses a landmark Cox regression design with HRS longitudinal data. A secondary predictive modeling analysis compares logistic regression, LASSO logistic regression, and random forest models for 10-year incident dementia prediction after the W3 landmark.

## Data Access
### Raw data
HRS data requires a free registration and data use agreement.
1. Download the following files (link: <https://hrsdata.isr.umich.edu/data-products/long-format-data-rand-hrs-longitudinal-file-2020>) and place them in `data/raw-data/`:
RAND HRS Longitudinal File 2020: `rlong_table.sas7bdat`, `hlong_table.sas7bdat`, `rwide_table.sas7bdat`.
2. Link to the Langa-Weir cognition data: <https://hrsdata.isr.umich.edu/data-products/langa-weir-classification-cognitive-function-1995-2022>.

## Software Requirements
The project was run in R. Main packages include `here`, `haven`, `tidyverse`, `survival`, `broom`, `gtsummary`, `flextable`, `gt`, `tidymodels`, `glmnet`, and `ranger`.

## Reproduction Instructions

**Step 1 — Data Processing:**
- Download raw HRS files and place them in `data/raw-data/` as described above. If you cannot download the raw files, start from Step 2 with the processed data already included in `data/processed-data/`.
- Run `code/processing-code/processingfile-hrs.qmd`

**Step 2 — Exploratory Analysis:**
- Run `code/eda-code/eda.qmd`

**Step 3 — Cox Regression Analysis (run in order):**
1. `code/analysis-code/statistical-analysis.R` (automatically sources `Result-Function.R`)
2. `code/analysis-code/Result-Table.R`
3. `code/analysis-code/Result-Figure.R`

**Step 4 — Predictive ML Model Comparison:**
- Run `code/analysis-code/ml-model-comparison.R`
- This script builds the W3 10-year prediction dataset, performs a 75/25 train-test split, uses 10-fold cross-validation, tunes the LASSO and random forest models, compares models against a null baseline, and saves ML tables and figures.

## Project Status
| Part | Status |
|------|--------|
| Part 1 — Data description and research question | Done |
| Part 2 — Data cleaning and EDA | Done |
| Part 3 — Initial analysis | Done |
| Part 4 — Full modeling | Done  |
| Part 5 — Final report for peer review | Done  |
| Part 6 — Final project revision | Done |
