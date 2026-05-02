# Association between BMI Variability and Incident of Dementia in Older Adults
**Author:** Liang Lai  
**Data:** Health and Retirement Study (HRS), 1996–2020
## Project Overview
This project examines whether long-term BMI variability (measured as coefficient of variation across repeated waves) is associated with incident dementia in older adults, using a landmark analysis design with the HRS longitudinal data.

## Data Access
### Raw data
HRS data requires a free registration and data use agreement.
1. Download the following files (link: <https://hrsdata.isr.umich.edu/data-products/long-format-data-rand-hrs-longitudinal-file-2020>) and place them in `data/raw-data/`:
RAND HRS Longitudinal File 2020: `rlong_table.sas7bdat`, `hlong_table.sas7bdat`, `rwide_table.sas7bdat`.
2. Link to the lange weir cognition data: <https://hrsdata.isr.umich.edu/data-products/langa-weir-classification-cognitive-function-1995-2022>.

## Reproduction Instructions

**Step 1 — Data Processing:**
- Download raw HRS files (see Data Access above. If you can't download, you can start from Step 2 with the processed data.)
- Run `code/processing-code/processingfile-hrs.qmd`

**Step 2 — Exploratory Analysis:**
- Run `code/eda-code/eda.qmd`

**Step 3 — Statistical Modeling (run in order):**
1. `code/analysis-code/statistical-analysis.R` (automatically sources `Result-Function.R`)
2. `code/analysis-code/Supplementary ML .R`
3. `code/analysis-code/Result-Table.R`
4. `code/analysis-code/Result-Figure.R`

**Step 4 — Render Manuscript:**
- Run `Product/manuscript/manuscript.qmd`

## Project Status
| Part | Status |
|------|--------|
| Part 1 — Data description and research question | Done |
| Part 2 — Data cleaning and EDA | Done |
| Part 3 — Initial analysis | Done |
| Part 4 — Full modeling | Done  |
| Part 5 — Final report for peer review | Done  |
| Part 6 — Final project revision | In progress |
