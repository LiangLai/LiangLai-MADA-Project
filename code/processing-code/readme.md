# Processing Code README
##  Important Note on Data Access
The raw data files used in this processing script are **not included in this repository** due to their large file sizes. To reproduce the full processing pipeline, you will need to download the following datasets manually:
### Required Downloads
**1. RAND HRS Longitudinal File 2020 (V2)**
- Files needed: `rlong_table.sas7bdat`, `hlong_table.sas7bdat`, `rwide_table.sas7bdat`
- Download from: [RAND HRS](https://hrsdata.isr.umich.edu/data-products/long-format-data-rand-hrs-longitudinal-file-2020)

**2. Langa-Weir Cognitive Classification Dataset (1995-2022)**
- File needed: `cogfinallimp_9522wide.dta`
- Download from: [HRS Website](https://hrsdata.isr.umich.edu/data-products/langa-weir-classification-cognitive-function-1995-2022)
Then put them into correct position

## Overview
This folder contains the data processing and cleaning scripts for the HRS cognition study.
## Files
- `processingfile-hrs.qmd` : Main data processing script (loading, cleaning, wrangling)

## Data Sources
- `rlong_table.sas7bdat` : RAND HRS Longitudinal File 2020 (respondent wave-varying variables)
- `hlong_table.sas7bdat` : RAND HRS Longitudinal File 2020 (household wave-varying variables)
- `rwide_table.sas7bdat` : RAND HRS Longitudinal File 2020 (time-invariant demographic variables)
- `cogfinallimp_9522wide.dta` : Langa-Weir Cognitive Classification Dataset (1995-2022)

## Processing Steps
### Step 1: Load Data
### Step 2: Reshape Langa-Weir Data
### Step 3: Merge HRS Data and Langa-Weir Data
### Step 4: Select Variables
- Keep only variables needed for analysis based on codebook
- Includes demographics, health conditions, insurance, cognitive outcomes

### Step 5: Data Cleaning and Create New Variables
- Drop variables with >50% missing (`R_GOVOT`, `R_PMBMI`, `RAWTSAMP`, `cogivewmode`)
- Recode chronic condition variables (0/1/3/4/5/6) to binary (0/1/NA)
- Create `CHRONIC_N`: count of chronic conditions (hypertension, diabetes, cancer, lung disease, heart disease, stroke)
- Create `HAS_INS`: any health insurance coverage (0=No, 1=Yes, NA=all missing)
- Create `AGE`: age at interview (`STUDYYR - RABYEAR`), filter to age 50+
- `RACE_ETH3`: race/ethnicity (Non-Hispanic White, Non-Hispanic Black, Hispanic or Other)
- `SMOKESTATUS`: smoking status (Never smoker, Ever/Current smoker)
- `EMP_STATUS`: employment status (Employed, Unemployed/Not in LF, Retired, Disabled)
- `BMI_CAT`: BMI categories (Thin, Normal, Overweight, Obese)
- `LN_WEALTH`: log household wealth using `asinh()` transformation
- `LN_HITOT`: log total income using `log(H_ITOT + 1)`
- `BASELINE_AGE`: age at baseline (1996 wave)
- `BASELINE_BMI`: BMI at baseline (1996 wave)
- Recode `RAEDUC`: merge categories 4 and 5 into "College or above"
- Convert categorical variables to factor with meaningful labels

### Step 6: Final Dataset
- Select only variables needed for regression analysis
- Filter to individuals with baseline data (1996 wave)
- Save as `hrs_cognition_final.rds`

## Output
- **10,324 unique individuals**
- **88,879 person-wave observations**
- **13 waves** (1996-2020, biennial)
- **Average 8.6 waves per person**