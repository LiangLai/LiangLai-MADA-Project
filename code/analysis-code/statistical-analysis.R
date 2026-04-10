###############################
# analysis script
#################################
#this script loads the processed, cleaned data, does a simple analysis and saves the results to the results folder
#load needed packages. make sure they are installed.
library(ggplot2) #for plotting
library(broom) #for cleaning up output from lm()
library(here) #for data loading/saving
library(tidyverse)
library(haven) 
library(stringr)
library(zoo)
library(survival)
library(forcats)
library(readr)
library(purrr)
library(rlang)

#note the use of the here() package and not absolute paths
data_location <- here::here("data","processed-data","hrs_final_long.rds")

#load data. 
hrs_final <- readRDS(data_location)

#### First model fit: Landmark analysis at W3 (2000)
# Define eligible risk set for W3
# 1) be observed in all three waves: 1996, 1998, and 2000
# 2) have no cognitive impairment (CIND or dementia) up to the landmark
eligible_w3 <- hrs_final %>%
  filter(STUDYYR %in% c(1996, 1998, 2000)) %>%
  group_by(HHIDPN) %>%
  summarise(
    n_waves_present   = n_distinct(STUDYYR),
    any_cog_imp_to_lm = any(COGFUNCTION %in% c("CIND", "Demented"),
                             na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(
    n_waves_present == 3,   # Present in all 3 waves (1996, 1998, 2000)
    !any_cog_imp_to_lm      # No CIND or Dementia up to landmark
  ) %>%
  pull(HHIDPN)

cat("W3 eligible individuals:", length(eligible_w3), "\n")

# Build long survival data after landmark (2000)
# Construct the post-landmark survival dataset for W3
surv_w3 <- hrs_final %>%
  filter(HHIDPN %in% eligible_w3,
         STUDYYR >= 2000) %>%
  arrange(HHIDPN, STUDYYR) %>%
  group_by(HHIDPN) %>%
  mutate(
    tstart         = lag(STUDYYR, default = 2000),
    tstop          = STUDYYR,
    event_dementia = as.integer(COGFUNCTION == "Demented")
  ) %>%
  filter(tstart < tstop) %>%
  # Truncate rows after first dementia event
  filter(cumsum(lag(event_dementia, default = 0)) == 0) %>%
  ungroup()

cat("W3 survival data dimensions:", dim(surv_w3), "\n")
## Describe the W3 survival dataset
surv_w3 %>%
  summarise(
    n_people   = n_distinct(HHIDPN),
    n_rows     = n(),
    n_events   = n_distinct(HHIDPN[event_dementia == 1]),
    event_rate = round(n_distinct(HHIDPN[event_dementia == 1]) / 
                       n_distinct(HHIDPN) * 100, 1)
  ) %>%
print()

# For modeling, we will use a complete-case approach for the covariates of interest.
surv_w3_model <- surv_w3 %>%
  filter(
    !is.na(BMI_CV_W3),
    !is.na(BASELINE_AGE),
    !is.na(RAGENDER),
    !is.na(RACE_ETH3),
    !is.na(RAEDUC),
    !is.na(LN_WEALTH),
    !is.na(CHRONIC_N),
    !is.na(HAS_INS),
    !is.na(EMP_STATUS),
    !is.na(SMOKESTATUS),
    !is.na(R_DEPRES)
  )

write_csv(surv_w3_model, here("results", "tables", "w3_survival_data_summary.csv"))

#  Fit model 
#  Model 1: Unadjusted Cox model
#  Model 2: Adjusted Cox model with covariates
#  Model 3: Adjusted Cox model with categorical BMI variability (tertiles)
fit_w3_unadjusted <- coxph(
  Surv(tstart, tstop, event_dementia) ~
    BMI_CV_W3 + cluster(HHIDPN),
  data = surv_w3_model
)
print(summary(fit_w3_unadjusted))

# Test proportional hazards assumption
ph_w3_unadjusted <- cox.zph(fit_w3_unadjusted)
print(ph_w3_unadjusted)

# Tidy model output and save
cox_w3_unadjusted_table <- tidy(
  fit_w3_unadjusted,
  exponentiate = TRUE,
  conf.int = TRUE
) %>%
  mutate(model = "Unadjusted")

write_csv(cox_w3_unadjusted_table, here("results", "tables", "cox_w3_unadjusted_results.csv"))


# Adjusted Cox model with covariates
# This model adjusts for baseline demographic, socioeconomic, and health-related covariates.
fit_w3_full <- coxph(
  Surv(tstart, tstop, event_dementia) ~
    BMI_CV_W3 +
    BASELINE_AGE + RAGENDER + RACE_ETH3 + RAEDUC +
    LN_WEALTH + CHRONIC_N + HAS_INS +
    EMP_STATUS + SMOKESTATUS + R_DEPRES +
    cluster(HHIDPN),
  data = surv_w3_model
)
summary(fit_w3_full)

# Test proportional hazards assumption
ph_w3_full <- cox.zph(fit_w3_full)
print(ph_w3_full)

# Tidy model output and save
cox_w3_adjusted_table <- tidy(
  fit_w3_full, exponentiate = TRUE, conf.int = TRUE
) %>%
  mutate(model = "Adjusted")

write_csv(cox_w3_adjusted_table,here("results", "tables", "cox_w3_adjusted_results.csv"))

# Model 3: Adjusted Cox model with categorical BMI variability (tertiles)
fit_w3_cat <- coxph(
  Surv(tstart, tstop, event_dementia) ~
    BMIV_CAT_W3 +
    BASELINE_AGE + RAGENDER + RACE_ETH3 + RAEDUC +
    LN_WEALTH + CHRONIC_N + HAS_INS +
    EMP_STATUS + SMOKESTATUS + R_DEPRES +
    cluster(HHIDPN),
  data = surv_w3_model
)
summary(fit_w3_cat)

# Test proportional hazards assumption
ph_w3_cat <- cox.zph(fit_w3_cat)
print(ph_w3_cat)

# Tidy model output and save
cox_w3_categorical_table <- tidy(
  fit_w3_cat,
  exponentiate = TRUE,
  conf.int = TRUE
) %>%
  mutate(model = "Categorical")

write_csv(cox_w3_categorical_table, here("results", "tables", "cox_w3_categorical_results.csv"))

bind_rows(
  cox_w3_unadjusted_table,
  cox_w3_adjusted_table,
  cox_w3_categorical_table
) %>%
  write_csv(here("results", "tables", "cox_all_models_combined.csv"))

# Forest plot for categorical BMI variability groups (Model 3)
# We will plot the hazard ratios for the medium and high variability groups compared to the low variability reference group.
forest_data_cat <- cox_w3_categorical_table %>%
  mutate(
    term_label = case_when(
      term == "BMIV_CAT_W3Medium variability"                      ~ "Medium vs Low variability",
      term == "BMIV_CAT_W3High variability"                        ~ "High vs Low variability",
      term == "BASELINE_AGE"                                        ~ "Age at baseline",
      term == "RAGENDERFemale"                                      ~ "Female (ref: Male)",
      term == "RACE_ETH3Non-Hispanic Black"                         ~ "Non-Hispanic Black (ref: White)",
      term == "RACE_ETH3Hispanic or Other"                          ~ "Hispanic or Other (ref: White)",
      term == "RAEDUCHigh school graduate"                          ~ "High school graduate (ref: <HS)",
      term == "RAEDUCCollege or above"                              ~ "College or above (ref: <HS)",
      term == "RAEDUCGED"                                           ~ "GED (ref: <HS)",
      term == "LN_WEALTH"                                           ~ "Log household wealth",
      term == "CHRONIC_N"                                           ~ "No. of chronic conditions",
      term == "HAS_INSHas insurance"                                ~ "Has insurance (ref: None)",
      term == "EMP_STATUSRetired"                                   ~ "Retired (ref: Employed)",
      term == "EMP_STATUSDisabled"                                  ~ "Disabled (ref: Employed)",
      term == "EMP_STATUSUnemployed/Not in labor force"             ~ "Not in labor force (ref: Employed)",
      term == "SMOKESTATUSEver/Current smoker"                      ~ "Ever/current smoker (ref: Never)",
      term == "R_DEPRESYes"                                         ~ "Depression (ref: No)",
      TRUE ~ term
    ),
    # Assign group for color coding
    group = case_when(
      str_detect(term, "BMIV_CAT")                              ~ "Exposure",
      term %in% c("BASELINE_AGE", "RAGENDERFemale",
                  "RACE_ETH3Non-Hispanic Black",
                  "RACE_ETH3Hispanic or Other")                 ~ "Demographics",
      term %in% c("RAEDUCHigh school graduate",
                  "RAEDUCCollege or above",
                  "RAEDUCGED", "LN_WEALTH")                    ~ "Socioeconomic",
      TRUE                                                      ~ "Health"
    )
  ) %>% 
  mutate(term_label = fct_inorder(term_label))

# Draw forest plot — bottom to top, exposure highlighted in red
forest_cat <- ggplot(forest_data_cat,
                     aes(x = estimate, y = term_label, color = group)) +
  # Reference line at HR = 1
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "gray50", linewidth = 0.5) +
  # Confidence intervals
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.25, linewidth = 0.6) +
  # Point estimates
  geom_point(size = 2.8) +
  # Exposure in red, all covariates in gray
  scale_color_manual(values = c(
    "Exposure"      = "#E24B4A",
    "Demographics"  = "#555555",
    "Socioeconomic" = "#555555",
    "Health"        = "#555555"
  )) +
  scale_x_log10(breaks = c(0.5, 0.75, 1, 1.25, 1.5, 2, 3)) +
  labs(
    title    = "Adjusted hazard ratios for incident dementia — categorical BMI variability (W3)",
    subtitle = "Cox model with robust SE (clustered by individual); low variability = reference",
    x        = "Hazard ratio (log scale, reference line at HR = 1)",
    y        = NULL
  ) +
  theme_minimal() +
  theme(
    legend.position  = "none",
    panel.grid.minor = element_blank(),
    axis.text.y      = element_text(size = 9)
  )

# Save figure
ggsave(
  here("results", "figures", "figureS1_forest_cat_full.png"),
  plot   = forest_cat,
  width  = 9,
  height = 7,
  dpi    = 300
)
