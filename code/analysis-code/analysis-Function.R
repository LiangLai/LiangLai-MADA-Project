###############################
# analysis script
#################################
#this script loads the processed, cleaned data, does a simple analysis and saves the results to the results folder
# Load required packages
library(ggplot2)
library(broom)
library(here)
library(tidyverse)
library(haven)
library(stringr)
library(zoo)
library(survival)
library(forcats)
library(readr)
library(purrr)
library(rlang)

# Define data location
data_location <- here::here("data","processed-data","hrs_final_long.rds")

# load data
hrs_final <- readRDS(data_location)


# Function to tidy PH test results
tidy_ph_test <- function(ph_obj, wave, model_label) {
  ph_df <- as.data.frame(ph_obj$table)
  ph_df$term <- rownames(ph_df)
  rownames(ph_df) <- NULL
  
  ph_df %>%
    as_tibble() %>%
    rename(
      chisq = chisq,
      p_value = p
    ) %>%
    mutate(
      wave = paste0("W", wave),
      model = model_label
    ) %>%
    select(wave, model, term, chisq, p_value, everything())
}

# Function to build landmark survival data for one wave
# For each landmark wave:
# 1) define the eligible risk set
# 2) build post-landmark survival data
# 3) apply complete-case filtering for model covariates

build_landmark_dataset <- function(data, wave) {
  
  # Convert wave number to landmark year
  # W3 = 2000, W4 = 2002, ..., W9 = 2012
  landmark_year <- 1994 + 2 * wave
  
  # Define all years that must be observed before/through the landmark
  years_needed <- seq(1996, landmark_year, by = 2)
  
  # Wave-specific exposure variable names
  bmi_cv_var   <- paste0("BMI_CV_W", wave)
  bmiv_cat_var <- paste0("BMIV_CAT_W", wave)
  
  # Step 1: Define eligible risk set
  # Eligible individuals must:
  # - be observed in every required wave up to the landmark
  # - have no CIND or dementia up to the landmark
  eligible_ids <- data %>%
    filter(STUDYYR %in% years_needed) %>%
    group_by(HHIDPN) %>%
    summarise(
      n_waves_present = n_distinct(STUDYYR),
      any_cog_imp_to_lm = any(COGFUNCTION %in% c("CIND", "Demented"), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(
      n_waves_present == length(years_needed),
      !any_cog_imp_to_lm
    ) %>%
    pull(HHIDPN)
  
  cat("Wave", wave, "- eligible individuals:", length(eligible_ids), "\n")
  
  # Step 2: Construct post-landmark survival data
  surv_data <- data %>%
    filter(
      HHIDPN %in% eligible_ids,
      STUDYYR >= landmark_year
    ) %>%
    arrange(HHIDPN, STUDYYR) %>%
    group_by(HHIDPN) %>%
    mutate(
      # Start time for each interval
      tstart = lag(STUDYYR, default = landmark_year),
      # End time for each interval
      tstop = STUDYYR,
      # Event indicator for dementia
      event_dementia = as.integer(COGFUNCTION == "Demented")
    ) %>%
    # Remove zero-length intervals
    filter(tstart < tstop) %>%
    # Keep observations only up to the first dementia event
    filter(cumsum(lag(event_dementia, default = 0)) == 0) %>%
    ungroup()
  
  cat("Wave", wave, "- survival data dimensions:", dim(surv_data), "\n")
  
  # Step 3: Create a summary of the survival dataset
  surv_summary <- surv_data %>%
    summarise(
      wave = paste0("W", wave),
      landmark_year = landmark_year,
      n_people = n_distinct(HHIDPN),
      n_rows = n(),
      n_events = n_distinct(HHIDPN[event_dementia == 1]),
      event_rate = round(
        n_distinct(HHIDPN[event_dementia == 1]) / n_distinct(HHIDPN) * 100,
        1
      )
    )
  
  print(surv_summary)
  
  # Step 4: Apply complete-case filtering for modeling
  surv_model <- surv_data %>%
    filter(
      !is.na(.data[[bmi_cv_var]]),
      !is.na(.data[[bmiv_cat_var]]),
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
  
  # Also create a summary after complete-case filtering
  model_summary <- surv_model %>%
    summarise(
      wave = paste0("W", wave),
      landmark_year = landmark_year,
      n_people_model = n_distinct(HHIDPN),
      n_rows_model = n(),
      n_events_model = n_distinct(HHIDPN[event_dementia == 1]),
      event_rate_model = round(
        n_distinct(HHIDPN[event_dementia == 1]) / n_distinct(HHIDPN) * 100,
        1
      )
    )
  
  list(
    eligible_ids = eligible_ids,
    surv_data = surv_data,
    surv_model = surv_model,
    surv_summary = surv_summary,
    model_summary = model_summary
  )
}

# Function to run Cox models for one wave
# This function fits:
# Model 1: Unadjusted continuous BMI variability
# Model 2: Adjusted continuous BMI variability
# Model 3: Adjusted categorical BMI variability

run_cox_models_by_wave <- function(data, wave) {
  
  bmi_cv_var   <- paste0("BMI_CV_W", wave)
  bmiv_cat_var <- paste0("BMIV_CAT_W", wave)
  
  # Model 1: Unadjusted continuous BMI variability
  form_unadj <- as.formula(
    paste0(
      "Surv(tstart, tstop, event_dementia) ~ ",
      bmi_cv_var,
      " + cluster(HHIDPN)"
    )
  )
  
  fit_unadj <- coxph(form_unadj, data = data)
  ph_unadj  <- cox.zph(fit_unadj)
  
  tbl_unadj <- tidy(
    fit_unadj,
    exponentiate = TRUE,
    conf.int = TRUE
  ) %>%
    mutate(
      wave = paste0("W", wave),
      model = "Unadjusted"
    )
  
  ph_tbl_unadj <- tidy_ph_test(ph_unadj, wave, "Unadjusted")
  
  # Model 2: Adjusted continuous BMI variability
  form_adj <- as.formula(
    paste0(
      "Surv(tstart, tstop, event_dementia) ~ ",
      bmi_cv_var,
      " + BASELINE_AGE + RAGENDER + RACE_ETH3 + RAEDUC + ",
      "LN_WEALTH + CHRONIC_N + HAS_INS + ",
      "EMP_STATUS + SMOKESTATUS + R_DEPRES + ",
      "cluster(HHIDPN)"
    )
  )
  
  fit_adj <- coxph(form_adj, data = data)
  ph_adj  <- cox.zph(fit_adj)
  
  tbl_adj <- tidy(
    fit_adj,
    exponentiate = TRUE,
    conf.int = TRUE
  ) %>%
    mutate(
      wave = paste0("W", wave),
      model = "Adjusted"
    )
  
  ph_tbl_adj <- tidy_ph_test(ph_adj, wave, "Adjusted")
  
  # Model 3: Adjusted categorical BMI variability
  form_cat <- as.formula(
    paste0(
      "Surv(tstart, tstop, event_dementia) ~ ",
      bmiv_cat_var,
      " + BASELINE_AGE + RAGENDER + RACE_ETH3 + RAEDUC + ",
      "LN_WEALTH + CHRONIC_N + HAS_INS + ",
      "EMP_STATUS + SMOKESTATUS + R_DEPRES + ",
      "cluster(HHIDPN)"
    )
  )
  
  fit_cat <- coxph(form_cat, data = data)
  ph_cat  <- cox.zph(fit_cat)
  
  tbl_cat <- tidy(
    fit_cat,
    exponentiate = TRUE,
    conf.int = TRUE
  ) %>%
    mutate(
      wave = paste0("W", wave),
      model = "Categorical"
    )
  
  ph_tbl_cat <- tidy_ph_test(ph_cat, wave, "Categorical")
  
  list(
    fit_unadj = fit_unadj,
    fit_adj = fit_adj,
    fit_cat = fit_cat,
    tbl_unadj = tbl_unadj,
    tbl_adj = tbl_adj,
    tbl_cat = tbl_cat,
    ph_unadj = ph_tbl_unadj,
    ph_adj = ph_tbl_adj,
    ph_cat = ph_tbl_cat
  )
}

# Build landmark datasets for W3-W9
waves <- 3:9

landmark_results <- map(waves, ~ build_landmark_dataset(hrs_final, .x))
names(landmark_results) <- paste0("W", waves)

# Combine and save survival data summaries
surv_summary_all <- map_dfr(landmark_results, "surv_summary")
model_summary_all <- map_dfr(landmark_results, "model_summary")

write_csv(
  surv_summary_all,
  here("results", "tables", "w3_w9_survival_data_summary.csv")
)

write_csv(
  model_summary_all,
  here("results", "tables", "w3_w9_model_data_summary.csv")
)

# 5. Run Cox models for W3-W9
all_model_results <- map(
  waves,
  function(w) {
    wave_name <- paste0("W", w)
    run_cox_models_by_wave(
      data = landmark_results[[wave_name]]$surv_model,
      wave = w
    )
  }
)
names(all_model_results) <- paste0("W", waves)

# 6. Combine model estimates across waves
cox_unadjusted_all <- map_dfr(all_model_results, "tbl_unadj")
cox_adjusted_all   <- map_dfr(all_model_results, "tbl_adj")
cox_categorical_all <- map_dfr(all_model_results, "tbl_cat")

cox_all_models <- bind_rows(
  cox_unadjusted_all,
  cox_adjusted_all,
  cox_categorical_all
)

# 7. Combine proportional hazards test results
ph_unadjusted_all <- map_dfr(all_model_results, "ph_unadj")
ph_adjusted_all   <- map_dfr(all_model_results, "ph_adj")
ph_categorical_all <- map_dfr(all_model_results, "ph_cat")

ph_all_models <- bind_rows(
  ph_unadjusted_all,
  ph_adjusted_all,
  ph_categorical_all
)


# 8. Save all model outputs
write_csv(
  cox_unadjusted_all,
  here("results", "tables", "cox_w3_w9_unadjusted_results.csv")
)

write_csv(
  cox_adjusted_all,
  here("results", "tables", "cox_w3_w9_adjusted_results.csv")
)

write_csv(
  cox_categorical_all,
  here("results", "tables", "cox_w3_w9_categorical_results.csv")
)

write_csv(
  cox_all_models,
  here("results", "tables", "cox_w3_w9_all_models_combined.csv")
)

write_csv(
  ph_all_models,
  here("results", "tables", "cox_w3_w9_ph_tests.csv")
)

