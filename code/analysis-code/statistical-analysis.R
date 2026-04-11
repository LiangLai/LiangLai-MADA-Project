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

# source the function to build landmark datasets and run Cox models 
source(here::here("code", "analysis-code", "Result-Function.R"))

# Run all three models for W3 and save results in a list
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

# Run Cox models for W3-W9
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

# Combine model estimates across waves
cox_unadjusted_all    <- map_dfr(all_model_results, "tbl_unadj")
cox_adjusted_all      <- map_dfr(all_model_results, "tbl_adj")
cox_cat_unadjusted_all <- map_dfr(all_model_results, "tbl_cat_unadj")
cox_cat_adjusted_all   <- map_dfr(all_model_results, "tbl_cat_adj")

cox_all_models <- bind_rows(
  cox_unadjusted_all,
  cox_adjusted_all,
  cox_cat_unadjusted_all,
  cox_cat_adjusted_all
)

# Combine proportional hazards test results
ph_unadjusted_all     <- map_dfr(all_model_results, "ph_unadj")
ph_adjusted_all       <- map_dfr(all_model_results, "ph_adj")
ph_cat_unadjusted_all <- map_dfr(all_model_results, "ph_cat_unadj")
ph_cat_adjusted_all   <- map_dfr(all_model_results, "ph_cat_adj")

ph_all_models <- bind_rows(
  ph_unadjusted_all,
  ph_adjusted_all,
  ph_cat_unadjusted_all,
  ph_cat_adjusted_all
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
  cox_cat_unadjusted_all,
  here("results", "tables", "cox_w3_w9_cat_unadjusted_results.csv")
)

write_csv(
  cox_cat_adjusted_all,
  here("results", "tables", "cox_w3_w9_cat_adjusted_results.csv")
)

write_csv(
  cox_all_models,
  here("results", "tables", "cox_w3_w9_all_models_combined.csv")
)

write_csv(
  ph_all_models,
  here("results", "tables", "cox_w3_w9_ph_tests.csv")
)
