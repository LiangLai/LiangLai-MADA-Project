###############################
# analysis script
#
#this script loads the processed, cleaned data, does a simple analysis
#and saves the results to the results folder

#load needed packages. make sure they are installed.
library(ggplot2) #for plotting
library(broom) #for cleaning up output from lm()
library(here) #for data loading/saving
library(tidyverse)
library(haven) 
library(stringr)
library(zoo)
#note the use of the here() package and not absolute paths
data_location <- here::here("data","processed-data","hrs_final_long.rds")

#load data. 
hrs_final <- readRDS(data_location)

#### First model fit
# fit linear model using height as outcome, weight as predictor
# Define eligible risk set for W3
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

surv_w3 %>%
  summarise(
    n_people   = n_distinct(HHIDPN),
    n_rows     = n(),
    n_events   = n_distinct(HHIDPN[event_dementia == 1]),
    event_rate = round(n_distinct(HHIDPN[event_dementia == 1]) / 
                       n_distinct(HHIDPN) * 100, 1)
  ) %>%
  print()

#  Cox model
library(survival)

fit_w3 <- coxph(
  Surv(tstart, tstop, event_dementia) ~
    BMI_CV_W3 +
    BASELINE_AGE + RAGENDER + RACE_ETH3 + RAEDUC +
    LN_WEALTH + CHRONIC_N + HAS_INS +
    EMP_STATUS + SMOKESTATUS + R_DEPRES +
    cluster(HHIDPN),
  data = surv_w3
)
summary(fit_w3)
cox.zph(fit_w3) 
cox_table <- tidy(fit_w3, exponentiate = TRUE, conf.int = TRUE)
write_csv(cox_table, here("results", "tables", "cox_w3_results.csv"))