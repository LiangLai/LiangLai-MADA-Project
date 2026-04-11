############################################
# Result-RSF.R
# Supplementary Random Survival Forest analysis
# W3 and W6 analytic samples
#
# Prerequisites: Run statistical-analysis.R first
#   so that landmark_results are in your environment.
############################################

library(tidyverse)
library(here)
library(readr)
library(randomForestSRC)


############################################
# PART 1: W3 Random Survival Forest
############################################

#-------------------------------------------
# 1. Build a person-level W3 dataset
#-------------------------------------------
# landmark_results$W3$surv_model is long-format survival data
# We convert it to one row per person for RSF.

w3_rsf_data <- landmark_results$W3$surv_model %>%
  group_by(HHIDPN) %>%
  summarise(
    # Follow-up time after the W3 landmark (year 2000)
    time_to_event = max(tstop) - 2000,
    
    # Event indicator: 1 if dementia occurred during follow-up
    event_dementia = max(event_dementia),
    
    # Exposure
    BMI_CV_W3 = first(BMI_CV_W3),
    
    # Baseline covariates
    BASELINE_AGE = first(BASELINE_AGE),
    RAGENDER     = first(RAGENDER),
    RACE_ETH3    = first(RACE_ETH3),
    RAEDUC       = first(RAEDUC),
    LN_WEALTH    = first(LN_WEALTH),
    CHRONIC_N    = first(CHRONIC_N),
    HAS_INS      = first(HAS_INS),
    EMP_STATUS   = first(EMP_STATUS),
    SMOKESTATUS  = first(SMOKESTATUS),
    R_DEPRES     = first(R_DEPRES),
    .groups = "drop"
  ) %>%
  mutate(
    RAGENDER    = factor(RAGENDER),
    RACE_ETH3   = factor(RACE_ETH3),
    RAEDUC      = factor(RAEDUC),
    HAS_INS     = factor(HAS_INS),
    EMP_STATUS  = factor(EMP_STATUS),
    SMOKESTATUS = factor(SMOKESTATUS),
    R_DEPRES    = factor(R_DEPRES)
  )

# Quick check
cat("W3 RSF data:", dim(w3_rsf_data)[1], "persons,",
    sum(w3_rsf_data$event_dementia), "events\n\n")

#-------------------------------------------
# 2. Fit W3 Random Survival Forest
#-------------------------------------------
set.seed(1234)

rsf_w3 <- rfsrc(
  Surv(time_to_event, event_dementia) ~
    BMI_CV_W3 +
    BASELINE_AGE + RAGENDER + RACE_ETH3 + RAEDUC +
    LN_WEALTH + CHRONIC_N + HAS_INS +
    EMP_STATUS + SMOKESTATUS + R_DEPRES,
  data = w3_rsf_data,
  ntree = 1000,
  importance = "permute"
)

#-------------------------------------------
# 3. W3 RSF summary
#-------------------------------------------
rsf_w3_summary <- tibble(
  analysis             = "W3 Random Survival Forest",
  n                    = nrow(w3_rsf_data),
  n_events             = sum(w3_rsf_data$event_dementia),
  ntree                = rsf_w3$ntree,
  oob_prediction_error = tail(rsf_w3$err.rate, 1)
)

print(rsf_w3_summary)

write_csv(
  rsf_w3_summary,
  here("results", "tables", "rsf_w3_summary.csv")
)

#-------------------------------------------
# 4. W3 Variable importance
#-------------------------------------------
vimp_w3 <- tibble(
  variable   = names(rsf_w3$importance),
  importance = as.numeric(rsf_w3$importance)
) %>%
  mutate(label = case_when(
    variable == "BMI_CV_W3"    ~ "BMI Variability (CV%)",
    variable == "BASELINE_AGE" ~ "Age at Baseline",
    variable == "RAGENDER"     ~ "Sex",
    variable == "RACE_ETH3"    ~ "Race/Ethnicity",
    variable == "RAEDUC"       ~ "Education",
    variable == "LN_WEALTH"    ~ "Log Wealth",
    variable == "CHRONIC_N"    ~ "Chronic Conditions",
    variable == "HAS_INS"      ~ "Insurance Status",
    variable == "EMP_STATUS"   ~ "Employment Status",
    variable == "SMOKESTATUS"  ~ "Smoking Status",
    variable == "R_DEPRES"     ~ "Depression",
    TRUE ~ variable
  )) %>%
  arrange(desc(importance))

print(vimp_w3)

write_csv(
  vimp_w3,
  here("results", "tables", "rsf_w3_variable_importance.csv")
)

#-------------------------------------------
# 5. W3 Variable importance plot
#-------------------------------------------
vimp_w3_plot <- vimp_w3 %>%
  mutate(label = forcats::fct_reorder(label, importance)) %>%
  ggplot(aes(x = importance, y = label)) +
  geom_col(fill = "#2c7bb6", width = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title = "Variable Importance from W3 Random Survival Forest",
    x = "Permutation Importance",
    y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "gray40", size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )

print(vimp_w3_plot)

ggsave(
  filename = here("results", "figures", "figureS_rsf_w3_variable_importance.png"),
  plot = vimp_w3_plot,
  width = 8,
  height = 5,
  dpi = 300
)

cat("\n--- W3 RSF analysis complete ---\n\n")


############################################
# PART 2: W6 Random Survival Forest
############################################

#-------------------------------------------
# 1. Build a person-level W6 dataset
#-------------------------------------------
# landmark_results$W6$surv_model is long-format survival data
# We convert it to one row per person for RSF.

w6_rsf_data <- landmark_results$W6$surv_model %>%
  group_by(HHIDPN) %>%
  summarise(
    # Follow-up time after the W6 landmark (year 2006)
    time_to_event = max(tstop) - 2006,
    
    # Event indicator: 1 if dementia occurred during follow-up
    event_dementia = max(event_dementia),
    
    # Exposure
    BMI_CV_W6 = first(BMI_CV_W6),
    
    # Baseline covariates
    BASELINE_AGE = first(BASELINE_AGE),
    RAGENDER     = first(RAGENDER),
    RACE_ETH3    = first(RACE_ETH3),
    RAEDUC       = first(RAEDUC),
    LN_WEALTH    = first(LN_WEALTH),
    CHRONIC_N    = first(CHRONIC_N),
    HAS_INS      = first(HAS_INS),
    EMP_STATUS   = first(EMP_STATUS),
    SMOKESTATUS  = first(SMOKESTATUS),
    R_DEPRES     = first(R_DEPRES),
    .groups = "drop"
  ) %>%
  mutate(
    RAGENDER    = factor(RAGENDER),
    RACE_ETH3   = factor(RACE_ETH3),
    RAEDUC      = factor(RAEDUC),
    HAS_INS     = factor(HAS_INS),
    EMP_STATUS  = factor(EMP_STATUS),
    SMOKESTATUS = factor(SMOKESTATUS),
    R_DEPRES    = factor(R_DEPRES)
  )

# Quick check
cat("W6 RSF data:", dim(w6_rsf_data)[1], "persons,",
    sum(w6_rsf_data$event_dementia), "events\n\n")

#-------------------------------------------
# 2. Fit W6 Random Survival Forest
#-------------------------------------------
set.seed(1234)

rsf_w6 <- rfsrc(
  Surv(time_to_event, event_dementia) ~
    BMI_CV_W6 +
    BASELINE_AGE + RAGENDER + RACE_ETH3 + RAEDUC +
    LN_WEALTH + CHRONIC_N + HAS_INS +
    EMP_STATUS + SMOKESTATUS + R_DEPRES,
  data = w6_rsf_data,
  ntree = 1000,
  importance = "permute"
)

#-------------------------------------------
# 3. W6 RSF summary
#-------------------------------------------
rsf_w6_summary <- tibble(
  analysis             = "W6 Random Survival Forest",
  n                    = nrow(w6_rsf_data),
  n_events             = sum(w6_rsf_data$event_dementia),
  ntree                = rsf_w6$ntree,
  oob_prediction_error = tail(rsf_w6$err.rate, 1)
)

print(rsf_w6_summary)

write_csv(
  rsf_w6_summary,
  here("results", "tables", "rsf_w6_summary.csv")
)

#-------------------------------------------
# 4. W6 Variable importance
#-------------------------------------------
vimp_w6 <- tibble(
  variable   = names(rsf_w6$importance),
  importance = as.numeric(rsf_w6$importance)
) %>%
  mutate(label = case_when(
    variable == "BMI_CV_W6"    ~ "BMI Variability (CV%)",
    variable == "BASELINE_AGE" ~ "Age at Baseline",
    variable == "RAGENDER"     ~ "Sex",
    variable == "RACE_ETH3"    ~ "Race/Ethnicity",
    variable == "RAEDUC"       ~ "Education",
    variable == "LN_WEALTH"    ~ "Log Wealth",
    variable == "CHRONIC_N"    ~ "Chronic Conditions",
    variable == "HAS_INS"      ~ "Insurance Status",
    variable == "EMP_STATUS"   ~ "Employment Status",
    variable == "SMOKESTATUS"  ~ "Smoking Status",
    variable == "R_DEPRES"     ~ "Depression",
    TRUE ~ variable
  )) %>%
  arrange(desc(importance))

print(vimp_w6)

write_csv(
  vimp_w6,
  here("results", "tables", "rsf_w6_variable_importance.csv")
)

#-------------------------------------------
# 5. W6 Variable importance plot
#-------------------------------------------
vimp_w6_plot <- vimp_w6 %>%
  mutate(label = forcats::fct_reorder(label, importance)) %>%
  ggplot(aes(x = importance, y = label)) +
  geom_col(fill = "#2c7bb6", width = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title = "Variable Importance from W6 Random Survival Forest",
    x = "Permutation Importance",
    y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "gray40", size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )

print(vimp_w6_plot)

ggsave(
  filename = here("results", "figures", "figureS_rsf_w6_variable_importance.png"),
  plot = vimp_w6_plot,
  width = 8,
  height = 5,
  dpi = 300
)