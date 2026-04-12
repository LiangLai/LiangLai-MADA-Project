############################################
# Result-Table.R
# This script creates:
#   Table 1: Baseline characteristics of the W3 analytic sample by BMI variability group
#   Table 2: Association between BMI variability and incident dementia (W3, W6, W9)
# Prerequisites: Run statistical-analysis.R first
############################################

library(tidyverse)
library(here)
library(gt)
library(gtsummary)
library(survival)

# ============================================================
# TABLE 1: Baseline Characteristics by BMI Variability Group
# ============================================================


# 1. Create baseline dataset for the W3 analytic sample


# Get one row per person from the W3 analytic sample with BMIV group
table1_ids <- landmark_results$W3$surv_model %>%
  distinct(HHIDPN, BMIV_CAT_W3)

# Get 1996 baseline records
baseline_1996 <- hrs_final %>%
  filter(STUDYYR == 1996) %>%
  distinct(HHIDPN, .keep_all = TRUE) %>%
  select(
    HHIDPN,
    BASELINE_AGE,
    BASELINE_BMI,
    RAGENDER,
    RACE_ETH3,
    RAEDUC,
    LN_WEALTH,
    CHRONIC_N,
    HAS_INS,
    EMP_STATUS,
    SMOKESTATUS,
    R_DEPRES
  )

# Merge baseline data with W3 BMIV group
table1_data <- baseline_1996 %>%
  inner_join(table1_ids, by = "HHIDPN") %>%
  mutate(
    BMIV_CAT_W3 = factor(
      BMIV_CAT_W3,
      levels = c("Low variability", "Medium variability", "High variability"),
      labels = c(
        "Low (≤3.9%)",
        "Mid (>3.9% to ≤6.5%)",
        "High (>6.5%)"
      )
    ),
    # Recode chronic conditions into 3 categories
    CHRONIC_N = case_when(
      CHRONIC_N == 0 ~ "0",
      CHRONIC_N == 1 ~ "1",
      CHRONIC_N >= 2 ~ "≥2"
    ),
    CHRONIC_N   = factor(CHRONIC_N, levels = c("0", "1", "≥2")),
    RAGENDER    = factor(RAGENDER),
    RACE_ETH3   = factor(RACE_ETH3),
    RAEDUC      = factor(RAEDUC),
    HAS_INS     = factor(HAS_INS),
    EMP_STATUS  = factor(EMP_STATUS),
    SMOKESTATUS = factor(SMOKESTATUS),
    R_DEPRES    = factor(R_DEPRES)
  )


# 2. Create Table 1 with gtsummary

table1_gtsummary <- table1_data %>%
  select(
    BMIV_CAT_W3,
    BASELINE_AGE,
    BASELINE_BMI,
    RAGENDER,
    RACE_ETH3,
    RAEDUC,
    LN_WEALTH,
    CHRONIC_N,
    HAS_INS,
    EMP_STATUS,
    SMOKESTATUS,
    R_DEPRES
  ) %>%
  tbl_summary(
    by = BMIV_CAT_W3,
    statistic = list(
      all_continuous()  ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      all_continuous() ~ 2
    ),
    missing = "no",
    label = list(
      BASELINE_AGE ~ "Age at baseline",
      BASELINE_BMI ~ "Baseline BMI",
      RAGENDER     ~ "Sex",
      RACE_ETH3    ~ "Race/ethnicity",
      RAEDUC       ~ "Education",
      LN_WEALTH    ~ "Log household wealth",
      CHRONIC_N    ~ "Number of chronic conditions",
      HAS_INS      ~ "Insurance status",
      EMP_STATUS   ~ "Employment status",
      SMOKESTATUS  ~ "Smoking status",
      R_DEPRES     ~ "Depression"
    )
  ) %>%
  add_overall() %>%
  add_p() %>%
  bold_labels()


# 3. Convert to gt and add title/footnote

table1_gt <- table1_gtsummary %>%
  as_gt() %>%
  tab_header(
    title = md("**Table 1. Baseline characteristics of the W3 analytic sample by BMI variability group**")
  ) %>%
  tab_source_note(
    source_note = md(
      "Baseline characteristics were measured in 1996.
      BMI variability groups were defined using W9 tertile cutpoints: Low ≤ 3.9%, Mid > 3.9% to ≤ 6.5%, and High > 6.5%."
    )
  )

# Print
table1_gt

table1_ft <- table1_gtsummary %>%
  as_flex_table() %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 9, part = "all") %>%
  bold(part = "header") %>%
  align(align = "center", part = "header") %>%
  align(j = 1, align = "left", part = "body") %>%
  padding(padding = 3, part = "all") %>%
  line_spacing(space = 1, part = "all") %>%
  autofit() %>%
  set_table_properties(layout = "autofit", width = 1)

saveRDS(table1_ft, here("results", "tables", "table1_ft.rds"))
# 4. Save Table 1

saveRDS(table1_gt, here("results", "tables", "table1_gt.rds"))
saveRDS(table1_gtsummary, here("results", "tables", "table1_gtsummary.rds"))
write_csv(
  table1_gtsummary %>% as_tibble(),
  here("results", "tables", "table1_baseline_by_bmiv_group.csv")
)



# ============================================================
# TABLE 2: Association Between BMI Variability and Incident
#           Dementia Across Landmark Waves (W3, W6, W9)
# ============================================================


# 1. Load saved model results

cox_results <- read_csv(
  here("results", "tables", "cox_w3_w9_all_models_combined.csv"),
  show_col_types = FALSE
)

model_summary <- read_csv(
  here("results", "tables", "w3_w9_model_data_summary.csv"),
  show_col_types = FALSE
)

target_waves <- c("W3", "W6", "W9")

cox_results <- cox_results %>%
  filter(wave %in% target_waves)

model_summary <- model_summary %>%
  filter(wave %in% target_waves)

# 2. Helper functions
format_hr_ci <- function(est, low, high, digits = 2) {
  if (any(is.na(c(est, low, high)))) return("")
  sprintf(
    paste0("%.", digits, "f (%.", digits, "f, %.", digits, "f)"),
    est, low, high
  )
}

format_p <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

extract_result <- function(data, wave_value, model_value, term_value) {
  out <- data %>%
    filter(wave == wave_value, model == model_value, term == term_value)
  
  if (nrow(out) == 0) {
    return(tibble(
      estimate = NA_real_, conf.low = NA_real_,
      conf.high = NA_real_, p.value = NA_real_
    ))
  }
  
  out %>% slice(1) %>% select(estimate, conf.low, conf.high, p.value)
}


# 3. Define model structure
model_structure <- tibble(
  model_label        = c("Model 1", "Model 2"),
  continuous_source  = c("Unadjusted", "Adjusted"),
  categorical_source = c("Categorical_Unadjusted", "Categorical_Adjusted")
)


# 4. Build Table 2 body

table2_data <- crossing(
  wave = target_waves,
  model_label = model_structure$model_label
) %>%
  left_join(model_structure, by = "model_label") %>%
  mutate(
    wave_label = case_when(
      wave == "W3" ~ "Wave 3",
      wave == "W6" ~ "Wave 6",
      wave == "W9" ~ "Wave 9",
      TRUE ~ wave
    )
  ) %>%
  left_join(
    model_summary %>%
      transmute(wave, events_at_risk = paste0(n_events_model, "/", n_people_model)),
    by = "wave"
  ) %>%
  rowwise() %>%
  mutate(
    # Continuous BMIV
    continuous_term = paste0("BMI_CV_", wave),
    continuous_res  = list(extract_result(cox_results, wave, continuous_source, continuous_term)),
    continuous_hr   = format_hr_ci(continuous_res$estimate[[1]], continuous_res$conf.low[[1]], continuous_res$conf.high[[1]]),
    p_for_trend     = format_p(continuous_res$p.value[[1]]),
    
    # Categorical BMIV
    mid_term  = paste0("BMIV_CAT_", wave, "Medium variability"),
    high_term = paste0("BMIV_CAT_", wave, "High variability"),
    
    mid_res  = list(extract_result(cox_results, wave, categorical_source, mid_term)),
    high_res = list(extract_result(cox_results, wave, categorical_source, high_term)),
    
    low_cat  = "1 (ref)",
    mid_cat  = format_hr_ci(mid_res$estimate[[1]], mid_res$conf.low[[1]], mid_res$conf.high[[1]]),
    high_cat = format_hr_ci(high_res$estimate[[1]], high_res$conf.low[[1]], high_res$conf.high[[1]])
  ) %>%
  ungroup() %>%
  mutate(
    model_label = factor(model_label, levels = c("Model 1", "Model 2")),
    wave_label  = factor(wave_label, levels = c("Wave 3", "Wave 6", "Wave 9"))
  ) %>%
  arrange(wave_label, model_label) %>%
  select(wave_label, model_label, events_at_risk,
         continuous_hr, low_cat, mid_cat, high_cat, p_for_trend)


# 5. Create Table 2 with gt

table2_gt <- table2_data %>%
  gt(
    rowname_col   = "model_label",
    groupname_col = "wave_label"
  ) %>%
  cols_label(
    events_at_risk = "Events / at risk",
    continuous_hr  = md("Continuous HR (95% CI)"),
    low_cat        = "Low",
    mid_cat        = "Mid",
    high_cat       = "High",
    p_for_trend    = md("*p* for trend")
  ) %>%
  tab_spanner(
    label   = md("**BMIV categories, HR (95% CI)**"),
    columns = c(low_cat, mid_cat, high_cat)
  ) %>%
  tab_header(
    title = md("**Table 2. Association of BMI variability and risk of all-cause dementia across W3, W6, and W9**")
  ) %>%
  cols_align(
    align   = "center",
    columns = c(events_at_risk, continuous_hr, low_cat, mid_cat, high_cat, p_for_trend)
  ) %>%
  tab_options(
    table.font.size         = px(16),
    heading.title.font.size = px(18),
    row_group.font.weight   = "bold",
    data_row.padding        = px(5)
  ) %>%
  tab_source_note(
    source_note = md(
      "Model 1 is the unadjusted model and Model 2 is the adjusted model.
      Continuous HR represents the change in dementia risk per 1-unit increase in BMIV (CV, %).
      For categorical analyses, BMIV groups were defined using W9 tertile cutpoints: Low ≤ 3.9%, Mid > 3.9% to ≤ 6.5%, and High > 6.5%."
    )
  )

table2_ft <- table2_data %>%
  rename(
    Wave = wave_label,
    Model = model_label,
    `Events / at risk` = events_at_risk,
    `Continuous HR (95% CI)` = continuous_hr,
    Low = low_cat,
    Mid = mid_cat,
    High = high_cat,
    `p for trend` = p_for_trend
  ) %>%
  mutate(
    Wave = factor(Wave, levels = c("Wave 3", "Wave 6", "Wave 9")),
    Model = factor(Model, levels = c("Model 1", "Model 2"))
  ) %>%
  arrange(Wave, Model) %>%
  mutate(
    Wave = as.character(Wave),
    Model = as.character(Model)
  ) %>%
  flextable() %>%
  add_header_row(
    values = c("", "", "", "", "BMIV categories, HR (95% CI)", ""),
    colwidths = c(1, 1, 1, 1, 3, 1)
  ) %>%
  merge_v(j = "Wave") %>%
  valign(j = "Wave", valign = "top") %>%
  bold(j = "Wave", bold = TRUE, part = "body") %>%
  bold(part = "header") %>%
  align(j = c("Wave", "Model"), align = "left", part = "body") %>%
  align(
    j = c("Events / at risk", "Continuous HR (95% CI)", "Low", "Mid", "High", "p for trend"),
    align = "center",
    part = "all"
  ) %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 9, part = "all") %>%
  padding(padding = 3, part = "all") %>%
  line_spacing(space = 1, part = "all") %>%
  autofit() %>%
  set_table_properties(layout = "autofit", width = 1) %>%
  add_footer_lines(
    values = c(
      "Model 1 is the unadjusted model and Model 2 is the adjusted model.",
      "Continuous HR represents the change in dementia risk per 1-unit increase in BMIV (CV, %).",
      "BMIV groups were defined using W9 tertile cutpoints: Low ≤ 3.9%, Mid > 3.9% to ≤ 6.5%, and High > 6.5%."
    )
  ) %>%
  fontsize(size = 8, part = "footer") %>%
  italic(part = "footer")

saveRDS(table2_ft, here("results", "tables", "table2_ft.rds"))

write_csv(
  table2_data,
  here("results", "tables", "table2_cox_results_w3_w6_w9.csv")
)
