############################################
# ml-model-comparison.R
# Predictive model comparison for incident dementia
# W3 landmark, 10-year prediction horizon
############################################

library(tidyverse)
library(tidymodels)
library(here)
library(readr)
library(ranger)
library(glmnet)

# Fixed seed so the train/test split, CV folds, and tuning grids are reproducible.
set.seed(20260502)

# Create output folders before writing any generated tables or figures.
dir.create(here("results", "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("results", "figures"), recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. Load processed longitudinal dataset
# ------------------------------------------------------------

hrs_final <- readRDS(here("data", "processed-data", "hrs_final_long.rds"))

# ------------------------------------------------------------
# 2. Build W3 landmark prediction dataset
# ------------------------------------------------------------

landmark_year <- 2000
prediction_horizon_year <- landmark_year + 10
years_needed <- c(1996, 1998, 2000)

# Eligible participants are observed through W3 and free of CIND/dementia through the landmark wave. This mirrors the landmark design used in the Cox analysis and ensures predictors are measured before the prediction period.
eligible_w3_ids <- hrs_final %>%
  filter(STUDYYR %in% years_needed) %>%
  group_by(HHIDPN) %>%
  summarise(
    n_waves_present = n_distinct(STUDYYR),
    any_cog_imp_to_landmark = any(COGFUNCTION %in% c("CIND", "Demented"), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_waves_present == length(years_needed), !any_cog_imp_to_landmark) %>%
  pull(HHIDPN)

# Follow-up is restricted to the fixed 10-year prediction horizon. This avoids labeling participants with longer follow-up as cases simply because they were observed for more time.
w3_surv_data <- hrs_final %>%
  filter(
    HHIDPN %in% eligible_w3_ids,
    STUDYYR >= landmark_year,
    STUDYYR <= prediction_horizon_year
  ) %>%
  arrange(HHIDPN, STUDYYR) %>%
  group_by(HHIDPN) %>%
  mutate(
    # Counting-process style intervals are used only to identify incident
    # dementia before the 10-year horizon.
    tstart = lag(STUDYYR, default = landmark_year),
    tstop = STUDYYR,
    event_dementia_interval = as.integer(COGFUNCTION == "Demented")
  ) %>%
  filter(tstart < tstop) %>%
  # Keep the dementia event interval but remove any later intervals.
  filter(cumsum(lag(event_dementia_interval, default = 0)) == 0) %>%
  ungroup() %>%
  # Complete-case filtering keeps the modeling comparison simple and ensures all models are trained on the same observations.
  filter(
    !is.na(BMI_CV_W3),
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

# Convert follow-up records to one row per person. Cases are people with dementia during the 10-year window. Non-events are retained only if they were observed through the horizon, so their binary outcome status is known.
w3_person_data <- w3_surv_data %>%
  group_by(HHIDPN) %>%
  summarise(
    event_dementia = max(event_dementia_interval, na.rm = TRUE),
    last_observed_year = max(tstop, na.rm = TRUE),
    BMI_CV_W3 = first(BMI_CV_W3),
    BASELINE_AGE = first(BASELINE_AGE),
    RAGENDER = first(RAGENDER),
    RACE_ETH3 = first(RACE_ETH3),
    RAEDUC = first(RAEDUC),
    LN_WEALTH = first(LN_WEALTH),
    CHRONIC_N = first(CHRONIC_N),
    HAS_INS = first(HAS_INS),
    EMP_STATUS = first(EMP_STATUS),
    SMOKESTATUS = first(SMOKESTATUS),
    R_DEPRES = first(R_DEPRES),
    .groups = "drop"
  ) %>%
  filter(event_dementia == 1 | last_observed_year >= prediction_horizon_year) %>%
  mutate(
    # Put the event class first so yardstick treats "dementia" as the positive event by default.
    event_dementia = factor(
      if_else(event_dementia == 1, "dementia", "no_dementia"),
      levels = c("dementia", "no_dementia")
    ),
    across(
      c(RAGENDER, RACE_ETH3, RAEDUC, HAS_INS, EMP_STATUS, SMOKESTATUS, R_DEPRES),
      as.factor
    )
  ) %>%
  select(-last_observed_year)

# Save the exact ML-ready dataset for reproducibility.
write_csv(w3_person_data, here("results", "tables", "ml_w3_prediction_dataset.csv"))

# Save basic sample counts
ml_sample_summary <- w3_person_data %>%
  summarise(
    landmark = "W3",
    landmark_year = landmark_year,
    prediction_horizon_year = prediction_horizon_year,
    n_people = n(),
    n_events = sum(event_dementia == "dementia"),
    event_rate = mean(event_dementia == "dementia")
  )

write_csv(ml_sample_summary, here("results", "tables", "ml_w3_sample_summary.csv"))

# ------------------------------------------------------------
# 3. Train/test split and cross-validation folds
# ------------------------------------------------------------

data_split <- initial_split(w3_person_data, prop = 0.75, strata = event_dementia)
train_data <- training(data_split)
test_data <- testing(data_split)

# Ten-fold CV is used within the training data for tuning and model comparison.
# The held-out test set is used only for final performance evaluation.
cv_folds <- vfold_cv(train_data, v = 10, strata = event_dementia)

# ------------------------------------------------------------
# 4. Preprocessing recipes
# ------------------------------------------------------------

# The base recipe handles categorical variables and removes zero-variance predictors.
base_recipe <- recipe(event_dementia ~ ., data = train_data) %>%
  update_role(HHIDPN, new_role = "id") %>%
  step_unknown(all_nominal_predictors()) %>%
  step_novel(all_nominal_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_zv(all_predictors())

# Linear and penalized models are normalized after dummy encoding.
linear_recipe <- base_recipe %>%
  step_normalize(all_numeric_predictors())

# Tree-based models do not require numeric normalization.
tree_recipe <- base_recipe

# Count encoded predictors so the random-forest mtry grid has a valid range.
n_encoded_predictors <- tree_recipe %>%
  prep(training = train_data) %>%
  bake(new_data = NULL) %>%
  select(-event_dementia, -HHIDPN) %>%
  ncol()

# ------------------------------------------------------------
# 5. Model specifications
# ------------------------------------------------------------

# Standard logistic regression is the full multivariable baseline model.
logistic_spec <- logistic_reg() %>%
  set_engine("glm") %>%
  set_mode("classification")

# LASSO logistic regression tunes the penalty and performs regularization.
lasso_spec <- logistic_reg(penalty = tune(), mixture = 1) %>%
  set_engine("glmnet") %>%
  set_mode("classification")

# Random forest is a nonlinear ensemble model. Permutation importance is requested for interpretation.
rf_spec <- rand_forest(trees = 500, mtry = tune(), min_n = tune()) %>%
  set_engine("ranger", importance = "permutation") %>%
  set_mode("classification")

# ------------------------------------------------------------
# 6. Workflows, metrics, and tuning grids
# ------------------------------------------------------------

logistic_wf <- workflow() %>% add_recipe(linear_recipe) %>% add_model(logistic_spec)
lasso_wf <- workflow() %>% add_recipe(linear_recipe) %>% add_model(lasso_spec)
rf_wf <- workflow() %>% add_recipe(tree_recipe) %>% add_model(rf_spec)

# The primary comparison metric is ROC AUC, but PR AUC and Brier score are also included for a more complete picture of model performance. For Brier score, lower is better; for ROC AUC and PR AUC, higher is better.
ml_metrics <- metric_set(roc_auc, pr_auc, brier_class)

# LASSO tunes over penalty values on the log10 scale.
lasso_grid <- grid_regular(penalty(range = c(-4, 0)), levels = 20)

# Random forest tuning grids are generated with a space-filling algorithm to efficiently explore the hyperparameter space. The mtry range is set from 2 to the total number of predictors, and the min_n range is set from 5 to 40 based on common defaults and the sample size.
rf_grid <- grid_space_filling(
  mtry(range = c(2L, min(20L, n_encoded_predictors))),
  min_n(range = c(5L, 40L)),
  size = 12
)

# ------------------------------------------------------------
# 7. Cross-validation and tuning
# ------------------------------------------------------------

# Logistic regression has no tuning parameters, so it is evaluated across the same CV folds used by the tuned models.
logistic_res <- fit_resamples(
  logistic_wf,
  resamples = cv_folds,
  metrics = ml_metrics,
  control = control_resamples(save_pred = TRUE)
)

# Tune the LASSO penalty by cross-validation.
lasso_res <- tune_grid(
  lasso_wf,
  resamples = cv_folds,
  grid = lasso_grid,
  metrics = ml_metrics,
  control = control_grid(save_pred = TRUE, save_workflow = TRUE)
)

# Tune random-forest hyperparameters by cross-validation.
rf_res <- tune_grid(
  rf_wf,
  resamples = cv_folds,
  grid = rf_grid,
  metrics = ml_metrics,
  control = control_grid(save_pred = TRUE, save_workflow = TRUE)
)

# Null model: predict the training event rate in each CV fold. For this constant-risk model, ROC AUC is 0.5, PR AUC equals the event prevalence in the assessment fold, and the Brier score is computed directly.
null_cv_performance <- map2_dfr(cv_folds$splits, cv_folds$id, function(split, fold_id) {
  fold_train <- analysis(split)
  fold_assess <- assessment(split)
  p_train <- mean(fold_train$event_dementia == "dementia")
  p_assess <- mean(fold_assess$event_dementia == "dementia")
  brier <- mean((as.integer(fold_assess$event_dementia == "dementia") - p_train)^2)

  tibble(
    model = "Null model",
    .metric = c("roc_auc", "pr_auc", "brier_class"),
    .estimate = c(0.5, p_assess, brier),
    fold_id = fold_id
  )
}) %>%
  group_by(model, .metric) %>%
  summarise(
    mean = mean(.estimate),
    std_err = sd(.estimate) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

# Choose final tuning parameters using CV ROC AUC as the primary metric.
best_lasso <- select_best(lasso_res, metric = "roc_auc")
best_rf <- select_best(rf_res, metric = "roc_auc")

# Save selected hyperparameters so the tuning decision is transparent.
write_csv(best_lasso, here("results", "tables", "ml_best_lasso_params.csv"))
write_csv(best_rf, here("results", "tables", "ml_best_rf_params.csv"))

# Pull the metrics corresponding to the selected tuning parameters.
collect_selected_metrics <- function(res, model_name, best_params = NULL) {
  metrics <- collect_metrics(res)

  if (!is.null(best_params)) {
    tuning_cols <- intersect(names(best_params), names(metrics))
    metrics <- metrics %>% inner_join(best_params, by = tuning_cols)
  }

  metrics %>%
    mutate(model = model_name, .before = 1) %>%
    select(model, .metric, mean, std_err, n)
}

# Combine null baseline and model CV performance into one table. For Brier score, lower is better; for ROC AUC and PR AUC, higher is better.
cv_performance <- bind_rows(
  null_cv_performance,
  collect_selected_metrics(logistic_res, "Logistic regression"),
  collect_selected_metrics(lasso_res, "LASSO logistic regression", best_lasso),
  collect_selected_metrics(rf_res, "Random forest", best_rf)
) %>%
  mutate(sort_value = if_else(.metric == "brier_class", -mean, mean)) %>%
  arrange(.metric, desc(sort_value)) %>%
  select(-sort_value)

write_csv(cv_performance, here("results", "tables", "ml_cv_performance.csv"))

# ------------------------------------------------------------
# 8. Final model fitting and test-set evaluation
# ------------------------------------------------------------

# Fit each final model on the full training set before evaluating on the test set.
final_logistic_fit <- fit(logistic_wf, data = train_data)
final_lasso_fit <- fit(finalize_workflow(lasso_wf, best_lasso), data = train_data)
final_rf_fit <- fit(finalize_workflow(rf_wf, best_rf), data = train_data)

# Return dementia probabilities for the test set. Class predictions are not needed because performance is evaluated with threshold-free metrics.
get_test_predictions <- function(fit_obj, model_name, new_data) {
  predict(fit_obj, new_data = new_data, type = "prob") %>%
    bind_cols(new_data %>% select(event_dementia)) %>%
    mutate(model = model_name, .before = 1)
}

test_predictions <- bind_rows(
  get_test_predictions(final_logistic_fit, "Logistic regression", test_data),
  get_test_predictions(final_lasso_fit, "LASSO logistic regression", test_data),
  get_test_predictions(final_rf_fit, "Random forest", test_data)
)

# The null test performance is computed analytically from the train event rate and test event prevalence.
test_event_rate <- mean(test_data$event_dementia == "dementia")
null_train_rate <- mean(train_data$event_dementia == "dementia")
null_test_brier <- mean((as.integer(test_data$event_dementia == "dementia") - null_train_rate)^2)

# Combine null and fitted-model test performance for direct comparison.
test_performance <- bind_rows(
  tibble(
    model = "Null model",
    .metric = c("roc_auc", "pr_auc", "brier_class"),
    .estimator = "binary",
    .estimate = c(0.5, test_event_rate, null_test_brier)
  ),
  test_predictions %>%
    group_by(model) %>%
    ml_metrics(truth = event_dementia, .pred_dementia) %>%
    ungroup()
) %>%
  mutate(sort_value = if_else(.metric == "brier_class", -.estimate, .estimate)) %>%
  arrange(.metric, desc(sort_value)) %>%
  select(-sort_value)

write_csv(test_predictions, here("results", "tables", "ml_test_predictions.csv"))
write_csv(test_performance, here("results", "tables", "ml_test_performance.csv"))

# ------------------------------------------------------------
# 9. Variable importance and figures
# ------------------------------------------------------------

# Random forest permutation importance is saved as an interpretability summary.
variable_importance <- final_rf_fit %>%
  extract_fit_engine() %>%
  pluck("variable.importance") %>%
  enframe(name = "feature", value = "importance") %>%
  arrange(desc(importance)) %>%
  slice_head(n = 20)

write_csv(variable_importance, here("results", "tables", "ml_rf_variable_importance.csv"))

model_order <- c("Null model", "Logistic regression", "LASSO logistic regression", "Random forest")

# Figure 1 compares CV performance against the null model baseline.
model_comparison_data <- cv_performance %>%
  mutate(
    metric_label = recode(
      .metric,
      roc_auc = "ROC AUC (higher is better)",
      pr_auc = "PR AUC (higher is better)",
      brier_class = "Brier score (lower is better)"
    ),
    model = factor(model, levels = model_order)
  )

# Use the null model as a red dashed reference line in each metric panel.
null_reference <- model_comparison_data %>%
  filter(model == "Null model") %>%
  transmute(metric_label, null_mean = mean)

model_comparison_plot <- model_comparison_data %>%
  filter(model != "Null model") %>%
  ggplot(aes(x = mean, y = model)) +
  geom_vline(
    data = null_reference,
    aes(xintercept = null_mean),
    linetype = "dashed",
    color = "firebrick",
    linewidth = 0.6,
    inherit.aes = FALSE
  ) +
  geom_errorbar(
    aes(xmin = mean - 1.96 * std_err, xmax = mean + 1.96 * std_err),
    width = 0.2,
    color = "gray40"
  ) +
  geom_point(size = 3, color = "#1f78b4") +
  facet_wrap(~ metric_label, scales = "free_x") +
  labs(
    title = "Cross-Validated Model Performance",
    subtitle = "Red dashed line = null model; error bars = mean +/- 1.96 SE across 10 CV folds",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "gray40", size = 10)
  )

ggsave(here("results", "figures", "ml_model_comparison.png"), model_comparison_plot, width = 10, height = 5, dpi = 300)

# Figure 2 shows discrimination on the independent test set.
roc_curve_plot <- test_predictions %>%
  group_by(model) %>%
  roc_curve(truth = event_dementia, .pred_dementia) %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity, color = model)) +
  geom_abline(linetype = "dashed", color = "gray60") +
  geom_path(linewidth = 1) +
  coord_equal() +
  labs(title = "Test Set ROC Curves", x = "1 - Specificity", y = "Sensitivity", color = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"), legend.position = "bottom")

ggsave(here("results", "figures", "ml_roc_curves.png"), roc_curve_plot, width = 7, height = 6, dpi = 300)

# Figure 3 shows precision-recall performance, which is useful for rare outcomes. The dashed line is the test-set event rate.
pr_curve_plot <- test_predictions %>%
  group_by(model) %>%
  pr_curve(truth = event_dementia, .pred_dementia) %>%
  ggplot(aes(x = recall, y = precision, color = model)) +
  geom_hline(yintercept = test_event_rate, linetype = "dashed", color = "gray60") +
  geom_path(linewidth = 1) +
  labs(
    title = "Test Set Precision-Recall Curves",
    subtitle = "Dashed line = test event rate (null model precision)",
    x = "Recall",
    y = "Precision",
    color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"), legend.position = "bottom")

ggsave(here("results", "figures", "ml_pr_curves.png"), pr_curve_plot, width = 7, height = 6, dpi = 300)

# Figure 4 shows the top random-forest predictors.
importance_plot <- variable_importance %>%
  slice_head(n = 10) %>%
  mutate(feature = fct_reorder(feature, importance)) %>%
  ggplot(aes(x = importance, y = feature)) +
  geom_col(fill = "#1f78b4") +
  labs(
    title = "Random Forest Permutation Variable Importance",
    subtitle = "Top 10 predictors",
    x = "Permutation importance",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

ggsave(here("results", "figures", "ml_rf_variable_importance.png"), importance_plot, width = 8, height = 6, dpi = 300)

print(ml_sample_summary)
print(cv_performance)
print(test_performance)

