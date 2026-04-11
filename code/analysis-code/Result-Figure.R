############################################
# Result-Figure.R
# Figure 1: Trend plot across landmark waves
# Panel A = Continuous BMIV
# Panel B = Categorical BMIV
############################################

library(tidyverse)
library(readr)
library(here)
library(ggplot2)
library(patchwork)


#  Read model results
cox_results <- read_csv(
  here("results", "tables", "cox_w3_w9_all_models_combined.csv"),
  show_col_types = FALSE
)

#  Keep only adjusted/full-model results
waves_to_plot <- 3:9

# Continuous adjusted BMIV results
continuous_data <- cox_results %>%
  filter(
    model == "Adjusted",
    term %in% paste0("BMI_CV_W", waves_to_plot)
  ) %>%
  mutate(
    wave_num = as.numeric(str_extract(term, "[3-9]$")),
    landmark_year = 1994 + 2 * wave_num,
    group = "Continuous"
  ) %>%
  select(wave_num, landmark_year, group, estimate, conf.low, conf.high)

# Categorical adjusted BMIV results: Mid
mid_data <- cox_results %>%
  filter(
    model == "Categorical_Adjusted",
    term %in% paste0("BMIV_CAT_W", waves_to_plot, "Medium variability")
  ) %>%
  mutate(
    wave_num = as.numeric(str_extract(term, "[3-9]")),
    landmark_year = 1994 + 2 * wave_num,
    group = "Mid"
  ) %>%
  select(wave_num, landmark_year, group, estimate, conf.low, conf.high)

# Categorical adjusted BMIV results: High
high_data <- cox_results %>%
  filter(
    model == "Categorical_Adjusted",
    term %in% paste0("BMIV_CAT_W", waves_to_plot, "High variability")
  ) %>%
  mutate(
    wave_num = as.numeric(str_extract(term, "[3-9]")),
    landmark_year = 1994 + 2 * wave_num,
    group = "High"
  ) %>%
  select(wave_num, landmark_year, group, estimate, conf.low, conf.high)

# Low is the reference group
low_data <- tibble(
  wave_num = waves_to_plot,
  landmark_year = 1994 + 2 * waves_to_plot,
  group = "Low",
  estimate = 1,
  conf.low = 1,
  conf.high = 1
)

# Combine categorical groups
categorical_data <- bind_rows(
  low_data,
  mid_data,
  high_data
) %>%
  mutate(
    group = factor(group, levels = c("Low", "Mid", "High"))
  )

# Add horizontal offsets for categorical plot
# This keeps Low, Mid, and High separated within the same year.

offset_map <- tibble(
  group = factor(c("Low", "Mid", "High"), levels = c("Low", "Mid", "High")),
  offset = c(-0.25, 0, 0.25)
)

categorical_data <- categorical_data %>%
  left_join(offset_map, by = "group") %>%
  mutate(
    x_pos = landmark_year + offset
  )

# Panel A: Continuous BMIV
panel_a <- ggplot(
  continuous_data,
  aes(x = landmark_year, y = estimate)
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "gray50",
    linewidth = 0.6
  ) +
  geom_errorbar(
    aes(ymin = conf.low, ymax = conf.high),
    width = 0.25,
    linewidth = 0.7
  ) +
  geom_point(size = 2.8) +
  scale_x_continuous(
    breaks = c(2000, 2002, 2004, 2006, 2008, 2010, 2012)
  ) +
  labs(
    title = "Panel A. Continuous BMIV",
    x = "Landmark year",
    y = "Hazard ratio (HR)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

# Panel B: Categorical BMIV
panel_b <- ggplot(
  categorical_data,
  aes(x = x_pos, y = estimate, color = group)
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "gray50",
    linewidth = 0.6
  ) +
  geom_errorbar(
    aes(ymin = conf.low, ymax = conf.high),
    width = 0.08,
    linewidth = 0.7
  ) +
  geom_point(size = 2.8) +
  scale_x_continuous(
    breaks = c(2000, 2002, 2004, 2006, 2008, 2010, 2012),
    labels = c(2000, 2002, 2004, 2006, 2008, 2010, 2012)
  ) +
  labs(
    title = "Panel B. Categorical BMIV",
    x = "Landmark year",
    y = "Hazard ratio (HR)",
    color = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "top"
  )

# Combine the two panels
figure1_trend <- panel_a / panel_b 

# Print in Viewer
figure1_trend

# Save outputs
ggsave(
  filename = here("results", "figures", "figure1_trend_two_panels.png"),
  plot = figure1_trend,
  width = 9,
  height = 10,
  dpi = 300
)
