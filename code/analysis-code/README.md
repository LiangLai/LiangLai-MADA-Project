# Analysis Code

This folder contains the scripts used to run the main Cox regression analysis, predictive ML model comparison, and manuscript-ready result tables and figures.

Run order:

1. `statistical-analysis.R`: loads the processed HRS dataset, builds W3-W9 landmark datasets, runs Cox proportional hazards models, and saves model summaries.
2. `Result-Table.R`: creates formatted baseline and Cox regression tables.
3. `Result-Figure.R`: creates the landmark trend figure.
4. `ml-model-comparison.R`: runs the W3 10-year predictive modeling analysis with logistic regression, LASSO logistic regression, random forest, and a null baseline. It saves cross-validation metrics, test-set metrics, best tuning parameters, ROC and precision-recall curves, and random forest variable importance.

`Result-Function.R` contains helper functions sourced by `statistical-analysis.R`.

`Supplementary ML .R` is an older exploratory Random Survival Forest script and is not required to reproduce the final manuscript.
