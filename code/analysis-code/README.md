# Analysis Code

This folder contains the scripts used to run the main statistical analyses and generate manuscript-ready result tables and figures.

Run order:

1. `statistical-analysis.R`: loads the processed HRS dataset, builds W3-W9 landmark datasets, runs Cox proportional hazards models, and saves model summaries.
2. `Supplementary ML .R`: runs supplementary Random Survival Forest analyses for W3 and W6 landmark samples.
3. `Result-Table.R`: creates formatted baseline and Cox regression tables.
4. `Result-Figure.R`: creates the landmark trend figure.

`Result-Function.R` contains helper functions sourced by `statistical-analysis.R`.
