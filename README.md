# Predicting Translation from RNA Profiles

Undergraduate Research with supervisor Dr. Can Cenik, UT Austin


This project builds machine learning models to understand how RNA expression patterns predict translation efficiency (TE) — how actively genes are converted into proteins — in human HEK293T cells.
The ultimate goal is to map how drugs or perturbations shift translation activity across the genome.

## Project Pipeline
### 1. Data Loading and Alignment

Input files:

RNA_HEK293T.csv — normalized RNA expression (samples × genes)

TE_HEK293T.csv — normalized translation efficiency (samples × genes)

The notebook aligns shared samples and gene IDs, producing consistent samples × genes matrices for modeling.

### 2. Baseline Model: Per-Gene Linear Regression

Predicts each gene’s TE from its own RNA level independently:
TE(gene) ~ RNA(gene)

Evaluated with:
Pearson correlation (r)
Mean Absolute Error (MAE)
Root Mean Squared Error (RMSE)

Result: Most genes show weak or negative correlation, as TE is defined relative to RNA.

### 3. Multivariate Model: Ridge Regression

Moves beyond per-gene models by using the full RNA profile to predict translation of each gene or cluster.

Regularization (ridge penalty) prevents overfitting by shrinking noisy coefficients.

Evaluated via 5-fold cross-validation.

### 4. Cluster-Level Lasso Regression (RPS/RPL genes)

Focused on ribosomal protein genes (RPS, RPL), which strongly influence translation.

Builds one model per cluster:

Predicts mean TE across all RPS or RPL genes.

Uses the full RNA expression profile as input.

Lasso regression selects only the most informative genes, improving interpretability.

Outputs include:

Model weights per gene

Top predictors and their direction of effect

Cross-validated error metrics



## Technical Stack

Python / Colab

pandas, numpy, scikit-learn, matplotlib, networkx

Models: LinearRegression, RidgeCV, LassoCV, GraphicalLassoCV
