# Predicting Translation from RNA Profiles

Undergraduate Research with supervisor Dr. Can Cenik, UT Austin


This project builds machine learning models to understand how RNA expression patterns predict translation efficiency (TE, how actively genes are converted into proteins).
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

### 5. MultiTask Lasso — Vector-Based Translation Model

This stage replaces the earlier cluster-mean approach with a vector-output model that jointly predicts translation efficiency across all ribosomal genes (RPS, RPL) simultaneously.

Input (X): RNA expression from non-ribosomal genes

Output (Y): Translation efficiency values for ribosomal genes

The model learns a shared sparse set of predictors - non-ribosomal genes whose expression patterns explain variation in ribosomal translation.

Joint L1 penalty encourages feature sharing across all outputs. Integration with CIGS Dataset: Each compound’s RNA profile from the CIGS HEK293T dataset is used as input. The trained model predicts the translational response vector (changes in RPS/RPL TE). Compounds are ranked by magnitude and direction of predicted translational impact.

### 6. Drug Perturbation Integration (CIGS Dataset)

Uses the Chemical-Induced Gene Signature (CIGS) dataset (https://cigs.iomicscloud.com/) to connect drug-induced RNA changes to predicted translational responses.

Each compound’s RNA signature (HEK293T) is input into the model, and the model predicts translation efficiency shifts across RPS/RPL genes

Compounds are ranked by magnitude and direction of predicted translational impact

Outcome: Provides a computational approach to link small molecules with translation-level changes 

## Technical Stack

Python / Colab

pandas, numpy, scikit-learn, matplotlib, networkx

Models: LinearRegression, RidgeCV, LassoCV, MultiTaskLassoCV, GraphicalLassoCV
