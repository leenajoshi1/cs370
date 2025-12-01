# Predicting Translation from RNA Profiles

Undergraduate research with supervisor Dr. Can Cenik (UT Austin)

This project documents a semester-long effort to answer a simple but motivating question: *Can we forecast ribosomal translation efficiency directly from transcriptome snapshots, and then use that model to reason about small-molecule perturbations?* The RNA to TE Prediction notebook walks through how the HEK293T RNA/TE matrices were curated, how interpretable regression models were built, and how drug signatures (CIGS) were projected into predicted ribosomal translation shifts. 

## Notebook Walkthrough

Each numbered section corresponds to a contiguous block of cells in the notebook. The text explains why the block exists, how it ties back to the research question, and what outputs you should see when things run correctly.

### Step 1 – Environment setup (Cells 1–2)
* Imports NumPy/Pandas/Matplotlib + sklearn, fixes the random seed, and defines reproducibility knobs (`RANDOM_SEED`, `N_SPLITS`, `USE_PCA`).
* Result: consistent logging/plotting behavior across reruns.

### Step 2 – Resolve data paths (Cells 3–4)
* Verifies the repository structure, locates the RNA/TE CSVs and CIGS workbooks under `data/`, and fails fast if anything is missing.
* Result: console prints the resolved `data/` directory plus explicit filenames for RNA, TE, counts, and metadata inputs.

### Step 3 – Reusable helpers (Cells 5–6)
* Cleans up CSV indices, infers ID columns, computes Pearson safely, maps transcripts to symbols, and defines utility masks (e.g., `rps_rpl_mask`).
* Result: helper functions (`drop_index_like`, `safe_pearsonr`, `clr`, etc.) that keep downstream cells concise.

### Step 4 – Align RNA and TE matrices (Cells 7–8)
* Reads `RNA_HEK293T.csv` and `TE_HEK293T.csv`, trims to shared samples/genes, renames ID columns, and logs resulting shapes.
* Result: harmonized RNA/TE tables with identical sample axes so the TE modeling never hits alignment errors.

### Step 5 – Build modeling matrices (Cells 9–10)
* Converts the aligned DataFrames into NumPy matrices (`X` for RNA predictors, `Y` for TE targets) and stores ordered gene metadata.
* Result: `X`, `Y`, `genes`, and `gene_symbols` ready for every model in later steps.

### Step 6 – Cluster-level Lasso sanity checks (Cells 11–13)
* Focuses on RPS/RPL clusters, runs cross-validated `LassoCV`, refits the best alpha on full data, and visualizes the highest-weight RNA regulators.
* Result highlights:
	* Printed metrics (`r`, MAE, RMSE) for RPS/RPL predictions.
	* Median alpha + non-zero feature counts to gauge sparsity.
	* Tables/plots of RNA drivers that most strongly shift each cluster.

### Step 7 – MultiTask Lasso for ribosomal TE (Cells 14–18)
* Separates non-ribosomal predictors from ribosomal targets, handles missing values via CLR-style filling, and trains a shared `MultiTaskLasso` (`mtl_model`).
* Result highlights:
	* Train RMSE + chosen alpha so you can confirm optimization and hyperparameters.
	* Model summary with sparsity, per-target RMSE, and a ranked list of multi-target regulators (the earlier noisy coefficient heatmap remains retired).
	* Mean/median RMSE across outputs for a quick model-health metric between retrains.

### Step 8 – CIGS ingestion and prediction setup (Cells 19–21)
* Rebuilds XLSX parsing with standard-library `zipfile` + streaming `ElementTree.iterparse`, strips note rows, and infers plate:well/sample IDs so metadata merges survive header drift.
* Filters metadata to HEK293T / 24 h / valid doses, applies the same CLR transform used in training, aligns features to `feature_names`, and feeds compounds through `mtl_model` to obtain `pred_TE` (compounds × ribosomal genes).
* Result highlights: console logs raw/filtered shapes, overlapping feature counts, predicted TE dimensions, and a top-10 compound table sorted by L2 magnitude alongside dataset counts (compounds, target genes, predictors).

### Step 9 – Visualize top compound responses (Cell 22)
* Confirms `pred_TE` and `ranked` exist, then renders a heatmap for the highest-ranking compounds.
* Result: figure showing the top translational perturbations across ribosomal genes to quickly spot expected positives/negatives.

### Step 10 – Rank compounds with combined scores (Cell 23)
* Collapses each predicted TE vector into z-scored mean and magnitude components, then builds `drug_summary` with a combined score.
* Result: sorted DataFrame that highlights the most extreme translational modulators whether directional or broad-spectrum.

### Step 11 – Baseline per-gene regression (Cell 24)
* Fits simple `LinearRegression` models (TE ~ RNA for the same gene) to provide a conservative baseline.
* Result: median R²/RMSE printouts so you know the floor before multitask sharing.

### Step 12 – Ridge baselines across ribosomal genes (Cells 25–26)
* Uses `RidgeCV` per ribosomal gene plus aggregate summaries to benchmark against the multitask approach.
* Result: counts of fitted ridge models and median/mean error metrics for apples-to-apples comparisons.

### Step 13 – Compound score diagnostics (Cell 27)
* Builds scatter + histogram views of the standardized mean-effect versus magnitude components to see which compounds are directional specialists versus broad-spectrum modulators.
* Result: bivariate plot with combined-score color/size encoding, a distribution view for all compounds, and a printed top-10 list for downstream reporting.

### Step 14 – Ribosomal difficulty profile (Cell 28)
* Uses the per-target RMSE array directly instead of a coarse heatmap to flag ribosomal genes that remain challenging for the multitask model.
* Result: paired bar charts for the hardest vs. easiest dozen targets plus summary tables that guide where to invest more features or experiments.

## Result Snapshot
* **Step 6 – Cluster Lasso** – Validates that ribosomal clusters carry interpretable RNA drivers (Cells 11–13).
* **Step 7 – MultiTask Lasso** – Captures the joint translation landscape; pay special attention to the sparsity and per-target error readouts (Cells 14–18).
* **Step 8 – CIGS projection** – Connects drug signatures to predicted translation shifts; the top-10 table matches what was shared during lab meetings (Cells 19–21).
* **Steps 9–10 – Compound figures** – Heatmap + combined z-score ranking make it easy to pitch candidate compounds (Cells 22–23).
* **Steps 11–12 – Baselines** – Linear + ridge models provide a sanity floor so we know the multitask model truly helps (Cells 24–26).
* **Step 13 – Diagnostics** – Scatter + histogram views of the combined scores make it easy to discuss whether a compound is directionally biased or a broad-spectrum ribosome modulator (Cell 27).
* **Step 14 – Difficulty readout** – Hard/easy ribosomal targets are summarized with bar charts instead of the earlier, inaccurate heatmap (Cell 28).

### Current Highlights (Nov 2025 run)
These numbers come from the most recent end-to-end execution. Treat them as targets when you rerun or extend the notebook.
* MultiTask Lasso train RMSE: ~0.19 log2 units with alpha 0.01 and ~83% sparsity.
* Per-target RMSE median: ~0.17; mean: ~0.21 across 87 ribosomal outputs.
* Top predicted translational amplifiers in latest CIGS pass: homoharringtonine, cycloheximide, and emetine (expected positive controls). Several HDAC inhibitors show strong negative mean effect, suggesting a translational repression axis that deserves wet-lab validation.
* Compound score diagnostics currently show ~30% of compounds dominated by directional shifts (|z_mean| > z_l2); the rest act as broad-spectrum ribosome modulators.
* Ribosomal difficulty profiling shows a median RMSE gap of ~0.11 between the hardest and easiest dozen targets, which now serves as the progress-tracking metric instead of the retired coefficient heatmap.
* Ridge baseline median R² ≈ 0.42 vs. the multitask sparsity/readouts, reinforcing why we invest in the multitask approach.


## Future Research Directions (Drug-Focused)
These examples emphasize translating model predictions into actionable drug insights—both computational follow-ups and suggested experiments.
1. **Compound clustering by TE signature** – Use the predicted `pred_TE` matrix to cluster drugs by shared translational fingerprints, then relate each cluster to known mechanisms (HDAC inhibition, ribosome poison, etc.).
2. **Repurposing exploration** – Cross-reference top translational modulators with FDA-approved indications; flag candidates whose predicted TE effects suggest new uses (e.g., HDAC inhibitors as translation dampers in hyperactive ribosome diseases).

## Technical Stack
* Python + Jupyter/VS Code
* pandas, numpy, matplotlib, seaborn
* scikit-learn models: `LinearRegression`, `RidgeCV`, `LassoCV`, `MultiTaskLasso`, plus custom CLR + XLSX utilities
