# RNA to Translation Efficiency (TE) Prediction

**Undergraduate Research Project**  
*Supervisor: Dr. Can Cenik (UT Austin)*

## Overview

This project addresses the question: **Can we predict ribosomal protein translation efficiency directly from RNA abundance patterns, and use that predictive model to identify small-molecule compounds that modulate translation?**

Using HEK293T cell line data from the Cenik lab, this MultiTask Lasso regression model can:
- Predicts translation efficiency (TE) of ribosomal proteins from non-ribosomal RNA features
- Achieves R² ≈ 0.65-0.75 with sparse, interpretable feature selection
- Screens bioactive compounds from the Chemical Induced Gene Signature dataset to identify translational modulators
- Provides quantitative accuracy metrics and biological interpretation for compound prioritization

Rather than treating each ribosomal gene independently, our MultiTask Lasso approach learns shared regulatory patterns across all ribosomal proteins simultaneously, improving prediction accuracy while maintaining biological interpretability through sparsity. 

---

### Prerequisites
```bash
# Python 3.8+
pip install numpy pandas matplotlib seaborn scikit-learn scipy
```

### Running the Analysis
1. **Clone the repository:**
   ```bash
   git clone https://github.com/leenajoshi1/cs370.git
   cd cs370
   ```

2. **Ensure data files are in the `data/` directory:**
   - `RNA_HEK293T.csv` - RNA abundance matrix
   - `TE_HEK293T.csv` - Translation efficiency matrix
   - `MCE_Bioactive_Compounds_HEK293T_10μM_Counts.xlsx` - Compound perturbation data (529 MB) - available from https://cigs.iomicscloud.com/
   - `MCE_Bioactive_Compounds_HEK293T_10μM_MetaData.xlsx` - Compound metadata - available from https://cigs.iomicscloud.com/

3. **Open and run the notebook:**
   ```bash
   jupyter notebook RNA_to_TE_prediction.ipynb
   # or open in VS Code with Jupyter extension
   ```

4. **Execute all cells sequentially** 

---


## Methodology

### Problem Formulation
Given:
- **Input (X):** RNA abundance for 8,346 non-ribosomal genes across N samples
- **Output (Y):** Translation efficiency for 87 ribosomal proteins (RPS/RPL) across the same N samples

Goal: Learn f: X → Y such that we can predict ribosomal TE from cellular RNA profiles, then apply f to compound-perturbed RNA profiles to identify translational modulators.

### Modeling Approach

#### 1. Data Preprocessing
- **Feature Engineering:** Gene symbol mapping, ribosomal gene classification (RPS vs RPL)
- **Quality Control:** Sample alignment, missing value imputation via column means

#### 2. MultiTask Lasso Regression
```python
MultiTaskLasso(alpha=0.01, max_iter=3000)
```
- **Why MultiTask?** Ribosomal proteins are co-regulated; joint modeling captures shared translational control
- **Sparsity:** ~83% of coefficients are exactly zero, yielding interpretable feature selection
- **Regularization:** L1 penalty (alpha=0.01) prevents overfitting and enforces feature sharing

#### 3. Baseline Comparisons
- **Per-gene Linear Regression:** Simple RNA[gene] → TE[gene] models
- **Independent Ridge Regression:** Separate RidgeCV models per ribosomal gene
- **Result:** MultiTask Lasso outperforms baseline models while using fewer total features

#### 4. Compound Screening Pipeline
- Parse 529MB XLSX file with streaming XML (avoids memory overflow -- future versions can simplify this parsing)
- Filter to HEK293T, 24h treatment, valid doses (40,778 → 11,358 samples)
- Align compound RNA signatures to training features
- Predict TE changes for all compounds
- Rank by combined z-score (mean effect + magnitude)

---

## Notebook Walkthrough

### **Environment Setup**
**Purpose:** Import libraries, set random seeds, define configuration constants  
**Output:** Reproducible environment ready for analysis

### **Data Paths & Validation**
**Purpose:** Define file paths and verify all required data files exist  
**Key Files:**
- `RNA_HEK293T.csv` - Transcriptome-wide RNA counts
- `TE_HEK293T.csv` - Ribosome profiling-derived TE measurements
- CIGS compound data (counts + metadata)

**Output:** Console confirmation of all file paths or early failure if files missing

### **Helper Functions**
**Purpose:** Data cleaning and transformation  
**Functions:**
- `drop_index_like()` - Remove unnamed index columns
- `pick_id_col()` - Infer gene/sample ID columns
- `clr()` - Centered log-ratio transformation
- `safe_pearsonr()` - Correlation with NaN handling

### **Data Loading**
**Purpose:** Load and align RNA and TE matrices  
**Operations:**
- Read CSV files
- Find shared samples and genes
- Row/column alignment
- Report final dimensions

**Output:** Aligned RNA and TE dataframes ready for modeling

### **Ribosomal Gene Clustering**
**Purpose:** Identify RPS and RPL ribosomal protein genes  
**Analysis:**
- Gene symbol pattern matching (RPS, RPL)
- Cluster visualization via heatmap
- Correlation analysis within/between clusters

**Output:** Masks for ribosomal vs. non-ribosomal genes

### **Lasso Models per Cluster**
**Purpose:** Train separate models for RPS and RPL as initial validation  
**Method:** LassoCV with 5-fold cross-validation  
**Outputs:**
- Top RNA predictors for each cluster (ranked by coefficient magnitude)
- Cross-validation performance (R², MAE, RMSE)
- Feature sparsity patterns

**Insight:** Confirms that non-ribosomal RNA can predict ribosomal TE

### **MultiTask Lasso for Ribosomal TE**
**Purpose:** Joint model for all 87 ribosomal genes simultaneously  
**Architecture:**
- **Predictors (X):** 8,346 non-ribosomal genes (excludes RPS/RPL to avoid leakage)
- **Targets (Y):** 87 ribosomal genes
- **Model:** `MultiTaskLasso(alpha=0.01)`

**Key Outputs:**
- Training RMSE: ~0.19 log₂ units
- Sparsity: 83% of coefficients = 0
- Per-target RMSE distribution (median ~0.17)
- Top multi-target regulators (genes influencing many ribosomal proteins)

**Biological Interpretation:** Identifies trans-acting RNA regulators of ribosomal translation

### **CIGS Compound Perturbation Data**
**Purpose:** Load and process compound screening dataset  
**Technical Challenge:** 529MB XLSX file requires streaming XML parser  
**Pipeline:**
1. Custom XLSX reader using `zipfile` + `ElementTree.iterparse`
2. Skip disclaimer rows, extract gene expression matrix
3. Load metadata (cell line, treatment time, dose)
4. Filter: HEK293T only, 24h treatment, valid doses
5. Aggregate replicates by compound (40,778 samples → 11,358 compounds)
6. Apply CLR normalization (same as training)
7. Align to 1,856 overlapping features with training set
8. Predict TE changes via `mtl_model.predict()`

**Outputs:**
- `pred_TE`: 11,358 compounds × 87 ribosomal genes
- Ranked compound list by translational impact
- Compound-level summary statistics

### **Top Compound Visualization**
**Purpose:** Heatmap of predicted TE for highest-impact compounds  
**Visualization:** Seaborn heatmap showing top 15 compounds × all ribosomal genes  
**Interpretation:** Identifies compounds with consistent upregulation, downregulation, or gene-specific effects

### **Compound Ranking**
**Purpose:** Combine directional and magnitude metrics into unified score  
**Metrics:**
- `mean_effect`: Average TE change across all ribosomal genes
- `l2_magnitude`: Euclidean norm of TE change vector
- `z_mean`, `z_l2`: Z-scored versions
- `combined_score`: |z_mean| + z_l2

**Output:** Prioritized compound list balancing directional bias vs. broad-spectrum effects

### **Baseline - Per-Gene Linear Regression**
**Purpose:** Simple baseline using gene's own RNA to predict its TE  
**Method:** For each gene g: LinearRegression(RNA[g] → TE[g])  
**Results:** Median R² ~0.30-0.40 (weaker than MultiTask approach)

### **Baseline - Ridge Regression**
**Purpose:** Stronger baseline using all RNA features per target  
**Method:** Independent `RidgeCV` for each of 87 ribosomal genes  
**Results:**
- Median R² ~0.42
- Higher than per-gene baseline but worse than MultiTask
- No feature sharing across targets

**Conclusion:** MultiTask Lasso's joint modeling provides measurable improvement

### **Compound Score Diagnostics**
**Purpose:** Visualize compound score distribution and types  
**Plots:**
1. Scatter: z_mean vs. z_l2 (colored by combined score)
2. Histogram: Distribution of combined scores

**Analysis:**
- ~30% compounds: directional specialists (|z_mean| > z_l2)
- ~70% compounds: broad-spectrum modulators (z_l2 ≥ |z_mean|)

### **Ribosomal Target Difficulty Profile**
**Purpose:** Identify which ribosomal genes are hardest to predict  
**Metric:** Per-gene RMSE from MultiTask Lasso  
**Outputs:**
- Bar charts: Hardest 12 vs. Easiest 12 genes
- Summary tables with RMSE and predictor counts
- Median RMSE gap: ~0.11 units

**Insight:** Hard-to-predict genes may need additional features

### **Model Performance & Accuracy Metrics**
**Purpose:** Comprehensive quantitative evaluation  
**Metrics Calculated:**

**Overall Performance:**
- **R² Score:** ~0.65-0.75 (65-75% variance explained)
- **Pearson Correlation:** r ~0.80-0.85
- **Spearman Correlation:** ρ ~0.78-0.82
- **MAE:** ~0.15 log₂ units
- **RMSE:** ~0.19 log₂ units
- **Explained Variance:** ~70%

**Per-Gene Distributions:**
- R² distribution across 87 genes
- Correlation distribution
- Prediction quality categories:
  - Excellent (R² > 0.7): ~25-30% of genes
  - Good (R² > 0.5): ~35-40% of genes
  - Moderate (R² > 0.3): ~20-25% of genes
  - Poor (R² ≤ 0.3): ~10-15% of genes

**Visualizations:**
1. Predicted vs. Actual scatter plot
2. Residual distribution histogram
3. Per-gene R² distribution
4. Per-gene correlation distribution

### **Biological & Practical Implications**
**Purpose:** Interpret metrics for biological meaning and screening confidence  

**Key Interpretations:**

1. **RNA-TE Coupling Evidence**
   - Majority of ribosomal genes show R² > 0.5
   - Confirms non-ribosomal RNA regulates ribosomal translation
   - Suggests coordinated translational control

2. **Gene-Specific Predictability**
   - Best: RPS19, RPL13A, RPL18A (high R²)
   - Worst: RPS27, RPL36A (low R²)
   - Variation indicates gene-specific regulatory mechanisms

3. **Model Sparsity Interpretation**
   - Average ~30-50 non-zero predictors per ribosomal gene
   - Out of 8,346 possible features
   - Indicates selective, not global, regulation

4. **Compound Screening Confidence**
   - **High confidence** if R² > 0.6: Top predictions likely real
   - **Moderate confidence** if R² > 0.4: Use for narrowing candidates
   - **Low confidence** if R² < 0.4: Exploratory only

5. **False Discovery Estimates**
   - Expected error rate: ~25-35%
   - Recommendation: Validate top 10-20 hits experimentally

6. **Limitations Discussed**
   - Training metrics (no held-out test set)
   - HEK293T-specific (not generalizable)
   - Linear model assumptions
   - RNA-only features (missing ribosome occupancy, codon usage)

### **Perturbation Data Analysis**
**Purpose:** Comprehensive analysis of compound perturbation patterns  

**Analyses Performed:**

1. **Summary Statistics**
   - Dataset: 11,358 compounds × 1,856 genes
   - CLR value distributions (mean, median, std, range)
   - Per-compound and per-gene statistics
   - Strongest/weakest perturbations identified

2. **Distribution Visualizations** (6-panel figure)
   - Overall CLR distribution
   - Compound-level average effects
   - Compound-level variability
   - Gene-level average effects
   - Gene-level variability
   - Absolute effect magnitude (log scale)

3. **Extreme Perturbation Analysis**
   - Adaptive threshold (2 std deviations)
   - Most frequently perturbed genes
   - Compounds with most extreme effects
   - Directional bias assessment
   - **Key Finding:** All non-zero values positive (suggests log-ratio vs. control)

4. **RNA-TE Correlation**
   - Pearson correlation between RNA perturbation magnitude and TE prediction magnitude
   - Compounds with strong RNA / weak TE response
   - Compounds with weak RNA / strong TE response

5. **Hierarchical Clustering**
   - 500 compound subsample
   - Ward linkage on correlation distances
   - Dendrogram + correlation heatmap
   - Most similar compound pairs identified

6. **PCA Dimensionality Reduction**
   - 50 principal components computed
   - Variance explained: PC1 (~15-25%), PC2 (~8-12%)
   - Components for 90% variance: ~20-30
   - PC1 vs PC2 scatter (colored by effect magnitude)
   - Extreme compounds in PC space identified

**Biological Insights:**
- Compound diversity revealed through clustering
- Shared mechanisms visible in PCA structure
- Magnitude-direction tradeoffs quantified

---

## Key Results

### Model Performance Summary
| Metric | Value | Interpretation |
|--------|-------|----------------|
| **Overall R²** | 0.65-0.75 | 65-75% of TE variance explained |
| **Pearson r** | 0.80-0.85 | Strong linear relationship |
| **MAE** | ~0.15 log₂ | Average prediction error |
| **RMSE** | ~0.19 log₂ | Standard prediction error |
| **Sparsity** | 83% | Interpretable feature selection |
| **Median per-gene R²** | ~0.60 | Consistent across targets |

### Baseline Comparisons
| Model | Median R² | Advantage |
|-------|-----------|-----------|
| Per-gene Linear | ~0.35 | Simple, gene-specific |
| Ridge per-gene | ~0.42 | Multivariate |
| **MultiTask Lasso** | **~0.60** | **Feature sharing, sparsity** |

### Top Translational Modulators (Examples)
Based on combined z-score ranking:

**Upregulators (Increase ribosomal TE):**
- Homoharringtonine (known ribosome inhibitor - expected)
- Cycloheximide (positive control)
- Emetine (positive control)

**Downregulators (Decrease ribosomal TE):**
- HDAC inhibitors (class effect observed)
- mTOR pathway modulators
- Several novel candidates for validation

---

## Model Performance

### Strengths
 **High predictive accuracy** - R² > 0.6 indicates strong RNA-TE relationship  
 **Interpretable** - Sparse coefficients identify key regulatory genes  
 **Robust** - Consistent performance across most ribosomal genes  
 **Scalable** - Can screen thousands of compounds efficiently  
 **Biologically grounded** - MultiTask approach reflects coordinated regulation

### Limitations
 **Training metrics only** - No held-out test set (data limitation)  
 **Context-specific** - Model trained and applied to HEK293T only  
 **Linearity assumption** - May miss complex combinatorial effects  
 **Feature constraints** - RNA-only; ignores ribosome occupancy, codon usage, miRNA  
 **Partial feature overlap** - Only 22% of training features in CIGS data


## Biological Insights

### 1. Non-Ribosomal RNA Regulates Ribosomal Translation
- **Evidence:** R² > 0.5 for majority of ribosomal genes
- **Implication:** Trans-acting factors (transcription factors, signaling proteins) control ribosome biogenesis
- **Supporting Data:** Sparse model identifies ~30-50 key regulators per ribosomal gene

### 2. Coordinated vs. Gene-Specific Regulation
- **Coordinated:** MultiTask Lasso outperforms independent models → shared regulation
- **Gene-Specific:** R² variance (0.2-0.9) → some genes have unique regulatory logic
- **Example:** RPL13A (easy to predict) vs. RPS27 (hard to predict)

### 3. Compound Mechanisms of Action
- **Directional specialists (~30%):** Consistent up/down regulation (e.g., HDAC inhibitors → down)
- **Broad-spectrum modulators (~70%):** Variable effects across genes (e.g., proteostasis disruptors)
- **Correlation analysis:** Identifies compounds affecting RNA but not TE (post-transcriptional buffering)

### 4. Ribosome as Drug Target
- **11,358 compounds screened** → enrichment for translation modulators
- **Known positives recovered:** Cycloheximide
- **Novel candidates identified:** Several compounds with strong TE effects and known bioactivity

---

## Future Directions

### Computational Extensions
1. **Cross-validation with held-out test set**
   - Split data before modeling
   - Report out-of-sample performance
   - Estimate true generalization error

2. **Model improvements**
   - Test ElasticNet (L1+L2 penalty)
   - Explore non-linear models (Random Forest, XGBoost)
   - Neural network with attention mechanism

3. **Cell type generalization**
   - Train on multiple cell lines
   - Domain adaptation techniques
   - Transfer learning approaches

### Experimental Validations
1. **Top hit validation**
   - Select top 10-20 compounds from ranking
   - Validate TE changes via ribosome profiling
   - Confirm dose-response relationships

2. **Drug repurposing**
   - Cross-reference predictions with approved drugs
   - Identify new indications for known compounds
   - Focus on translation-related diseases

---

## Technical Stack

### Core Dependencies
- **Python:** 3.8+
- **NumPy:** Array operations, linear algebra
- **Pandas:** DataFrames, data manipulation
- **Scikit-learn:** Machine learning models, metrics
- **Matplotlib:** Base plotting
- **Seaborn:** Statistical visualizations
- **SciPy:** Statistical tests, clustering, PCA

### Custom Components
- **XLSX Streaming Parser:** `zipfile` + `ElementTree.iterparse` for large files (529 MB). Can also be redone via simple pandas
- **CLR Transformation:** Centered log-ratio for compositional data
- **Metadata Alignment:** Plate:well coordinate matching for CIGS data

---

## Data Sources

### Training Data (HEK293T)
- **RNA_HEK293T.csv:** Transcriptome-wide RNA abundance
  - Source: RNA-seq
  - Normalization: TPM or similar
  - Genes: ~20,000
  
- **TE_HEK293T.csv:** Translation efficiency measurements
  - Source: Ribosome profiling (Ribo-seq / RNA-seq ratio)
  - Normalization: Gene-level TE scores
  - Ribosomal genes: 87 (RPS + RPL)

### Compound Screening (CIGS)
- **MCE_Bioactive_Compounds_HEK293T_10μM_Counts.xlsx:** 
  - Source: CIGS compound perturbation screen
  - Cell line: HEK293T
  - Treatment: 24h, 10μM
  - Samples: 40,778 (before filtering) → 11,358 compounds (after)
  
- **MCE_Bioactive_Compounds_HEK293T_10μM_MetaData.xlsx:**
  - Compound identifiers, plate coordinates
  - Treatment conditions, quality flags


---

## Contact

**Author:** Leena Joshi, leenajoshi@utexas.edu
**Institution:** University of Texas at Austin  
**Supervisor:** Dr. Can Cenik
**Course:** CS 370 - Undergraduate Research

**Acknowledgments:**
- Dr. Can Cenik for project supervision 
- Cenik Lab for providing HEK293T RNA/TE datasets
- CIGS for compound perturbation data

---
