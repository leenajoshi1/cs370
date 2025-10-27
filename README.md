Branch for Graphical Lasso Network Inference

Models the gene–gene relationships in translation.

Graphical Lasso learns a sparse precision matrix (inverse covariance), identifying conditional dependencies between genes.

Produces:

Partial correlation network (partial_correlations.csv)

Thresholded adjacency matrix (adjacency_thresholded.csv)

Visual graph showing major translational modules

## Technical Stack

Python / Colab

pandas, numpy, scikit-learn, matplotlib, networkx

Models: LinearRegression, RidgeCV, LassoCV, GraphicalLassoCV
