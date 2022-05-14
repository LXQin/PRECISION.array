# PRECISION.array v1.0.0

We propose a new R package PRECISION.array that provides a pipeline for resampling-based data simulation, data normalization, and classifier development for microRNA microarrays. The package offers two study designs for simulating data via a process dubbed ‘virtual rehybridization’ (confounding design and stratified design), three methods for normalizing the training data (median normalization, quantile normalization, and variance stabilizing normalization), seven methods for normalizing the test data (median normalization, frozen median normalization, quantile normalization, frozen quantile normalization, pooled quantile normalization, variance stabilizing normalization, and frozen variance stabilizing normalization), and seven methods for building a sample classifier (LASSO, PAM, ClaNC, DLDA, kNN, SVM, and Random Forest). This package can be used to study the performance of various study designs for microarray data generation and various methods for data normalization and sample classification, in connection with each other. And the full pipelines of normalization and classification are provided by the wrapper functions precision.simulate(), precision.simulate.class(), and precision.simulate.multiclass().

The full package documentation can be found [here](https://lxqin.github.io/PRECISION.array/).

The package can be installed in R:
```ruby
devtools::install_github("LXQin/PRECISION.array")
```

### REFERENCE

[1]: Qin LX, Huang HC, Begg CB. Cautionary note on cross validation in molecular classification. Journal of Clinical Oncology. 2016, http://ascopubs.org/doi/abs/10.1200/JCO.2016.68.1031.

[2]: Huang HC, Qin LX. PRECISION: an R package for assessing study design and data normalization for molecular classification studies using a unique pair of microarray datasets (2016). GitHub repository, https://github.com/LXQin/PRECISION.array.

[3]: Wu YL, Huang HC, Qin LX. Making External Validation Valid for Molecular Classifier Development. JCO Precision Oncology. 2021, 
https://ascopubs.org/doi/full/10.1200/PO.21.00103
