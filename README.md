# WGS-Analysis

This repository contains R scripts for risk score analysis and annotation clustering in An *et al.* (2018). 
The raw data are not included in this repository. 
For data generation, please refer to [https://github.com/sanderslab/WGS-pipeline](https://github.com/sanderslab/WGS-pipeline).

+ Risk score analysis
  + Compute the predictive R2 of rare annotations: ``risk_score/predictive_r2.R``
  + Compute the p-value: ``risk_score/significance_r2.R``
+ Annotation clustering
  + Spectral clustering of annotation categories: ``clustering/clustering.R``
