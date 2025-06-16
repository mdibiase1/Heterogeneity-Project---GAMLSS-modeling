# GAMLSS-Based Normative Modeling for Biological Phenotypes

This repository provides code and documentation for fitting GAMLSS (Generalized Additive Models for Location, Scale and Shape) models to biological phenotypes using the Sinh-Arcsinh (SHASH) distribution. This method supports robust modeling of non-Gaussian data and computes individualized z-scores for use in downstream normative modeling and stratification analyses.

## Overview

The pipeline includes:
  
- Model training using the SHASH distribution on healthy control data
- Z-score calculation for all subjects based on GAMLSS-derived distributional parameters
- Output writing to .mat files for MATLAB compatibility

## Application

The following R packages are required:
  
- gamlss
- R.matlab
- pracma

## Input 
- .mat file: Preprocessed data containing phenotype, age, sex, and diagnosis (a matrix of subjects on rows and columns = groups (healthy control group = column one and coluns 2-x are diagnostic groups) and 1 given to subject with a diagnosis)

## Output

- Output .mat file: Contains z-scores, centile estimates and fitted parameters

## Notes

- Z-scores computed reflect individual deviation from normative predictions after adjusting for age and sex using SHASH distribution parameters.

- The model is trained using healthy controls only

- The SHASH distribution enables modeling of skewness and kurtosis, improving fit to real-world phenotype distributions


## Citation

If using this method, please cite the original reference for SHASH:
  
  Jones, M. C., & Pewsey, A. (2009). Sinh-arcsinh distributions. Biometrika, 96(4), 761–780.
https://doi.org/10.1093/biomet/asp053
