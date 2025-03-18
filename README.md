# README

## Introduction

This repository contains the code and data for inferring RSV infections from serological data using the RSV PreF protein.

## Installation

To install the `serojump` package, use the following command in R:

```{r}

devtools::install_github("seroanalytics/serojump")

```

## R Files Description

### 1. `main.R`

Load the relevant libraries for analysis.

### 2. `clean_fudan.R`

This script preprocesses the raw serological data (not included in this repo, please email me for the data), including cleaning, normalization, and transformation steps necessary for analysis with `serojump`.

### 3. `model_base.R/model_base_hier_2.R/model_base_hier_3.R`

This script defined the model structure required for input into `serojump`. One is `model_base_hier2.R` which has hierarchical effects on age for antibody kinetics and the COP.

### 4. `model_run.R`

This script fits the statistical models to the preprocessed data to infer RSV infections using `serojump`.

### 5. `manu_fudan_burden.R/manu_fudan_abkinetics.R/manu_fudan_COP.R`

This script analyzes the results from the model fitting, generating the figures required for the paper.

### 6. `utils.R`

This script contains various helper functions used throughout the analysis, such as data manipulation, plotting functions, and utility functions.

## Contact

For any questions or issues, please contact the repository maintainer `david_hodgson@lshtm.ac.uk`.
