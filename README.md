# WFDA

## Weighted Functional Data Analysis for the Calibration of a Ground Motion Model in Italy

This repository contains the implementations of the methodologies of Weighted Functional Data Analysis discussed in the paper "Weighted Functional Data Analysis for the Calibration of a Ground Motion Model in Italy". The repository contains all the necessary data and scripts to replicate the analysis conducted in the simulation study and in the case study. The scripts lead to the estimation of the functional regression coefficients of a Ground Motion Model calibrated for Italy, together with the quantification of the uncertainty associated to the estimates.

## Weighted Functional Data Analysis
Weighted Functional Data Analysis adapts the classical techniques of smoothing and function-on-function linear regression for the analysis of reconstructed partially observed functional data.

## Structure of the repository
The repository is composed of three folders:
* `code`: All scripts for replicating the analyses from the simulation study and the case study are reported here.
* `data`: The folder contains the flat file ITA18_SA_flatfile.csv, containing peak, duration and energy parameters, as well as the spectral acceleration
ordinates SA calculated assuming 5% damping in the range 0.01-10s and associated metadata of ITA18 dataset. 
* `output`: The folder contains all results and figures which can be replicated with the scripts in the `code` folder. The scripts automatically save in this folder the intermediate and final results and all figures.

### Code
* `simulation-set1.R`: The first simulation set discussed in Section 4.2 of the paper. The estimates from the weighted and the unweighted methodology are compared over different reconstruction methods.
* `simulation-set2.R`: The second simulation set discussed in Section 4.2 of the paper. The estimates from the weighted methodology are compared over multiple alternative definitions of the functional weights.
* `simulation-set3.R`: The second simulation set discussed in Section 2.3 of the supplementary material. The performance of the weighted methodology is assessed in scenarios with a varying fraction of partially observed data.
* `preprocessing-and-coll-analysis.R`: The code reproduces the preprocessing of the seismological data (loaded from the flat file in saved in the data folder), defines and saves the functional covariates as discussed in Section 2 of the manuscript, and performs the collinearity analysis presented in Section 3.1 of Supplement A.
* `application-calibration.R`: Calibration of the methodology for the case study, as described in Section 5.1 of the Manuscript. Specifically, the script allows to identify the best penalization parameters for the functional regression coefficients, and to select the optimal reconstruction strategy and functional weights for the seismic data.
* `application-main.R`: The script fits the final ground motion model, evaluates a bootstrap sample of the functional regression coefficients and performs a comparison with the scalar benchmark, producing the results presented in Section 5.2 of the Manuscript.
* `functions`: Folder containing all functions which have been implemented and are used to simulate partially observed functional data, perform the weighted analysis, estimate prediction errors via cross-validation, and plot results.

All analytical results and images are automatically saved in the output folder.

## Installation

Run file `install_packages.R` to have an automatic installation of the required R packages.
