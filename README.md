# WFDA

## Weighted Functional Data Analysis for the Calibration of a Ground Motion Model in Italy

This repository contains the implementations of the methodologies of Weighted Functional Data Analysis discussed in the paper "Weighted Functional Data Analysis for the Calibration of a Ground Motion Model in Italy". The repository contains all the necessary data and scripts to replicate the analysis conducted in the simulation study and in the case study. The scripts lead to the estimation of the functional regression coefficients of a Ground Motion Model calibrated for Italy, together with the quantification of the uncertainty associated to the estimates.

## Weighted Functional Data Analysis
Weighted Functional Data Analysis adapts the classical techniques of smoothing and function-on-function linear regression for the analysis of reconstructed partially observed functional data.

## Structure of the repository
The repository is composed of:
* `preprocessing.R`: R script that operates the preprocessing of data loaded from scratch, namely from the files in .mat extension. The code loads data in the .mat files, extract the useful information and creates and saves in directiory `DATA` all .RData files that are going to be directly loaded in the analysis. In particular, the function covariates are generated from the original seismic parameters via a smoothing technique.
* `main.R`: R script that performs the WFDA analysis, from curves reconstruction via linear extrapolation to the weighted function-on-function regression and diagnostic.
* `main-unwgt.R`: R script that performs the unweighted FDA analysis, from curves reconstruction via linear extrapolation to the function-on-function regression and diagnostic.
* `main-bootstrap.R`: R script that performs bootstrap sampling from the distribution of the regression functional coefficients.
* `blist-selection.R`: R script that performs a preliminary analysis on the most appropriate choice the q penalization parameters introduced in the objective function of the penalized weighted least squares criterion for regression. The analysis consists in a 5-fold cross-validation that identifies the optimal parameters as those minimizing the Mean Squared Error.
* `reconstruction-methods-comparison.R`: R script that performs a cross-validation analysis between some reconstruction methods present in the literature and the extrapolation method associated to the weighted analysis. The methods are compared by means of the Mean Squared Error.
* `DATA`: Folder containing the data of the case study. In particular, `curves.RData` is the file containing the values of the registered curves at the sampling instants. It is a matrix of dimension N x n, where N is the number of sampling instant and n is the numerosity of the curves. `data.RData` contains the values of the original seismic parameters, i.e. Joyner-Boore distance (dJB), magnitude (MAG), style-of-faulting (SoF), shear-wave velocity (VS30) and frequencies of registration U_hp and V_hp. `xlist.RData` contains a 9-element list of the functional covariates, after a smoothing procedure was performed in a preprocessing step. `events.RData` contains information about the event associated to each curve, namely the event identification number (event.id), the longitude and latitude (event.long, event.lat).
* `methods`: Folder containing all functions used in the analysis. Each function in the folder contains a brief description of its usage and of its input and return parameters.  `methods/plots` contains the functions used to create the plots and to do diagnostic on the results of the analysis. `methods/Regression` contains the functions used to perform the function-on-function linear regression.
* `Simulation`: Folder containing all codes and results related to the simulation study discussed in the paper.

## Installation

Run file `install_packages.R` to have an automatic installation of the required R packages.
