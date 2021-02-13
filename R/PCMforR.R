# Package documentation
#' PCMforR: Pattern Component Modelling for R
#'
#' This package allows you to easily conduct a Pattern
#' Component Modelling (PCM) analysis on Representational Similarity
#' Analysis (RSA) data (e.g. from functional brain imaging data) in
#' R. PCM allows you to take a set of hypothetical structures - Independent
#' Pattern Components, or IPCs - by which data might be represented in a given
#' brain region - and compare these structures to your data to determine which
#' hypothesis best fits your data using a Bayesian Information Criterion (BIC)
#' metric. For more details on PCM analysis, see
#' [Kryklywy, Ehlers, Beukers, et al., 2020](https://doi.org/10.1101/2020.09.24.310383).
#'
#' @section Functions:
#' * `run_PCM()`: The main runnable function of the package. Given input representational similarity analysis (RSA)
#' data from different regions of interest (ROIs) in the brain along with hypothetical Information Pattern Component (IPC)
#' models, outputs statistics about the IPC models that best fit your data.
#' * `get_data_by_ROI()`: Extracts data for each ROI, then runs Bayesian Information Criterion (BIC) analyses and
#' returns a summary table. This function is run in a loop by `run_PCM()` for each ROI.
#' * `run_BIC()`: Returns Bayesian Information Criterion (BIC) information comparing participant data for that ROI to the input
#' IPCs, in order to determine which IPCs best fit the data. This function is run in a loop by `get_data_by_ROI()`.
#' * `compare_IPC()`: Tracks, identifies, and records the best-matching IPCs, along with associated statistical information,
#' for output.
#' * `check_levels()`: Runs the loop that checks whether additional paths are needed to find the best IPC fit. This function is
#' run in a loop by `compare_IPC()`.
#' * `follow_paths()`: Checks individual levels of the IPC search path, and determines whether the model has improved or not.
#'
#' @section Acknowledgments:
#' The code underpinning this package was conceptualized and written by Dr. James Kryklywy (james.kryklywy[at]gmail.com). It was
#' reformatted into an R package by Brandon Forys (brandon.forys[at]alumni.ubc.ca). The project was supervised by
#' Dr. Rebecca Todd (becket.todd[at]psych.ubc.ca) at the University of British Columbia in Vancouver, BC, Canada.
#'
#' @references
#' If you use this package in your own work, please cite it as follows:
#'
#' Kryklywy, J.H., Forys, B.J., & Todd, R.M., (2021). Pattern Component Modelling for R (PCMforR), R package, version `r packageVersion("PCMforR")`.
#'
#' BibTeX:
#' ```
#' [at]Manual{,
#' title = {PCMforR: Pattern Component Modelling for R},
#' author = {James H. Krylywy and Brandon J. Forys and Rebecca M. Todd},
#' year = {2021},
#' note = {R package version 0.1.0},
#' }
#' ````
#'
#' @docType package
#' @name PCMforR
NULL
# > NULL
