# run_PCM
#' A function that accepts input to run a Pattern Component Modelling analysis of Representational Similarity
#' Analysis (RSA) data.

#' Accepts input information about the structure of your data.
#' @param input_dir A string to the folder containing an input file with a name defined in `POI_file`, as well
#' as RSA data files for each *region of interest (ROI)* (which must be structured in the form `[sample size x number of comparisons]`).
#' @param output_dir A string to the folder where you would like all outputs to be written.
#' @param POI_file A string defining the path to a CSV file containing one row for each  *patterns of interest (POI)* which you want
#' to fit to your data, and one column for each comparison made, indicating the pattern of expected correlations for that POI. The first column
#' of this file should be the names of the POIs.  See (Kryklywy, Ehlers, Beukers, et al., 2020) for a discussion of the setup of -
#' and rationale behind - POIs.
#' @param analysis_type Possible values:  `'one_sample'`, `'two_sample'`, `'individual'`, `'cross_val'`. A string indicating the
#' type of analysis that you would like to carry out on the data.
#' - `one_sample`: Train a model on one sample and test it on the same sample.
#' - `two_sample`: Train a model on one sample and test it on a different sample.
#' - `individual`: Train a model on one sample and test it on one or more individuals' data.
#' - `cross_val`: Cross-validation - iteratively train a model on one sample and test it on a held-out sample.
#' @param holdout (default = 2) If `analysis_type` is `cross_val`, is an integer defining the number of values in your data to hold out on each iteration.
#' @param num_iters (default = 1000) If `analysis_type` is `cross_val`, is an integer defining the number CV iterations to run.
#' @returns a summary of PCM data indicating the POIs that best fit the data.
#' @export
run_PCM <- function(input_dir, output_dir, POI_file, analysis_type, holdout = 2, num_iters = 1000) {
  # 1. Read in data from each file in input path
  input_files <- grepl(list.files(input_dir), '*_vert.csv')
  data_inputs <- lapply(input_files, process_inputs, analysis_type)
  # Extract a list of data for all ROIs
  data_to_use <- data_inputs[1]
  # Extract ROI names
  ROIs <- data_inputs[2]

  # 2.Loop through iterations of training and testing
  # 2.1. Split data
  # split_data.R
  # 2.2. Calculate and update BIC scores
  # calc_BIC.R
}
