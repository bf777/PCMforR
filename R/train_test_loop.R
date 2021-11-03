# train_test_loop
#' A function that handles each iteration of training and testing a PCM model.
#'
#' @param input_filename A file containing *representational similarity analysis (RSA)* data for a given *region of interest (ROI)*,
#' in the form `[sample size x number of comparisons]`. This file's name should be of the format `ROI_vert.csv`, where `ROI` is the
#' name of a brain region from which the data was obtained (e.g. `ACC`).
#' @param ROI The name of the current *region of interest (ROI)* within which to run the model.
#' #' @param analysis_type Possible values:  `'one_sample'`, `'two_sample'`, `'individual'`, `'cross_val'`. A string indicating the
#' type of analysis that you would like to carry out on the data.
#' - `one_sample`: Train a model on one sample and test it on the same sample.
#' - `two_sample`: Train a model on one sample and test it on a different sample.
#' - `individual`: Train a model on one sample and test it on one or more individuals' data.
#' - `cross_val`: Cross-validation - iteratively train a model on one sample and test it on a held-out sample.
#' @param holdout (default = 2) If `analysis_type` is `cross_val`, is an integer defining the number of values in your data to hold out on each iteration.
#' @param num_iters (default = 1000) If `analysis_type` is `cross_val`, is an integer defining the number CV iterations to run.
train_test_loop <- function(input_filename, ROI, POI_names, POIs, analysis_type, holdout, best_BICs) {
  # Split data into train and test
  split_data_outputs <- split_data(input_filename, analysis_type, holdout)

  # Get split train and test data
  train_data <- split_data_outputs[1]
  test_data <- split_data_outputs[2]

  # Initialize horizontal and vertical levels
  horiz_level <- vector()
  vert_level <- vector()

  # Initialize horizontal and vertical level indices
  horiz_level_idx <- 1
  vert_level_idx <- 1

  # Initialize POIs dataframe for lm
  POI_placeholder <- rep(data.frame(matrix(1, nrow = nrow(POIs[[1]]),
                                           ncol = ncol(POIs[[1]]))), length(POIs))
  POIs_df <- data.frame(POI_placeholder)
  colnames(POIs_df) <- POI_names

  # Train the model at the current level by finding combination of paths with lowest BIC
  # run_BIC_at_level.R
  best_POIs_BICs <- lapply(train_data, run_BIC_at_level, POI_names, POIs, POIs_df, horiz_level, vert_level, best_BICs)
}
