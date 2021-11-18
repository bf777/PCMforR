# train_test_loop
#' A function that handles each iteration of training and testing a PCM model.
#'
#' @param input_filename A file containing *representational similarity analysis (RSA)* data for a given *region of interest (ROI)*,
#' in the form `[sample size x number of comparisons]`. This file's name should be of the format `ROI_vert.csv`, where `ROI` is the
#' name of a brain region from which the data was obtained (e.g. `ACC`).
#' @param ROI The name of the current *region of interest (ROI)* within which to run the model.
#' @param analysis_type Possible values:  `'one_sample'`, `'two_sample'`, `'individual'`, `'cross_val'`. A string indicating the
#' type of analysis that you would like to carry out on the data.
#' - `one_sample`: Train a model on one sample and test it on the same sample.
#' - `two_sample`: Train a model on one sample and test it on a different sample.
#' - `individual`: Train a model on one sample and test it on one or more individuals' data.
#' - `cross_val`: Cross-validation - iteratively train a model on one sample and test it on a held-out sample.
#' @param holdout (default = 2) If `analysis_type` is `cross_val`, is an integer defining the number of values in your data to hold out on each iteration.
#' @param num_iters (default = 1000) If `analysis_type` is `cross_val`, is an integer defining the number CV iterations to run.
train_test_loop <- function(data_for_ROI, ROI, POI_names, POIs_list, analysis_type,
                            holdout, best_BICs, output_dir) {

  print(paste('ROI:', ROI), quote = FALSE)

  # Split data into train and test
  split_data_outputs <- split_data(data_for_ROI, analysis_type, holdout)

  # Get split train and test data
  train_data <- split_data_outputs[1]
  test_data <- split_data_outputs[2]

  # Initialize POIs placeholder dataframe for lm
  POIs_to_use <- vector()

  # Train the model at the current level by finding combination of paths with lowest BIC
  # run_BIC_at_level.R
  weighted_POIs_at_iter <- lapply(train_data, run_BIC_at_level, POIs_list, POI_names,
                                  POIs_to_use, analysis_type, ROI, output_dir)
}
