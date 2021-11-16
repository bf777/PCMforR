# split_data.R
#' A function that splits data into training and testing data, depending on the analysis type.
#' @param input_filename A file containing *representational similarity analysis (RSA)* data for a given *region of interest (ROI)*,
#' in the form `[sample size x number of comparisons]`. This file's name should be of the format `ROI_vert.csv`, where `ROI` is the
#' name of a brain region from which the data was obtained (e.g. `ACC`).
#' @param analysis_type Possible values:  `'one_sample'`, `'two_sample'`, `'individual'`, `'cross_val'`. A string indicating the
#' type of analysis that you would like to carry out on the data.
#' - `one_sample`: Train a model on one sample and test it on the same sample.
#' - `two_sample`: Train a model on one sample and test it on a different sample.
#' - `individual`: Train a model on one sample and test it on one or more individuals' data.
#' - `cross_val`: Cross-validation - iteratively train a model on one sample and test it on a held-out sample.
#' @param holdout (default = 2) If `analysis_type` is `cross_val`, is an integer defining the number of values in your data to hold out on each iteration.
split_data <- function(input_filename, analysis_type, holdout) {
  # If analysis_type == 'two_sample', check which sample to use for training (contains `train` in name) and which to use for
  # testing (contains `test` in name)
  if(analysis_type == 'two_sample') {
    if ('train' %in% input_filename) {
      held_out_idx <- 0
      train_data <- input_df
    } else if ('test' %in% input_filename) {
      held_out_idx <- 0
      test_data <- input_df
    }
  }
  # If analysis_type == 'cross_val', randomly sample from the data according to the holdout size
  # If analysis_type == 'individual', use one row as data for testing
  if (analysis_type == 'cross_val') {
    held_out_idx <- sample.int(nrows(input_df), size = holdout)
    train_data <- input_df[-(held_out_idx)]
    test_data <- input_df[held_out_idx]
  } else if (analysis_type == 'individual') {
    held_out_idx <- 0
    train_data <- input_df
    test_data <- input_df
  } else if (analysis_type == 'one_sample') {
    held_out_idx <- 0
    train_data <- input_df
    test_data <- input_df
  }
  return(c(train_data, test_data))
}
