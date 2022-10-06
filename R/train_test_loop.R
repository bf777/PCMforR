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
#' @param holdout (default = 2) If `analysis_type` is `cross_val`, is an integer defining the number of values
#' in your data to hold out on each iteration.
#' @param output_dir A string to the folder where you would like all outputs to be written.
#' @param iter The current iteration (relevant for multi-iteration analyses such as cross-validation and individual analyses).
#' @param verbose (default = FALSE) If `TRUE`, output detailed information during computational runs.
train_test_loop <- function(data_for_ROI, ROI, ROI_idx, POI_names, POIs_list, analysis_type,
                            holdout, output_dir, CV_log, IND_log, summary_df, iter, verbose) {
  cat(paste('iter:', iter, '\n'))

  # TRAINING

  # Split data into train and test
  split_data_outputs <- split_data(data_for_ROI, analysis_type, holdout, iter)

  # Get split train and test data
  train_data <- split_data_outputs[[1]]
  test_data <- split_data_outputs[[2]]
  held_out_idx <- split_data_outputs[[3]]

  # If `analysis_type` == `cross_val`, apply holdout to all POI data as well
  # Create variable for each POI
  for (i in seq_along(POIs_list)) {
    data_for_POI <- POIs_list[i][[1]]
    POI_name <- POI_names[i]
    if (analysis_type == 'cross_val') {
      data_for_POI_train <- data.frame(data_for_POI[-held_out_idx,])
      data_for_POI_train <- unlist(data_for_POI_train)
      assign(POI_name, data_for_POI_train, envir = .GlobalEnv)
    } else if (analysis_type == 'individual') {
      data_for_POI_train <- data.frame(data_for_POI[iter,])
      data_for_POI_train <- unlist(data_for_POI_train)
      assign(POI_name, data_for_POI_train, envir = .GlobalEnv)
    } else {
      data_for_POI <- unlist(data.frame(data_for_POI))
      assign(POI_name, data_for_POI, envir = .GlobalEnv)
    }
  }

  # Initialize POIs placeholder dataframe for lm
  POIs_to_use <- vector()

  # Train the model at the current level by finding combination of paths with lowest BIC
  # run_BIC_at_level.R
  weighted_POIs_at_iter <- run_BIC_at_level(train_data, test_data, POIs_list, POI_names,
                                            POIs_to_use, analysis_type, ROI, ROI_idx,
                                            output_dir, held_out_idx, CV_log, IND_log,
                                            summary_df, iter, verbose)

  # TESTING
  if (analysis_type == 'cross_val') {
    ROI_reconstruction <- unlist(list(data.frame(weighted_POIs_at_iter[[4]])))
  } else {
    ROI_reconstruction <- unlist(list(weighted_POIs_at_iter[[4]]))
  }

  test_data <- unlist(test_data)

  lm_test <- lm(test_data ~ ROI_reconstruction)
  lm_test_summary <- summary(lm_test)

  CV_log <- weighted_POIs_at_iter[[2]]
  IND_log <- weighted_POIs_at_iter[[3]]

  if (analysis_type == 'cross_val') {
    # t statistics
    CV_log[iter, ((2 * (length(POI_names)) + 2) + 1):((2 * (length(POI_names)) + 2) + 4)] <- lm_test_summary$coefficients[2,1:4]

    # Overall R^2
    CV_log[iter, (2 * (length(POI_names)) + 2) + 5] <- lm_test_summary$adj.r.squared

    # F statistics
    CV_log[iter, ((2 * (length(POI_names)) + 2) + 6):((2 * (length(POI_names)) + 2) + 8)] <- lm_test_summary$f[1:3]

    data_logger(CV_log, 'CV_log', analysis_type, ROI, output_dir, n_path)
  } else if (analysis_type == 'individual') {
    # t statistics
    IND_log[iter, ((2 * (length(POI_names)) + 2) + 1):((2 * (length(POI_names)) + 2) + 4)] <- lm_test_summary$coefficients[2,1:4]

    # Overall R^2
    IND_log[iter, (2 * (length(POI_names)) + 2) + 5] <- lm_test_summary$adj.r.squared

    # F statistics
    IND_log[iter, ((2 * (length(POI_names)) + 2) + 6):((2 * (length(POI_names)) + 2) + 8)] <- lm_test_summary$f[1:3]
    data_logger(IND_log, 'IND_log', analysis_type, ROI, output_dir, n_path)
  }

  return(list(weighted_POIs_at_iter[[1]], CV_log, IND_log, ROI_reconstruction, weighted_POIs_at_iter[[5]]))
}
