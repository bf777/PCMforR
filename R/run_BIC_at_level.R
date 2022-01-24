# run_BIC_at_level
#
#' A function that calculates and stores a *Bayesian Information Criterion (BIC)*
#' value at the current horizontal and vertical *pattern of interest (POI)*
#' level.
#'
#' @param train_data The data on which the *pattern component modelling (PCM)*
#' model will be trained.
#' @param POIs_list A list of *patterns of interest (POI)* which you want
#' to fit to your data, with the shape `[sample size x number of comparisons]`.
#' @param POI_names The names of the *patterns of interest (POI)* which you want
#' to fit to your data.
#' @param POIs_to_use A dataframe containing each previously identified POI for use
#' in the BIC calculations.
#' @param analysis_type Possible values:  `'one_sample'`, `'two_sample'`,
#' `'individual'`, `'cross_val'`. A string indicating the type of analysis that
#' you would like to carry out on the data.
#' @param ROI The name of the current *region of interest (ROI)* within which to run the model.
#' @param held_out_idx If `analysis_type` is `cross_val`,indicates the held-out values for this iteration.
#' @param iter The current iteration (relevant for multi-iteration analyses such as cross-validation and individual analyses).
run_BIC_at_level <- function(train_data, POIs_list, POI_names, POIs_to_use, analysis_type,
                             ROI, output_dir, held_out_idx, CV_log, iter) {

  # Unlist train_data for lm
  train_data <- unlist(train_data)

  # For max number of levels to look at when training model
  max_vert_level <- length(POIs_list)

  # Initialize list of all best BICs
  best_BICs <- vector()

  # Initialize horizontal and vertical levels
  horiz_level_idx <- 1
  vert_level_idx <- 1

  # Initialize BIC log
  BIC_log <- data.frame(matrix(0, nrow = length(POIs_list),
                               ncol = length(POIs_list) + 1))
  colnames(BIC_log) <- c(POI_names, 'best')

  # Initialize list of BIC logs for each vertical level
  BIC_log_list <- vector()

  # Initialize list of horizontal levels for each vertical level
  horiz_level_list <- vector()

  # Initialize BIC paths dataframe
  BIC_paths <- data.frame(matrix(NA, nrow = length(POIs_list),
                                 ncol = length(POIs_list)))

  # Initialize min_BIC at -1
  min_BIC <- -1

  # Initialize path number at 1
  n_path <- 1

  while(vert_level_idx < max_vert_level) {
    if (analysis_type != 'cross_val') {
      print(paste('Checking level', vert_level_idx), quote = FALSE)
    }

    # Calculate all BIC scores at level
    BICs <- unlist(lapply(POI_names, calc_lm_BIC, POIs_to_use, train_data))

    # find which BIC is the best, save it in best_BICs, and index its position
    min_BIC <- min(BICs)
    best_BICs <- c(best_BICs, min_BIC)
    min_BIC_idx <- match(min_BIC, BICs)

    # horiz_level is list of equivalent BIC scores at horizontal level
    horiz_level <- vector()

    # add POI with new best BIC to POIs_to_use so that it can be included in
    # next level of BIC
    # Record all BIC scores at vertical level in BIC log
    BIC_log[vert_level_idx, 1:length(POIs_list)] <- BICs
    POIs_to_use <- c(POIs_to_use, POI_names[min_BIC_idx])
    BIC_log[vert_level_idx, length(POIs_list) + 1] <- min_BIC_idx
    BIC_paths[horiz_level_idx, vert_level_idx] <- min_BIC_idx

    # If BIC improved vs. previous best BIC more than criterion, check if BIC is better than other BICS +/- 2

    # If we're not on the first vertical level (we have a prior BIC to compare to),
    # compare BICs
    if (vert_level_idx >= 3) {
      criterion <- 2
      previous_min_BIC <- best_BICs[length(best_BICs) - 1]
      min_BIC_two_back <- best_BICs[length(best_BICs) - 2]
      if (abs(min_BIC - previous_min_BIC) > criterion) {
        # If best BIC already found in horizontal level, check if it's the same
        # horizontal level
        if (length(horiz_level) == 0) {
          if (length(duplicated(BICs))) {
            equivalent_BICs <- which(duplicated(BICs) == BICs)
            # If we have equivalency and we're not already in a horizontal level,
            # record indices of equivalent BICs in horiz_level
            horiz_level <- c(horiz_level, equivalent_BICs)
            horiz_level_list <- c(horiz_level_list, horiz_level)
            # Continue at next vertical level
          }
        }
        # If not, continue directly to next vertical level and save current BIC
        # log table

        } else {
          # If BIC not improved vs criterion:
          # Check if min BIC < best BIC
          # print(min_BIC)
          # print(previous_min_BIC)
          if (min_BIC < previous_min_BIC) {
            data_logger(BIC_log, 'BIC_log', analysis_type, ROI, output_dir, 'best')
          } else if (min_BIC == min_BIC_two_back) {
            data_logger(BIC_log, 'BIC_log', analysis_type, ROI, output_dir, 'best')
            break
          }

        # Check if horizontal level is empty, so long as we're not on the first
        # vertical level
        while (length(horiz_level_list[vert_level_idx]) == 0) {
          if (vert_level_idx == 1) {
            # If on first vertical level, continue to testing
            break
          } else {
            # If not, check horizontal level at previous vertical level
            vert_level_idx <- vert_level_idx - 1
            print(paste("Returning to level", vert_level_idx, quote = FALSE))
          }
          }
        }
      # If the horizontal level at the current vertical level is not empty,
      # continue BIC
      horiz_level_idx <- horiz_level_idx + 1
      n_path <- n_path + 1

      # Export BIC logs for path
      data_logger(BIC_log, 'BIC_log', analysis_type, ROI, output_dir, n_path)
    }

  # If we're on the first level, go to the next level without recording BIC
  vert_level_idx <- vert_level_idx + 1
  }
  # If BIC not improved vs. previous and empty level found, finish training

  # Weight POIs_list
  POIs_to_add_weighting <- paste(POIs_to_use, collapse = " + ")
  formula_to_use_weighting <- paste0('train_data ~ ', POIs_to_add_weighting)
  POI_weights_lm <- lm(as.formula(formula_to_use_weighting))

  # Get betas
  POI_weights <- summary(POI_weights_lm)$coefficients[,1]
  POI_weights_for_log <- rep(0, length(POI_names))
  POI_weights_for_log[match(unique(POIs_to_use), POI_names)] <- POI_weights

  # Get overall coefficient (intercept)
  POI_overall_fit <- summary(POI_weights_lm)$coefficients[1,1]

  # Get overall R^2
  POI_R_squared <- summary(POI_weights_lm)$adj.r.squared

  # Get F statistic (df1, df2, f)
  POI_F_stats <- summary(POI_weights_lm)$f[1:3]

  # If we're doing individual or cross-validation, record POI_weights, POI_R_squared, POI_F_stats for iteration
  if (analysis_type == 'cross_val') {
    # Column 1: ROI
    CV_log[iter, 1] <- ROI

    # Column 2: analysis type
    CV_log[iter, 2] <- analysis_type

    # Identified POIs
    identified_POIs <- rep(0, length(POI_names))
    identified_POIs[match(unique(POIs_to_use), POI_names)] <- 1
    CV_log[iter, 3:length(POI_names)] <- identified_POIs

    # POI weights
    CV_log[iter, (length(POI_names) + 3):(2 * (length(POI_names)) + 2)] <- POI_weights_for_log

    # Overall fit
    CV_log[iter, (2 * (length(POI_names)) + 2) + 1] <- POI_overall_fit

    # Overall R^2
    CV_log[iter, (2 * (length(POI_names)) + 2) + 5] <- POI_R_squared

    # F statistics
    CV_log[iter, ((2 * (length(POI_names)) + 2) + 6):((2 * (length(POI_names)) + 2) + 8)] <- POI_F_stats

    data_logger(CV_log, 'CV_log', analysis_type, ROI, output_dir, n_path)
  }

  weighted_POIs_list <- vector()

  for (POI_to_weight_idx in seq(1,length(POIs_list))) {
    POI_to_weight <- unlist(POIs_list[POI_to_weight_idx])
    POI_to_weight_name <- POI_names[POI_to_weight_idx]
    if(POI_to_weight_name %in% POIs_to_use) {
      POI_to_weight <- POI_to_weight * POI_weights[POI_to_weight_name]
      weighted_POIs_list <- c(weighted_POIs_list, POI_to_weight)
    }
  }

  # Export all BIC paths
  data_logger(BIC_paths, 'BIC_paths', analysis_type, ROI, output_dir, n_path)

  # Continue to testing
  return(list(weighted_POIs_list, CV_log))
}

calc_lm_BIC <- function(POI, POIs_to_use, train_data) {
  # set up lm with train data ~ current POI + all identified POIs

  # If more than one POI has been identified, add identified POIs to lm
  if (length(POIs_to_use) != 0) {
    POIs_to_add <- paste(POIs_to_use, collapse = " + ")
    formula_to_use <- paste0('train_data ~ ', POI, ' + ', POIs_to_add)
  } else {
    formula_to_use <- paste0('train_data ~ ', POI)
  }

  # Run lm
  lm_at_level <- lm(as.formula(formula_to_use))

  # Calculate BIC on lm
  BIC_at_level <- BIC(lm_at_level)
  return(BIC_at_level)
}
