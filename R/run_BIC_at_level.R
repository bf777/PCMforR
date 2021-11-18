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
run_BIC_at_level <- function(train_data, POIs_list, POI_names, POIs_to_use, analysis_type,
                             ROI, output_dir) {

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

  # Initialize min_BIC at -1
  min_BIC <- -1

  while(vert_level_idx < max_vert_level) {
    print(paste('Checking level', vert_level_idx), quote = FALSE)

    # Calculate all BIC scores at level
    BICs <- unlist(lapply(POI_names, calc_lm_BIC, POIs_to_use, train_data))

    # Record all BIC scores at vertical level in BIC log
    BIC_log[vert_level_idx, 1:length(POIs_list)] <- BICs

    # find which BIC is the best, save it in best_BICs, and index its position
    min_BIC <- min(BICs)
    best_BICs <- c(best_BICs, min_BIC)
    min_BIC_idx <- match(min_BIC, BICs)

    # horiz_level is list of equivalent BIC scores at horizontal level
    horiz_level <- vector()

    # Make a vector of all BICs except for the min one (for comparison)
    other_BICs <- BICs[-min_BIC_idx]

    # TODO: figure out how to best compare best BIC to other BICs

    # If BIC improved vs. criterion, check if BIC is better than other BICS +/- 2
    criterion <- 2
    for (other_BIC in other_BICs) {
      if (abs(min_BIC - other_BIC) > criterion) {
        # If best BIC already found in horizontal level, check if it's the same
        # horizontal level
        if (length(horiz_level) > 0) {
          if (length(duplicated(BICs))) {
            equivalent_BICs <- which(duplicated(BICS) == BICs)
            # If we have equivalency, record indices of equivalent BICs in horiz_level
            horiz_level <- c(horiz_level, equivalent_BICs)
            horiz_level_list <- c(horiz_level_list, horiz_level)
            # Continue at next vertical level
          }
        }
          # If not, continue directly to next vertical level and save current BIC
          # log table

          # add POI with new best BIC to POIs_to_use so that it can be included in
          # next level of BIC
          POIs_to_use <- c(POIs_to_use, POI_names[min_BIC_idx])
          BIC_log[vert_level_idx, length(POIs_list) + 1] <- min_BIC_idx
          BIC_log_list <- c(BIC_log_list, BIC_log)
          vert_level_idx <- vert_level_idx + 1
        } else {
          # If BIC not improved vs criterion:
          # Check if min BIC < best BIC
          if (min_BIC < best_BICs[length(best_BICs) - 1]) {
            best_BIC <- min_BIC
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
    }
    # If the horizontal level at the current vertical level is not empty,
    # continue BIC
    horiz_level_idx <- horiz_level_idx + 1

    # If BIC not improved vs. previous and empty level found, finish training

    # Weight POIs_list
    # POIs_to_add_weighting <- sapply(POIs_to_use, paste, collapse = " + ")
    # formula_to_use_weighting <- paste0('train_data ~ ', POIs_to_add_weighting)
    # weighted_POIs_lm <- lm(formula_to_use_weighting)
    # weighted_POIs <- summary(weighted_POIs_lm)$coefficients[1,1]

    # Export BIC logs
    # data_logger(BIC_log_list, 'BIC_log', analysis_type, ROI, output_dir)

    # Continue to testing
    # return(weighted_POIs)
  }
}

calc_lm_BIC <- function(POI, POIs_to_use, train_data) {
  # set up lm with train data ~ current POI + all identified POIs
  POIs_to_add <- sapply(POIs_to_use, paste, collapse = " + ")

  # If more than one POI has been identified, add identified POIs to lm
  if (length(POIs_to_add) != 0) {
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
