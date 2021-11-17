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

  # For max number of levels to look at when training model
  max_vert_level <- length(POIs_list)

  # Initialize list of all best BICs
  best_BICs <- vector()

  # Initialize horizontal and vertical levels
  horiz_level_idx <- 1
  vert_level_idx <- 1

  # Initialize BIC log
  BIC_log <- data.frame(matrix(0, nrow = length(POIs_list),
                               ncol = length(POIs_list)))
  colnames(BIC_log) <- POI_names

  # Initialize list of BIC logs for each vertical level
  BIC_log_list <- vector()

  # Initialize min_BIC at -1
  min_BIC <- -1

  while(vert_level_idx < max_vert_level) {
    noquote(paste('Checking level', vert_level_idx))

    # Calculate all BIC scores at level
    BICs <- unlist(lapply(POI_names, calc_lm_BIC, POIs_to_use, unlist(train_data)))

    # Record all BIC scores at vertical level in BIC log
    BIC_log[vert_level_idx,] <- BICs
    print(BIC_log)

    # find which BIC is the best, save it in best_BICs, and index its position
    min_BIC <- min(BICs)
    best_BICs <- c(best_BICs, min_BIC)
    min_BIC_idx <- match(min_BIC, BICs)

    # horiz_level is current vertical level in BIC_log
    horiz_level <- unlist(BIC_log[vert_level_idx,])

    # If BIC improved vs. previous, save index of min BIC to current horiz_level
    if (length(best_BICs) > 1) {
    if (min_BIC < best_BICs[length(best_BICs) - 1]) {

      # If best BIC already found in horizontal level, check if it's the same
      # horizontal level
      if (min_BIC_idx %in% horiz_level) {
        # If at same horizontal level, record index of best BIC
        BIC_log[vert_level_idx, length(POIs_list) + 1] <- min_BIC_idx
        vert_level_idx <- vert_level_idx + 1

        # add POI with new best BIC to POIs_to_use so that it can be included in
        # next level of BIC
        POIs_to_use[min_BIC_idx] <- POI_names[min_BIC_idx]

        # Continue at next vertical level
      } else {
        # If not, continue directly to next vertical level and save current BIC
        # log table
        vert_level_idx <- vert_level_idx + 1
        BIC_log_list <- c(BIC_log_list, BIC_log)
      }
    } else {
      # If BIC not improved vs previous:

      # TODO: add extra BIC check/path recording

      # Check if horizontal level is empty, so long as we're not on the first
      # vertical level
      while (vert_level_idx != 1) {
        if (length(vert_levels[vert_level_idx]) == 0) {
          } else {
            # If not, check horizontal level at previous vertical level
            vert_level_idx <- vert_level_idx - 1
          }
        }
      # If the horizontal level at the current vertical level is not empty,
      # continue BIC
      horiz_level_idx <- horiz_level_idx + 1
      break
    }
    }
  }

  # If BIC not improved vs. previous and empty level found, finish training

  # Weight POIs_list
  # lm(train_data ~ POIs_df[best_BICs])


  # Export BIC logs
  data_logger(BIC_log_list, 'BIC_log', analysis_type, ROI, output_dir)

  # Continue to testing
  }

calc_lm_BIC <- function(POI, POIs_to_use, train_data) {
  # set up lm with train data ~ current POI + all identified POIs
  POIs_to_add <- sapply(POIs_to_use, paste, collapse = " + ")

  # If moe than one POI has been identified, add identified POIs to lm
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
