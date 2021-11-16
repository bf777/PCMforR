# run_BIC_at_level
#
#' A function that calculates and stores a *Bayesian Information Criterion (BIC)*
#' value at the current horizontal and vertical *pattern of interest (POI)*
#' level.
#'
#' @param train_data The data on which the *pattern component modelling (PCM)*
#' model will be trained.
#' @param POIs A list of *patterns of interest (POI)* which you want
#' to fit to your data, and one column for each comparison made,
#' indicating the pattern of expected correlations for that POI.
#' #' @param POIs The names of the *patterns of interest (POI)* which you want
#' to fit to your data.
#' @param POIs_df A dataframe containing each previously identified POI for use
#' in the BIC calculations.
#' @param analysis_type Possible values:  `'one_sample'`, `'two_sample'`,
#' `'individual'`, `'cross_val'`. A string indicating the type of analysis that
#' you would like to carry out on the data.
#' @param ROI The name of the current *region of interest (ROI)* within which to run the model.
run_BIC_at_level <- function(train_data, POIs, POI_names, POIs_df, analysis_type,
                             ROI, output_dir) {

  # For max number of levels to look at, train model
  max_vert_level <- length(POIs)

  # Initialize list of all best BICs
  best_BICs

  # Initialize horizontal and vertical levels
  horiz_level_idx <- 1
  vert_level_idx <- 1

  # Initialize BIC log
  BIC_log <- data.frame(matrix(0, nrows = length(POIs),
                               ncols = length(POIs) + 1))
  colnames(BIC_log) <- c(POI_names, 'best_BIC_idx')

  # Initialize list of BIC logs for each vertical level
  BIC_log_list <- vector()


  # Initialize min_BIC at -1
  min_BIC <- -1

  while(vert_level_idx < max_vert_level) {
    print(paste('Checking level', vert_level_idx))

    # Calculate all BIC scores at level
    BICs <- unlist(lapply(POIs, calc_lm_BIC, POIs_df, train_data))

    # Record all BIC scores at vertical level in BIC log
    BIC_log[vert_level_idx,] <- BICs

    # find which BIC is the best, and index its position
    min_BIC <- min(BICs)
    min_BIC_idx <- match(min_BIC, BICs)

    # horiz_level is current vertical level in BIC_log
    horiz_level <- unlist(BIC_log[vert_level_idx,])

    # If BIC improved vs. previous, save index of min BIC to current horiz_level
    if (min_BIC < best_BICs[length(best_BICs) - 1]) {

      # If best BIC already found in horizontal level, check if it's the same
      # horizontal level
      if (min_BIC_idx %in% horiz_level) {
        # If at same horizontal level, record index of best BIC
        BIC_log[vert_level_idx, length(POIs) + 1] <- min_BIC_idx
        vert_level_idx <- vert_level_idx + 1

        # add POI with new best BIC to POIs_df so that it can be included in
        # next level of BIC
        POIs_df[,min_BIC_idx] <- POIs[[min_BIC_idx]]

        # Continue at next vertical level
      } else {
        # If not, continue directly to next vertical level and save current BIC
        # log table
        vert_level_idx <- vert_level_idx + 1
        BIC_log_list.append(BIC_log)
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
        } else {
          # If the horizontal level at the current vertical level is not empty,
          # continue BIC
          horiz_level_idx <- horiz_level_idx + 1
          break
        }
    }
  }

  # If BIC not improved vs. previous and empty level found, finish training

  # Weight POIs
  # lm(train_data ~ POIs_df[best_BICs])


  # Export BIC logs
  data_logger(BIC_log_list, 'BIC_log', analysis_type)

  # Continue to training
  }

calc_lm_BIC <- function(POI, POIs_df, train_data) {
  POIs_df_for_lm <- data.frame(cbind(train_data, POI, POIs_df))
  lm_at_level <- lm(train_data ~ ., data = POIs_df_for_lm)
  BIC_at_level <- BIC(lm_at_level)
  return(BIC_at_level)
}
