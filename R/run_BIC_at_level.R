# run_BIC_at_level
#
#' A function that calculates and stores a *Bayesian Information Criterion (BIC)*
#' value at the current horizontal and vertical *pattern of interest (POI)* level.
#'
#' @param train_data The data on which the *pattern component modelling (PCM)* model will be trained.
#' @param horiz_level The *horizontal level* to explore when searching for an optimal combination of POIs.
#' @param vert_level The *vertical level* to explore when searching for an optimal combination of POIs.
run_BIC_at_level <- function(train_data, POI_names, POIs, POIs_df, horiz_level, vert_level, best_BICs) {
  # calculate BIC for each POI at horizontal level
  BICs <- unlist(lapply(POIs, calc_lm_BIC, POIs_df, train_data))
  # find which BIC is the best, and index its position
  min_BIC <- min(BICs)
  min_BIC_idx <- match(min_BIC, BICs)
  # append min BIC to list of min BIC scores, and update the list
  best_BICs.append(min_BIC_idx)
  assign('best_BICs', best_BICs, envir = globalenv())

  # add POI with new best BIC to POIs_df so that it can be included in next level of BIC
  POIs_df[,min_BIC_idx] <- POIs[[min_BIC_idx]]

  # If BIC improved vs. previous, save index of min BIC to horiz_level
  if (min_BIC < best_BICs[length(best_BICs) - 1]) {
    horiz_level.append(min_BIC_idx)
    assign('horiz_level', horiz_level, envir = globalenv())
    # TODO: Ask about equivalent horiz_level
  } else {
    # If BIC not improved vs previous:
    # Check if horizontal level is empty
    if (length(horiz_level) == 0) {
      # Check if we're on the first vertical level
      if (vert_level == 1) {
        # Weight POIs
        # lm(train_data ~ POIs_df[best_BICs])
        # Continue to training
      } else {
        # Check previous vertical level and return to BIC calculations
        assign('vert_level', vert_level - 1, envir = globalenv())
        # Finished here
      }
    }
  }
}

calc_lm_BIC <- function(POI, POIs_df, train_data) {
  POIs_df_for_lm <- data.frame(cbind(train_data, POI, POIs_df))
  lm_at_level <- lm(train_data ~ ., data = POIs_df_for_lm)
  BIC_at_level <- BIC(lm_at_level)
  return(BIC_at_level)
}
