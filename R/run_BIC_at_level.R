# run_BIC_at_level
#
#' A function that calculates and stores a *Bayesian Information Criterion (BIC)*
#' value at the current horizontal and vertical *pattern of interest (POI)* level.
#'
#' @param train_data The data on which the *pattern component modelling (PCM)* model will be trained.
#' @param horiz_level The *horizontal level* to explore when searching for an optimal combination of POIs.
#' @param vert_level The *vertical level* to explore when searching for an optimal combination of POIs.
run_BIC_at_level <- function(train_data, POI_names, POIs, POIs_df, horiz_level, vert_level) {
  # calculate BIC for each POI at horizontal level
  BICs <- unlist(lapply(POIs, calc_lm_BIC, POIs_df, train_level, horiz_level, vert_level))
  # find which BIC is the best, and index its position
  min_BIC <- min(BICs)
  min_BIC_idx <- match(min_BIC, BICs)
  # save index of min BIC as new horiz_level
  horiz_level <<- min_BIC_idx
}

calc_lm_BIC <- function(POI, POIs_df, train_data, horiz_level, vert_level) {
  POIs_df_for_lm <- data.frame(cbind(train_data, POI, POIs_df))
  lm_at_level <- lm(train_data ~ ., data = POIs_df_for_lm)
  BIC_at_level <- BIC(lm_at_level)
  return(BIC_at_level)
}
