# run_BIC_at_level
#
#' A function that calculates and stores a *Bayesian Information Criterion (BIC)*
#' value at the current horizontal and vertical *pattern of interest (POI)* level.
#'
#' @param train_data The data on which the *pattern component modelling (PCM)* model will be trained.
#' @param horiz_level The *horizontal level* to explore when searching for an optimal combination of POIs.
#' @param vert_level The *vertical level* to explore when searching for an optimal combination of POIs.
run_BIC_at_level <- function(train_data, POIs, horiz_level, vert_level) {
  # calculate BIC for all POIs at horizontal level
  POIs_lms <- lapply(POIs, calc_lm_BIC, train_level, horiz_level, vert_level)

}

calc_lm_BIC <- function(POI, train_data, horiz_level, vert_level) {
  lm_at_level <- lm(train_data ~ POI)
  BIC_at_level <- BIC(lm_at_level)
}
