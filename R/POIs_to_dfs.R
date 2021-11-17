# POIs_to_dfs
#' A function that takes individual rows of a dataframe defining all POIs and
#' converts it into its own dataframe for comparison to input data.
POIs_to_dfs <- function(POI_row, sample_size) {
  POI_row_df <- data.frame(t(POI_row))
  POI_df <- do.call("rbind", replicate(sample_size, POI_row_df, simplify = FALSE))
  return(POI_df)
}
