# process_inputs
#' A function that processes inputs to run a Pattern Component Modelling analysis of Representational Similarity
#' Analysis (RSA) data.
#'
#' @param input_filename A file containing *representational similarity analysis (RSA)* data for a given *region of interest (ROI)*,
#' in the form `[sample size x number of comparisons]`. This file's name should be of the format `ROI_vert.csv`, where `ROI` is the
#' name of a brain region from which the data was obtained (e.g. `ACC`).
process_inputs <- function(input_filename) {
  # Read in the data
  input_df <- read.csv(input_filename)
  # Extract the ROI name from the data
  ROI_name <- strsplit(input_filename, '_')[1]
  # Return outputs
  return(c(input_df, ROI_name))
}
