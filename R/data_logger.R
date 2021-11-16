# data_logger
#
#' A function that handles recording of data output
#'
#' @param input_data The data to be recorded to a .csv file.
#' @param data_type The type of data to be recorded.
#' @param POIs_df A dataframe containing each previously identified POI for use
#' in the BIC calculations.
#' @param analysis_type Possible values:  `'one_sample'`, `'two_sample'`,
#' `'individual'`, `'cross_val'`. A string indicating the type of analysis that
#' you would like to carry out on the data.
#' @param ROI The name of the current *region of interest (ROI)* within which to run the model.
#' @param output_dir A string to the folder where you would like all outputs to be written.
data_logger <- function(input_data, data_type, POIs_df, analysis_type, ROI,
                        output_dir) {
  # Interim logging of BIC scores to .csv (if not doing cross-validation)
  if (data_type == 'BIC_log' & analysis_type != 'cross_val') {
    # Record BICs for each level
    for (BIC_log_idx in seq(1, length(BIC_log_list))) {
    write.csv(input_data[BIC_log_idx], file.path(output_dir,
                                                 paste('BIC_log_level_', BIC_log_idx,
                                                       '_ROI_', ROI,
                                                       '_analysis', analysis_type,
                                                       '.csv', sep = '')))
    }
    # Best BICs
    write.csv(input_data[length(input_data)], file.path(output_dir,
                                                        paste('BIC_log_ROI_', ROI,
                                                              '_analysis_', analysis_type,
                                                              '.csv', sep = '')))
  }
}
