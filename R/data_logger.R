# data_logger
#
#' A function that handles recording of data output
#'
#' @param input_data The data to be recorded to a .csv file.
#' @param data_type The type of data to be recorded.
#' @param analysis_type Possible values:  `'one_sample'`, `'two_sample'`,
#' `'individual'`, `'cross_val'`. A string indicating the type of analysis that
#' you would like to carry out on the data.
#' @param ROI The name of the current *region of interest (ROI)* within which to run the model.
#' @param output_dir A string to the folder where you would like all outputs to be written.
data_logger <- function(input_data, data_type, analysis_type, ROI, output_dir, n_path) {
  # Interim logging of BIC scores to .csv (if not doing cross-validation)
  if (data_type == 'BIC_log' & analysis_type != 'cross_val') {
    # Best BICs
    if(n_path == 'best') {
      write.csv(input_data, file.path(output_dir, paste('BIC_log_ROI_', ROI,
                                                        '_analysis_', analysis_type,
                                                        '.csv', sep = '')), row.names = FALSE)
    } else {
      # BICs for current path
      write.csv(input_data, file.path(output_dir, paste('BIC_log_ROI_', ROI,
                                                        '_analysis_', analysis_type,
                                                        '_path_', n_path,
                                                        '.csv', sep = '')), row.names = FALSE)
    }
  } else if (data_type == 'BIC_paths') {
    write.table(input_data, file.path(output_dir, paste('BIC_paths_ROI_', ROI,
                                                      '_analysis_', analysis_type,
                                                      '.csv', sep = '')), row.names = FALSE, col.names=FALSE, sep = ',')
  } else if (data_type == 'CV_log') {
    # Record CV log
    write.csv(input_data, file.path(output_dir, paste('CV_log_ROI_', ROI,
                                                      '_analysis_', analysis_type,
                                                      '.csv', sep = '')))
  } else if (data_type == 'summary') {
    # Record summary data
  }
}
