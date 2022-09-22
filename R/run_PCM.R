# run_PCM
#' A function that accepts input to run a Pattern Component Modelling analysis of Representational Similarity
#' Analysis (RSA) data.

#' Accepts input information about the structure of your data.
#' @param input_dir A string to the folder containing an input file with a name defined in `POI_file`, as well
#' as RSA data files in CSV format for each *region of interest (ROI)* (whose names must end with the string `_vert.csv` and
#' which must be structured in the form `[sample size x number of comparisons]`).
#' @param output_dir A string to the folder where you would like all outputs to be written.
#' @param POI_file A string defining the path to a CSV file containing one row for each  *pattern of interest (POI)* which you want
#' to fit to your data, and one column for each comparison made, indicating the pattern of expected correlations for that POI. The first column
#' of this file should be the names of the POIs.  See (Kryklywy, Ehlers, Beukers, et al., 2020) for a discussion of the setup of -
#' and rationale behind - POIs.
#' @param analysis_type Possible values:  `'one_sample'`, `'two_sample'`, `'individual'`, `'cross_val'`. A string indicating the
#' type of analysis that you would like to carry out on the data.
#' - `one_sample`: Train a model on one sample and test it on the same sample.
#' - `two_sample`: Train a model on one sample and test it on a different sample.
#' - `individual`: Train a model on one sample and test it on one or more individuals' data.
#' - `cross_val`: Cross-validation - iteratively train a model on one sample and test it on a held-out sample.
#' @param holdout (default = 2) If `analysis_type` is `cross_val`, is an integer defining the number of values in your data to hold out on each iteration.
#' @param num_iters (default = 1000) If `analysis_type` is `cross_val`, is an integer defining the number CV iterations to run.
#' @returns a summary of PCM data indicating the POIs that best fit the data.
#' @export
run_PCM <- function(input_dir, output_dir, POI_file, analysis_type, holdout = 2, num_iters = 1000) {
  # 1. Read in data from each file in input path
  input_files <- list.files(input_dir, full.names = TRUE)
  ROI_files <- input_files[grepl('*_vert.csv', input_files)]

  # Extract a list of data for all ROIs
  data_to_use <- lapply(ROI_files, function(filename) read.csv(filename,
                                                               header = FALSE))

  # Extract ROI names
  ROIs <- lapply(ROI_files, function(filename) strsplit(basename(filename), '_')[[1]][1])

  # Get POIs
  POI_df <- read.csv(POI_file, header = FALSE)
  POI_names <- unlist(POI_df[1])
  POIs <- POI_df[2:ncol(POI_df)]

  # 2. Loop through iterations of training and testing
  # train_test_loop.R
  # Calculate sample size from number of rows in first data input file
  sample_size <- nrow(data_to_use[[1]])
  if (analysis_type == 'one_sample' | analysis_type == 'two_sample') {
    num_iters <- 1
  } else if (analysis_type == 'individual') {
    num_iters <- sample_size
  }
  # Initialize dataframe of identified POIs
  POIs_found_df <- data.frame(matrix(0, nrow = num_iters,
                                     ncol = length(POI_names)))

  # Convert each row of POIs into its own dataframe for lm
  POIs_list <- apply(POIs, 1, POIs_to_dfs, sample_size)

  # Initialize summary file
  summary_df  <- data.frame(matrix(NA, nrow = length(data_to_use),
                                   ncol = (length(POIs_list) * 2) + 11))
  colnames(summary_df) <- c('ROI', 'analysis_type', POI_names,
                        paste(POI_names, '-?', sep = ''),
                        'Recon-?', 'St.Err', 't', 'p', 'Adj.r', 'F', 'df1',
                        'df2', 'Branches')

  cat(paste('Analysis type:', analysis_type, '\n'))

  # else if (analysis_type == 'cross_val'), num_iters = number of iterations specified in num_iters (default = 1000)
  for (i in seq_along(data_to_use)) {
    data_for_ROI <- data_to_use[i]
    ROI <- ROIs[i]
    cat('\n')
    cat(paste('ROI:', ROI, '\n'))

    # Initialize CV log (will remain blank if not using analysis_type = `cross_val`)
    CV_log <- data.frame(matrix(0, nrow = num_iters,
                                ncol = (length(POIs_list) * 2) + 11 + holdout))
    colnames(CV_log) <- c('ROI', 'analysis_type', POI_names,
                          paste(POI_names, '-?', sep = ''),
                          'Recon-?', 'St.Err', 't', 'p', 'Adj.r', 'F', 'df1',
                          'df2', 'Branches', paste(rep('HO_', holdout),
                                                   seq(1, holdout),
                                                   sep = ''))

    # Initialize inidvidual log (will remain blank if not using analysis_type = `individual`)
    IND_log <- data.frame(matrix(0, nrow = num_iters,
                                ncol = (length(POIs_list) * 2) + 11))
    colnames(IND_log) <- c('ROI', 'analysis_type', POI_names,
                          paste(POI_names, '-?', sep = ''),
                          'Recon-?', 'St.Err', 't', 'p', 'Adj.r', 'F', 'df1',
                          'df2', 'Branches')

    for (iter in seq(1, num_iters)) {
      train_test_loop_output_at_iter <- train_test_loop(data_for_ROI, ROI, i,
                                                        POI_names, POIs_list,
                                                        analysis_type, holdout,
                                                        output_dir, CV_log, IND_log,
                                                        summary_df, iter)
      CV_log <- train_test_loop_output_at_iter[[2]]
      IND_log <- train_test_loop_output_at_iter[[3]]
      summary_df <- train_test_loop_output_at_iter[[5]]
    }

    if (analysis_type == 'cross_val') {
      summary_df[i, 3:((2 * (length(POI_names)) + 2) + 8)] <- colMeans(CV_log[,3:ncol(CV_log)])
    } else if (analysis_type == 'individual') {
      summary_df[i, 3:((2 * (length(POI_names)) + 2) + 8)] <- colMeans(IND_log[,3:ncol(IND_log)])
    }
  }

  data_logger(summary_df, 'summary', analysis_type, '', output_dir, n_path)

  cat('\n')
  cat(paste('PCMforR has finished running! You can find the outputs in', output_dir, '\n'))
}
