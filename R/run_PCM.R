# run_PCM
#' A function that accepts input to run a Pattern Component Modelling analysis of Representational Similarity
#' Analysis (RSA) data.

#' Accepts input information about the structure of your data.
#' @param input_path A string to the folder containing an input file with a name defined in `IPC_file`, as well
#' as RSA data files for each ROI (which must be structured in the form `[sample size x number of comparisons]`).
#' @param output_path A string to the folder where you would like all outputs to be written.
#' @param ROIs A vector of strings naming *regions of interest (ROIs)* in your dataset. PCMs will be run on
#' data from each ROI in the order named here. Make sure that your RSA data is split by ROI into individual
#' .csv files.
#' @param IPC_names A vector of strings naming *information pattern components (IPCs)* to which you want
#' to fit your data. Each IPC will be represented as a row in the file defined in `IPC_file`, equal in length
#' to the number of comparisons in your dataset (which is equal to the cumulative sum of the number of
#' conditions in your dataset). Each cell of this row should contain the expected correlation
#' of that comparison. See (Kryklywy, Ehlers, Beukers, et al., 2020) for a discussion of the setup of -
#' and rationale behind - IPCs.
#' @param IPC_file A string defining the path to a CSV file containing one row for each IPC (see `IPC_names`),
#' and one column for each comparison made, indicating the pattern of expected correlations for that IPC.
#' @param IPC_prior_level A vector containing integers for each level of priors that you would like to check.
#' @param DataType A string indicating the type of data structure that you would like to use. Currently,
#' PCMforR only supports `'US'` (Unconditioned stimulus) as a data input.
#' @param Analysis Possible values:  `'Group'`, `'CV'`, `'Individual'`, `'Mixed'`. A string indicating the
#' type of analysis that you would like to carry out on the data. `'Group'` will analyze data for all
#' participants as a group. `'CV'` will conduct a cross-validation analysis on the group data,re-running the
#' IPC comparison to the data iteratively in order to achieve the best possible fit on the data.
#' `'Individual'` and `'Mixed'` are currently not fully implemented.
#' @param SampleSize An integer indicating the number of participants in your data.
#' @param num_Conditions An integer indicating the number of trial conditions in your dataset. This will be
#' used to determine the number of comparisons in your dataset.
#' @param Condition_names A vector of strings naming the conditions in your dataset (NOTE: not currently used).
#' @param Hold If using cross-validation (CV), is an integer defining the number of values to hold out on each iteration.
#' @param cv_iterations If using cross-validation (CV), is an integer defining the number CV iterations to run.
#' @returns Does not return a specific output, but outputs a summary of PCM data indicating the IPCs that best fit the data.
#' @export
run_PCM <- function(input_path, output_path, ROIs, IPC_names, IPC_file, IPC_prior_level, DataType, Analysis, SampleSize,
                    num_Conditions, Condition_names, Hold, cv_iterations = 1000) {
  # Determine other data characteristics from input
  num_Comparisons <- sum(1:num_Conditions)
  num_IPC <- length(IPC_names)

  # Dataset size defaults
  if (Analysis == "Group") {
    repetitions <- 1
    train_n <- SampleSize
    test_n <- SampleSize
  } else if (Analysis == "Mixed") {
    repetitions <- SampleSize
    train_n <- SampleSize
    test_n <- 1
  } else if (Analysis == "Individual") {
    repetitions <- SampleSize
    train_n <- 1
    test_n <- 1
  } else if (Analysis == "CV") {
    repetitions <- cv_iterations
    train_n <- (SampleSize - Hold)
    test_n <- Hold
  }

  # Set data log defaults
  data_file <- c("")
  data_name <- c("PCM")

  # Import IPCs and size them appropriately
  All_IPCs <- read.table(file.path(input_path, IPC_file), header = FALSE, sep = ",")
  testIPCs <- unlist(list(matrix(, nrow = 1, ncol = num_Comparisons * test_n)))
  trainIPCs <- unlist(list(matrix(, nrow = 1, ncol = num_Comparisons * train_n)))

  IPC_list <- seq(1, num_IPC)

  # Restructure IPCs into column-wise form for easier comparison
  IPCs_types_list <- assign_IPCs(num_IPC, IPC_names, All_IPCs, test_n, train_n, testIPCs, trainIPCs, num_Comparisons)
  testIPCs <- IPCs_types_list[[1]]
  trainIPCs <- IPCs_types_list[[2]]

  # Prepare master summary
  if (Analysis == "Mixed") {
    Summary_Headers <- c(
      "Analysis", "ROI", IPC_names, paste("g", IPC_names, "-?", sep = ""),
      "Recon-?", "St.Err", "T", "p", "Adj.r", "f", "df1", "df2", paste("i", IPC_names, "-?", sep = "")
    )
    Summary <- matrix(, nrow = (length(ROIs) + 1), ncol = (3 * length(ROIs) + 9))
  } else {
    Summary_Headers <- c(
      "Analysis", "ROI", IPC_names, paste(IPC_names, "-?", sep = ""),
      "Recon-?", "St.Err", "T", "p", "Adj.r", "f", "df1", "df2", "Branches"
    )
    Summary <- matrix(, nrow = (length(ROIs) + 1), ncol = length(Summary_Headers))
  }
  colnames(Summary) <- Summary_Headers
  ROI_nums <- seq(1, length(ROIs))

  # Conduct PCMs for each ROI
  Summary_output <- lapply(
    ROI_nums, get_data_by_ROI, ROIs, Analysis, data_file, Summary, train_n, test_n, data_name,
    repetitions, num_IPC, num_Comparisons, testIPCs, trainIPCs, SampleSize, input_path, output_path
  )
  for (summary_row in 1:length(ROIs)) {
    Summary[summary_row, 1:length(Summary_Headers)] <- Summary_output[[summary_row]][summary_row + 1, 1:length(Summary_Headers)]
  }
  Summary <- Summary[1:length(ROIs), 1:length(Summary_Headers)]

  # Finally, output summary table
  write.table(Summary, file.path(
    output_path,
    paste(Analysis, "_", DataType, "_summary.csv", sep = "")
  ), row.names = FALSE, sep = ",")
  print(paste("PCM complete! It was run on", SampleSize, "participants in", length(ROIs) ,"ROIs. You can view the outputs in the output folder:", output_path))
}
