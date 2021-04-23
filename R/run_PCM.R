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
