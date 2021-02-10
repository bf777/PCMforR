# prep_calc_BIC
#' A series of functions for extracting data for each ROI and conducting a Bayesian Information Criterion (BIC)
#' analysis on data from each ROI.

#' Extracts data for each ROI, then runs Bayesian Information Criterion (BIC) analyses and returns a summary table. This
#' function is run in a loop by `run_PCM()` for each ROI.
#' @returns A Summary matrix, containing summary information on the best-matching IPCSs for output to the final data file.
#' @export
get_data_by_ROI <- function(ROI_num, ROIs, Analysis, data_file, Summary, train_n, test_n,
                            data_name, repetitions, num_IPC, num_Comparisons, testIPCs, trainIPCs, SampleSize, input_path,
                            output_path) {
  ROI <- ROIs[ROI_num]
  print(paste("Running PCM on ROI:", ROI))
  ROI_Data <- read.table(file.path(
    input_path,
    paste(ROI, "_US_vert.csv", sep = "")
  ), header = FALSE, sep = ",")

  for (timing in 1:(length(data_file))) {
    file_name <- data_file[timing]
    dataset_name <- data_name[timing]
    if (Analysis == "Mixed") {
      DataType2 <- "CS"
    } else {
      DataType2 <- DataType
    } # CS
    assign(paste(dataset_name, "ROI_Data", sep = ""), read.table(file.path(input_path, paste(ROI, "_", DataType2, "_", file_name, "vert.csv", sep = "")),
      header = FALSE, sep = ","
    ))
  }

  Summary[(ROI_num + 1), 1] <- Analysis
  Summary[(ROI_num + 1), 2] <- ROI

  if (Analysis == "CV" | Analysis == "Group") {
    Headers <- c(IPC_names, paste(IPC_names, "-?", sep = ""), "Recon-?", "St.Err", "T", "p", "Adj.r", "f", "df1", "df2", "Branches")
    CV_Log <- matrix(0, nrow = (repetitions + 1), ncol = length(Headers))
    colnames(CV_Log) <- Headers
  } else if (Analysis == "Mixed") {
    Mixed_headers <- c(
      IPC_names, paste(IPC_names, "-?", sep = ""), "Recon-?", "St.Err", "T", "p", "Adj.r", "f", "df1", "df2",
      paste("i", IPC_names, "-?", sep = "")
    )
    CV_Log <- matrix(0, nrow = (repetitions + 1), ncol = length(Mixed_headers))
    colnames(CV_Log) <- Mixed_headers
  }

  mci_reps <- seq(1, repetitions)
  for (mci in mci_reps) {
    BIC_output_list <- run_BIC(
      mci, num_IPC, IPC_names, train_n, test_n, num_Comparisons, ROI_Data, CV_Log, ROI,
      testIPCs, trainIPCs, DataType2, Analysis, ROI_num, Summary, SampleSize, output_path
    )
    Summary <- BIC_output_list[[1]]
    Best_Output <- BIC_output_list[[2]]
    All_Paths <- BIC_output_list[[3]]
    CV_Log <- BIC_output_list[[4]]
  }

  # update Summary table
  if (Analysis == "CV") {
    for (cols in 1:(2 * length(ROIs) + 9)) {
      Summary[(ROI_num), (cols + 2)] <- mean(as.numeric(CV_Log[1:(repetitions), (cols)]))
    }
    write.table(CV_Log, file.path(output_path, paste(ROI, Analysis, "_", DataType, "_log.csv", sep = "")),
      row.names = FALSE, sep = ","
    )
  } else if (Analysis == "Individual") {
    for (cols in 1:(2 * length(ROIs) + 8)) {
      Summary[(ROI_num + 1), (cols + 1)] <- mean(as.numeric(IND_Log[2:(repetitions + 1), (cols)]))
    }
    write.table(IND_Log, file.path(output_path, paste(ROI, Analysis, "_", DataType, "_byINDV.csv", sep = "")),
      row.names = FALSE, col.names = FALSE, sep = ","
    )
  } else if (Analysis == "Mixed") {
    for (cols in 1:47) {
      Summary[(ROI_num + 1), (cols + 1)] <- mean(as.numeric(IND_Log[2:(repetitions + 1), (cols)]))
    }
    write.table(IND_Log, file.path(output_path, paste(ROI, Analysis, "_", DataType, "_byINDV.csv", sep = "")),
      row.names = FALSE, col.names = FALSE, sep = ","
    )
  } else {
    BIC_Chain <- matrix(, nrow = 12, ncol = num_IPC)
    Chain_Header <- IPC_names
    for (heads in 1:num_IPC) {
      BIC_Chain[1, heads] <- Chain_Header[heads]
      BIC_Chain[2:12, heads] <- Best_Output[1:11, heads]
    }
    write.table(BIC_Chain, file.path(
      output_path,
      paste(DataType2, "BIC_CHAIN_", ROI, Analysis, ".csv", sep = "")
    ), row.names = FALSE, col.names = FALSE, sep = ",")
    write.table(All_Paths, file.path(
      output_path,
      paste(DataType2, "All_Paths", ROI, Analysis, ".csv", sep = "")
    ), row.names = FALSE, col.names = FALSE, sep = ",")
  }
  # close ROI loop
  print("ROI Complete")
  return(Summary)
}

#' Returns Bayesian Information Criterion (BIC) information comparing participant data for that ROI to the input IPCs,
#' in order to determine which IPCs best fit the data. This function is run in a loop by `get_data_by_ROI()`.
#' @return A list of outputs containing information about the IPC fits in the dataset, as follows:
#' * `Summary`: A summary table of the BIC outputs, to be used in the output tables.
#' * `Best_Output`: The IPC that best fits the data, according to BIC.
#' * `All_Paths`: A list of the paths that were searched in order to find the IPC that best fits the data.
#' * `CV_Log`: If cross-validation was conducted, an output of each iteration of cross-validation in which the IPC that best
#' fits the data was found.
#' @export
run_BIC <- function(mci, num_IPC, IPC_names, train_n, test_n, num_Comparisons, ROI_Data, CV_Log, ROI,
                    testIPCs, trainIPCs, DataType2, Analysis, ROI_num, Summary, SampleSize, output_path) {
  BICS_Output1 <- matrix(0, nrow = num_IPC, ncol = (num_IPC + 1))
  num_chains <- 0
  All_Paths <- matrix(0, nrow = 1000, ncol = num_IPC)

  usROI_Train <- matrix(, nrow = train_n, ncol = num_Comparisons)
  usROI_Test <- matrix(, nrow = test_n, ncol = num_Comparisons)

  # assign sampling order
  if (Analysis == "CV") {
    subset_sample <- sample(1:SampleSize, SampleSize, replace = FALSE)
  } else {
    subset_sample <- matrix(1:SampleSize)
  }

  # divide and prep data for train and test
  if (Analysis == "Group") {
    usROI_Train <- unlist(list(ROI_Data))
    usROI_Test <- unlist(list(ROI_Data))
  } else if (Analysis == "Mixed") {
    usROI_Train <- unlist(list(ROI_Data))
    usROI_Test <- unlist(list(ROI_Data[mci, 1:num_Comparisons]))
  } else if (Analysis == "Individual") {
    usROI_Train <- unlist(list(ROI_Data[mci, 1:num_Comparisons]))
    usROI_Test <- unlist(list(ROI_Data[mci, 1:num_Comparisons]))
  } else if (Analysis == "CV") {
    TestIndex <- (train_n + 1)
    for (rowi in TestIndex:SampleSize) {
      if (DataType == "US") {
        testsubj <- subset_sample[rowi]
        Trow <- (rowi - train_n)
        LeftOut <- ((2 * length(ROIs) + 9) + Trow)
        CV_Log[(mci), LeftOut] <- testsubj
        for (coli in 1:num_Comparisons) {
          usROI_Test[Trow, coli] <- ROI_Data[testsubj, coli]
        }
      }

      for (rowi in 1:train_n) {
        trainsubj <- subset_sample[rowi]
        for (cols in 1:num_Comparisons) {
          usROI_Train[rowi, cols] <- ROI_Data[trainsubj, cols]
        }
      }
    }
    usROI_Test <- unlist(list(usROI_Test))
  }

  # set train datasets
  if (DataType != "CS") {
    TrainSet <- matrix(usROI_Train, nrow = (train_n * num_Comparisons), ncol = 1)
    BIC_Type <- 1
  } else {
    TrainSet <- matrix(, nrow = (train_n * num_Comparisons), ncol = 3)
    TrainSet[, 1] <- eROI_Train
    TrainSet[, 2] <- mROI_Train
    TrainSet[, 3] <- lROI_Train
    BIC_Type <- 3
  }

  # set number of BIC identification datasets
  for (BIC_ID in 1:BIC_Type) {
    ROI_Train <- TrainSet[, BIC_ID]

    # level 1 BIC-GFBS
    BIC_list <- vector(mode = "list", length = num_IPC)
    for (IPC in 1:num_IPC) {
      BIC_list[[IPC]] <- BIC(lm(ROI_Train ~ unlist(trainIPCs[[IPC]])))
    }

    BIC_list <- unlist(BIC_list)
    # track and identify best 1st score, note it in BICS_Output1 and as variable "Best1"
    for (IPC in 1:num_IPC) {
      BICS_Output1[1, 1:num_IPC] <- BIC_list
    }

    BICS_Output1[1, (num_IPC + 1)] <- which.min(BICS_Output1[1, 1:(num_IPC)])
    Best1 <- BICS_Output1[1, (num_IPC + 1)]

    Currentlevel <- 0
    loop_end <- 0

    # Initially set all levels as 1
    Lv_list <- vector(mode = "list", length = num_IPC)
    test_Lv_list <- vector(mode = "list", length = num_IPC)

    for (level in 1:num_IPC) {
      Lv_list[[level]] <- rep(1, length(trainIPCs[[level]]))
      test_Lv_list[[level]] <- rep(1, length(testIPCs[[level]]))
    }

    # Begin loop of component comparisons: LV2 and >
    for (idx_IPC in 2:num_IPC) {
      if (loop_end == 0) {
        output_list <- compare_IPC(
          idx_IPC, num_IPC, Currentlevel, IPC_names, Best1, ROI_Train, DataType2, loop_end, trainIPCs, testIPCs,
          Lv_list, test_Lv_list, BICS_Output1, All_Paths, ROI, Analysis, num_chains, ThisPath
        )
        loop_end <- output_list[[1]]
        Lv_list <- output_list[[2]]
        test_Lv_list <- output_list[[3]]
        BICS_Output1 <- output_list[[4]]
        Best1 <- output_list[[5]]
        All_Paths <- output_list[[6]]
        num_chains <- output_list[[7]]
        ThisPath <- output_list[[8]]
        BestBic <- output_list[[9]]
      }
    }

    Best_Output <- BICS_Output1

    # Initialize variables
    checking_levels <- seq(1, 6)
    Splits <- 0
    Prior_Best_IPCs <- matrix(0, ncol = num_IPC, nrow = 1)
    Do_Check <- 1
    Level_Paths <- matrix(0, nrow = num_IPC, ncol = 1)
    model_improved <- 1

    # Check levels
    num_chains <- check_levels(
      checking_levels, BICS_Output1, Prior_Best_IPCs, Do_Check, Level_Paths, num_IPC, IPC_names, IPC_prior_level,
      ROI_Train, trainIPCs, Lv_list, DataType2, model_improved, num_chains, All_Paths, ThisPath, BestBic, output_path
    )

    # reset level to final IPCs
    count <- 1
    checkduplicate <- count
    while (Best_Output[checkduplicate + 1] != 0) {
      Bestprior <- Best_Output[count, (num_IPC + 1)]
      IPC_prior_level[Bestprior] <- trainIPCs[[count]]
      if (Analysis == "Mixed") {
        IPC_prior_level[Bestprior] <- testIPCs[[count]]
      }
      count <- count + 1
      checkduplicate <- count
      if ((Best_Output[count, (Best_Output[count, (num_IPC + 1)])]) > (Best_Output[(count - 1), (Best_Output[(count - 1), (num_IPC + 1)])] - 2)) {
        count <- (count - 1)
        break
      }
    }

    # run regression on identified scores (Component and Lvs correspondence noted in Best_Output table)
    if (DataType == "US") {
      ROI_Train_reg <- run_ROI_lm(ROI_Train, Lv_list)
      if (Analysis == "Mixed") {
        ROI_Test_reg <- run_ROI_lm(usROI_Test, test_Lv_list)
      }
      # else if (DataType == "Mixed") {
      #   ROI_Train_reg <- run_ROI_lm(ROI_Train, Lv_list)
      #   eROI_Train_reg <- run_ROI_lm(eROI_Test, Lv_list)
      #   mROI_Train_reg <- run_ROI_lm(mROI_Test, Lv_list)
      #   iROI_Train_reg <- run_ROI_lm(iROI_Test, Lv_list)
      # }
    }
  }
  TE_Ones <- (Lv_list[[1]])[1:length(usROI_Test)]
  # set up reconstruction
  ROI_reconstruction <- TE_Ones * summary(ROI_Train_reg)$coefficients[1, 1]

  # begin building the reconstructed matrix by level
  for (models in 1:count) {

    # number of useful levels = total_BIC_iteration -1 (final iteration is a redundant one; iBIC = iteration of BIC)
    iBIC <- models
    models <- models + 1

    Contrib_Model <- unlist(testIPCs[Best_Output[iBIC, (num_IPC + 1)]])


    if (Analysis == "CV") {
      CV_Log[(mci), Best_Output[iBIC, (num_IPC + 1)]] <- 1
      CV_Log[(mci), (Best_Output[iBIC, (num_IPC + 1)] + num_IPC)] <- summary(ROI_Train_reg)$coefficients[models, 1]
      CV_Log[(mci), (2 * length(ROIs) + 9)] <- num_chains
    } else if (Analysis == "Group") {
      Summary[(ROI_num + 1), (Best_Output[iBIC, (num_IPC + 1)] + 2)] <- 1
      Summary[(ROI_num + 1), (Best_Output[iBIC, (num_IPC + 1)] + (num_IPC + 2))] <- summary(ROI_Train_reg)$coefficients[models, 1]
      Summary[(mci + 1), (2 * length(ROIs) + 10)] <- num_chains
    } else if (Analysis != "Group") {
      IND_Log[(mci + 1), Best_Output[iBIC, (num_IPC + 1)]] <- 1
      IND_Log[(mci + 1), (Best_Output[iBIC, (num_IPC + 1)] + num_IPC)] <- summary(ROI_Train_reg)$coefficients[models, 1]
      if (Analysis == "Mixed") {
        IND_Log[(mci + 1), (Best_Output[iBIC, (num_IPC + 1)] + (2 * length(ROIs) + 8))] <- summary(ROI_Test_reg)$coefficients[models, 1]
      }
    }

    # update the reconstruction to include each identified model multiplied by its beta coefficient
    # (this is the end of the building loop)
    ROI_reconstruction <- ROI_reconstruction + ((Contrib_Model) * summary(ROI_Train_reg)$coefficients[models, 1])
  }

  # format reconstructed data for regression
  ROI_reconstruction <- unlist(list(ROI_reconstruction))
  ROI_Test_reg <- (lm(usROI_Test ~ ROI_reconstruction))

  if (Analysis == "CV") {
    CV_Log[(mci), (2 * length(IPC_names) + 1):(2 * length(IPC_names) + 4)] <- summary(ROI_Test_reg)$coefficients[2, 1:4]
    CV_Log[(mci), (2 * length(IPC_names) + 5)] <- summary(ROI_Test_reg)$adj.r.squared
    CV_Log[(mci), (2 * length(IPC_names) + 6):(2 * length(IPC_names) + 8)] <- summary(ROI_Test_reg)$f[1:3]
    CV_Log[(mci), (2 * length(IPC_names) + 9)] <- num_chains
  } else if (Analysis != "Group") {
    IND_Log[(mci + 1), (2 * length(ROIs) + 1):(2 * length(ROIs) + 4)] <- summary(ROI_Test_reg)$coefficients[2, 1:4]
    IND_Log[(mci + 1), (2 * length(ROIs) + 5)] <- summary(ROI_Test_reg)$adj.r.squared
    IND_Log[(mci + 1), (2 * length(ROIs) + 6):(2 * length(ROIs) + 8)] <- summary(ROI_Test_reg)$f[1:3]
  } else {
    Summary[(ROI_num + 1), (2 * length(IPC_names) + 3):(2 * length(IPC_names) + 6)] <- summary(ROI_Train_reg)$coefficients[2, 1:4]
    Summary[(ROI_num + 1), (2 * length(IPC_names) + 7)] <- summary(ROI_Train_reg)$adj.r.squared
    Summary[(ROI_num + 1), (2 * length(IPC_names) + 8):(2 * length(IPC_names) + 10)] <- summary(ROI_Train_reg)$f[1:3]
    Summary[(ROI_num + 1), (2 * length(IPC_names) + 11)] <- num_chains
  }

  # write all paths for iteration if group
  All_Paths[All_Paths == 0] <- NA
  writepaths <- All_Paths[complete.cases(All_Paths[, 1])]

  # end repetition loop
  return(list(Summary, Best_Output, All_Paths, CV_Log))
}
