# eval BIC paths
#' A series of functions to compare IPCs and determine which IPCs provide the best fit to the data. All of these
#' functions are used within the loop of `run_BIC()`.

#' Tracks, identifies, and records the best-matching IPCs, along with associated statistical information, for output.
#' @returns A list of outputs to facilitate the updating of IPCs and associated BIC scores, as follows:
#' * `loop_end`: Set to 1 if the current best BIC score matches the previous best (i.e. no more improvement in the model);
#' defaults to 0 during run of loop.
#' * `Lv_list`: The list of best-matching IPCs for training data for a given ROI, updated every loop. IPCs that have not been
#' identified in the ROI are assigned a matrix of 1s in this list.
#' * `test_Lv_list`: Same as `Lv_list`, except for test data (especially for use in cross-validation).
#' * `BICS_Output_1`: The row of BIC scores for each IPCs, for output to the Summary table.
#' * `Best1`: The index of the best-matching IPC for a given ROI.
#' * `All_Paths`: A list of the paths that were searched in order to find the IPC that best fits the data.
#' * `num_chains`: The number of paths that were searched in order to find the IPC that best fits the data. The lower this
#' number, the fewer paths it took to reach the best-matching IPC (which may suggest that the IPC contributed more strongly
#' to stimulus representations in that ROI).
#' * `This_Path`: The current path that was searched in order to find the IPC that best fits the data.
#' * `BestBic`: The BIC score for the best-matching IPC for a given ROI.
#' @export
compare_IPC <- function(i, num_IPC, CurrentLevel, IPC_names, Best1, ROI_Train, DataType2, loop_end, trainIPCs, testIPCs,
                        Lv_list, test_Lv_list, BICS_Output1, All_Paths, ROI, Analysis, num_chains, ThisPath) {
  # identify previous best component and set it as Lv (the rest, if not set on an earlier iteration, will remain as 1s)
  Lv_list[[i]] <- unlist(trainIPCs[Best1])
  test_Lv_list[[i]] <- unlist(testIPCs[Best1])

  # Run BIC component identification will all previously identified component +
  # all component individually (doubled ones do not change the BIC score)

  # track and identify best 1st score, note it in BICS_Output1 at correct iteration and as variable "Best1"
  BICS_Output1[i, 1:num_IPC] <- unlist(lapply(trainIPCs, run_ipc_lm, ROI_Train, Lv_list))
  BICS_Output1[i, (num_IPC + 1)] <- which.min(BICS_Output1[i, 1:num_IPC])
  Best1 <- BICS_Output1[i, (num_IPC + 1)]

  # Compare current previous best. if matching, exit function.
  check <- i - 1
  BestPath <- matrix(0, nrow = 1, ncol = num_IPC)
  ThisPath <- matrix(0, nrow = 1, ncol = num_IPC)
  Chain_Header <- c(IPC_names, "Best")

  # Initialize BestBic value in case it's not defined
  BestBic <- 0
  if ((BICS_Output1[check, BICS_Output1[check, (num_IPC + 1)]] - 2) < BICS_Output1[i, Best1]) {
    num_chains <- 1

    for (Num_IPCs in 1:(i - 1)) {
      BestPath[Num_IPCs] <- BICS_Output1[Num_IPCs, (num_IPC + 1)]
    }
    All_Paths[num_chains, 1:num_IPC] <- BestPath
    BestBic <- BICS_Output1[(Num_IPCs), BICS_Output1[(Num_IPCs), (num_IPC + 1)]]
    Best_Output <- BICS_Output1

    # write bic table
    if (Analysis == "Group") {
      BIC_Chain <- matrix(, nrow = 12, ncol = (num_IPC + 1))
      for (heads in 1:(num_IPC + 1)) {
        BIC_Chain[1, heads] <- Chain_Header[heads]
        BIC_Chain[2:12, heads] <- BICS_Output1[1:11, heads]
      }

      write.table(BIC_Chain, file.path(output_path, paste(DataType2, num_chains, "_BIC_CHAIN_", ROI, Analysis, ".csv", sep = "")),
        row.names = FALSE, col.names = FALSE, sep = ","
      )
    }

    # leave BIC level loop
    loop_end <- 1
  }
  return(list(loop_end, Lv_list, test_Lv_list, BICS_Output1, Best1, All_Paths, num_chains, ThisPath, BestBic, Chain_Header))
}

#' Runs the loop that checks whether additional paths are needed to find the best IPC fit. This function is run in a loop
#' by `compare_IPC()`.
#' @return `num_chains`: The number of paths that were searched in order to find the IPC that best fits the data. The lower this
#' number, the fewer paths it took to reach the best-matching IPC (which may suggest that the IPC contributed more strongly
#' to stimulus representations in that ROI).
#' @export
check_levels <- function(checking_levels, BICS_Output1, Prior_Best_IPCs, Do_Check, Level_Paths, num_IPC, IPC_names, IPC_prior_level,
                         ROI_Train, trainIPCs, Lv_list, DataType2, model_improved, num_chains, All_Paths, ThisPath, BestBic, output_path,
                         Chain_Header, ROI, Analysis) {
  for (checking_level in checking_levels) {
    # begin checking
    BICS_Output_alternate <- matrix(0, nrow = num_IPC, ncol = (num_IPC + 1))
    Splits <- 0
    if (Analysis == "Group") {
      print(paste("Checking level", checking_level, sep = " "))
    }
    Level2Check <- BICS_Output1[checking_level, 1:(num_IPC + 1)]
    for (previous_bests in 1:checking_level) {
      Prior_Best_IPCs[previous_bests] <- BICS_Output1[previous_bests, (num_IPC + 1)]
    }
    if (BICS_Output1[checking_level, (num_IPC + 1)] == 0) {
      break
    }
    follow_paths_output <- follow_paths(
      checking_level, Level2Check, Splits, BICS_Output1, BICS_Output_alternate, Prior_Best_IPCs,
      Do_Check, Level_Paths, num_IPC, IPC_names, IPC_prior_level,
      ROI_Train, trainIPCs, Lv_list, DataType2, model_improved, num_chains, All_Paths, ThisPath,
      BestBic, output_path, Chain_Header, ROI, Analysis
    )
    model_improved <- follow_paths_output[[1]]
    num_chains <- follow_paths_output[[2]]
    if (model_improved == 0) {
      print(paste("Checking level ", checking_level, ": Model did not improve!", sep = ''))
      return(num_chains)
      break
    }
  }
  return(num_chains)
}

#' Checks individual levels of the IPC search path, and determines whether the model has improved or not.
#' @returns A list of outputs indicating the results of the level check, as follows:
#' * `model_improved`: Has a value of 1 if the model moved after the new path search (which continues the loop). If
#' the model did not improve after the last search, this has a value of 0, which breaks the search loop and outputs
#' the results of the search to a series of tables.
#' * `num_chains`: The number of paths that were searched in order to find the IPC that best fits the data. The lower this
#' number, the fewer paths it took to reach the best-matching IPC (which may suggest that the IPC contributed more strongly
#' to stimulus representations in that ROI).
#' @export
follow_paths <- function(checking_level, Level2Check, Splits, BICS_Output1, BICS_Output_alternate, Prior_Best_IPCs,
                         Do_Check, Level_Paths, num_IPC, IPC_names, IPC_prior_level,
                         ROI_Train, trainIPCs, Lv_list, DataType2, model_improved, num_chains, All_Paths, ThisPath,
                         BestBic, output_path, Chain_Header, ROI, Analysis) {
  # begin specific level check
  # check if additional paths needed, and fill a list of each
  Level2Check[Level2Check[(num_IPC + 1)]] <- Level2Check[Level2Check[(num_IPC + 1)]] + 2
  if (which.min(Level2Check[1:num_IPC]) %in% Level_Paths | which.min(Level2Check[1:num_IPC]) %in% Prior_Best_IPCs) {
    checking_level <- checking_level + 1
    Level_Paths <- matrix(0, nrow = 1, ncol = num_IPC)
  } else {
    More_Paths <- 1
    Level_Paths <- Prior_Best_IPCs
    while (More_Paths == 1) {
      Splits <- Splits + 1
      Level_Paths[Splits + checking_level] <- which.min(Level2Check[1:num_IPC])
      Level2Check[Level2Check[(num_IPC + 1)]] <- Level2Check[Level2Check[(num_IPC + 1)]] + 2
      Level2Check[(num_IPC + 1)] <- which.min(Level2Check[1:num_IPC])
      if (which.min(Level2Check[1:num_IPC]) %in% Level_Paths) {
        More_Paths <- 0
      }
    }
    # follow each new path
    sC <- 1
    while (Level_Paths[checking_level + sC] != 0) {
      if (Analysis == "Group") {
        print(paste("Checking level", checking_level, "Subpath ", sC, sep = " "))
      }
      check_path <- (sC + checking_level)
      sC <- sC + 1
      subSplits <- 0
      subLevel_Paths <- matrix(0, nrow = 1, ncol = num_IPC)

      # Set up priors as empty
      for (Levels in (checking_level):num_IPC) {
        assign(paste("check_LV", Levels, sep = ""), "TR_Ones")
      }

      # Fill priors with known data
      Prior_Best_IPCs[checking_level] <- Level_Paths[(check_path)]
      BICS_Output_alternate[checking_level, 1:num_IPC] <- BICS_Output1[checking_level, 1:num_IPC]
      BICS_Output_alternate[checking_level, (num_IPC + 1)] <- Level_Paths[check_path]
      if (checking_level > 1) {
        BICS_Output_alternate[1:(checking_level - 1), 1:(num_IPC + 1)] <- BICS_Output1[1:(checking_level - 1), 1:(num_IPC + 1)]
      }

      for (buildlevel in 1:(checking_level)) {
        Bestprior <- BICS_Output_alternate[buildlevel, (num_IPC + 1)]
        IPC_prior_level[Bestprior] <- paste("check_LV", buildlevel, sep = "")
      }
      BICS_Output_alternate[(checking_level + 1), 1:(num_IPC)] <- unlist(lapply(trainIPCs, run_ipc_lm, ROI_Train, Lv_list))
      BICS_Output_alternate[(checking_level + 1), (num_IPC + 1)] <- which.min(BICS_Output_alternate[(checking_level + 1), 1:num_IPC])

      Level2Check <- BICS_Output_alternate[(checking_level + 1), 1:(num_IPC + 1)]

      # determine if model improved
      if ((BICS_Output_alternate[
        checking_level,
        BICS_Output_alternate[checking_level, (num_IPC + 1)]
      ] - 2) < (BICS_Output_alternate[(checking_level + 1), BICS_Output_alternate[(checking_level + 1), (num_IPC + 1)]])) {
        model_improved <- 0
        num_chains <- num_chains + 1
        for (Num_IPCs in 1:(checking_level + 1)) {
          ThisPath[Num_IPCs] <- BICS_Output_alternate[Num_IPCs, (num_IPC + 1)]
        }
        All_Paths[num_chains, 1:num_IPC] <- ThisPath
        if (BestBic > BICS_Output_alternate[checking_level, BICS_Output_alternate[checking_level, (num_IPC + 1)]]) {
          BestBic <- BICS_Output_alternate[(checking_level), BICS_Output_alternate[checking_level, (num_IPC + 1)]]
          BestPath <- ThisPath
          Best_Output <- BICS_Output_alternate
        }
        if (Analysis == "Group") {
          BIC_Chain <- matrix(, nrow = 12, ncol = (num_IPC + 1))
          for (heads in 1:(num_IPC + 1)) {
            BIC_Chain[1, heads] <- Chain_Header[heads]
            BIC_Chain[2:12, heads] <- BICS_Output_alternate[1:11, heads]
          }
          if (DataType != "CS") {
            write.table(BIC_Chain, paste(DataType2, num_chains, "_BIC_CHAIN_", ROI, Analysis, ".csv", sep = ""), row.names = FALSE, col.names = FALSE, sep = ",")
          } else {
            if (BIC_ID == 1) {
              write.table(BIC_Chain, file.path(
                output_path,
                paste("e", DataType2, num_chains, "_BIC_CHAIN_", ROI, Analysis, ".csv", sep = "")
              ), row.names = FALSE, col.names = FALSE, sep = ",")
            } else if (BIC_ID == 2) {
              write.table(BIC_Chain, file.path(
                output_path,
                paste("m", DataType2, num_chains, "_BIC_CHAIN_", ROI, Analysis, ".csv", sep = "")
              ), row.names = FALSE, col.names = FALSE, sep = ",")
            } else if (BIC_ID == 3) {
              write.table(BIC_Chain, file.path(
                output_path,
                paste("l", DataType2, num_chains, "_BIC_CHAIN_", ROI, Analysis, ".csv", sep = "")
              ), row.names = FALSE, col.names = FALSE, sep = ",")
            }
          }
        }
      } else {
        model_improved <- 1
      }
    }
  }
  return(list(model_improved, num_chains))
}
