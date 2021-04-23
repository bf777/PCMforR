assign_IPCs <- function(num_IPC, IPC_names, All_IPCs, test_n, train_n, testIPCs, trainIPCs, num_Comparisons) {
  # Initialize empty lists for the output
  testIPCs_list <- vector(mode = "list", length = num_IPC)
  trainIPCs_list <- vector(mode = "list", length = num_IPC)
  # For each IPC, restructure the IPC such that it is stacked element-wise (i.e. repeat each value
  # by the number of comparisons).
  for (get_IPC in 1:num_IPC) {
    whichIPC <- IPC_names[get_IPC]
    testIPCs <- rep(All_IPCs[get_IPC, 1:num_Comparisons], each = test_n)
    trainIPCs <- rep(All_IPCs[get_IPC, 1:num_Comparisons], each = train_n)
    testIPCs_list[[get_IPC]] <- testIPCs
    trainIPCs_list[[get_IPC]] <- trainIPCs
  }
  # Returns a list containing a list of test IPCs and a list of train IPCs.
  return(list(testIPCs_list, trainIPCs_list))
}

run_ipc_lm <- function(IPC, ROI_Train, Lv_list) {
  IPC_df <- data.frame(ROI_Train = unlist(ROI_Train), IPC = unlist(IPC))
  IPC_df_colnames <- c("ROI_Train", "IPC")
  for (Lv_val in 1:length(Lv_list)) {
    IPC_df <- data.frame(cbind(IPC_df, unlist(Lv_list[[Lv_val]])))
    IPC_df_colnames <- append(IPC_df_colnames, paste("Lv", toString(Lv_val - 1), sep = ""))
  }
  colnames(IPC_df) <- IPC_df_colnames
  IPC_lm <- BIC(lm(ROI_Train ~ ., data = IPC_df))
  return(IPC_lm)
}

run_ROI_lm <- function(ROI, Lv_list) {
  IPC_df <- data.frame(ROI = unlist(ROI))
  ROI_df_colnames <- c("ROI")
  for (Lv_val in 1:length(Lv_list)) {
    IPC_df <- data.frame(cbind(IPC_df, unlist(Lv_list[[Lv_val]])))
    ROI_df_colnames <- append(ROI_df_colnames, paste("Lv", toString(Lv_val - 1), sep = ""))
  }
  colnames(IPC_df) <- ROI_df_colnames
  IPC_lm <- lm(ROI ~ ., data = IPC_df)
  return(IPC_lm)
}
