library(logger)

behavioral_main <- c("outcome", "PE", "absPE", "expectation")
behavioral_support <- c("trialInBlock", "trialInMiniBlock", "learning_miniblock", "half", "choiceApproach", "abs_expectation", "RT")
cardiac_main <- c("phaseFeedb", "sysDiaFeedb", "changeIBI4_3", "changeIBI5_4", "changeIBI6_5", "changeIBI7_6")
cardiac_support <- c("IBI4", "IBI5", "IBI6", "IBI7", "msDistToFeedb") # strict cardiac clean, removed "IBI3" from here, moved to more_vars2
cardiac_support_more_vars <- c("RestRMSSD", "RestAveIBI") # previously ran with others here (for other R, but here i would like to redo it with HRV outlier removal only but also without SDNN, only RMSSD) "IBI1", "IBI2", "IBI8", "msDistToFeedb", "changeIBI2_1", "changeIBI3_2", "changeIBI8_7", "RestSDNN". These are the settings until 24.02 # also ! included in very strict cardiac clean, along with the previous one
cardiac_support_more_vars2 <- c("IBI3","beforeFeedbAveIBI", "beforeFeedbRMSSD", "beforeFeedbSDNN", "beforeStimAveIBI", "beforeStimRMSSD", "beforeStimSDNN", "afterStimAveIBI", "afterStimRMSSD", "afterStimSDNN", "afterFeedbAveIBI", "afterFeedbRMSSD", "afterFeedbSDNN") # also moved "RestAveIBI" to cardiac_support_more_vars from here
neural_main <- c("Pz__p3b_Win", "Pz__p3b_Lose") # could add p3a here if we include it , "Pz__p3b"
neural_support <- c("Pz__p3b__withBase", "Pz__p3b_Win__withBase", "Pz__p3b_Lose__withBase", "Pz__baseline", "FCz__frn", "FCz__p3a", "FCz__baseline", "EKG1__p3b", "EKG1__p3a", "EKG1__p3b_Win", "EKG1__p3b_Lose", "EKG1__baseline", "EOGV__p3b_Win", "EOGV__p3b_Lose", "EOGH__p3b_Win", "EOGH__p3b_Lose", "EOGV__p3b", "EOGH__p3b", "FCz__frn__withBase", "FCz__p3a__withBase")
factors <- c("learning_miniblock", "half", "choiceApproach", "outcome", "sysDiaFeedb", "subID", "absPE_2groups","expectation_3groups", "first_next","first_last", "trueCorrect_nextTrial", "switchStay_nextTrialAction", "switchStay_sameCueAction", "switchStay_nextTrialStrategy")
tfr_main <- c("Pz__high_alpha__baseline__TFR", "Pz__low_alpha__baseline__TFR")
tfr_support <- c("Pz__high_alpha__p3b_Win__TFR", "Pz__high_alpha__p3b_Lose__TFR", "Pz__low_alpha__p3b_Win__TFR", "Pz__low_alpha__p3b_Lose__TFR", "Pz__high_alpha__baseline_tfr__TFR", "Pz__low_alpha__baseline_tfr__TFR", "Pz__high_alpha__p3b_Win__TFR__withBase", "Pz__high_alpha__p3b_Lose__TFR__withBase", "Pz__low_alpha__p3b_Win__TFR__withBase", "Pz__low_alpha__p3b_Lose__TFR__withBase", "Pz__high_alpha__p3b__TFR", "Pz__low_alpha__p3b__TFR", "Pz__high_alpha__p3b__TFR__withBase", "Pz__low_alpha__p3b__TFR__withBase", "Pz__high_alpha__baseline__TFR", "Pz__low_alpha__baseline__TFR", "Pz__high_alpha__baseline__TFR__withBase", "Pz__low_alpha__baseline__TFR__withBase")
output_vars <- c("changeIBI4_3", "changeIBI5_4", "changeIBI6_5", "changeIBI7_6","Pz__p3b_Win", "Pz__p3b_Lose", "Pz__p3b")
next_trial_vars <- c('Pz__baseline_next_time', 'FCz__baseline_next_time', 'Pz__high_alpha__baseline__TFR_next_time', 
                     'FCz__high_alpha__baseline__TFR_next_time', 'Pz__high_alpha__baseline__TFR__withBase_next_time', 
                     'FCz__high_alpha__baseline__TFR__withBase_next_time',
                     'Pz__p3b_Win_next_time','Pz__p3b_Lose_next_time', 'FCz__p3a_next_time')

scale_all_vars <- c(behavioral_main, behavioral_support, cardiac_main, cardiac_support, cardiac_support_more_vars, neural_main, neural_support, tfr_main, tfr_support, cardiac_support_more_vars2)

scale_vars <- function(vars2scale, df) {
  scale_df <- data.frame(matrix(nrow = nrow(df), ncol = 0))  
  for (i in vars2scale) {
    if (is.numeric(df[[i]])) {  # Check if the column is numeric
      scale_df[[paste0("scale_", i)]] <- c(scale(df[[i]])) 
    }
  }
  return(scale_df)  
}

remove_outliers <- function(df, outlier_vars, stats_vars = NULL, output_values = FALSE, 
                            output_drop_results = FALSE, output_drop_results_per_subject = FALSE) {
  # Ensure the selected columns exist and are numeric
  outlier_cols <- outlier_vars[sapply(df[outlier_vars], is.numeric)]
  scaled_cols <- paste0("scale_", outlier_cols)  # Use scaled columns for outlier detection
  
  # Default stats_vars to outlier_vars if not provided
  if (is.null(stats_vars)) {
    stats_vars <- outlier_vars
  }
  
  stats_cols <- stats_vars[sapply(df[stats_vars], is.numeric)]  # Ensure stats_vars are numeric
  
  # Identify rows without outliers (values within ±3 SD)
  outlier_mask <- rowSums(abs(df[, scaled_cols]) > 3) == 0
  df_filtered <- df[outlier_mask, ]
  
  # Initialize output lists
  drop_results <- NULL
  drop_results_per_subject <- NULL
  
  # Calculate number of outliers removed
  if (output_drop_results) {
    total_rows <- nrow(df)
    removed_rows <- total_rows - sum(outlier_mask)  # Number of rows removed
    
    drop_results <- data.frame(
      variable = c(outlier_cols, "overall"),  # Use outlier variable names without "scale_"
      num_removed = c(colSums(abs(df[, scaled_cols]) > 3), removed_rows),
      percent_removed = c(colSums(abs(df[, scaled_cols]) > 3) / total_rows * 100, removed_rows / total_rows * 100)
    )
    
    rownames(drop_results) <- NULL  # Remove row names to prevent confusion
  }
  
  # Calculate per-subject outlier counts
  if (output_drop_results_per_subject && ("subID" %in% colnames(df))) {
    unique_subjects <- unique(df$subID)  # Get unique subject IDs
    
    # Initialize an empty matrix to store per-subject outlier counts
    drop_results_per_subject <- matrix(0, nrow = length(outlier_cols) + 1, 
                                       ncol = length(unique_subjects) + 1)
    rownames(drop_results_per_subject) <- c(outlier_cols, "overall")
    colnames(drop_results_per_subject) <- c(as.character(unique_subjects), "overall")
    
    for (sub in unique_subjects) {
      sub_rows <- df$subID == sub  # Mask for subject's rows
      sub_outlier_mask <- !outlier_mask[sub_rows]  # Find dropped rows per subject
      
      # Count per variable
      for (i in seq_along(outlier_cols)) {
        drop_results_per_subject[i, sub] <- sum(abs(df[sub_rows, scaled_cols[i]]) > 3, na.rm = TRUE)
      }
      
      # Count total dropped rows for this subject
      drop_results_per_subject["overall", sub] <- sum(sub_outlier_mask, na.rm = TRUE)
    }
    
    # Compute total per variable across subjects
    drop_results_per_subject[, "overall"] <- rowSums(drop_results_per_subject, na.rm = TRUE)
    
    # Convert matrix to dataframe
    drop_results_per_subject <- as.data.frame(drop_results_per_subject)
  }
  
  # Always check if any subject has fewer than 150 remaining observations
  if ("subID" %in% colnames(df)) {
    remaining_counts <- table(df_filtered$subID)  # Count remaining rows per subject
    subjects_with_few_observations <- names(remaining_counts[remaining_counts < 150])
    
    if (length(subjects_with_few_observations) > 0) {
      warning_message <- paste("Warning: Some subjects have fewer than 150 remaining observations:",
                               paste(subjects_with_few_observations, collapse = ", "))
      log_warn(warning_message)  # Log the warning
      warning(warning_message)   # Print the warning to the console
    }
  }
  
  # Compute old and new mean & SD for stats_cols
  stats_df <- NULL
  if (output_values) {
    stats_df <- data.frame(
      variable = stats_cols,  # Store variable names correctly
      old_mean = sapply(stats_cols, function(col) mean(df[[col]], na.rm = TRUE)),
      old_sd = sapply(stats_cols, function(col) sd(df[[col]], na.rm = TRUE)),
      new_mean = sapply(stats_cols, function(col) mean(df_filtered[[col]], na.rm = TRUE)),
      new_sd = sapply(stats_cols, function(col) sd(df_filtered[[col]], na.rm = TRUE))
    )
    stats_df$thresh_lower <- stats_df$old_mean - 3 * stats_df$old_sd
    stats_df$thresh_upper <- stats_df$old_mean + 3 * stats_df$old_sd
    rownames(stats_df) <- NULL
  }
  
  # Return results based on requested outputs
  result_list <- list(df_filtered = df_filtered)
  if (output_values) result_list$stats_df <- stats_df
  if (output_drop_results) result_list$drop_results <- drop_results
  if (output_drop_results_per_subject) result_list$drop_results_per_subject <- drop_results_per_subject
  
  return(result_list)
}

preprocess_df <- function(df, strict_cardiac_clean = FALSE, very_strict_cardiac_clean = FALSE, clean_neural = FALSE, clean_tfr = FALSE, output_cleaning_results = FALSE, rescale_after_outliers = TRUE, next_trial = FALSE) {
  
  df$choiceApproach <- ifelse(df$choiceApproach == 1, "Approach", "Avoid")
  
  df$Pz__p3b <- ifelse(df$outcome == "Win", df$Pz__p3b_Win, df$Pz__p3b_Lose)
  df$Pz__p3b__withBase <- ifelse(df$outcome == "Win", df$Pz__p3b_Win__withBase, df$Pz__p3b_Lose__withBase)
  df$EKG1__p3b <- ifelse(df$outcome == "Win", df$EKG1__p3b_Win, df$EKG1__p3b_Lose)
  df$EOGH__p3b <- ifelse(df$outcome == "Win", df$EOGH__p3b_Win, df$EOGH__p3b_Lose)
  df$EOGV__p3b <- ifelse(df$outcome == "Win", df$EOGV__p3b_Win, df$EOGV__p3b_Lose)
  df$Pz__high_alpha__p3b__TFR <- ifelse(df$outcome == "Win", 
                                        df$Pz__high_alpha__p3b_Win__TFR, 
                                        df$Pz__high_alpha__p3b_Lose__TFR)
  df$Pz__low_alpha__p3b__TFR <- ifelse(df$outcome == "Win", 
                                       df$Pz__low_alpha__p3b_Win__TFR, 
                                       df$Pz__low_alpha__p3b_Lose__TFR)
  df$Pz__high_alpha__p3b__TFR__withBase <- ifelse(df$outcome == "Win", 
                                                  df$Pz__high_alpha__p3b_Win__TFR__withBase, 
                                                  df$Pz__high_alpha__p3b_Lose__TFR__withBase)
  df$Pz__low_alpha__p3b__TFR__withBase <- ifelse(df$outcome == "Win", 
                                                 df$Pz__low_alpha__p3b_Win__TFR__withBase, 
                                                 df$Pz__low_alpha__p3b_Lose__TFR__withBase)
  
  df$abs_expectation <- abs(df$expectation)
  
  df[factors] <- lapply(df[factors], as.factor)
  df$sysDiaFeedb <- factor(df$sysDiaFeedb, levels = c("Systole", "Diastole"))
  df$first_next <- factor(df$first_next, levels = c("first", "next"))
  df$first_last <- factor(df$first_last, levels = c("first", "last"))
  df$outcome <- factor(df$outcome, levels = c("Win", "Lose"))
  df$learning_miniblock <- factor(df$learning_miniblock, levels = c("first", "next", "reversal", "last"), labels = c("First", "Next", "Reversal", "Last"))
  df$expectation_3groups <- factor(df$expectation_3groups, levels = c("high", "medium", "low"), labels = c("High", "Medium", "Low"))
  df$absPE_2groups <- factor(df$absPE_2groups, levels = c("high", "low"), labels = c("High", "Low"))
  
  df$changeIBI4_3_2groups <- factor(df$changeIBI4_3_2groups, levels = c("high", "low"), labels = c("High", "Low"))
  df$changeIBI5_4_2groups <- factor(df$changeIBI5_4_2groups, levels = c("high", "low"), labels = c("High", "Low"))
  df$changeIBI6_5_2groups <- factor(df$changeIBI6_5_2groups, levels = c("high", "low"), labels = c("High", "Low"))
  df$changeIBI7_6_2groups <- factor(df$changeIBI7_6_2groups, levels = c("high", "low"), labels = c("High", "Low"))
  
  perf_by_miniblock <- c("switchStay_nextTrialAction", 
                         "switchNeg_nextTrialAction",
                         "stayPos_nextTrialAction",
                         "switchStay_nextTrialStrategy", 
                         "switchNeg_nextTrialStrategy",
                         "stayPos_nextTrialStrategy",
                         "switchStay_sameCueAction", 
                         "switchNeg_sameCueAction",
                         "stayPos_sameCueAction",
                         "trueCorrect")
  # Step 1: Ensure all are factors with levels 0 and 1
  df <- df %>%
    mutate(across(all_of(perf_by_miniblock), ~ factor(.x, levels = c(0, 1))))
  # Step 2: Reshape to long format
  df_long <- df %>%
    select(subID, learning_miniblock, all_of(perf_by_miniblock)) %>%
    pivot_longer(cols = all_of(perf_by_miniblock), names_to = "variable", values_to = "value")
  # Step 3: Compute proportion of 1s per subID, miniblock, variable
  df_prop <- df_long %>%
    group_by(subID, learning_miniblock, variable) %>%
    summarise(proportion = mean(value == "1", na.rm = TRUE), .groups = "drop")
  # Step 4: Create new variable names
  df_prop <- df_prop %>%
    mutate(varname = paste0("percent_", variable, "_", learning_miniblock)) %>%
    select(subID, varname, proportion)
  # Step 5: Reshape to wide format
  df_wide <- df_prop %>%
    pivot_wider(names_from = varname, values_from = proportion)
  # Step 6: Join back to original df (if desired)
  df <- left_join(df, df_wide, by = "subID")
  
  scale_df <- scale_vars(scale_all_vars, df)
  df <- cbind(df, scale_df)
  
  if (next_trial){
    scale_df <- scale_vars(next_trial_vars, df)
    df <- cbind(df, scale_df)
  }
  
  check_vars <- cardiac_main
  
  if (strict_cardiac_clean) {
    check_vars <- c(check_vars, cardiac_support)
    if (very_strict_cardiac_clean) {
      check_vars <- c(check_vars, cardiac_support_more_vars)
    }
  } else {
    if (very_strict_cardiac_clean) {
      check_vars <- c(check_vars, cardiac_support, cardiac_support_more_vars)
    }
  }
  
  if (clean_neural) {
    check_vars <- c(check_vars, neural_main)
  }
  
  if (clean_tfr) {
    check_vars <- c(check_vars, tfr_main)
  }
  
  if (clean_tfr|clean_neural) {
    df <- df[df$bad_epoch!=1, ]
  }
  
  if (!output_cleaning_results) {
    df_clean <- remove_outliers(df, check_vars)
  } else {
    df_clean <- remove_outliers(df, check_vars, scale_all_vars, TRUE, TRUE, TRUE)
  }
  
  if (rescale_after_outliers) {
    scale_df <- scale_vars(scale_all_vars, df_clean$df_filtered)
    df_clean$df_filtered[names(scale_df)] <- scale_df
  }
  
  return(df_clean)
}

compare_drop_results_per_subject <- function(df, check_vars_list) {
  
  results_list <- list()
  differences_list <- list()
  
  for (i in seq_along(check_vars_list)) {
    check_vars <- check_vars_list[[i]]
    res <- remove_outliers(df, check_vars, 
                           output_values = FALSE, 
                           output_drop_results = FALSE, 
                           output_drop_results_per_subject = TRUE)
    
    if (!is.null(res$drop_results_per_subject) && nrow(res$drop_results_per_subject) > 0) {
      results_list[[paste0("check_", i)]] <- res$drop_results_per_subject
    } else {
      warning(paste("Cleaning", i, "did not produce a valid drop_results_per_subject table."))
    }
  }
  
  if (length(results_list) < 2) {
    stop("At least two valid check_vars sets are required for comparison.")
  }
  
  for (i in 1:(length(results_list) - 1)) {
    for (j in (i + 1):length(results_list)) {
      pair_name <- paste0("check_", i, "_check_", j)
      
      overall_i <- results_list[[paste0("check_", i)]]
      overall_j <- results_list[[paste0("check_", j)]]
      
      if (!("overall" %in% rownames(overall_i)) || !("overall" %in% rownames(overall_j))) {
        warning(paste("One of the cleaning results in pair", pair_name, "does not have an 'overall' row."))
        next
      }
      
      overall_row_i <- overall_i["overall", , drop = FALSE]
      overall_row_j <- overall_j["overall", , drop = FALSE]
      
      overall_i_numeric <- as.numeric(as.character(overall_row_i))
      overall_j_numeric <- as.numeric(as.character(overall_row_j))
      
      diff_numeric <- overall_i_numeric - overall_j_numeric
      
      diff_df <- as.data.frame(t(diff_numeric))
      rownames(diff_df) <- "overall"
      
      differences_list[[pair_name]] <- diff_df
    }
  }
  
  output <- list(
    results = results_list,
    differences = differences_list
  )
  
  return(output)
}