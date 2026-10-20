library(logger)
library(ggplot2)
library(moments)       # For skewness & kurtosis
library(nortest)       # For Anderson-Darling test
library(fitdistrplus)  # For distribution fitting
library(sn)            # For skew-normal distribution
library(goftest)  # For Cramér-von Mises test
library(diptest) # Hartigan's Dip test
library(reshape2)
library(RColorBrewer)

plot_distributions <- function(df, plot_vars, vis_path, plot_factors = FALSE, plot_width = 8, plot_height = 6, 
                               show_plots = FALSE, save_plots = TRUE, plot_format = "svg",
                               fill_color = "gray", show_title = TRUE, title_string = "", units_string = "", dpi = 150) {
  
  # Ensure only existing columns are used
  plot_vars <- plot_vars[plot_vars %in% names(df)]
  
  # Sort variables into numeric and factor types
  numeric_vars <- plot_vars[sapply(df[plot_vars], is.numeric)]
  factor_vars <- plot_vars[sapply(df[plot_vars], is.factor)]
  
  # Loop through numeric variables and plot histograms
  for (var in numeric_vars) {
    # Determine title and units
    final_title <- if (title_string == "") map_variable_name(var) else title_string
    final_x_label <- map_variable_name(var)
    final_y_label <- "Count"
    final_units <- if (units_string == "") map_variable_units(var) else units_string
    
    if (!is.na(final_units) && nzchar(final_units)) {
      if (grepl("mu", final_units)) {
        final_x_label <- as.expression(bquote(.(final_x_label) ~ "(" * .(parse(text = final_units)[[1]]) * ")"))
      } else {
        final_x_label <- paste0(final_x_label, " (", final_units, ")")
      }
    }
    
    # Create plot
    p <- ggplot(df, aes(.data[[var]])) +
      geom_histogram(binwidth = (max(df[[var]], na.rm = TRUE) - min(df[[var]], na.rm = TRUE)) / 30, 
                     fill = fill_color, color = "white") +
      labs(title = if (show_title) final_title else NULL, 
           x = final_x_label, y = final_y_label) 
    
    # Save or show plot
    if (save_plots) {
      file_name <- file.path(vis_path, paste0("distribution_", var, ".", plot_format))
      ggsave(file_name, p, width = plot_width, height = plot_height)
    }
    if (show_plots) {
      print(p)
    }
  }
  
  # Loop through factor variables if plot_factors is TRUE
  if (plot_factors) {
    for (var in factor_vars) {
      final_title <- if (title_string == "") map_variable_name(var) else title_string
      final_x_label <- map_variable_name(var)
      final_y_label <- "Count"
      
      p <- ggplot(df, aes(.data[[var]])) +
        geom_bar(fill = fill_color, color = "white") +
        labs(title = if (show_title) final_title else NULL, 
             x = final_x_label, y = final_y_label) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      # Save or show plot
      if (save_plots) {
        file_name <- file.path(vis_path, paste0("distribution_", var, ".", plot_format))
        ggsave(file_name, p, width = plot_width, height = plot_height, dpi = dpi)
      }
      if (show_plots) {
        print(p)
      }
      return(p)
    }
  }
  return(p)
}

map_variable_units <- function(variable_name) {
  # Default: Use map_variable_name(variable_name)
  units <- ""
  
  # If "scale_" is in the variable name, set units to "A.u."
  if (grepl("scale_", variable_name, fixed = TRUE)) {
    return("A.u.")
  }
  
  # If "msDistToFeedb" is in the variable name, set units to "ms"
  if (grepl("msDistToFeedb", variable_name, fixed = TRUE)) {
    return("ms")
  }
  
  # If "IBI" is in the variable name, set units to "ms"
  if (grepl("IBI", variable_name, fixed = TRUE)) {
    return("ms")
  }
  
  # List of EEG/EKG/EOG channels
  eeg_channels <- c("Pz", "FCz", "EKG1", "EOGV", "EOGH")
  
  # If an EEG channel, "TFR", AND "__withBase" are in the variable name, set units to "A.u."
  if (any(grepl(paste(eeg_channels, collapse = "|"), variable_name)) &&
      grepl("TFR", variable_name, fixed = TRUE) &&
      grepl("__withBase", variable_name, fixed = TRUE)) {
    return("A.u.")
  }
  
  # If an EEG channel AND "TFR" are in the variable name, set units to µV²
  if (any(grepl(paste(eeg_channels, collapse = "|"), variable_name)) &&
      grepl("TFR", variable_name, fixed = TRUE)) {
    #return(expression(mu ~ "V"^2))  # µV²
    return("mu * 'V'^2")
  }
  
  # If any of these channels are in the variable name, set units to "µV" (microvolts)
  if (any(grepl(paste(eeg_channels, collapse = "|"), variable_name))) {
    #return(expression(mu ~ "V"))  # µV
    return("mu * 'V'")
  }
  return(units)
}

suggest_brms_family <- function(variable, 
                                output_results = FALSE, 
                                limit_student_df = FALSE, 
                                df_threshold = 30, 
                                default_family = "student", 
                                verbose = FALSE) {
  
  # Custom Student-t distribution functions
  dstudent_t <- function(x, mu, sigma, df) {
    dt((x - mu) / sigma, df = df) / sigma
  }
  pstudent_t <- function(q, mu, sigma, df) {
    pt((q - mu) / sigma, df = df)
  }
  qstudent_t <- function(p, mu, sigma, df) {
    mu + sigma * qt(p, df = df)
  }
  rstudent_t <- function(n, mu, sigma, df) {
    mu + sigma * rt(n, df = df)
  }
  
  variable <- na.omit(variable)
  results <- list()
  
  # 1. If variable is a factor
  if (is.factor(variable)) {
    n_levels <- length(levels(variable))
    log_info(paste("A variable is a factor with", n_levels, "levels"))
    if (n_levels == 2) {
      suggested_family <- "bernoulli"
      if (verbose){
        conclusion <- sprintf("Based on the number of levels (%d == 2), the suggested family is %s", n_levels, suggested_family)
        log_info(conclusion)
        print(conclusion)
      }
    } else {
      suggested_family <- "categorical"
      if (verbose){
        conclusion <- sprintf("Based on the number of levels (%d != 2), the suggested family is %s", n_levels, suggested_family)
        log_info(conclusion)
        print(conclusion)
      }
    }
    if (output_results) {
      results <- list(
        "Variable Type" = "Factor",
        "Number of Levels" = n_levels,
        "Suggested Family" = suggested_family
      )
      return(results)
    } else {
      return(suggested_family)
    }
  }
  
  # 2. If a variable is binary (numeric 0/1)
  if (all(variable %in% c(0, 1))) {
    suggested_family <- "binomial"
    log_info(paste("A variable is numeric (0/1)"))
    if (verbose){
      conclusion <- paste("Based on the numeric (0/1) structure, the suggested family is", suggested_family)
      log_info(conclusion)
      print(conclusion)
    }
    if (output_results) {
      results <- list(
        "Variable Type" = "Binary (Numeric 0/1)",
        "Suggested Family" = suggested_family
      )
      return(results)
    } else {
      return(suggested_family)
    }
  }
  
  attributes(variable) <- NULL  
  
  if (!is.numeric(variable)) {
    warning("Variable is not numeric. Attempting to convert to numeric.")
    variable <- as.numeric(as.character(variable))  # Convert factor/character to numeric
    variable <- na.omit(variable)
  }
  if (length(variable) <= 1) {
    stop("Error: Variable must have at least two numeric values for distribution fitting.")
  }
  
  # 3. Detect possible count data (non-negative integers)
  suggested_family_count <- NULL
  if (all(variable == floor(variable)) && all(variable >= 0)) {
    warn_message <- "Variable appears to be count data (non-negative integers)"
    log_warn(warn_message)
    warning(warn_message)
    var_mean <- mean(variable)
    var_var <- var(variable)
    var_ratio <- var_var / var_mean
    if (var_var > var_mean) {
      suggested_family_count <- "negbinomial"
      if (verbose){
        conclusion <- sprintf("If a variable is indeed a count variable, then the suggested family is %s based on variance-to-mean ratio %.3f, where variance (%.3f) is larger than mean (%.3f).",
          suggested_family_count, var_ratio, var_var, var_mean)
        log_info(conclusion)
        print(conclusion)
      }
    } else {
      suggested_family_count <- "poisson"
      if (verbose){
        conclusion <- sprintf("If a variable is indeed a count variable, then the suggested family is %s based on variance-to-mean ratio %.3f, where variance (%.3f) is not larger than mean (%.3f).",
          suggested_family_count, var_ratio, var_var, var_mean)
        log_info(conclusion)
        print(conclusion)
      }
    }
  }

  # 4. Uniformity check
  suggested_family_uniform <- NULL
  skew_val <- skewness(variable)
  kurt <- kurtosis(variable)
  ks_unif <- ks.test(variable, "punif", min(variable), max(variable))
  cvm_unif <- cvm.test((variable - min(variable)) / (max(variable) - min(variable)), "punif")
  # CVM Similar to the KS test but considers differences over the entire distribution rather than just the largest deviation. 
  # KS is a general goodness of fit test and is more sensitive to deviations in the middle of the distribution
  # KS measures Largest deviation between CDFs
  # CVM is about Overall difference between CDFs
  
  is_uniform <- (ks_unif$p.value > 0.05) && (cvm_unif$p.value > 0.05) &&
    (abs(skew_val) < 0.1) && (abs(kurt - 1.8) < 0.3)
  
  if (is_uniform) {
    warn_message <- sprintf("Variable appears to be approximately uniform: skewness (%.3f), kurtosis (%.3f), Kolmogorov-Smirnov test (p = %.5f), Cramér-von Mises test (p = %.5f)", 
                            skew_val, kurt, ks_unif$p.value, cvm_unif$p.value) # doesnt get displayed in console..?
    log_warn(warn_message)
    warning(warn_message)
    if (verbose){
      conclusion <- sprintf("The variable was classified as approximately uniform because (1) the Kolmogorov-Smirnov test (p = %.5f) did not detect significant deviation from a uniform distribution, (2) the Cramér-von Mises test (p = %.5f) also did not detect deviations from uniformity, (3) the sample skewness (%.3f) is close to zero, indicating symmetry, and (4) the sample kurtosis (%.3f) is near the expected value of 1.8 for a uniform distribution.", 
            ks_unif$p.value, cvm_unif$p.value, skew_val, kurt) # while this does..
      log_info(conclusion)
      print(conclusion)
    }
    if (all(variable >= 0) && all(variable <= 1)) {
      suggested_family_uniform <- "beta"
      if (verbose){
        conclusion <- paste("If a variable is indeed close to uniform, there is no directly fitting family in brms, but because the variable is between 0 and 1, the suggested family is", suggested_family_uniform)
        log_info(conclusion)
        print(conclusion)
      }
    } else {
      suggested_family_uniform <- default_family
      if (verbose){
        conclusion <- paste("If a variable is indeed close to uniform, there is no directly fitting family in brms, and because the variable is NOT between 0 and 1, the suggested family is default family", suggested_family_uniform)
        log_info(conclusion)
        print(conclusion)
      }
    }
  } else {
    if (verbose){
      conclusion <- sprintf("The variable was not classified as approximately uniform because (1) the Kolmogorov-Smirnov test (p = %.5f) detected significant deviation from a uniform distribution, or (2) the Cramér-von Mises test (p = %.5f) also detected deviations from uniformity, (3) or the sample skewness (%.3f) was not close to zero, indicating asymmetry, or (4) the sample kurtosis (%.3f) was not near the expected value of 1.8 for a uniform distribution.", 
                            ks_unif$p.value, cvm_unif$p.value, skew_val, kurt)
      log_info(conclusion)
      print(conclusion)
    }
  }
  
  # 5. For continuous data:
  sample_size <- length(variable)
  
  # Normality test: Shapiro-Wilk if n <= 5000; otherwise, Anderson-Darling.
  if (sample_size <= 5000) {
    nw_test <- shapiro.test(variable)
    normality_p <- nw_test$p.value
    normality_stat <- nw_test$statistic
    normality_test_name <- "Shapiro-Wilk"
  } else {
    ad_test <- ad.test(variable)
    normality_p <- ad_test$p.value
    normality_stat <- ad_test$statistic
    normality_test_name <- "Anderson-Darling"
  }
  
  # KS is more sensitive to deviations in the middle of the distribution, whereas AD is more sensitive to deviations in the tails. 
  # this is important to consider for tests on large samples
  # KS measures Largest deviation between CDFs
  # while AD is Weighted deviation emphasizing tails
  ks_test <- ks.test(variable, "pnorm", mean = mean(variable), sd = sd(variable))
  
  # Fit normal and custom Student-t distributions for AIC comparison.
  fit_gaussian <- fitdist(variable, "norm")
  fit_student <- fitdist(variable, "student_t", 
                         start = list(mu = mean(variable), sigma = sd(variable), df = 10))
  gof_results <- gofstat(list(fit_gaussian, fit_student))
  
  mean_AIC_1_perc <- 0.01 *mean(c(gof_results$aic[1], gof_results$aic[2]))
  abs_dif <- abs(gof_results$aic[1] - gof_results$aic[2])
  if (abs_dif < mean_AIC_1_perc) {
    warn_message <- sprintf("Gaussian (%.2f) and Student-t (%.2f) AIC values are within 1%% of their mean (%.2f < %.2f). Consider performing prior predictive checks for both likelihood families where relevant", 
                            gof_results$aic[1], gof_results$aic[2], abs_dif, mean_AIC_1_perc)
    log_warn(warn_message)
    warning(warn_message)
  }
  
  if (abs(skew_val) > 0.5) {
    warn_message <- sprintf("Absolute skewness (%.3f) is > 0.5. Consider performing prior predictive checks for skew-normal in any case", 
                            skew_val)
    log_warn(warn_message)
    warning(warn_message)
  }
  
  if (skew_val > 1) {
    warn_message <- sprintf("The variable is strongly right-skewed because its sample skewness (%.3f) exceeds 1. Consider performing prior predictive checks for exgaussian in any case",
                            skew_val)
    log_warn(warn_message)
    warning(warn_message)
  }
  
  dip_test <- dip.test(variable)
  if (dip_test$p.value < 0.05) {
    warn_message <- sprintf("Variable appears to be multimodal (Hartigan's Dip Test p-value: %.5f). Consider gp", dip_test$p.value)
    log_warn(warn_message)
    warning(warn_message)
  }
  
  # Extract the estimated degrees of freedom from the Student-t fit.
  df_est <- fit_student$estimate["df"]
  
  # Check if normal
  if (normality_p > 0.05 && ks_test$p.value > 0.05) {
    # If both normality tests pass, choose Gaussian.
    suggested_family <- "gaussian"
    if (verbose){
      conclusion <- sprintf("Neither AD or SW normality test (p = %.5f) nor KS test (p = %.5f) detected deviations from normality. Suggested family is %s", 
                            normality_p, ks_test$p.value, suggested_family)
      log_info(conclusion)
      print(conclusion)
    }
  } else if (kurt > 3 && gof_results$aic[2] < gof_results$aic[1]) {
    # Check if heavy-tailed
    if (verbose){
      conclusion <- sprintf("Either AD or SW normality test (p = %.5f) or KS test (p = %.5f) detected deviations from normality. The variable was classified as heavy-tailed because its kurtosis (%.3f) exceeds 3 and the Student-t model (AIC = %.2f) provides a better fit than the Gaussian model (AIC = %.2f)", 
                            normality_p, ks_test$p.value, kurt, gof_results$aic[2], gof_results$aic[1])
      log_info(conclusion)
      print(conclusion)
    }
    # Heavy-tailed: decide based on limit_student_df setting.
    if (limit_student_df) {
      if (df_est <= df_threshold) {
        suggested_family <- "student"
        if (verbose){
          conclusion <- sprintf("The variable was classified as heavy-tailed, and the estimated degrees of freedom (%.2f) are below the provided threshold (%.2f) suggesting that the estimated Student-t is meaningfully different from Gaussian, so the suggested family is %s", 
                                df_est, df_threshold, suggested_family)
          log_info(conclusion)
          print(conclusion)
        }
      } else {
        suggested_family <- "gaussian"
        if (verbose){
          conclusion <- sprintf("The variable was classified as heavy-tailed, but the estimated degrees of freedom (%.2f) are higher than the provided threshold (%.2f). This indicates that the distribution is almost indistinguishable from Gaussian, so the suggested family is %s. However, consider both Student-t and Gaussian for prior predictive checks", 
                                df_est, df_threshold, suggested_family)
          log_info(conclusion)
          print(conclusion)
        }
      }
    } else {
      suggested_family <- "student"
      if (verbose){
        conclusion <- sprintf("The variable was classified as heavy-tailed, the suggested family is %s", suggested_family)
        log_info(conclusion)
        print(conclusion)
      }
    }
  } else {
    # Check if skewed and how
    if (verbose){
      conclusion <- sprintf("Either AD or SW normality test (p = %.5f) or KS test (p = %.5f) detected deviations from normality but no heavy tails were observed based on kurtosis (%.3f) or comparing fit of Student-t (AIC = %.2f) and Gaussian models (AIC = %.2f)", 
                            normality_p, ks_test$p.value, kurt, gof_results$aic[2], gof_results$aic[1])
      log_info(conclusion)
      print(conclusion)
    }
    
    skew_model <- selm(variable ~ 1, family = "SN")
    cp_params <- skew_model@param$cp  
    gamma1_est <- cp_params["gamma1"]
    cp_var <- skew_model@param.var$cp
    if ("gamma1" %in% rownames(cp_var) && "gamma1" %in% colnames(cp_var)) {
      gamma1_se <- sqrt(cp_var["gamma1", "gamma1"])
    } else {
      gamma1_se <- NA
      warning("Variance of gamma1 (skewness parameter) not available. The skewness test may not be reliable.")
    }
    
    if (!is.na(gamma1_se) && gamma1_se > 0) {
      z_value <- gamma1_est / gamma1_se
      skew_p <- 2 * (1 - pnorm(abs(z_value)))
      if (is.na(skew_p)) skew_p <- 0
    } else {
      z_value <- NA
      skew_p <- NA
    }
    
    if (!is.na(skew_p) && abs(skew_val) > 0.5 && skew_p < 0.05) {
      suggested_family <- "skew_normal"
      if (verbose){
        conclusion <- sprintf("The variable was classified as significantly skewed because its sample skewness (%.3f) exceeds 0.5. The skewness test yielded a p-value of %.5f, indicating that the skewness is statistically significant. Therefore, the suggested family is %s", 
                              skew_val, skew_p, suggested_family)      
        log_info(conclusion)
        print(conclusion)
      }
    } else if (skew_val > 1) {
      suggested_family <- "exgaussian"
      if (verbose){
        conclusion <- sprintf("The variable is classified as strongly right-skewed because its sample skewness (%.3f) exceeds 1. Even if the skewness test is non-significant, the high magnitude of skew suggests that an ex-Gaussian distribution is more appropriate. Therefore, the suggested family is %s", 
                              skew_val, suggested_family)
        log_info(conclusion)
        print(conclusion)
      }
    } else {
      suggested_family <- default_family
      if (verbose){
        conclusion <- sprintf("The variable does not exhibit substantial skewness. The sample skewness (%.3f) is within the acceptable range, and/or the skewness test yielded a p-value of %.5f, indicating no statistically significant skew. The default family %s is selected.", 
              skew_val, skew_p, default_family)
        log_info(conclusion)
        print(conclusion)
      }
    }
  }
  
  # Prepare detailed results if requested.
  if (output_results) {
    results <- list(
      "Variable Type" = "Continuous",
      "Sample Size" = sample_size,
      # Normality Tests
      "Normality Test" = normality_test_name,
      "Normality Test Statistic" = normality_stat,
      "Normality Test p-value" = normality_p,
      "KS Test Statistic" = ks_test$statistic,
      "KS Test p-value" = ks_test$p.value,
      # Skewness & Kurtosis
      "Skewness (sample)" = skew_val,
      "Kurtosis" = kurt,
      # AIC Values
      "AIC (Normal)" = gof_results$aic[1],
      "AIC (Student)" = gof_results$aic[2],
      # Student-t df
      "Estimated df from Student-t" = df_est,
      # Skewness Tests
      "Gamma1 (CP) Estimate" = if (exists("gamma1_est")) gamma1_est else NA,
      "Gamma1 Standard Error" = if (exists("gamma1_se")) gamma1_se else NA,
      "Gamma1 Z-value" = if (exists("z_value")) z_value else NA,
      "Gamma1 p-value" = if (exists("skew_p")) skew_p else NA,
      # Count & Uniformity Checks
      "Mean" = if (exists("var_mean")) var_mean else NA,
      "Variance" = if (exists("var_var")) var_var else NA,
      "Variance-to-Mean Ratio" = if (exists("var_ratio")) var_ratio else NA,
      "KS Test p-value (Uniformity Check)" = ks_unif$p.value,
      "CvM Test p-value (Uniformity Check)" = cvm_unif$p.value,
      "Hartigan's Dip Test p-value" = dip_test$p.value,
      # Suggested Families
      "Default Family" = default_family,
      "Suggested Count Family" = suggested_family_count,
      "Suggested Uniform Family" = suggested_family_uniform,
      "Suggested Family (Continuous)" = suggested_family
    )
    # If limit_student_df is enabled, include df details.
    if (limit_student_df && normality_p <= 0.05) {
      results[["df_threshold"]] <- df_threshold
    }
    return(results)
  } else {
    return(suggested_family)
  }
}

description_for_rr <- "The selection of the likelihood family was based on a series of statistical tests examining the distributional properties of the variable. First, factor variables were classified as binomial (if binary) or categorical (if multi-level). Numeric binary variables (0/1) were also classified as binomial. For variables appearing to be count data (non-negative integers), we flagged them for potential modeling using Poisson or negative binomial families, based on whether the variance exceeded the mean. Next, we performed uniformity checks using the Kolmogorov-Smirnov (KS) test, Cramér-von Mises (CvM) test, skewness, and kurtosis. If the variable met uniformity criteria, we suggested beta (if bounded between [0,1]) or defaulted to the predefined default_family otherwise.
For continuous variables, we conducted normality tests using Shapiro-Wilk (if n ≤ 5000) or Anderson-Darling (otherwise), alongside a Kolmogorov-Smirnov test for normality. If both tests failed to reject normality, the Gaussian likelihood was selected. If normality was rejected and kurtosis exceeded 3, we compared Gaussian and Student-t models via Akaike Information Criterion (AIC). If Student-t provided a better fit (lower AIC), we classified the variable as heavy-tailed and suggested a Student-t distribution. If heavy tails were absent, we performed skewness testing using a skew-normal fit. If skewness was statistically significant (p < 0.05) and absolute skewness exceeded 0.5, the skew-normal likelihood was selected. In cases where we observed extreme right-skew (skewness > 1) even in the absence of significant skeweness, an ex-Gaussian likelihood was preferred. Otherwise, the default_family was retained."
additional_description_with_limited_df_for_student <- "If limit_student_df = TRUE, an additional constraint was applied when selecting the Student-t likelihood. In this case, the estimated degrees of freedom (df_est) from the Student-t fit were compared against a threshold (df_threshold). If df_est < df_threshold, Student-t was retained. However, if df_est ≥ df_threshold, the distribution was considered effectively Gaussian, and a Gaussian likelihood was selected instead. This accounts for scenarios where the Student-t model is nearly identical to Gaussian but would still be selected due to small AIC differences."

corr_and_plot <- function(df, cols, method = 'spearman') {
  
  for(col in cols) {
    if(is.factor(df[[col]]) && length(levels(df[[col]])) == 2) {
      # Convert factor to numeric (0, 1)
      df[[col]] <- as.numeric(df[[col]]) - 1
    }
  }
  
  # Subset the dataframe to include only the desired columns
  sub_df <- df[, cols, drop = FALSE]
  
  # Compute the correlation matrix
  cor_mat <- cor(sub_df, use = "complete.obs", method = method)
  
  mapped_names <- sapply(rownames(cor_mat), map_variable_name)
  rownames(cor_mat) <- mapped_names
  colnames(cor_mat) <- mapped_names
  
  # Convert the matrix to a data frame for output (optional)
  cor_df <- as.data.frame(cor_mat)
  
  # Melt the correlation matrix into long format for ggplot
  cor_melt <- melt(cor_mat)
  cor_melt$Var2 <- factor(cor_melt$Var2, levels = rev(levels(cor_melt$Var2)))
  
  cor_melt$Correlation <- cor_melt$value
  
  # Create the ggplot heatmap
  p <- ggplot(cor_melt, aes(x = Var1, y = Var2, fill = Correlation)) +
    geom_tile() +
    scale_fill_gradientn(colors = rev(brewer.pal(11, "PuOr"))) +
    geom_text(aes(label = round(value, 2)), color = "black", size = 4) +
    theme_minimal() +
    labs(title = "Correlation Matrix of Covariates", x = "", y = "") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Return both the data frame and the ggplot object in a list
  return(list(correlation_matrix = cor_df, plot = p))
}