library(dplyr)
library(ggplot2)

plot_progression <- function(df, var_col, ci = FALSE, add_lines = FALSE, 
                            average_across_blocks = FALSE, colors = "darkgrey") {

  # Convert the target variable to numeric.
  # For a two-level factor, code the first level as 1 and the second as 0.
  if (is.factor(df[[var_col]]) && nlevels(df[[var_col]]) == 2) {
    df$num_value <- as.numeric(df[[var_col]] == levels(df[[var_col]])[1])
  } else {
    df$num_value <- as.numeric(df[[var_col]])
  }
  
  if (!average_across_blocks) {
    # Use the full 300 trials (trialCount) as x-axis.
    # Compute summary statistics for each trialCount.
    if (ci) {
      mean_df <- df %>%
        group_by(trialCount) %>%
        summarize(
          mean_value = mean(num_value, na.rm = TRUE),
          sd_value   = sd(num_value, na.rm = TRUE),
          n          = sum(!is.na(num_value)),
          se         = sd_value / sqrt(n),
          lower      = mean_value - 1.96 * se,
          upper      = mean_value + 1.96 * se
        ) %>%
        ungroup()
    } else {
      mean_df <- df %>%
        group_by(trialCount) %>%
        summarize(mean_value = mean(num_value, na.rm = TRUE)) %>%
        ungroup()
    }
    
    # Prepare data for individual subject lines (only if not using CI).
    if (!ci) {
      df_processed <- df %>%
        arrange(subID, trialCount) %>%
        group_by(subID) %>%
        # Create a segment variable: if the difference between successive trialCount is not 1,
        # then increment the segment (this prevents connecting across gaps).
        mutate(segment = cumsum(if_else(trialCount - lag(trialCount, default = first(trialCount)) == 1, 0, 1))) %>%
        ungroup()
    }
    
    # Start building the plot.
    p <- ggplot()
    
    # Add vertical dashed lines first so that they appear in the background.
    if (add_lines) {
      max_trial <- max(df$trialCount, na.rm = TRUE)
      # Compute breaks every 50 trials.
      breaks <- seq(50, max_trial, by = 50)
      # When plotting all 300 trials, remove a line at the final trial if it is exactly max_trial.
      breaks <- breaks[breaks != max_trial]
      p <- p + geom_vline(xintercept = breaks, 
                          linetype   = "dashed", 
                          color      = "lightgrey", 
                          size       = 0.5)
    }
    
    # Set axis labels (using map_variable_name for y).
    p <- p + labs(x = "Trial number", y = map_variable_name(var_col))
    
    # Add the data layers.
    if (ci) {
      p <- p + geom_ribbon(data = mean_df, 
                           aes(x = trialCount, ymin = lower, ymax = upper),
                           fill  = colors[1], 
                           alpha = 0.5)
    } else {
      p <- p + geom_line(data = df_processed,
                         aes(x = trialCount, y = num_value, group = interaction(subID, segment)),
                         color = "lightgrey", alpha = 0.3)
    }
    
    # Overlay the mean progression line.
    p <- p + geom_line(data = mean_df,
                       aes(x = trialCount, y = mean_value),
                       color = colors[1], size = 1.2)
    
  } else {
    # average_across_blocks is TRUE.
    # Average across blocks so that each subject has one value per trialInBlock (100 trials).
    subj_avg <- df %>%
      group_by(subID, trialInBlock) %>%
      summarize(num_value = mean(num_value, na.rm = TRUE)) %>%
      ungroup()
    
    # For individual subject lines, compute segments based on trialInBlock.
    if (!ci) {
      subj_avg <- subj_avg %>%
        arrange(subID, trialInBlock) %>%
        group_by(subID) %>%
        mutate(segment = cumsum(if_else(trialInBlock - lag(trialInBlock, default = first(trialInBlock)) == 1, 0, 1))) %>%
        ungroup()
    }
    
    # Compute summary statistics by trialInBlock across subjects.
    if (ci) {
      mean_df <- subj_avg %>%
        group_by(trialInBlock) %>%
        summarize(
          mean_value = mean(num_value, na.rm = TRUE),
          sd_value   = sd(num_value, na.rm = TRUE),
          n          = sum(!is.na(num_value)),
          se         = sd_value / sqrt(n),
          lower      = mean_value - 1.96 * se,
          upper      = mean_value + 1.96 * se
        ) %>%
        ungroup()
    } else {
      mean_df <- subj_avg %>%
        group_by(trialInBlock) %>%
        summarize(mean_value = mean(num_value, na.rm = TRUE)) %>%
        ungroup()
    }
    
    # Build the plot.
    p <- ggplot()
    # Add a vertical dashed line at trial 50.
    if (add_lines) {
      p <- p + geom_vline(xintercept = 50, 
                          linetype   = "dashed", 
                          color      = "lightgray", 
                          size       = 0.5)
    }
    
    p <- p + labs(x = "Trial number", y = map_variable_name(var_col))
    
    if (ci) {
      p <- p + geom_ribbon(data = mean_df,
                           aes(x = trialInBlock, ymin = lower, ymax = upper),
                           fill  = colors[1],
                           alpha = 0.5)
    } else {
      p <- p + geom_line(data = subj_avg,
                         aes(x = trialInBlock, y = num_value, group = interaction(subID, segment)),
                         color = "lightgrey", alpha = 0.3)
    }
    
    # Overlay the mean progression line.
    p <- p + geom_line(data = mean_df,
                       aes(x = trialInBlock, y = mean_value),
                       color = colors[1], size = 1.2)
  }
  
  return(p)
}

plot_dot <- function(vec1, vec2, label1 = "Group1", label2 = "Group2",
                     xlab = "Group", ylab = "Value", title = "Dot Plot",
                     dot_size = 0.5, binwidth = 0.5, colors = c("blue", "red")) {
  # Combine the two vectors into a single data frame with a grouping variable.
  df <- data.frame(
    value = c(vec1, vec2),
    group = factor(c(rep(label1, length(vec1)), rep(label2, length(vec2))),levels = c(label1,label2))
  )
  
  # Create the dot plot using geom_dotplot.
  p <- ggplot(df, aes(x = group, y = value, fill = group)) +
    geom_dotplot(binaxis = "y", stackdir = "center",color = NA) +
    labs(x = xlab, y = ylab, title = title) +
    scale_fill_manual(values = colors) +
    theme(legend.position = "none")
  
  return(p)
}
