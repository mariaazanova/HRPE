library(cowplot)
library(RColorBrewer)
library(ggplot2)
library(ggsci)

plot_theme_default <- theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        axis.title.x = element_text(size = 12, margin = margin(t = 4)),
        axis.title.y = element_text(size = 12))
theme_set(plot_theme_default)
update_geom_defaults("text", list(size = 4))

plot_theme_default <- theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        axis.title.x = element_text(size = 14, margin = margin(t = 6)),
        axis.title.y = element_text(size = 14))
theme_set(plot_theme_default)
update_geom_defaults("text", list(size = 6))

plot_theme_word <- theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        axis.title.x = element_text(size = 10, margin = margin(t = 4)),
        axis.title.y = element_text(size = 10))

get_colors <- function(palette = "PuOr") {
  palette_colors <- brewer.pal(n = 9, name = palette)  
  colors <- palette_colors[c(2, 8, 1, 9, 3, 7, 4, 6, 5)] 
  return(colors)
}

get_colors_2 <- function(palette = "Dark2") {
  palette_colors <- brewer.pal(n = 8, name = palette)  
  colors <- palette_colors[c(1,4,3,2, 8, 7, 6, 5)] 
  return(colors)
}

get_colors_3 <- function(palette = "PRGn") {
  palette_colors <- brewer.pal(n = 9, name = palette)  
  colors <- palette_colors[c(2, 8, 1, 9, 3, 7, 4, 6, 5)] 
  return(colors)
}

get_colors_purples <- function(palette = "Purples") {
  palette_colors <- brewer.pal(n = 8, name = palette)  
  colors <- palette_colors[c(4,5,6,7,8)] 
  return(colors)
}

get_colors_greens <- function(palette = "Greens") {
  palette_colors <- brewer.pal(n = 8, name = palette)  
  colors <- palette_colors[c(4,5,6,7,8)] 
  return(colors)
}



npg_pal <- pal_npg(palette=c("nrc"),alpha = 1)
colors3 <- npg_pal(10)

save_figure <- function(fig, aspect_ratio, save_path, filename, ext) {
  ggsave(
    filename = file.path(save_path, paste0(filename, ".", ext)),
    plot = fig,
    width = figure_setup$width,
    height = aspect_ratio * figure_setup$width,
    units = figure_setup$units,
    dpi = figure_setup$dpi,
    device = ext,
    bg = "white"
  )
}

figure_setup <- list(
  width = 190,
  units = "mm",
  dpi = 300
)

font_setup <- list(
  title = 10,
  subtitle = 8,
  axis.title = 8,
  axis.text = 7,
  legend.title = 9,
  legend.text = 8
)

plot_theme_default <-  theme_minimal() +
  theme(panel.grid = element_blank(), 
        plot.title = element_text(size = font_setup$title), 
        plot.subtitle = element_text(size = font_setup$subtitle),
        axis.title.x = element_text(size = font_setup$axis.title),
        axis.text.x = element_text(size = font_setup$axis.text),
        axis.title.y = element_text(size = font_setup$axis.title),
        axis.text.y = element_text(size = font_setup$axis.text),
        legend.title = element_text(size = font_setup$legend.title),
        legend.text = element_text(size = font_setup$legend.text))

theme_set(plot_theme_default)

# can add map_var_names to titles..
plot_subject_summary <- function(df, value_var, group_var, subject_var = "subID", facet_var = NULL,
                                 type = c("bar", "mean", "box", "violin", "dot", "line", "raincloud",
                                          "dotplot", "dotplot_violin", "dotplot_mean_ci", "dotplot_boxplot", "dotplot_raincloud"),
                                 remove_outliers = TRUE, lower_q = 0.025, upper_q = 0.975,
                                 fill_color = "skyblue", dotplot_binwidth = 5, dotplot_size = 1, dot_move = 0, shared_ylim = FALSE,
                                 x_title = NULL, y_title = NULL, plot_title = NULL) {
  type <- match.arg(type)
  stopifnot(all(c(value_var, group_var, subject_var) %in% names(df)))
  
  dot_nudge <- dot_move 
  
  # Drop rows with NA in any grouping variable
  drop_vars <- c(value_var, group_var, subject_var)
  if (!is.null(facet_var)) drop_vars <- c(drop_vars, facet_var)
  df <- df %>% filter(if_all(all_of(drop_vars), ~ !is.na(.)))
  
  if (!is.null(facet_var)) {
    facet_levels <- unique(df[[facet_var]])
    if (!identical(shared_ylim, FALSE) && is.null(shared_ylim)) {
      shared_ylim <- df %>%
        pull(.data[[value_var]]) %>%
        quantile(c(0.025, 0.975), na.rm = TRUE)
      # Optional padding
      ylim_range <- diff(shared_ylim)
      shared_ylim <- shared_ylim + c(-0.05, 0.05) * ylim_range
    }
    plot_list <- lapply(facet_levels, function(facet_level) {
      df_sub <- df[df[[facet_var]] == facet_level, , drop = FALSE]
      single_plot <- plot_subject_summary(
        df = df_sub,
        value_var = value_var,
        group_var = group_var,
        subject_var = subject_var,
        type = type,
        facet_var = NULL,  # prevent recursive loop
        shared_ylim = shared_ylim, 
        remove_outliers = remove_outliers,
        lower_q = lower_q,
        upper_q = upper_q,
        fill_color = fill_color,
        x_title = x_title,
        y_title = y_title,
        plot_title = paste0(plot_title, " – ", facet_var, ": ", facet_level),
        dotplot_size = dotplot_size,
        dot_move = dot_move,
        dotplot_binwidth = dotplot_binwidth
      )
      return(single_plot)
    })
    return(cowplot::plot_grid(plotlist = plot_list, ncol = 2))
  }
  
  if (remove_outliers) {
    df <- df %>%
      group_by(.data[[group_var]]) %>%
      mutate(
        lower = quantile(.data[[value_var]], lower_q, na.rm = TRUE),
        upper = quantile(.data[[value_var]], upper_q, na.rm = TRUE)
      ) %>%
      filter(.data[[value_var]] >= lower, .data[[value_var]] <= upper) %>%
      select(-lower, -upper) %>%
      ungroup()
  }
  
  # Check if group_var is constant per subject
  group_check <- df %>%
    group_by(.data[[subject_var]]) %>%
    summarise(n_unique = n_distinct(.data[[group_var]]), .groups = "drop")
  
  group_var_is_between <- all(group_check$n_unique == 1)
  
  if (group_var_is_between) {
    # Between-subject grouping (e.g., control vs treatment)
    subject_groups <- df %>%
      distinct(.data[[subject_var]], .data[[group_var]])
    
    subj_avg <- df %>%
      group_by(.data[[subject_var]]) %>%
      summarise(mean_value = mean(.data[[value_var]], na.rm = TRUE), .groups = "drop") %>%
      left_join(subject_groups, by = subject_var)
    
  } else {
    # Within-subject condition (e.g., condition A vs B)
    subj_avg <- df %>%
      group_by(.data[[subject_var]], .data[[group_var]]) %>%
      summarise(mean_value = mean(.data[[value_var]], na.rm = TRUE), .groups = "drop")
  }
  
  # Compute summary stats
  summary_df <- subj_avg %>%
    group_by(.data[[group_var]]) %>%
    summarise(
      mean = mean(mean_value),
      se = sd(mean_value) / sqrt(n()),
      ci = 1.96 * se,
      .groups = "drop"
    )
  
  # Base plot
  p <- ggplot()
  
  if (type == "bar") {
    p <- p +
      geom_col(data = summary_df, aes(x = .data[[group_var]], y = mean), fill = fill_color, width = 0.6) +
      geom_errorbar(data = summary_df,
                    aes(x = .data[[group_var]], ymin = mean - ci, ymax = mean + ci), width = 0.2)
  } else if (type == "mean") {
    p <- p +
      geom_point(data = summary_df, aes(x = .data[[group_var]], y = mean), size = dotplot_size, color = 'black') +
      geom_errorbar(data = summary_df,
                    aes(x = .data[[group_var]], ymin = mean - ci, ymax = mean + ci), width = 0.1)
  } else if (type == "box") {
    p <- ggplot(df, aes(x = .data[[group_var]], y = .data[[value_var]])) +
      geom_boxplot(fill = fill_color, outlier.shape = NA)
  } else if (type == "violin") {
    p <- ggplot(df, aes(x = .data[[group_var]], y = .data[[value_var]])) +
      geom_violin(fill = fill_color, trim = FALSE)
  } else if (type == "dot") {
    p <- ggplot(subj_avg, aes(x = .data[[group_var]], y = mean_value)) +
      geom_point(position = position_jitter(width = 0.15), size = dotplot_size, color = fill_color)
  } else if (type == "line") {
    p <- ggplot(subj_avg, aes(x = .data[[group_var]], y = mean_value, group = .data[[subject_var]])) +
      geom_line(color = "gray50", alpha = 0.4) +
      geom_point(color = fill_color, size = dotplot_size)
  } else if (type == "dotplot") {
    p <- ggplot(subj_avg, aes(x = .data[[group_var]], y = mean_value, fill = .data[[group_var]])) +
      geom_dotplot(binaxis = "y", stackdir = "center", binwidth = dotplot_binwidth, dotsize = dotplot_size, color = fill_color) +
      scale_fill_manual(values = rep(fill_color, length(unique(subj_avg[[group_var]])))) +
      guides(fill = "none")
  } else if (type == "raincloud") {
    if (!requireNamespace("ggdist", quietly = TRUE)) stop("Install 'ggdist' for raincloud")
    p <- ggplot(df, aes(x = .data[[group_var]], y = .data[[value_var]])) +
      ggdist::stat_halfeye(adjust = 0.5, width = 0.6, .width = 0.95,
                           slab_color = "gray40", fill = fill_color)
  }
  
  ## ---- COMBINED TYPES ----
  
  if (type == "dotplot_violin") {
    p <- ggplot(df, aes(x = .data[[group_var]], y = .data[[value_var]])) +
      geom_violin(fill = fill_color, trim = FALSE, alpha = 0.5) +
      geom_dotplot(data = subj_avg, aes(y = mean_value), binaxis = "y", stackdir = "center",
                   binwidth = dotplot_binwidth, dotsize = dotplot_size, color = fill_color, fill = fill_color,
                   position = position_nudge(x = dot_nudge))
  } else if (type == "dotplot_mean_ci") {
    p <- ggplot(summary_df, aes(x = .data[[group_var]], y = mean)) +
      geom_point(size = dotplot_size, color = 'black') +
      geom_errorbar(aes(ymin = mean - ci, ymax = mean + ci), width = 0.15) +
      geom_dotplot(data = subj_avg, aes(x = .data[[group_var]], y = mean_value),
                   binaxis = "y", stackdir = "center", binwidth = dotplot_binwidth, dotsize = dotplot_size,
                   color = fill_color, fill = fill_color,
                   position = position_nudge(x = dot_nudge))
  } else if (type == "dotplot_boxplot") {
    p <- ggplot(df, aes(x = .data[[group_var]], y = .data[[value_var]])) +
      geom_boxplot(fill = fill_color, outlier.shape = NA, alpha = 0.5) +
      geom_dotplot(data = subj_avg, aes(y = mean_value), binaxis = "y", stackdir = "center",
                   binwidth = dotplot_binwidth, dotsize = dotplot_size, color = fill_color, fill = fill_color,
                   position = position_nudge(x = dot_nudge))
  } else if (type == "dotplot_raincloud") {
    if (!requireNamespace("ggdist", quietly = TRUE)) stop("Install 'ggdist' for raincloud")
    p <- ggplot(df, aes(x = .data[[group_var]], y = .data[[value_var]])) +
      ggdist::stat_halfeye(adjust = 0.5, width = 0.6, .width = 0.95,
                           slab_color = "gray40", fill = fill_color) +
      geom_dotplot(data = subj_avg, aes(y = mean_value), binaxis = "y", stackdir = "center",
                   binwidth = dotplot_binwidth, dotsize = dotplot_size, color = fill_color, fill = fill_color,
                   position = position_nudge(x = dot_nudge))
  }
  
  ## ---- COMMON STYLING ----
  
  p <- p +
    labs(
      x = x_title %||% group_var,
      y = y_title %||% value_var,
      title = plot_title
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank()
    ) 
  if (!is.null(shared_ylim) && !identical(shared_ylim, FALSE)) {
    p <- p + coord_cartesian(ylim = shared_ylim)
  }
  
  return(p)
}


save_plot <- function(p, vis_path, file_name, plot_format = "svg", plot_width = 6, plot_height = 6, dpi = 150) {
  full_path <- file.path(vis_path, paste0(file_name, ".", plot_format))
  ggsave(filename = full_path, plot = p, width = plot_width, height = plot_height, dpi = dpi)
}

combine_plots <- function(plots, col_ = NULL, lettered = FALSE, rel_widths = NULL) {
  # Filter out any NULL or non-ggplot objects from the list.
  plots <- Filter(function(p) !is.null(p) && inherits(p, "ggplot"), plots)
  
  # Wrap each plot with ggdraw() to force rendering as a grob.
  plots <- lapply(plots, cowplot::ggdraw)
  
  # Ensure that at least two valid plots are provided.
  if (length(plots) < 2) {
    stop("At least two valid plots are required to combine.")
  }
  
  # Determine the number of columns: use col_ if provided, otherwise default logic.
  ncol_val <- if (!is.null(col_)) col_ else if (length(plots) == 2) 1 else 2
  
  # Generate labels if requested.
  labels <- if (lettered) LETTERS[seq_along(plots)] else NULL
  
  # Combine the plots using cowplot::plot_grid.
  if (is.null(rel_widths)) {
    combined <- cowplot::plot_grid(plotlist = plots, ncol = ncol_val, labels = labels, vjust = 1)
  } else {
    combined <- cowplot::plot_grid(plotlist = plots, ncol = ncol_val, labels = labels, vjust = 1, rel_widths = rel_widths)
  }
  
  return(combined)
}

map_variable_name <- function(var_name, say_ibi = F) {
  
  # for coefficient names
  var_name <- gsub("b_", "", var_name, fixed = TRUE)
  var_name <- gsub("scale_", "", var_name, fixed = TRUE)
  var_name <- gsub(":", " : ", var_name, fixed = TRUE)
  
  # Specific replacements
  var_name <- gsub("absPE", "Absolute PE", var_name, fixed = TRUE)
  var_name <- gsub("abs_expectation", "Absolute expectation", var_name, fixed = TRUE)
  var_name <- gsub("expectation", "Expectation", var_name, fixed = TRUE)
  var_name <- gsub("learning_miniblock", "Miniblock", var_name, fixed = TRUE)
  var_name <- gsub("phaseFeedb", "Cardiac phase", var_name, fixed = TRUE)
  var_name <- gsub("sysDiaFeedbDiastole", "Cardiac phase Diastole", var_name, fixed = TRUE)
  var_name <- gsub("msDistToFeedb", "Ms from R-peak -1 to Feedback", var_name, fixed = TRUE)
  var_name <- gsub("choiceApproach", "Approach", var_name, fixed = TRUE)
  var_name <- gsub("frn", "FRN", var_name, fixed = TRUE)
  var_name <- gsub("_TFR", "", var_name, fixed = TRUE)
  var_name <- gsub("tfr", "for TFR", var_name, fixed = TRUE)
  var_name <- gsub("_2groups", "", var_name, fixed = TRUE)
  var_name <- gsub("_3groups", "", var_name, fixed = TRUE)
  var_name <- gsub("trueCorrect", "Correct choice", var_name, fixed = TRUE)
  var_name <- gsub("outcomeLose", "Outcome Lose", var_name, fixed = TRUE)
  var_name <- gsub("sysDiaFeedb", "Cardiac phase", var_name, fixed = TRUE)
  var_name <- gsub("phaseDiastole", "phase Diastole", var_name, fixed = TRUE)
  var_name <- gsub("trialInBlock", "Trial in Block", var_name, fixed = TRUE)
  var_name <- gsub("trialInMiniBlock", "Trial in Miniblock", var_name, fixed = TRUE)
  var_name <- gsub("beforeFeedbAveIBI", "Average IBI before feedback", var_name, fixed = TRUE)
  var_name <- gsub("beforeFeedbSDNN", "SDNN before feedback", var_name, fixed = TRUE)
  var_name <- gsub("beforeFeedbRMSSD", "RMSSD before feedback", var_name, fixed = TRUE)
  var_name <- gsub("EKG1", "ECG", var_name, fixed = TRUE)
  var_name <- gsub("RtoF_", "Distance from feedback to IBI", var_name, fixed = TRUE)
  var_name <- gsub("Pz_p3b_Lose", "P300", var_name, fixed = TRUE)
  var_name <- gsub("Pz_p3b_Win", "P300", var_name, fixed = TRUE)
  
  # Handle "changeIBI" pattern
  var_name <- gsub("changeIBI(\\d+)_\\d+", "Change of IBI \\1", var_name, perl = TRUE)
  
  # Handle IBI numbering
  var_name <- gsub("IBI(\\d+)", "IBI \\1", var_name, perl = TRUE)
  
  # Handle "__" replacement with ", " and "withBase" with "baseline-corrected"
  var_name <- gsub("__", ", ", var_name, fixed = TRUE)
  var_name <- gsub("withBase", "baseline-corrected", var_name, fixed = TRUE)
  
  # IBI number mapping
  number_mapping <- c(
    "1" = "-4", "2" = "-3", "3" = "-2", "4" = "-1",
    "5" = "1", "6" = "2", "7" = "3", "8" = "4"
  )
  
  matches <- gregexpr("IBI (\\d+)", var_name, perl = TRUE)
  extracted <- regmatches(var_name, matches)
  
  replaced <- lapply(extracted, function(matches) {
    sapply(matches, function(match) {
      num <- sub(".* ", "", match)
      if (num %in% names(number_mapping)) {
        paste0("IBI ", number_mapping[[num]])  
      } else {
        match 
      }
    })
  })
  
  regmatches(var_name, matches) <- replaced
  
  if (!say_ibi){
    var_name <- gsub("IBI", "RR interval", var_name, fixed = TRUE)
  }
  
  # Start with a capital letter
  var_name <- sub("^([a-z])", "\\U\\1", var_name, perl = TRUE)
  
  # Ensure proper spacing and trimming
  var_name <- gsub("_", " ", var_name)
  var_name <- gsub("[[:space:]]+", " ", var_name)
  var_name <- trimws(var_name)
  
  return(var_name)
}
