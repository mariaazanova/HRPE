library(ggplot2)
library(reshape2)
library(eegUtils)
library(scales) 
library(RColorBrewer)
library(digest)

prefixes <- c("pseudo_", "baseline_", "long_", "tfr_")

plot_evoked_responses <- function(keywords, input_dir, output_dir, 
                                  labels = NULL, legend_title = "Condition",
                                  show_title = TRUE, title_string = "",
                                  plot_ci = FALSE, chosen_channels = c(),
                                  show_plots = FALSE, save_plots = TRUE,
                                  plot_format = "svg", palette = "Dark2", colors = c(),
                                  remove_grid = TRUE, add_0 = FALSE,
                                  add_shaded_areas = FALSE, pick_roi = "all",
                                  latencies_file = NULL,
                                  plot_width = 8, plot_height = 6,
                                  legend_position = "bottom", 
                                  tmin = -0.3, tmax = 0.7
                                  ) {
  
  if (!is.vector(keywords)) {
    stop("Keywords should be provided as a vector.")
  }
  if (!is.null(labels) && length(labels) != length(keywords)) {
    stop("The length of 'labels' must match the length of 'keywords'.")
  }
  
  final_latencies <- NULL
  if (!is.null(latencies_file)) {
    raw_latencies <- read.csv(latencies_file, header = TRUE)
    final_latencies <- as.data.frame(t(raw_latencies))  
    colnames(final_latencies) <- "Time"
    final_latencies$Name <- rownames(final_latencies) 
    rownames(final_latencies) <- NULL  
    final_latencies$Time <- as.numeric(final_latencies$Time)
  }
  
  mapping_file <- file.path(input_dir, "filename_mapping.csv")
  maps <- load_filename_mapping(mapping_file)
  orig2short  <- maps$orig2short

  extension <- "_Channels_CI.csv"
  virtual_files <- paste0(keywords, extension)
  actual_files <- ifelse(virtual_files %in% names(orig2short),
                         orig2short[virtual_files],
                         virtual_files)
  actual_files_fullpaths <- file.path(input_dir, actual_files)
  
  missing_files <- actual_files_fullpaths[!file.exists(actual_files_fullpaths)]
  existing_files <- actual_files_fullpaths[file.exists(actual_files_fullpaths)]
  missing_keywords <- keywords[!file.exists(actual_files_fullpaths)]
  existing_keywords <- keywords[file.exists(actual_files_fullpaths)]
  
  if (length(missing_files) > 0) {
    
    if (plot_ci) {
      log_warn("The following files were not found:\n", paste(missing_files, collapse = "\n"))
      log_warn("Expected evokeds with CI but files were missing. Trying to plot data without CI.")
      plot_ci = FALSE
    }
    
    modified_filenames <- sub("_CI\\.csv$", ".csv", missing_files)
    missing_files_modified <- modified_filenames[!file.exists(modified_filenames)]
    existing_files_modified <- modified_filenames[file.exists(modified_filenames)]
    existing_keywords_without_ci <- missing_keywords[file.exists(modified_filenames)]
    
    if (length(missing_files_modified) > 0) {
      log_warn("The following files were not found:\n", paste(missing_files_modified, collapse = "\n"))
      log_error("Not possible to plot:\n", paste(keywords, collapse = "\n"))
      stop("Not possible to plot because not all files were found.")
      existing_files <- c()
    } else {
      existing_files <- c(existing_files, existing_files_modified)
      existing_keywords <- c(existing_keywords, existing_keywords_without_ci)
    }
  }
  
  all_data <- lapply(existing_files, function(file) {
    data <- read.csv(file)
    file_label <- existing_keywords[existing_files==file]
    data$condition <- file_label
    return(data)
  })
  common_columns <- Reduce(intersect, lapply(all_data, colnames))
  all_data_common <- lapply(all_data, function(df) df[, common_columns, drop = FALSE])
  combined_data <- do.call(rbind, all_data_common)
  
  channels <- unique(sapply(strsplit(colnames(combined_data), "_"), `[`, 1))
  channels <- channels[!grepl("^Time", channels)]
  channels <- channels[!grepl("^condition", channels)]
  
  if (length(intersect(chosen_channels, channels)) > 0) {
    channels <- intersect(chosen_channels, channels)
    log_info("Only plotting chosen channel(s): {paste(chosen_channels, collapse = ", ")}.")
  }
  
  for (channel in channels) {
    
    channel_data <- data.frame(
      Time = combined_data$Time..s.,
      Mean = combined_data[[paste0(channel, "_mean..uV.")]],
      Condition = combined_data$condition
    )
    
    if (plot_ci) {
      channel_data$Lower <- combined_data[[paste0(channel, "_CI_lower..uV.")]]
      channel_data$Upper <- combined_data[[paste0(channel, "_CI_upper..uV.")]]
    }
    
    if (!is.null(labels)) {
      label_map <- setNames(labels, keywords)
      channel_data$Condition <- factor(label_map[channel_data$Condition], levels = labels)
    } else {
      channel_data$Condition <- factor(channel_data$Condition, levels = unique(channel_data$Condition))
    }
    
    if (legend_position=="right") {
      right_space <- 5
    } else {
      right_space <- 15
    }
    
    if (title_string != "" && !grepl(" at$", title_string)) {
      title_string <- paste0(title_string, " at")
    }
    
    p <- ggplot(channel_data, aes(x = Time, y = Mean, color = Condition, group = Condition)) 
    
    if (add_0) {
      p <- p + geom_vline(xintercept = 0, linetype = "dashed", color = "lightgray", linewidth = 0.5)
    }
    
    if (add_shaded_areas && !is.null(final_latencies)) {
      latencies_to_plot <- final_latencies
      
      if (pick_roi == "all") {
        if (channel == "FCz") {
          latencies_to_plot <- subset(final_latencies, !grepl("^p3b", final_latencies$Name))
        } else if (channel == "Pz") {
          latencies_to_plot <- subset(final_latencies, grepl("^p3b", final_latencies$Name))
        } else {
          latencies_to_plot <- data.frame()
        }
      } else {
        latencies_to_plot <- subset(final_latencies, Name %in% pick_roi)
      }
      
      for (i in 1:nrow(latencies_to_plot)) {
        latency <- latencies_to_plot$Time[i]
        p <- p + geom_rect(
          xmin = latency - 0.025, xmax = latency + 0.025,
          ymin = -Inf, ymax = Inf,
          alpha = 0.01, fill = "lightgray", inherit.aes = FALSE
        )
      }
    }
    
    p <- p + geom_line(linewidth = 1) +
      labs(
        title = if (show_title) paste(title_string, channel) else NULL,
        x = "Time (s)", y = expression("Amplitude (" * mu * "V)"), # "Amplitude (µV)"
        color = legend_title
      ) + 
      #theme_minimal() +
      theme(
        panel.grid = if (remove_grid) element_blank() else element_line(),
        legend.position = legend_position,
        plot.margin = margin(t = 5, r = right_space, b = 5, l = 5)
      ) +
      scale_x_continuous(
        limits = c(tmin, tmax), expand = c(0, 0),
        #breaks = seq(tmin, tmax, by = 0.1),
        breaks = seq(tmin, tmax, by = 0.2),
        labels = label_number(accuracy = 0.01)
      ) 
    
    if (length(colors) >= length(unique(channel_data$Condition))) {
      p <- p + 
        scale_color_manual(values = colors[1:length(unique(channel_data$Condition))]) +
        scale_fill_manual(values = colors[1:length(unique(channel_data$Condition))])
    } else {
      p <- p + 
        scale_color_brewer(palette = palette) +
        scale_fill_brewer(palette = palette) 
    }
    
    if (plot_ci) {
      p <- p +
        geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Condition), alpha = 0.2, color = NA) +
        guides(fill = "none") 
    }
    
    if (save_plots) {
      plot_path <- file.path(output_dir, paste0(paste(c(keywords, channel, "evo."), collapse = "_"), plot_format))
      plot_path <- handle_long_filename(plot_path)
      ggsave(plot_path, plot = p, width = plot_width, height = plot_height)
    }
    if (show_plots) {
      print(p)
    }
    return(p)
  }
}

plot_topographies <- function(keyword, input_dir, output_dir, palette = "PuOr",
                              show_plots = FALSE, save_plots = TRUE, plot_format = "svg",
                              chosen_roi = c(), fill_title = NULL, title = NULL,
                              plot_width = 8, plot_height = 6, adjust_colorscale = FALSE) {
  
  mapping_file <- file.path(input_dir, "filename_mapping.csv")
  maps <- load_filename_mapping(mapping_file)
  orig2short  <- maps$orig2short
  virtual_file <- paste0(keyword, "_topo.csv")
  virtual_filepath <- file.path(input_dir, virtual_file)
  actual_file <- ifelse(virtual_file %in% names(orig2short),
                         orig2short[virtual_file],
                         virtual_file)
  actual_filepath <- file.path(input_dir, actual_file)
  
  if (!file.exists(actual_filepath)){
    stop(sprintf("%s does not exist", actual_filepath))
  }
  
  data <- read.csv(actual_filepath)
  
  if (!"Channel" %in% colnames(data)) {
    stop("The data file must contain a 'Channel' column")
  }
  
  data <- data[!data$Channel %in% c("EOGH", "EOGV", "EKG1"), ]
  
  valid_columns <- setdiff(colnames(data), "Channel")
  if (length(chosen_roi) > 0) {
    if (!all(chosen_roi %in% valid_columns)) {
      warning("Some elements in chosen_roi are not valid column names in the data. Continuing with valid columns.")
      chosen_roi <- intersect(chosen_roi, valid_columns)
    }
  } else {
    chosen_roi <- valid_columns
  }
  
  if (is.null(fill_title)) {
    if (grepl("tfr_", keyword) && grepl("baseline_", keyword)) {
      fill_title <- "A.u."
    } else if (grepl("tfr_", keyword)) {
      fill_title <- expression("Power (" * mu * "V"^2 * ")")
    } else {
      fill_title <- expression("Amplitude (" * mu * "V)")
    }
  }
  
  for (roi in chosen_roi) {
    eeg_data <- data.frame(
      electrode = data$Channel,
      amplitude = data[[roi]]
    )
    
    p <- topoplot(data = eeg_data, palette = palette, fill_title = fill_title) 
    
    if (adjust_colorscale) {
      scale_fill <- adjust_color_scale(eeg_data$amplitude, palette)
      p <- p + scale_fill
    }
    
    if (!is.null(title)) {
      p <- p + ggtitle(title)
    }
    
    if (show_plots) {
      print(p)
    }
    
    # Save plot if requested
    if (save_plots) {
      plot_path <- file.path(output_dir, paste0(keyword, "_", roi, "_topo.", plot_format))
      plot_path <- handle_long_filename(plot_path)
      ggsave(plot_path, plot = p, width = plot_width, height = plot_height)
    }
    return(p)
  }
}

adjust_color_scale <- function(data, palette = "PuOr", oob = scales::squish) {
  power_min <- min(data, na.rm = TRUE)
  power_max <- max(data, na.rm = TRUE)
  
  if (power_min < 0 && power_max > 0) {
    # Mixed positive and negative values: symmetric limits around 0
    max_abs <- max(abs(power_min), abs(power_max))
    symmetric_limits <- c(-max_abs, max_abs)
    return(scale_fill_distiller(
      palette = palette, direction = -1, 
      limits = symmetric_limits,
      oob = oob
    ))
  } else if (power_min < 0 && power_max <= 0) {
    # Negative values only (including 0)
    return(scale_fill_gradientn(
      colors = RColorBrewer::brewer.pal(11, palette)[6:11],
      limits = c(power_min, power_max),
      oob = oob
    ))
  } else if (power_min >= 0 && power_max > 0) {
    # Positive values only (including 0)
    return(scale_fill_gradientn(
      colors = rev(RColorBrewer::brewer.pal(11, palette)[1:6]),
      limits = c(power_min, power_max),
      oob = oob
    ))
  } else {
    stop("Unexpected data range.")
  }
}

plot_tfr <- function(keyword, input_dir, output_dir, palette = "PuOr",
                     show_plots = FALSE, save_plots = TRUE, plot_format = "svg",
                     fill_title = NULL, title = NULL,
                     plot_width = 8, plot_height = 6,
                     tmin = -0.3, tmax = 0.7,
                     adjust_colorscale = FALSE){
  
  mapping_file <- file.path(input_dir, "filename_mapping.csv")
  maps <- load_filename_mapping(mapping_file)
  orig2short  <- maps$orig2short
  virtual_file <- paste0(keyword, ".csv")
  virtual_filepath <- file.path(input_dir, virtual_file)
  actual_file <- ifelse(virtual_file %in% names(orig2short),
                        orig2short[virtual_file],
                        virtual_file)
  actual_filepath <- file.path(input_dir, actual_file)
  
  if (!file.exists(actual_filepath)){
    stop(sprintf("%s does not exist", actual_filepath))
  }
  
  data <- read.csv(actual_filepath, row.names = 1)
  tfr_long <- melt(as.matrix(data), varnames = c("Frequency", "Time"), value.name = "Power")
  if (all(tfr_long$Power < 1e-5)) {
    tfr_long$Power <- tfr_long$Power * 1e12
  }
  tfr_long$Time <- sub("^X\\.", "-", tfr_long$Time)  
  tfr_long$Time <- sub("^X", "", tfr_long$Time)     
  tfr_long$Time <- as.numeric(tfr_long$Time)  
  
  if (is.null(fill_title)) {
    if (grepl("baseline_", keyword)) {
      fill_title <- "A.u."
    } else {
      fill_title <- expression("Power (" * "µV"^2 * ")")
      #fill_title = expression("Power (" * mu * "V"^2 * ")")
    } 
  }
  
  if (is.null(title)) {
    suffix <- process_suffix(keyword)$suffix
    title <- paste("TFR at", suffix)
  }
  
  if (adjust_colorscale) {
    scale_fill <- adjust_color_scale(tfr_long$Power, palette)
  } else {
    scale_fill <- scale_fill_distiller(palette = palette, direction = -1)
  }
  
  p <- ggplot(tfr_long, aes(x = Time, y = Frequency, fill = Power)) +
    geom_tile() +
    scale_fill +
    labs(
      title = title,
      x = "Time (s)", 
      y = "Frequency (Hz)", 
      
      fill = fill_title
    ) +
    #theme_minimal() +
    theme(
      panel.grid = element_blank(),                  # Remove the grid
      panel.background = element_rect(fill = "white", color = NA),  # White background
      axis.title.x = element_text(margin = margin(t = 10)),  # Adjust x-axis title position
      axis.title.y = element_text(margin = margin(r = 10)),  # Adjust y-axis title position
      axis.text.x = element_text(margin = margin(t = -5)),   # Move x-axis tick labels closer
      legend.position = "right",                          # Position legend to the right
      legend.margin = margin(l = 0)    
    ) +
    scale_x_continuous(
      limits = c(tmin, tmax), expand = c(0, 0), 
      breaks = seq(tmin, tmax, by = 0.1),                 # Set custom tick marks every 0.1
      labels = label_number(accuracy = 0.01)             # Format ticks to two decimal places
    )
  
  if (show_plots) {
    print(p)
  }
  if (save_plots) {
    plot_path <- file.path(output_dir, paste0(keyword, ".", plot_format))
    plot_path <- handle_long_filename(plot_path)
    ggsave(plot_path, plot = p, width = plot_width, height = plot_height)
  }
  return(p)
}

get_evoked_conditions <- function(GA_path){
  #files1 <- list.files(GA_path, pattern = "_Channels\\.csv$", full.names = TRUE)
  #files2 <- list.files(GA_path, pattern = "_Channels_CI\\.csv$", full.names = TRUE)
  mapping_file <- file.path(GA_path, "filename_mapping.csv")
  files_df <- list_mapped_files(GA_path, db_path = mapping_file, extension = "csv")
  vf <- files_df$virtual
  files1 <- vf[grepl("_Channels\\.csv$",    vf)]
  files2 <- vf[grepl("_Channels_CI\\.csv$", vf)]
  remove_extensions <- function(file_names) {
    gsub("(_Channels\\.csv|_Channels_CI\\.csv)$", "", file_names)
  }
  cleaned_list1 <- remove_extensions(files1)
  cleaned_list2 <- remove_extensions(files2)
  unique_files <- unique(c(cleaned_list1, cleaned_list2))
  conditions <- basename(unique_files)
  conditions <- conditions[!grepl("dif_", conditions)]
  return(conditions)
}

compare_first_word <- function(str1, str2, prefixes) {
  prefix_pattern <- paste(prefixes, collapse = "|")
  str1 <- sub(paste0("^(dif_|", prefix_pattern, ")*"), "", str1)
  str2 <- sub(paste0("^(dif_|", prefix_pattern, ")*"), "", str2)
  first_word1 <- strsplit(str1, "_")[[1]][1]
  first_word2 <- strsplit(str2, "_")[[1]][1]
  return(first_word1 == first_word2)
}

extract_prefix <- function(file_name, prefixes) {
  prefix_pattern <- paste(prefixes, collapse = "|")
  matched_prefix <- regmatches(file_name, regexpr(paste0("^(", prefix_pattern, ")*"), file_name))
  remaining <- sub(paste0("^(", prefix_pattern, ")*"), "", file_name)
  list(prefix = matched_prefix, remaining = remaining)
}

process_dif_file <- function(dif_file, GA_path, prefixes) {
  
  conditions <- get_evoked_conditions(GA_path)
  
  # Extract prefix and remaining part of dif_file
  parts <- extract_prefix(dif_file, prefixes)
  prefix <- parts$prefix
  file_base <- sub("^dif_", "", parts$remaining)  # Remove dif_ from file_base
  
  # Initialize result variables
  condition1 <- NA
  condition2 <- NA
  
  # Define candidate conditions based on prefix
  if (prefix != "") {
    candidate_conditions <- grep(paste0("^", prefix), conditions, value = TRUE)
  } else {
    # Filter conditions that do not start with any prefix
    prefix_pattern <- paste(prefixes, collapse = "|")
    candidate_conditions <- grep(paste0("^(?!", prefix_pattern, ")"), conditions, value = TRUE, perl = TRUE)
  }
  
  # Check for double underscore in file_base
  if (grepl("__", file_base)) {
    # Split by double underscore
    conditions_pair <- strsplit(file_base, "__")[[1]]
    condition1 <- conditions_pair[1]
    condition2 <- conditions_pair[2]
    
    # Ensure condition1 is relevant by matching all words in condition1 with file_base
    condition1_base <- sub(paste0("^(dif_|", paste(prefixes, collapse = "|"), ")*"), "", condition1)
    if (!all(sapply(strsplit(condition1_base, "_")[[1]], function(word) grepl(word, file_base)))) {
      condition1 <- NA
      condition2 <- NA
    }
    
    # Ensure condition1 and condition2 are compatible (same first word)
    if (!is.na(condition1) && !is.na(condition2)) {
      if (!compare_first_word(condition1, condition2, prefixes)) {
        condition1 <- NA
        condition2 <- NA
      }
    }
  } else {
    # Match by the first word
    for (condition_key1 in candidate_conditions) {
      condition_key1_base <- sub(paste0("^(dif_|", paste(prefixes, collapse = "|"), ")*"), "", condition_key1)
      
      if (!all(sapply(strsplit(condition_key1_base, "_")[[1]], function(word) grepl(word, file_base)))) {
        next
      }
      
      if (!grepl(condition_key1_base, file_base, fixed = TRUE)) {
        next
      }
      
      # Attempt to determine condition2
      condition_key2_candidate <- substring(file_base, nchar(condition_key1_base) + 2)  # Skip underscore
      
      # Add prefix back to key base before checking in candidate conditions
      if (prefix != "") {
        condition_key2_candidate <- paste0(prefix, condition_key2_candidate)
      }
      
      if (condition_key2_candidate %in% candidate_conditions) {
        if (compare_first_word(condition_key1, condition_key2_candidate, prefixes)) {
          condition1 <- condition_key1
          condition2 <- condition_key2_candidate
          
          # If condition1 and condition2 are the same, continue to next match
          if (condition1 == condition2) {
            next
          }
          
          break
        }
      }
    }
    
    if (is.na(condition1) || is.na(condition2)) {
      warning(paste("Could not resolve condition pair for", dif_file))
      log_warn(paste("Could not resolve condition pair for", dif_file))
    }
  }
  
  # Extract labels for the difference if conditions are resolved
  label_condition1 <- NA
  label_condition2 <- NA
  if (!is.na(condition1) && !is.na(condition2)) {
    condition1_words <- strsplit(condition1, "_")[[1]]
    condition2_words <- strsplit(condition2, "_")[[1]]

    # Compare words by order and remove duplicates
    unique_condition1 <- unique(mapply(
      function(w1, w2) if (w1 != w2 && w1 != "") w1 else "", 
      condition1_words,
      condition2_words
    ))
    
    unique_condition2 <- unique(mapply(
      function(w1, w2) if (w1 != w2 && w2 != "") w2 else "", 
      condition1_words,
      condition2_words
    ))
    
    label_condition1 <- paste(unique_condition1[!is.na(unique_condition1)], collapse = " ")
    label_condition2 <- paste(unique_condition2[!is.na(unique_condition2)], collapse = " ")
    
    label_condition1 <- trimws(label_condition1, which = "both")
    label_condition2 <- trimws(label_condition2, which = "both")
  }
  
  # Return results
  list(
    dif_file = dif_file,
    prefix = prefix,
    condition1 = condition1,
    condition2 = condition2,
    label_condition1 = label_condition1,
    label_condition2 = label_condition2
  )
}

decode_file_name <- function(file_name, add_suffixes = TRUE, Rpeak_in_name = FALSE) {
  
  process_main_info <- function(main_part) {
    main_part <- gsub("_2groups|_3groups|sysDiaFeedb|outcome|learning", "", main_part)

    main_part <- gsub("changeIBI(\\d+)_\\d+", "changeIBI\\1", main_part, perl = TRUE)
    
    # Handle phrases like expectation_high, moving high/low/medium before the associated word
    main_part <- gsub("([[:alnum:]]+)(_high|_low|_medium)", "\\2 \\1", main_part)
    
    if (Rpeak_in_name) {
      main_part <- gsub("changeIBI(\\d+)", "change of R-peak \\1", main_part, perl = TRUE) # Format changeIBI
      main_part <- gsub("IBI(\\d+)", "IBI of R-Peak \\1", main_part, perl = TRUE) # Format IBI
    } else {
      main_part <- gsub("IBI(\\d+)", "IBI \\1", main_part, perl = TRUE) 
      main_part <- gsub("changeIBI", "change of IBI", main_part, perl = TRUE) 
    }
    
    # Replace underscores with spaces and clean up extra spaces
    main_part <- gsub("_", " ", main_part)
    main_part <- gsub("[[:space:]]+", " ", main_part)
    main_part <- trimws(main_part)
    
    main_part
  }
  
  if (grepl("dif", file_name)) {
    log_info("Not decoding difference strings, use process_dif_file first")
    return(NULL)
  }
  
  if (grepl("^long_|^tfr_", file_name)) {
    file_name <- gsub("^long_|^tfr_", "", file_name)
  }
  
  pseudo_prefix <- ""
  if (grepl("pseudo_", file_name)) {
    pseudo_prefix <- "Pseudotrials for "
    file_name <- gsub("pseudo_", "", file_name)
  }
  
  baseline_suffix <- ""
  if (grepl("baseline_", file_name)) {
    baseline_suffix <- ", baseline-corrected"
    file_name <- gsub("baseline_", "", file_name)
  }
  
  corrected_suffix <- ""
  if (grepl("_corrected_", file_name)) {
    corrected_suffix <- ", pseudotrial-corrected"
    file_name <- sub("_corrected_.*", "", file_name)
  }
  
  systole_diastole_suffix <- ""
  if (grepl("(?i)systole|diastole", file_name, perl = TRUE)) {
    if (grepl("(?i)systole", file_name, perl = TRUE)) {
      file_name <- sub("(?i)systole", "", file_name, perl = TRUE) # Remove matched systole
      systole_diastole_suffix <- " in systole"
    } else if (grepl("(?i)diastole", file_name, perl = TRUE)) {
      file_name <- sub("(?i)diastole", "", file_name, perl = TRUE) # Remove matched diastole
      systole_diastole_suffix <- " in diastole"
    }
  }
  
  miniblock_suffix <- ""
  miniblock_pattern <- "(?i)(miniblock_reversal|miniblock_next|miniblock_first|miniblock_last)"
  if (grepl(miniblock_pattern, file_name, perl = TRUE)) {
    match <- regmatches(file_name, regexpr(miniblock_pattern, file_name, perl = TRUE))
    file_name <- sub(miniblock_pattern, "", file_name, perl = TRUE) 
    match <- gsub("miniblock_([[:alnum:]]+)", "\\1 miniblock", match)
    match <- gsub("_", " ", match)
    miniblock_suffix <- paste0(", ", tolower(match)) 
  }
  
  main_info <- process_main_info(file_name)
  
  # Ensure systole/diastole phrases and miniblock suffixes are at the end
  if (systole_diastole_suffix != "") {
    main_info <- paste0(main_info, systole_diastole_suffix)
  }
  if (miniblock_suffix != "") {
    main_info <- paste0(main_info, miniblock_suffix)
  }
  
  result <- paste0(pseudo_prefix, main_info)
  if (add_suffixes) {
    result <- paste0(result, corrected_suffix, baseline_suffix)
  }
  
  result <- gsub("[[:space:]]+", " ", result) # Remove extra spaces
  
  if (grepl("^ in systole", result)) {
    result <- sub("^ in systole", "Systole", result)
  } else if (grepl("^ in diastole", result)) {
    result <- sub("^ in diastole", "Diastole", result)
  }

  # Capitalize only the first word of the result
  result <- sub("^(.)", "\\U\\1", result, perl = TRUE)
  result <- gsub("r-peak", "R-peak", result, fixed = TRUE) 
  result <- gsub("abspe", "absPE", result, fixed = TRUE) 
  result <- gsub("Abspe", "AbsPE", result, fixed = TRUE) 
  result <- gsub(" pe ", " PE ", result, fixed = TRUE) 
  result <- gsub("Pe ", "PE ", result, fixed = TRUE) 
  result <- gsub("ChoiceApproach 0", "Avoid", result, fixed = TRUE)
  result <- gsub("ChoiceApproach 1", "Approach", result, fixed = TRUE)
  result <- gsub("ChoiceApproach", "Approach", result, fixed = TRUE)
  
  keywords <- c("Win", "Lose", "Approach", "Avoid", "expectation", "absPE", "miniblock")
  pattern <- paste0("(?i)\\b(", paste(keywords, collapse = "|"), ")\\b(?! in|, |$)")
  result <- gsub(pattern, "\\1 and", result, perl = TRUE)
  
  result <- gsub("\\bfor in\\b", "for", result, perl = TRUE)
  result <- gsub("\\bfor ,\\s*", "for ", result, perl = TRUE)
  result <- gsub("( and| for )$", "", result)
  result <- gsub("^Pseudo(?!trials)", "Pseudotrials", result, perl = TRUE)
  result <- gsub("afterFeedbAveIBI", "average IBI after feedback", result, perl = TRUE)
  result <- gsub("beforeFeedbAveIBI", "average IBI before feedback", result, perl = TRUE)
  result <- gsub("^Corrected", "Pseudotrial-corrected", result)
  result <- gsub("^, pseudotrial_corrected", "Pseudotrial-corrected", result)
  result <- gsub("and corrected$", ", pseudotrial-corrected", result)
  
  number_mapping <- c(
    "1" = "-4", "2" = "-3", "3" = "-2", "4" = "-1",
    "5" = "1", "6" = "2", "7" = "3", "8" = "4"
  )
  
  matches <- gregexpr("(R-peak|IBI) (\\d)", result, perl = TRUE)
  extracted <- regmatches(result, matches)  
  
  replaced <- lapply(extracted, function(matches) {
    sapply(matches, function(match) {
      prefix <- sub(" (\\d+)$", "", match)
      num <- sub(".* ", "", match)
      if (num %in% names(number_mapping)) {
        paste0(prefix, " ", number_mapping[[num]])  
      } else {
        match 
      }
    })
  })
  
  regmatches(result, matches) <- replaced
  
  return(result)
}

clean_condition <- function(condition1, condition2, label1, label2, condition_name = FALSE, remove_condition_name = FALSE, Rpeak_in_name = FALSE) {
  
  condition1 <- gsub("_2groups|_3groups", "", condition1) 
  condition2 <- gsub("_2groups|_3groups", "", condition2) 
  condition1 <- gsub("changeIBI(\\d+)_\\d+", "changeIBI\\1", condition1, perl = TRUE)
  condition2 <- gsub("changeIBI(\\d+)_\\d+", "changeIBI\\1", condition2, perl = TRUE)
  
  condition1_words <- unlist(strsplit(condition1, "_"))
  condition2_words <- unlist(strsplit(condition2, "_"))
  label1_words <- unlist(strsplit(label1, " "))
  label2_words <- unlist(strsplit(label2, " "))
  
  positions1 <- which(condition1_words %in% label1_words)
  positions2 <- which(condition2_words %in% label2_words)
  overlapping_positions <- intersect(positions1, positions2)
  
  if (length(overlapping_positions) == 0) {
    stop("No overlapping positions found for label words in conditions.")
  }
  
  condition_names <- c()
  cleaned_condition1 <- condition1_words
  cleaned_condition2 <- condition2_words
  
  for (pos in overlapping_positions) {
    # Extract condition name from the word at -1 position
    if (pos > 1) {
      condition_name1 <- condition1_words[pos - 1]
      condition_name2 <- condition2_words[pos - 1]
      
      if (condition_name1 != condition_name2) {
        stop("Condition names do not match at overlapping positions.")
      }
      
      condition_names <- c(condition_names, condition_name1)
    }
    
    cleaned_condition1[pos] <- ""
    cleaned_condition2[pos] <- ""
  }
  
  # Ensure cleaned conditions are identical
  cleaned_condition1 <- cleaned_condition1[cleaned_condition1 != ""]
  cleaned_condition2 <- cleaned_condition2[cleaned_condition2 != ""]
  
  if (!identical(cleaned_condition1, cleaned_condition2)) {
    stop("Cleaned conditions are not identical.")
  }
  
  cleaned_condition <- paste(cleaned_condition1, collapse = "_")
  condition_name_result <- paste(unique(condition_names), collapse = " and ")
  
  if (remove_condition_name && !is.null(condition_name_result)) {
    condition_name_words <- unlist(strsplit(condition_name_result, " and "))
    for (word in condition_name_words) {
      cleaned_condition <- gsub(word, "", cleaned_condition)
    }
    cleaned_condition <- gsub("_+", "_", cleaned_condition)  # Clean multiple underscores
    cleaned_condition <- gsub("^_|_$", "", cleaned_condition)  # Remove leading/trailing underscores
  }
  
  condition_name_result <- if (condition_name) condition_name_result else NULL
  
  if (!is.null(condition_name_result)) {
    
    if (Rpeak_in_name) {
      condition_name_result <- gsub("changeIBI(\\d+)", "change of R-peak \\1", condition_name_result, perl = TRUE) # Format changeIBI
      condition_name_result <- gsub("IBI(\\d+)", "IBI of R-Peak \\1", condition_name_result, perl = TRUE) # Format IBI
    } else {
      condition_name_result <- gsub("IBI(\\d+)", "IBI \\1", condition_name_result, perl = TRUE) 
      condition_name_result <- gsub("changeIBI", "change of IBI", condition_name_result, perl = TRUE) 
    }
    
    number_mapping <- c(
      "1" = "-4", "2" = "-3", "3" = "-2", "4" = "-1",
      "5" = "1", "6" = "2", "7" = "3", "8" = "4"
    )
    
    matches <- gregexpr("(R-peak|IBI) (\\d)", condition_name_result, perl = TRUE)
    extracted <- regmatches(condition_name_result, matches)  
    
    replaced <- lapply(extracted, function(matches) {
      sapply(matches, function(match) {
        prefix <- sub(" (\\d+)$", "", match)
        num <- sub(".* ", "", match)
        if (num %in% names(number_mapping)) {
          paste0(prefix, " ", number_mapping[[num]])  
        } else {
          match 
        }
      })
    })
    
    regmatches(condition_name_result, matches) <- replaced
    
    condition_name_result <- gsub("choiceApproach", "Approach", condition_name_result, fixed = TRUE)
    condition_name_result <- gsub("sysDiaFeedb", "Cardiac Phase", condition_name_result, fixed = TRUE)
    condition_name_result <- gsub("afterFeedbAveIBI", "average IBI after feedback", condition_name_result, perl = TRUE)
    condition_name_result <- gsub("beforeFeedbAveIBI", "average IBI before feedback", condition_name_result, perl = TRUE)
    
    condition_name_result <- sub("^(.)", "\\U\\1", condition_name_result, perl = TRUE)
  }
  
  if (condition_name) {
    return(list(cleaned_condition = cleaned_condition, condition_name = condition_name_result))
  }else{
    return(list(cleaned_condition = cleaned_condition))
  }
}

process_suffix <- function(file_name) {
  if (grepl("__", file_name)) {
    parts <- strsplit(file_name, "__")[[1]]
    suffix <- tail(parts, n = 1)
    main_part <- paste(head(parts, -1), collapse = "__")
    return(list(main_part = main_part, suffix = suffix))
  } else {
    return(list(main_part = file_name, suffix = NA))
  }
}

handle_long_filename <- function(file_path,
                                 max_length = 200,
                                 db_path    = NULL) {
  # (1) split directory & filename
  dir_name  <- dirname(file_path)
  file_name <- basename(file_path)
  
  # (2) only shorten if it’s too long
  if (nchar(file_name) <= max_length) {
    return(file_path)
  }
  
  # (3) compute an 8-char MD5 of the stem
  stem <- tools::file_path_sans_ext(file_name)
  ext  <- tools::file_ext(file_name)
  hash8 <- substr(digest::digest(stem, algo = "md5"), 1, 8)
  new_name <- if (nzchar(ext)) {
    paste0(hash8, ".", ext)
  } else {
    hash8
  }
  safe_path <- file.path(dir_name, new_name)
  
  # (4) record mapping
  map_file <- if (!is.null(db_path)) db_path else file.path(dir_name, "filename_mapping.csv")
  
  # Check if the mapping file exists and already contains the filename
  if (file.exists(map_file)) {
    existing_map <- read.csv(map_file, stringsAsFactors = FALSE)
    match_idx <- match(file_name, existing_map$original)
    if (!is.na(match_idx)) {
      # Already mapped: return the corresponding shortened name
      return(file.path(dir_name, existing_map$shortened[match_idx]))
    }
  }
  
  new_row  <- data.frame(original = file_name,
                         shortened = new_name,
                         stringsAsFactors = FALSE)
  if (!file.exists(map_file)) {
    # write header + first row
    write.csv(new_row, map_file, row.names = FALSE)
  } else {
    # append without header
    write.table(new_row, map_file, append = TRUE, sep = ",",
                row.names = FALSE, col.names = FALSE)
  }
  
  log_info(sprintf("Filename '%s' too long (%d chars) → using '%s'",
                  file_name, nchar(file_name), new_name))
  return(safe_path)
}

load_filename_mapping <- function(db_path) {
  if (!file.exists(db_path)) {
    return(list(orig2short = character(0), short2orig = character(0)))
  }
  df <- read.csv(db_path, stringsAsFactors = FALSE)
  orig2short <- setNames(df$shortened, df$original)
  short2orig <- setNames(df$original, df$shortened)
  list(orig2short = orig2short, short2orig = short2orig)
}

list_mapped_files <- function(dir,
                              db_path = NULL,
                              extension = NULL) {
  mapping_file <- if (!is.null(db_path)) db_path else file.path(dir, "filename_mapping.csv")
  maps <- load_filename_mapping(mapping_file)
  short2orig <- maps$short2orig
  
  # get filenames on disk (optionally filter by extension)
  pattern <- if (!is.null(extension)) paste0("\\.", extension, "$") else NULL
  actual_files <- list.files(dir, pattern = pattern, full.names = FALSE)
  
  # map back to original where needed
  virtual_files <- ifelse(actual_files %in% names(short2orig),
                          short2orig[actual_files],
                          actual_files)
  
  data.frame(actual  = actual_files,
             virtual = virtual_files,
             stringsAsFactors = FALSE)
}

process_filenames <- function(GA_path, use_case) {

  remove_extensions <- function(file_names, pattern) {
    cleaned_files <- gsub(pattern, "", file_names)
    basename(cleaned_files)
  }
  
  mapping_file <- file.path(GA_path, "filename_mapping.csv")
  files_df <- list_mapped_files(GA_path, db_path = mapping_file, extension = "csv")
  vf <- files_df$virtual
  
  # Initialize lists for results
  dif_files <- list()
  other_files <- list()
  suffix_list <- list()
  dif_suffix_list <- list()
  
  # Perform operations based on use_case
  if (use_case == "evoked") {
    files1 <- vf[grepl("_Channels\\.csv$",    vf)]
    files2 <- vf[grepl("_Channels_CI\\.csv$", vf)]
    #files1 <- list.files(GA_path, pattern = "_Channels\\.csv$", full.names = TRUE)
    #files2 <- list.files(GA_path, pattern = "_Channels_CI\\.csv$", full.names = TRUE)
    unique_files <- unique(c(remove_extensions(files1, "_Channels\\.csv$"), 
                             remove_extensions(files2, "_Channels_CI\\.csv$")))
    dif_files <- unique_files[grepl("dif_", unique_files)]
    other_files <- unique_files[!grepl("dif_", unique_files)]
    
  } else if (use_case == "topo") {
    #files <- list.files(GA_path, pattern = "_topo\\.csv$", full.names = TRUE)
    files <- vf[grepl("_topo\\.csv$", vf)]
    files <- remove_extensions(files, "_topo\\.csv$")
    files <- files[!grepl("tfr_", files)]
    dif_files <- files[grepl("dif_", files)]
    other_files <- files[!grepl("dif_", files)]
    
  } else if (use_case == "tfr_topo") {
    #files <- list.files(GA_path, pattern = "_topo\\.csv$", full.names = TRUE)
    files <- vf[grepl("_topo\\.csv$", vf)]
    files <- remove_extensions(files, "_topo\\.csv$")
    files <- files[grepl("tfr_", files)]
    dif_files <- files[grepl("dif_", files)]
    other_files <- files[!grepl("dif_", files)]
    processed_other_files <- lapply(other_files, process_suffix)
    other_files <- sapply(processed_other_files, `[[`, "main_part")
    suffix_list <- sapply(processed_other_files, `[[`, "suffix")
    processed_dif_files <- lapply(dif_files, process_suffix)
    dif_files <- sapply(processed_dif_files, `[[`, "main_part")
    dif_suffix_list <- sapply(processed_dif_files, `[[`, "suffix")
    
  } else if (use_case == "tfr_channel") {
    #files <- list.files(GA_path, pattern = "^tfr_.*\\.csv$", full.names = TRUE)
    files <- vf[grepl("^tfr_.*\\.csv$", vf)]
    files <- files[!grepl("_topo\\.csv$", files)]
    files <- remove_extensions(files, "\\.csv$")
    dif_files <- files[grepl("dif_", files)]
    other_files <- files[!grepl("dif_", files)]
    processed_other_files <- lapply(other_files, process_suffix)
    other_files <- sapply(processed_other_files, `[[`, "main_part")
    suffix_list <- sapply(processed_other_files, `[[`, "suffix")
    processed_dif_files <- lapply(dif_files, process_suffix)
    dif_files <- sapply(processed_dif_files, `[[`, "main_part")
    dif_suffix_list <- sapply(processed_dif_files, `[[`, "suffix")
  } else {
    stop("Invalid use_case. Use one of: 'evoked', 'topo', 'tfr_channel', 'tfr_topo'.")
  }
  
  # Return the results as a list
  return(list(dif_files = dif_files, 
              dif_suffix_list = dif_suffix_list, 
              other_files = other_files, 
              suffix_list = suffix_list))
}

generate_difference_titles <- function(condition_names, label1, label2, add_title = TRUE, add_suffix = FALSE, titles = NULL, dif_suffix = NULL) {
  difference_title <- vector("character", length(condition_names))
  
  for (i in seq_along(condition_names)) {
    condition_name <- condition_names[i]
    
    if (condition_name %in% c("Cardiac Phase", "Outcome")) {
      condition_name <- ""
    } else if (condition_name == "Miniblock") {
      condition_name <- "miniblocks"
    }
    
    if (nzchar(condition_name)) { 
      condition_name <- paste0(" ", tolower(substr(condition_name, 1, 1)), substr(condition_name, 2, nchar(condition_name)))
    }
    
    label1_low <- paste0(tolower(substr(label1[i], 1, 1)), substr(label1[i], 2, nchar(label1[i])))
    label2_low <- paste0(tolower(substr(label2[i], 1, 1)), substr(label2[i], 2, nchar(label2[i])))
    
    if (add_title && !is.null(titles) && nzchar(titles[i])) {
      prefix <- paste0(titles[i], ". ")
    } else {
      prefix <- ""
    }
    
    if (add_suffix && !is.null(dif_suffix)) {
      suffix <- paste0(condition_name, ". ", dif_suffix[i])
    } else {
      suffix <- condition_name
    }
    
    # Special case for "approach"
    if (condition_name == " approach") {
      difference_title[i] <- paste0(prefix, "Difference of approach and avoid", suffix)
    } else {
      difference_title[i] <- paste0(prefix, "Difference of ", label1_low, " and ", label2_low, suffix)
    }
  }
  return(difference_title)
}

