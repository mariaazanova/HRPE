library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(reshape2)

plot_radial_line <- function(df, 
                             angle_var = "phaseFeedb",    # Angle variable in degrees (0-360)
                             measure_var = "IBI5",        # Continuous measure variable (e.g., IBI5)
                             bin_size = 10,               # Bin size in degrees
                             aggregator = mean,           # Aggregation function (e.g. mean, median)
                             normalize = FALSE,           # Normalize aggregated measure to [0,1]?
                             smooth = FALSE,              # Apply LOESS smoothing?
                             add_ci = FALSE,              # Option to add CI shading (default FALSE)
                             say_ibi = FALSE,             # Replace "RR interval" with "IBI" in titles?
                             groups = NULL,               # Optional grouping variable(s)
                             title_string = "",           # Custom plot title
                             legend_title = "",           # Custom legend title
                             show_title = TRUE,
                             show_legend_title = TRUE,
                             show_legend = TRUE,
                             vis_path,                    # Path to save figure
                             plot_width = 8,
                             plot_height = 6,
                             show_plots = FALSE,
                             save_plots = TRUE,
                             plot_format = "svg",
                             dpi = 150,
                             colors = NULL,
                             expand_center = 0.25, 
                             move_border = -0.1,
                             add_grid = FALSE,
                             output_as_one = FALSE,
                             scale = NULL,                #  if provided, a two-element vector (min, max) for tick labels
                             num_ticks = NULL) {           # desired number of tick labels
  
  # 1. Validate required columns
  req_cols <- c(angle_var, measure_var)
  missing_cols <- setdiff(req_cols, names(df))
  if(length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse=", "))
  }
  
  # 2. Process grouping variables: convert to factors and drop NA/empty levels
  if(!is.null(groups)) {
    for(g in groups) {
      df[[g]] <- as.factor(df[[g]])
      df <- df %>% filter(!is.na(.data[[g]]) & .data[[g]] != "")
      df[[g]] <- droplevels(df[[g]])
    }
  }
  if(!is.null(groups) && length(groups) > 0) {
    main_group <- groups[1]
    extra_groups <- if(length(groups) > 1) groups[-1] else NULL
  } else {
    main_group <- NULL
    extra_groups <- NULL
  }
  
  final_legend_title <- if(legend_title != "") {
    legend_title
  } else if(!is.null(main_group)) {
    map_variable_name(main_group)
  } else {
    ""
  }
  if(!show_legend_title) final_legend_title <- NULL
  
  # 3. Bin the angle variable
  angle_breaks <- seq(0, 360, by = bin_size)
  if(tail(angle_breaks,1) < 360) angle_breaks <- c(angle_breaks, 360)
  
  df_binned <- df %>% 
    mutate(angle_bin = cut(.data[[angle_var]], breaks = angle_breaks, 
                           include.lowest = TRUE, right = FALSE))
  
  # 4. Aggregate the measure within each bin and for all groups (if provided)
  grouping_vars <- c("angle_bin")
  if(!is.null(groups) && length(groups) > 0) {
    grouping_vars <- c(grouping_vars, groups)
  }
  
  if(add_ci) {
    df_summary <- df_binned %>% 
      group_by(across(all_of(grouping_vars))) %>%
      summarise(n = n(),
                mean_val = mean(.data[[measure_var]], na.rm = TRUE),
                sd_val = sd(.data[[measure_var]], na.rm = TRUE),
                lower_val = mean_val - qt(0.975, df = n - 1)*sd_val/sqrt(n),
                upper_val = mean_val + qt(0.975, df = n - 1)*sd_val/sqrt(n),
                .groups = "drop") %>%
      mutate(measure_value = mean_val)
  } else {
    df_summary <- df_binned %>% 
      group_by(across(all_of(grouping_vars))) %>%
      summarise(measure_value = aggregator(.data[[measure_var]], na.rm = TRUE),
                .groups = "drop")
  }
  
  # Convert bin factor to numeric center (e.g., "[0,10)" -> 5)
  get_bin_center <- function(bin_label) {
    low_bound <- as.numeric(sub("[\\(\\[](.*),.*", "\\1", bin_label))
    low_bound + bin_size/2
  }
  df_summary <- df_summary %>% 
    mutate(angle_center = sapply(as.character(angle_bin), get_bin_center))
  
  # 5. Normalize the aggregated measure if requested
  if(normalize) {
    rng <- range(df_summary$measure_value, na.rm = TRUE)
    if(diff(rng) == 0) {
      df_summary$measure_value <- 0
      if(add_ci) {
        df_summary$lower_val <- 0
        df_summary$upper_val <- 0
      }
    } else {
      df_summary$measure_value <- (df_summary$measure_value - rng[1])/(rng[2]-rng[1])
      if(add_ci) {
        df_summary$lower_val <- (df_summary$lower_val - rng[1])/(rng[2]-rng[1])
        df_summary$upper_val <- (df_summary$upper_val - rng[1])/(rng[2]-rng[1])
      }
    }
  }
  
  # 6. Set plot value to the aggregated measure (and CI if applicable)
  df_summary <- df_summary %>% 
    mutate(plot_value = measure_value)
  if(add_ci) {
    df_summary <- df_summary %>% 
      mutate(plot_lower = lower_val,
             plot_upper = upper_val)
  }
  
  # 7. Apply smoothing if requested.
  # Updated: group by all grouping variables (i.e., groups) so that extra groups are preserved.
  if(smooth) {
    if(!is.null(groups) && length(groups) > 0) {
      if(add_ci) {
        df_smoothed <- df_summary %>% 
          group_by(across(all_of(groups))) %>% 
          do({
            mod_val <- loess(plot_value ~ angle_center, data = ., span = 0.3)
            mod_lower <- loess(plot_lower ~ angle_center, data = ., span = 0.3)
            mod_upper <- loess(plot_upper ~ angle_center, data = ., span = 0.3)
            new_angles <- seq(0, 360, length.out = 200)
            data.frame(angle_center = new_angles,
                       plot_value = predict(mod_val, newdata = data.frame(angle_center = new_angles)),
                       plot_lower = predict(mod_lower, newdata = data.frame(angle_center = new_angles)),
                       plot_upper = predict(mod_upper, newdata = data.frame(angle_center = new_angles)))
          }) %>% ungroup()
      } else {
        df_smoothed <- df_summary %>% 
          group_by(across(all_of(groups))) %>% 
          do({
            mod_val <- loess(plot_value ~ angle_center, data = ., span = 0.3)
            new_angles <- seq(0, 360, length.out = 200)
            data.frame(angle_center = new_angles,
                       plot_value = predict(mod_val, newdata = data.frame(angle_center = new_angles)))
          }) %>% ungroup()
      }
    } else {  # No grouping variables provided.
      if(add_ci) {
        mod_val <- loess(plot_value ~ angle_center, data = df_summary, span = 0.3)
        mod_lower <- loess(plot_lower ~ angle_center, data = df_summary, span = 0.3)
        mod_upper <- loess(plot_upper ~ angle_center, data = df_summary, span = 0.3)
        new_angles <- seq(0, 360, length.out = 200)
        df_smoothed <- data.frame(angle_center = new_angles,
                                  plot_value = predict(mod_val, newdata = data.frame(angle_center=new_angles)),
                                  plot_lower = predict(mod_lower, newdata = data.frame(angle_center=new_angles)),
                                  plot_upper = predict(mod_upper, newdata = data.frame(angle_center=new_angles)))
      } else {
        mod_val <- loess(plot_value ~ angle_center, data = df_summary, span = 0.3)
        new_angles <- seq(0, 360, length.out = 200)
        df_smoothed <- data.frame(angle_center = new_angles,
                                  plot_value = predict(mod_val, newdata = data.frame(angle_center=new_angles)))
      }
    }
    df_plot <- df_smoothed
  } else {
    df_plot <- df_summary
  }
  
  # 8. Convert polar coordinates with rotation (so that 360 is at the top)
  df_plot <- df_plot %>%
    mutate(bin_radians = (angle_center + 90) * pi / 180,
           x = plot_value * cos(bin_radians),
           y = plot_value * sin(bin_radians))
  
  # Helper function to create the radial plot for a given dataset and title
  create_radial_plot <- function(plot_data, plot_title) {
    p <- ggplot(plot_data, aes(x = angle_center, y = plot_value))
    if(add_ci) {
      if(!is.null(main_group)) {
        p <- p + geom_ribbon(aes_string(fill = main_group, group = main_group, ymin = "plot_lower", ymax = "plot_upper"), 
                             alpha = 0.2, color = NA)
      } else {
        p <- p + geom_ribbon(aes(ymin = plot_lower, ymax = plot_upper), 
                             fill = "gray", alpha = 0.2, color = NA)
      }
    }
    if(!is.null(main_group)) {
      p <- p + geom_line(aes(color = .data[[main_group]], group = .data[[main_group]]), size = 1) +
        geom_point(aes(color = .data[[main_group]]), size = 2)
    } else {
      default_color <- if(is.null(colors)) "gray" else colors[1]
      p <- p + geom_line(size = 1, color = default_color) +
        geom_point(size = 2, color = default_color)
    }
    
    p <- p + coord_polar(theta = "x") +
      scale_x_continuous(breaks = c(0,90,180,270), limits = c(0,360), 
                         labels = c("0°", "90°", "180°", "270°"))
    
    # --- Amplitude Tick Label Annotation ---
    # Determine tick range: use provided scale if available, otherwise use range from data.
    tick_range <- if(is.null(scale)) range(plot_data$plot_value, na.rm = TRUE) else scale
    if(is.null(num_ticks)){
      tick_effects <- pretty(tick_range)
    } else {
      tick_effects <- seq(from = tick_range[1], to = tick_range[2], length.out = num_ticks)
    }
    # Round tick labels: first try rounding to the nearest integer.
    tick_labels <- round(tick_effects, 0)
    if(length(unique(tick_labels)) == 1) {
      tick_labels <- round(tick_effects, 2)
      if(length(unique(tick_labels)) == 1) {
        stop("Error: Tick labels cannot be properly rounded; they result in a single value.")
      }
    }
    tick_positions <- tick_effects
    tick_data <- data.frame(angle = 360, y = tick_positions, label = tick_labels)
    p <- p + geom_text(data = tick_data, aes(x = angle, y = y, label = label))
    # --- End Amplitude Tick Label Annotation ---
    
    p <- p + scale_y_continuous(labels = NULL, 
                                expand = expansion(mult = c(expand_center, move_border)))
    
    # --- Title formatting block update ---
    final_plot_title <- tolower(plot_title)
    if(say_ibi){
      final_plot_title <- gsub("rr interval", "ibi", final_plot_title, ignore.case = TRUE)
    } else {
      final_plot_title <- gsub("rr interval", "rr interval", final_plot_title, ignore.case = TRUE)
    }
    final_plot_title <- gsub(" pe ", " PE ", final_plot_title, ignore.case = TRUE)
    final_plot_title <- paste0(toupper(substr(final_plot_title, 1, 1)), substr(final_plot_title, 2, nchar(final_plot_title)))
    final_plot_title <- gsub("\\bibi\\b", "IBI", final_plot_title, ignore.case = TRUE)
    final_plot_title <- gsub("\\brr\\b", "RR", final_plot_title, ignore.case = TRUE)
    # Append " (ms)" at the end of the title
    final_plot_title <- paste0(final_plot_title, " (ms)")
    # --- End of title formatting block ---
    
    p <- p + labs(title = if(show_title) final_plot_title else NULL,
                  x = NULL, y = NULL,
                  color = final_legend_title, fill = final_legend_title) +
      theme(legend.position = if(show_legend) "right" else "none",
            plot.title = element_text(hjust = 0.5))
    
    if(add_grid){
      p <- p + theme(panel.grid.major = element_line(color = "gray90"),
                     panel.grid.minor = element_line(color = "gray95"))
    }
    
    if(!is.null(main_group)) {
      lvls <- unique(plot_data[[main_group]])
      if(length(lvls) == 1) {
        chosen_colors <- if(is.null(colors)) "gray" else colors[1]
      } else {
        chosen_colors <- if(is.null(colors)) NULL else colors[seq_along(lvls)]
      }
      if(!is.null(chosen_colors)) {
        p <- p + scale_color_manual(values = setNames(chosen_colors, lvls)) +
          scale_fill_manual(values = setNames(chosen_colors, lvls), guide = "none")
      }
    }
    return(p)
  }
  
  # 9. Build plot(s) with appropriate title formatting and splitting by extra groups if needed
  if(is.null(extra_groups)) {
    base_title <- if(title_string != "") title_string else {
      if(!is.null(main_group)) paste(map_variable_name(measure_var, say_ibi), "for", map_variable_name(main_group, say_ibi))
      else map_variable_name(measure_var, say_ibi)
    }
  } else {
    base_title <- if(title_string != "") title_string else paste(map_variable_name(measure_var, say_ibi), "for")
  }
  
  # Helper to create extra info for title from extra groups
  make_extra_info <- function(subdat) {
    extra_info <- sapply(extra_groups, function(var) {
      paste(map_variable_name(var, say_ibi), as.character(unique(subdat[[var]]))[1])
    })
    paste(extra_info, collapse = if(length(extra_groups) > 1) ", " else " ")
  }
  
  plot_list <- list()
  if(is.null(extra_groups)) {
    plot_list[[1]] <- create_radial_plot(df_plot, base_title)
    plot_list_names <- if(!is.null(main_group)) map_variable_name(main_group, say_ibi) else gsub("[^A-Za-z0-9]+", "", measure_var)
    names(plot_list) <- plot_list_names
  } else {
    split_data <- split(df_plot, df_plot[, extra_groups, drop = FALSE])
    plot_list <- lapply(split_data, function(subdat) {
      extra_info <- make_extra_info(subdat)
      full_title <- paste(base_title, extra_info)
      create_radial_plot(subdat, full_title)
    })
    plot_list_names <- lapply(split_data, function(subdat) {
      file_parts <- if(!is.null(main_group)) c(map_variable_name(main_group, say_ibi)) else c()
      for(var in extra_groups) {
        file_parts <- c(file_parts, paste(var, as.character(unique(subdat[[var]]))[1], sep = "_"))
      }
      paste(file_parts, collapse = "_")
    })
    names(plot_list) <- unlist(plot_list_names)
  }
  
  # 10. File naming and saving
  file_prefix <- paste0("plot_radial_", measure_var, "_")
  
  if(output_as_one && length(plot_list) > 1) {
    library(cowplot)
    ncol_val <- if(length(plot_list) == 2) 1 else 2
    combined_plot <- cowplot::plot_grid(plotlist = plot_list, ncol = ncol_val)
    combined_name <- if(!is.null(groups)) paste(groups, collapse = "_") else gsub("[^A-Za-z0-9]+", "", measure_var)
    combined_file_name <- file.path(vis_path, paste0(file_prefix, "combined_", combined_name, ".", plot_format))
    if(save_plots) {
      ggsave(combined_file_name, combined_plot, width = plot_width, height = plot_height, dpi = dpi)
    }
    if(show_plots) print(combined_plot)
    return(combined_plot)
  } else {
    if(save_plots) {
      for(i in seq_along(plot_list)) {
        file_suffix <- names(plot_list)[i]
        file_name <- file.path(vis_path, paste0(file_prefix, file_suffix, ".", plot_format))
        ggsave(file_name, plot_list[[i]], width = plot_width, height = plot_height, dpi = dpi)
      }
    }
    if(show_plots) {
      for(p in plot_list) print(p)
    }
    return(plot_list)
  }
}

plot_radial_line_list <- function(df,
                                  angle_var = "phaseFeedb",  # Angle variable in degrees (0-360)
                                  pure_IBI = FALSE,          # Determines which variables to use: IBI vs changeIBI
                                  show_fewer_Rpeaks = TRUE,  # If TRUE, use only a subset of the columns
                                  bin_size = 10,             # Bin size in degrees
                                  aggregator = mean,         # Aggregation function (e.g. mean, median)
                                  normalize = FALSE,         # Normalize aggregated measure to [0,1]?
                                  smooth = FALSE,            # Apply LOESS smoothing?
                                  add_ci = FALSE,            # Option to add CI shading (default FALSE)
                                  say_ibi = FALSE,           # Replace "RR interval" with "IBI" in titles?
                                  groups = NULL,             # Optional grouping variable(s)
                                  title_string = "",         # Custom plot title
                                  legend_title = "",         # Custom legend title
                                  show_title = TRUE,
                                  show_legend_title = TRUE,
                                  show_legend = TRUE,
                                  colors = NULL,
                                  expand_center = 0.25,
                                  move_border = -0.1,
                                  add_grid = FALSE,
                                  scale = NULL,              # if provided, a two-element vector (min, max) for tick labels
                                  num_ticks = NULL) {        # desired number of tick labels
  
  # Determine prefix and select pivot columns based on pure_IBI.
  prefix <- if (pure_IBI) "IBI" else "changeIBI"
  pivot_cols <- names(df)[startsWith(names(df), prefix) & !grepl("_2groups", names(df))]
  
  # Subset columns if show_fewer_Rpeaks is TRUE.
  if (show_fewer_Rpeaks) {
    if (pure_IBI) {
      subset_idx <- 4:7
    } else {
      subset_idx <- 3:6
    }
    pivot_cols <- pivot_cols[subset_idx]
  }
  
  # For each selected measure variable, call plot_radial_line with saving disabled.
  raw_list <- lapply(pivot_cols, function(mvar) {
    plot_radial_line(
      df = df,
      angle_var = angle_var,
      measure_var = mvar,
      bin_size = bin_size,
      aggregator = aggregator,
      normalize = normalize,
      smooth = smooth,
      add_ci = add_ci,
      say_ibi = say_ibi,
      groups = groups,
      title_string = title_string,
      legend_title = legend_title,
      show_title = show_title,
      show_legend_title = show_legend_title,
      show_legend = show_legend,
      vis_path = "",      # Saving disabled
      plot_width = NA,    # Not used
      plot_height = NA,   # Not used
      show_plots = FALSE, # Do not display immediately
      save_plots = FALSE, # Do not save plots to file
      plot_format = "svg",
      dpi = NA,
      colors = colors,
      expand_center = expand_center,
      move_border = move_border,
      add_grid = add_grid,
      output_as_one = FALSE,  # Always output a list from each call
      scale = scale,
      num_ticks = num_ticks
    )
  })
  
  # Recursive function to flatten a nested list and return only ggplot objects.
  flatten_plots <- function(x) {
    out <- list()
    for (item in x) {
      if (inherits(item, "ggplot")) {
        out <- c(out, list(item))
      } else if (is.list(item)) {
        out <- c(out, flatten_plots(item))
      }
    }
    return(out)
  }
  
  flat_plots <- flatten_plots(raw_list)
  names(flat_plots) <- NULL
  return(flat_plots)
}

plot_line_less_smoothed <- function(df,
                      x_var = "phaseFeedb",
                      y_var = "IBI5",
                      bin_size = 10,
                      aggregator = mean,
                      normalize = FALSE,
                      smooth = FALSE,
                      add_ci = FALSE,
                      groups = NULL,
                      clean_x_var = FALSE,
                      title_str = "",
                      legend_title = "",
                      show_title = TRUE,
                      show_legend = TRUE,
                      colors = NULL,
                      num_ticks = NULL) {
  
  ## 1. Validate columns
  req_cols <- c(x_var, y_var)
  missing <- setdiff(req_cols, names(df))
  if (length(missing) > 0) {
    stop("Missing columns: ", paste(missing, collapse = ", "))
  }
  
  ## 2. Handle grouping
  if (!is.null(groups)) {
    for (g in groups) {
      df[[g]] <- as.factor(df[[g]])
      df <- df %>%
        filter(!is.na(.data[[g]]) & .data[[g]] != "") %>%
        droplevels()
    }
    main_group <- groups[1]
  } else {
    main_group <- NULL
  }
  
  if (clean_x_var) {
    x_all <- df[[x_var]]
    mu <- mean(x_all, na.rm = TRUE)
    sig <- sd(x_all, na.rm = TRUE)
    df <- df %>% filter(.data[[x_var]] >= mu - 3*sig,
                        .data[[x_var]] <= mu + 3*sig)
  }
  
  ## 3. Create bins
  x_vals <- df[[x_var]]
  breaks <- seq(min(x_vals, na.rm = TRUE), max(x_vals, na.rm = TRUE), by = bin_size)
  if (tail(breaks,1) < max(x_vals, na.rm=TRUE)) {
    breaks <- c(breaks, max(x_vals, na.rm=TRUE))
  }
  df_b <- df %>%
    mutate(x_bin = cut(.data[[x_var]], breaks = breaks, include.lowest = TRUE, right = FALSE))
  
  ## 4. Summarize within bins
  grp_vars <- c("x_bin", groups)
  if (add_ci) {
    df_sum <- df_b %>%
      group_by(across(all_of(grp_vars))) %>%
      summarise(
        n      = n(),
        mean_y = mean(.data[[y_var]], na.rm = TRUE),
        sd_y   = sd(.data[[y_var]], na.rm = TRUE),
        lower  = mean_y - qt(0.975, df = n - 1) * sd_y / sqrt(n),
        upper  = mean_y + qt(0.975, df = n - 1) * sd_y / sqrt(n),
        .groups = "drop"
      ) %>%
      rename(plot_value = mean_y)
  } else {
    df_sum <- df_b %>%
      group_by(across(all_of(grp_vars))) %>%
      summarise(
        plot_value = aggregator(.data[[y_var]], na.rm = TRUE),
        .groups = "drop"
      )
  }
  
  ## 5. Compute bin centers
  center_fun <- function(lbl) {
    nums <- as.numeric(unlist(strsplit(gsub("\\[|\\]|\\(|\\)", "", lbl), ",")))
    mean(nums)
  }
  df_sum <- df_sum %>%
    mutate(x_center = sapply(as.character(x_bin), center_fun))
  
  ## 6. Normalize if requested
  if (normalize) {
    rng <- range(df_sum$plot_value, na.rm = TRUE)
    if (diff(rng) != 0) {
      df_sum <- df_sum %>%
        mutate(
          plot_value = (plot_value - rng[1]) / diff(rng),
          lower      = if (add_ci) (lower - rng[1]) / diff(rng) else NULL,
          upper      = if (add_ci) (upper - rng[1]) / diff(rng) else NULL
        )
    }
  }
  
  ## 7. Smooth if requested
  if (smooth) {
    if (!is.null(main_group)) {
      df_plot <- df_sum %>%
        group_by(across(all_of(main_group))) %>%
        group_modify(~ {
          dat <- .
          new_x <- seq(min(dat$x_center), max(dat$x_center), length.out = 200)
          mod_plot <- loess(plot_value ~ x_center, data = dat, span = 0.3)
          out <- data.frame(x_center = new_x,
                            plot_value = predict(mod_plot, newdata = data.frame(x_center = new_x)))
          if (add_ci) {
            mod_l <- loess(lower ~ x_center, data = dat, span = 0.3)
            mod_u <- loess(upper ~ x_center, data = dat, span = 0.3)
            out$lower <- predict(mod_l, newdata = data.frame(x_center = new_x))
            out$upper <- predict(mod_u, newdata = data.frame(x_center = new_x))
          }
          out[[main_group]] <- unique(dat[[main_group]])
          out
        }) %>%
        ungroup()
    } else {
      new_x <- seq(min(df_sum$x_center), max(df_sum$x_center), length.out = 200)
      mod_plot <- loess(plot_value ~ x_center, data = df_sum, span = 0.3)
      df_plot <- data.frame(
        x_center   = new_x,
        plot_value = predict(mod_plot, newdata = data.frame(x_center = new_x))
      )
      if (add_ci) {
        mod_l <- loess(lower ~ x_center, data = df_sum, span = 0.3)
        mod_u <- loess(upper ~ x_center, data = df_sum, span = 0.3)
        df_plot$lower <- predict(mod_l, newdata = data.frame(x_center = new_x))
        df_plot$upper <- predict(mod_u, newdata = data.frame(x_center = new_x))
      }
    }
  } else {
    df_plot <- df_sum
  }
  
  ## 8. Build plot
  p <- ggplot(df_plot, aes(x = x_center, y = plot_value))
  
  if (add_ci) {
    if (!is.null(main_group)) {
      p <- p + geom_ribbon(aes_string(fill = main_group, group = main_group,
                                      ymin = "lower", ymax = "upper"), alpha = 0.2)
    } else {
      p <- p + geom_ribbon(aes(ymin = lower, ymax = upper), fill = "gray80", alpha = 0.2)
    }
  }
  
  if (!is.null(main_group)) {
    p <- p +
      geom_line(aes_string(color = main_group, group = main_group), size = 1) +
      geom_point(aes_string(color = main_group), size = 2)
  } else {
    p <- p + geom_line(size = 1) + geom_point(size = 2)
  }
  
  p <- p +
    labs(
      title = if (show_title) title_str else NULL,
      x     = map_variable_name(x_var),
      y     = if (normalize) paste0(map_variable_name(y_var), " (normalized)") else map_variable_name(y_var),
      color = if (show_legend) legend_title else NULL,
      fill  = if (show_legend) legend_title else NULL
    ) +
    theme(
      legend.position  = if (show_legend) "right" else "none",
      plot.title       = element_text(hjust = 0.5)
    )
  
  if (!is.null(main_group) && !is.null(colors)) {
    lvls <- levels(df[[main_group]])
    p <- p + scale_color_manual(values = setNames(colors[seq_along(lvls)], lvls)) +
      scale_fill_manual(values = setNames(colors[seq_along(lvls)], lvls))
  }
  
  print(p)
  invisible(p)
}

plot_line <- function(df,
                      x_var       = "phaseFeedb",    # x-axis variable
                      y_var       = "IBI5",          # y-axis variable
                      bin_size    = 10,                # bin size for x_var
                      aggregator  = mean,              # aggregation function
                      normalize   = FALSE,             # normalize y to [0,1]
                      smooth      = FALSE,             # use geom_smooth
                      add_ci      = FALSE,             # add confidence interval ribbon
                      groups      = NULL,              # grouping variable(s)
                      clean_x_var = FALSE,             # remove ±3 SD outliers on x_var
                      clean_sd = 3,
                      title_str   = "",              # custom title
                      legend_title= "",              # custom legend title
                      show_title  = TRUE,
                      show_legend = TRUE,
                      colors      = NULL) {
  
  # 1. Validate columns
  req_cols <- c(x_var, y_var)
  missing <- setdiff(req_cols, names(df))
  if (length(missing) > 0) stop("Missing columns: ", paste(missing, collapse = ", "))
  
  # 2. Handle grouping
  if (!is.null(groups)) {
    for (g in groups) {
      df[[g]] <- as.factor(df[[g]])
      df <- df %>% filter(!is.na(.data[[g]]) & .data[[g]] != "") %>% droplevels()
    }
    main_group <- groups[1]
  } else {
    main_group <- NULL
  }
  
  # 3. Clean x_var outliers
  if (clean_x_var) {
    x_values <- df[[x_var]]
    mu <- mean(x_values, na.rm = TRUE)
    sigma <- sd(x_values, na.rm = TRUE)
    df <- df %>% filter(.data[[x_var]] >= mu - clean_sd*sigma,
                        .data[[x_var]] <= mu + clean_sd*sigma)
  }
  
  # 4. Create bins and drop NA bins
  x_vals <- df[[x_var]]
  breaks <- seq(min(x_vals, na.rm = TRUE), max(x_vals, na.rm = TRUE), by = bin_size)
  if (tail(breaks, 1) < max(x_vals, na.rm = TRUE)) breaks <- c(breaks, max(x_vals, na.rm = TRUE))
  df_b <- df %>%
    mutate(x_bin = cut(.data[[x_var]], breaks = breaks, include.lowest = TRUE, right = FALSE)) %>%
    filter(!is.na(x_bin))
  
  # 5. Summarise
  grp_vars <- c("x_bin", groups)
  if (add_ci) {
    df_sum <- df_b %>% group_by(across(all_of(grp_vars))) %>%
      summarise(
        n      = n(),
        mean_y = mean(.data[[y_var]], na.rm = TRUE),
        sd_y   = sd(.data[[y_var]], na.rm = TRUE),
        lower  = ifelse(n > 1, mean_y - qt(0.975, df = n - 1) * sd_y / sqrt(n), NA),
        upper  = ifelse(n > 1, mean_y + qt(0.975, df = n - 1) * sd_y / sqrt(n), NA),
        .groups = "drop"
      ) %>% rename(plot_value = mean_y)
  } else {
    df_sum <- df_b %>% group_by(across(all_of(grp_vars))) %>%
      summarise(plot_value = aggregator(.data[[y_var]], na.rm = TRUE), .groups = "drop")
  }
  
  # 6. Ensure grouping factor levels
  if (!is.null(main_group)) {
    df_sum[[main_group]] <- factor(df_sum[[main_group]], levels = levels(df[[main_group]]))
  }
  
  # 7. Compute bin centers
  center_fun <- function(lbl) {
    clean_lbl <- gsub("[^0-9.,-]", "", lbl)
    nums <- as.numeric(strsplit(clean_lbl, ",")[[1]])
    mean(nums, na.rm = TRUE)
  }
  df_sum <- df_sum %>% mutate(x_center = sapply(as.character(x_bin), center_fun))
  
  # 8. Normalize
  if (normalize) {
    rng <- range(df_sum$plot_value, na.rm = TRUE)
    if (diff(rng) != 0) {
      df_sum <- df_sum %>% mutate(
        plot_value = (plot_value - rng[1]) / diff(rng),
        lower      = if (add_ci) (lower - rng[1]) / diff(rng) else NA,
        upper      = if (add_ci) (upper - rng[1]) / diff(rng) else NA
      )
    }
  }
  
  # 9. Build ggplot
  p <- ggplot(df_sum, aes(x = x_center, y = plot_value))
  if (smooth) {
    if (!is.null(main_group)) {
      p <- p + geom_smooth(aes(color = .data[[main_group]], fill = .data[[main_group]], group = .data[[main_group]]),
                           method = "loess", se = add_ci, span = 0.3)
    } else {
      p <- p + geom_smooth(method = "loess", se = add_ci, span = 0.3)
    }
  } else {
    if (add_ci) {
      if (!is.null(main_group)) {
        p <- p + geom_ribbon(aes(ymin = lower, ymax = upper, fill = .data[[main_group]], group = .data[[main_group]]), alpha = 0.2)
      } else {
        p <- p + geom_ribbon(aes(ymin = lower, ymax = upper), fill = "gray80", alpha = 0.2)
      }
    }
    if (!is.null(main_group)) {
      p <- p + geom_line(aes(color = .data[[main_group]], group = .data[[main_group]]), size = 1) +
        geom_point(aes(color = .data[[main_group]]), size = 2)
    } else {
      p <- p + geom_line(size = 1) + geom_point(size = 2)
    }
  }
  
  # 10. Labels & theme
  p <- p + labs(
    x     = map_variable_name(x_var),
    y     = if (normalize) paste0(map_variable_name(y_var), " (normalized)") else map_variable_name(y_var),
    title = if (show_title) title_str else NULL,
    color = if (show_legend) legend_title else NULL,
    fill  = if (show_legend) legend_title else NULL
  ) +
    theme(
      legend.position = if (show_legend) "right" else "none",
      plot.title      = element_text(hjust = 0.5)
    )
  
  # 11. Manual colors
  if (!is.null(main_group) && !is.null(colors)) {
    if (is.null(names(colors))) {
      names(colors) <- levels(df[[main_group]])
    }
    p <- p + scale_color_manual(values = colors) +
      scale_fill_manual(values = colors)
  }
  
  print(p)
  invisible(p)
}

plot_line_list <- function(df,
                           x_var           = 'phaseFeedb',
                           pure_IBI      = FALSE,
                           show_fewer_Rpeaks     = TRUE,
                           bin_size        = 10,
                           aggregator      = mean,
                           normalize       = FALSE,
                           smooth          = FALSE,
                           add_ci          = FALSE,
                           groups          = NULL,
                           title_str       = '',
                           legend_title    = '',
                           show_title      = TRUE,
                           show_legend     = TRUE,
                           colors          = NULL,
                           less_smoothed = FALSE) {
  # select columns
  prefix <- if(pure_IBI) 'IBI' else 'changeIBI'
  vars <- names(df)[startsWith(names(df), prefix) & !grepl("_2groups", names(df))]
  
  # Subset columns if show_fewer_Rpeaks is TRUE.
  if (show_fewer_Rpeaks) {
    if (pure_IBI) {
      subset_idx <- 4:7
    } else {
      subset_idx <- 3:6
    }
    vars <- vars[subset_idx]
  }
  
  if (less_smoothed){
    plots <- lapply(vars, function(yv) {
      plot_line_less_smoothed(df = df,
                x_var = x_var,
                y_var = yv,
                bin_size = bin_size,
                aggregator = aggregator,
                normalize = normalize,
                smooth = smooth,
                add_ci = add_ci,
                groups = groups,
                title_str = title_str,
                legend_title = legend_title,
                show_title = show_title,
                show_legend = show_legend,
                colors = colors)
    })
  } else {
  plots <- lapply(vars, function(yv) {
    plot_line(df = df,
              x_var = x_var,
              y_var = yv,
              bin_size = bin_size,
              aggregator = aggregator,
              normalize = normalize,
              smooth = smooth,
              add_ci = add_ci,
              groups = groups,
              title_str = title_str,
              legend_title = legend_title,
              show_title = show_title,
              show_legend = show_legend,
              colors = colors)
  })
  }
  return(plots)
}

plot_smoothed_ibi <- function(df, groups = NULL, pure_IBI = FALSE, bin_size = 50, 
                              min_ms = NA, max_ms = NA, colors = NULL, 
                              add_line = TRUE, plot_ci = TRUE, show_title = TRUE, 
                              title_string = "", legend_title = "", show_legend_title = TRUE, 
                              show_legend = TRUE, vis_path, plot_width = 8, plot_height = 6, 
                              show_plots = FALSE, save_plots = TRUE, plot_format = "svg", 
                              dpi = 150, output_as_one = FALSE, col_ = NULL, 
                              smooth_further = FALSE, say_ibi = FALSE) {
  
  # Add an id column to keep track of trials
  df$id <- 1:nrow(df)
  
  # Choose columns based on pure_IBI flag.
  if (pure_IBI) {
    rtoF_cols <- paste0("RtoF_", 1:8)
    measure_cols <- paste0("IBI", 1:8)
  } else {
    rtoF_cols <- paste0("RtoF_", 2:8)  # drop RtoF_1
    measure_cols <- paste0("changeIBI", 2:8, "_", 1:7)
  }
  
  # Create long-format data.
  n <- nrow(df)
  k <- length(rtoF_cols)
  data_long <- data.frame(
    id = rep(df$id, times = k),
    rtoF = as.vector(as.matrix(df[, rtoF_cols])),
    measure = as.vector(as.matrix(df[, measure_cols]))
  )
  if (!is.null(groups)) {
    for (g in groups) {
      data_long[[g]] <- rep(df[[g]], times = k)
    }
  }
  
  if (is.na(min_ms)) {
    min_ms <- min(data_long$rtoF, na.rm = TRUE)
  }
  if (is.na(max_ms)) {
    max_ms <- max(data_long$rtoF, na.rm = TRUE)
  }
  data_long <- subset(data_long, rtoF >= min_ms & rtoF <= max_ms)
  
  # Remove NA/empty levels for grouping variables.
  if (!is.null(groups)) {
    for (var in groups) {
      #data_long[[var]] <- as.factor(data_long[[var]])
      data_long <- data_long %>% filter(!is.na(.data[[var]]) & .data[[var]] != "")
      data_long[[var]] <- droplevels(data_long[[var]])
      data_long[[var]] <- factor(
        data_long[[var]],
        levels = levels(df[[var]])     
      )
    }
  }
  
  # Define bin breaks and centers.
  breaks <- seq(min_ms, max_ms, by = bin_size)
  if (max(breaks) < max_ms) {
    breaks <- c(breaks, max_ms)
  }
  bin_centers <- breaks[-length(breaks)] + diff(breaks) / 2
  data_long$bin <- cut(data_long$rtoF, breaks = breaks, labels = bin_centers,
                       include.lowest = TRUE, right = FALSE)
  data_long$bin <- as.numeric(as.character(data_long$bin))
  
  # Compute summary statistics.
  library(dplyr)
  if (!is.null(groups)) {
    summary_stats <- data_long %>%
      group_by(bin, across(all_of(groups))) %>%
      summarise(n = n(),
                Mean = mean(measure, na.rm = TRUE),
                SD = sd(measure, na.rm = TRUE),
                Lower = Mean - qt(0.975, df = n - 1) * SD / sqrt(n),
                Upper = Mean + qt(0.975, df = n - 1) * SD / sqrt(n),
                .groups = "drop")
  } else {
    summary_stats <- data_long %>%
      group_by(bin) %>%
      summarise(n = n(),
                Mean = mean(measure, na.rm = TRUE),
                SD = sd(measure, na.rm = TRUE),
                Lower = Mean - qt(0.975, df = n - 1) * SD / sqrt(n),
                Upper = Mean + qt(0.975, df = n - 1) * SD / sqrt(n),
                .groups = "drop")
  }
  
  # Determine main and extra groups.
  if (!is.null(groups) && length(groups) > 0) {
    main_group <- groups[1]
    extra_groups <- if (length(groups) > 1) groups[-1] else NULL
  } else {
    main_group <- NULL
    extra_groups <- NULL
  }
  
  ### new
  if (!is.null(main_group)) { 
  summary_stats[[main_group]] <- 
    factor(summary_stats[[main_group]],
           levels = levels(df[[main_group]]))
  }
  
  # Build legend title.
  final_legend_title <- if (legend_title != "") legend_title else if (!is.null(main_group)) map_variable_name(main_group) else ""
  if (!show_legend_title) final_legend_title <- NULL
  
  # Build base title.
  if (is.null(extra_groups)) {
    base_title <- if (title_string != "") title_string else {
      if (pure_IBI) paste("RR interval for", map_variable_name(main_group))
      else paste("RR interval change for", map_variable_name(main_group))
    }
  } else {
    base_title <- if (title_string != "") title_string else {
      if (pure_IBI) "RR interval for" else "RR interval change for"
    }
    make_extra_info <- function(subdat) {
      extra_info <- sapply(extra_groups, function(var) {
        paste(map_variable_name(var), as.character(unique(subdat[[var]]))[1])
      })
      paste(extra_info, collapse = ", ")
    }
    base_title <- paste(base_title, make_extra_info(summary_stats))
  }
  if (say_ibi) {
    base_title <- gsub("RR interval", "IBI", base_title, ignore.case = TRUE)
  }
  
  # Set y-axis label.
  default_y_label <- if (pure_IBI) "RR interval (ms)" else "RR interval change (ms)"
  if (say_ibi) {
    default_y_label <- gsub("RR interval", "IBI", default_y_label, ignore.case = TRUE)
  }
  
  # Helper plotting function.
  plot_fun <- function(dat, default_title) {
    dat[[main_group]] <- 
      factor(dat[[main_group]], 
             levels = levels(df[[main_group]]))
    
    # 2) optionally sort by bin so lines draw in order
    dat <- dat[order(dat$bin), ]
    
    # 3) if colors is an unnamed vector, give it the
    #    names of your factor levels (so we map by name)
    if (!is.null(colors) && is.null(names(colors))) {
      names(colors) <- levels(dat[[main_group]])
    }
    
    #p <- ggplot(dat, aes(x = bin, y = Mean))
    p <- ggplot(dat, aes(
      x     = bin,
      y     = Mean,
      color = .data[[main_group]],
      group = .data[[main_group]]
    ))
    if (add_line) {
      p <- p + geom_vline(xintercept = 0, linetype = "dashed", 
                          color = "lightgray", size = 0.5)
    }
    if (smooth_further) {
      if (!is.null(main_group)) {
        # First, plot the smoothed CI with no legend...
        p <- p + geom_smooth(aes(fill = .data[[main_group]]),
                             method = "loess", se = TRUE, span = 0.3, show.legend = FALSE)
        # ...then add the smoothed line with legend mapping.
        p <- p + geom_smooth(aes(color = .data[[main_group]]),
                             method = "loess", se = FALSE, span = 0.3, show.legend = TRUE)
      } else {
        p <- p + geom_smooth(method = "loess", se = TRUE, span = 0.3,
                             color = if (is.null(colors)) "black" else colors[1],
                             fill = "gray")
      }
    } else {
      if (plot_ci) {
        if (!is.null(main_group)) {
          p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = .data[[main_group]]), 
                               alpha = 0.2, color = NA) +
            guides(fill = "none")
        } else {
          p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), 
                               fill = "gray", alpha = 0.2, color = NA)
        }
      }
      if (!is.null(main_group)) {
        p <- p + geom_line(aes(color = .data[[main_group]]), size = 1) +
          geom_point(aes(color = .data[[main_group]]))
      } else {
        default_color <- if (is.null(colors)) "black" else colors[1]
        p <- p + geom_line(size = 1, color = default_color) +
          geom_point(color = default_color)
      }
    }
    
    # Format title: convert to lowercase then replace as required.
    final_title <- tolower(default_title)
    final_title <- gsub("ibi", "IBI", final_title, ignore.case = TRUE)
    final_title <- gsub("rr", "RR", final_title, ignore.case = TRUE)
    final_title <- gsub(" pe", " PE", final_title, ignore.case = TRUE)
    final_title <- paste0(toupper(substr(final_title, 1, 1)), substr(final_title, 2, nchar(final_title)))
    
    p <- p + labs(x = "Time from feedback (ms)",
                  y = default_y_label,
                  title = if (show_title) final_title else NULL,
                  color = final_legend_title)
    
    #if (!is.null(main_group) && !is.null(colors)) {
      #lvls <- unique(dat[[main_group]])
     # lvls <- levels(df[[main_group]])
      #chosen_colors <- colors[seq_along(lvls)]
      #p <- p + #scale_color_manual(values = setNames(chosen_colors, lvls)) +
        #scale_fill_manual(values = setNames(chosen_colors, lvls), guide = "none")
       # scale_color_manual(
      #    breaks = levels(df[[main_group]]),            # e.g. c("High","Low")
      #    values = setNames(colors, levels(df[[main_group]]))
      #  ) +
       # scale_fill_manual(
      #    breaks = levels(df[[main_group]]),
      #    values = setNames(colors, levels(df[[main_group]])),
      #    guide  = "none"
      #  )
    #}
    
    if (!is.null(main_group) && !is.null(colors)) {
      p <- p +
        scale_color_manual(
          breaks = levels(df[[main_group]]),          #  ⇦ force legend key order
          values = setNames(colors, levels(df[[main_group]]))
        ) +
        scale_fill_manual(
          breaks = levels(df[[main_group]]),
          values = setNames(colors, levels(df[[main_group]])),
          guide  = "none"
        )
    }
    
    p <- p + theme_minimal() +
      theme(panel.grid = element_blank(),
            legend.position = if (show_legend) "right" else "none")
    return(p)
  }
  
  # Create plots and assign file names.
  plot_list <- list()
  if (!is.null(extra_groups)) {
    split_data <- split(summary_stats, summary_stats[, extra_groups, drop = FALSE])
    plot_list <- lapply(split_data, function(subdat) {
      plot_fun(subdat, base_title)
    })
    plot_list_names <- lapply(split_data, function(subdat) {
      file_parts <- c()  # Note: do not include main_group in file name when extra groups exist.
      for (var in extra_groups) {
        file_parts <- c(file_parts, paste(var, as.character(unique(subdat[[var]]))[1], sep = "_"))
      }
      paste(file_parts, collapse = "_")
    })
    names(plot_list) <- unlist(plot_list_names)
  } else {
    plot_list[[1]] <- plot_fun(summary_stats, base_title)
    names(plot_list) <- if (!is.null(main_group)) map_variable_name(main_group) else "smoothed_ibi"
  }
  
  # File prefix.
  file_prefix <- if (pure_IBI) "plot_smoothed_pureibi_" else "plot_smoothed_ibi_"
  
  if (output_as_one && length(plot_list) > 1) {
    ncol_val <- if (!is.null(col_)) col_ else if (length(plot_list) == 2) 1 else 2
    combined_plot <- cowplot::plot_grid(plotlist = plot_list, ncol = ncol_val)
    combined_name <- paste(groups, collapse = "_")
    combined_file_name <- file.path(vis_path, paste0(file_prefix, "combined_", combined_name, ".", plot_format))
    if (save_plots) {
      ggsave(combined_file_name, combined_plot, width = plot_width, height = plot_height, dpi = dpi)
    }
    if (show_plots) print(combined_plot)
    return(combined_plot)
  } else {
    if (save_plots) {
      for (i in seq_along(plot_list)) {
        file_suffix <- names(plot_list)[i]
        file_name <- file.path(vis_path, paste0(file_prefix, file_suffix, ".", plot_format))
        ggsave(file_name, plot_list[[i]], width = plot_width, height = plot_height, dpi = dpi)
      }
    }
    if (show_plots) {
      for (p in plot_list) print(p)
    }
    return(plot_list)
  }
}

compute_distances <- function(df) {
  # Compute R-to-feedback distances based on the assumption that:
  # - Feedback occurs between Rpeak 4 and Rpeak 5.
  # - msdistToFeedb is the distance (in ms) from Rpeak 4 to feedback.
  # - IBI2, IBI3, IBI4, IBI5, IBI6, IBI7, IBI8 correspond to the interbeat intervals between successive R-peaks.
  
  # For Rpeak 4: negative because feedback happens after it.
  df$RtoF_4 <- -df$msDistToFeedb
  
  # For Rpeak 5: subtract msdistToFeedb from IBI5 (which spans feedback).
  df$RtoF_5 <- df$IBI5 - df$msDistToFeedb
  
  # For earlier R-peaks: subtract successive IBI intervals.
  df$RtoF_3 <- df$RtoF_4 - df$IBI4
  df$RtoF_2 <- df$RtoF_3 - df$IBI3
  df$RtoF_1 <- df$RtoF_2 - df$IBI2
  
  # For later R-peaks: add successive IBI intervals.
  df$RtoF_6 <- df$RtoF_5 + df$IBI6
  df$RtoF_7 <- df$RtoF_6 + df$IBI7
  df$RtoF_8 <- df$RtoF_7 + df$IBI8
  
  # Return a mini data frame with only the RtoF variables.
  mini_df <- df[, paste0("RtoF_", 1:8)]
  return(mini_df)
}

plot_rpeak_dist <- function(df, color_rpeaks = FALSE, fill_colors = 'gray',
                            show_fewer_Rpeaks = TRUE, group_var = NULL, 
                            bin_frame = TRUE, binwidth = 50, opacity = 1,
                            show_title = TRUE, title_string = "",
                            legend_title = "",
                            show_legend_title = TRUE, show_legend = TRUE,
                            vis_path, plot_width = 8, plot_height = 6,
                            show_plots = FALSE, save_plots = TRUE, plot_format = "svg",
                            dpi = 150, add_line = TRUE) {
  
  # Define the RtoF columns and reshape the data.
  rtoF_vars <- paste0("RtoF_", 1:8)
  long_df <- melt(df, measure.vars = rtoF_vars,
                  variable.name = "Rpeak", value.name = "time_from_feedback")
  
  # Full mapping: RtoF_1 -> -4, RtoF_2 -> -3, RtoF_3 -> -2, RtoF_4 -> -1,
  #              RtoF_5 -> 1,  RtoF_6 -> 2,  RtoF_7 -> 3,  RtoF_8 -> 4.
  rename_mapping_full <- c("RtoF_1" = "-4", "RtoF_2" = "-3", "RtoF_3" = "-2", "RtoF_4" = "-1",
                           "RtoF_5" = "1",  "RtoF_6" = "2",  "RtoF_7" = "3",  "RtoF_8" = "4")
  
  # Optionally keep only a subset of R-peaks.
  if (show_fewer_Rpeaks) {
    long_df <- long_df[long_df$Rpeak %in% c("RtoF_4", "RtoF_5", "RtoF_6", "RtoF_7"), ]
    rename_mapping <- rename_mapping_full[c("RtoF_4", "RtoF_5", "RtoF_6", "RtoF_7")]
  } else {
    rename_mapping <- rename_mapping_full
  }
  
  # Determine the bin border color.
  bin_border <- ifelse(bin_frame, "white", NA)
  
  # Build the plot based on whether a grouping variable is provided.
  if (!is.null(group_var)) {
    if (!group_var %in% names(long_df)) {
      stop("group_var not found in the data frame.")
    }
    if (!is.factor(long_df[[group_var]])) {
      long_df[[group_var]] <- as.factor(long_df[[group_var]])
    }
    # Remove rows with NA in the grouping variable.
    long_df <- long_df[!is.na(long_df[[group_var]]), ]
    
    if (color_rpeaks) {
      # Create a combined grouping variable to separate histogram bins by Rpeak
      # while keeping the fill determined only by group_var.
      long_df$group_comb <- interaction(long_df$Rpeak, long_df[[group_var]])
      p <- ggplot(long_df, aes_string(x = "time_from_feedback", fill = group_var)) 
      if (add_line) {
        p <- p + geom_vline(xintercept = 0, linetype = "dashed", 
                            color = "lightgray", size = 0.5)
      }
      p <- p + geom_histogram(aes(group = group_comb), binwidth = binwidth, color = bin_border,
                       position = "identity", alpha = opacity) +
        scale_fill_manual(values = fill_colors) +
        labs(title = "Distribution of R-peaks Relative to Feedback",
             x = "Time from feedback (ms)",
             y = "Count",
             fill = group_var)
    } else {
      # When not coloring by rpeak, just aggregate by the group.
      p <- ggplot(long_df, aes_string(x = "time_from_feedback", fill = group_var)) 
      if (add_line) {
        p <- p + geom_vline(xintercept = 0, linetype = "dashed", 
                            color = "lightgray", size = 0.5)
      }
      p <- p + geom_histogram(binwidth = binwidth, color = bin_border,
                       position = "identity", alpha = opacity) +
        scale_fill_manual(values = fill_colors) +
        labs(title = "Distribution of R-peaks Relative to Feedback",
             x = "Time from feedback (ms)",
             y = "Count",
             fill = group_var)
    }
    
  } else {
    # When no grouping variable is provided.
    if (color_rpeaks) {
      # Ensure fill_colors is a vector of the appropriate length.
      if (length(fill_colors) == 1) {
        fill_colors <- rep(fill_colors, length(rename_mapping))
      } else if (length(fill_colors) < length(rename_mapping)) {
        fill_colors <- rep(fill_colors, length.out = length(rename_mapping))
      } else {
        fill_colors <- fill_colors[1:length(rename_mapping)]
      }
      
      p <- ggplot(long_df, aes(x = time_from_feedback, fill = Rpeak)) 
        if (add_line) {
          p <- p + geom_vline(xintercept = 0, linetype = "dashed", 
                              color = "lightgray", size = 0.5)
        }
      p <- p + geom_histogram(binwidth = binwidth, color = bin_border,
                       position = "identity", alpha = opacity) +
        scale_fill_manual(values = fill_colors,
                          labels = rename_mapping,
                          name = "R-peak") +
        labs(title = "Distribution of R-peaks Relative to Feedback",
             x = "Time from feedback (ms)",
             y = "Count")
      
    } else {
      p <- ggplot(long_df, aes(x = time_from_feedback)) 
      if (add_line) {
        p <- p + geom_vline(xintercept = 0, linetype = "dashed", 
                            color = "lightgray", size = 0.5)
      }
      p <- p + geom_histogram(binwidth = binwidth, fill = fill_colors[1], color = bin_border,
                       position = "identity", alpha = opacity) +
        labs(title = "Distribution of R-peaks Relative to Feedback",
             x = "Time from feedback (ms)",
             y = "Count")
    }
  }
  
  # Update plot title if needed.
  if (show_title) {
    if (title_string != "") {
      p <- p + labs(title = title_string)
    }
  } else {
    p <- p + labs(title = NULL)
  }
  
  # Determine the legend title to use.
  # If legend_title is provided, use it; otherwise, if group_var is provided,
  # use map_variable_name(group_var), and if not, default to "Rpeak" when coloring rpeaks.
  if (legend_title == "") {
    if (!is.null(group_var)) {
      legend_title_used <- map_variable_name(group_var)
    } else if (color_rpeaks) {
      legend_title_used <- "Rpeak"
    } else {
      legend_title_used <- NULL
    }
  } else {
    legend_title_used <- legend_title
  }
  
  # Update the legend title if a legend is used.
  if (!is.null(legend_title_used) && ( !is.null(group_var) || color_rpeaks )) {
    p <- p + labs(fill = if (show_legend_title) legend_title_used else "")
  }
  
  # Optionally hide the legend.
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }
  
  # Saving the plot.
  # Construct the filename based on group_var and color_rpeaks.
  file_name <- "rpeak_distribution"
  if (color_rpeaks) {
    file_name <- paste0(file_name, "_separate")
  }
  if (!is.null(group_var)) {
    file_name <- paste0(file_name, "_", group_var)
  }
  file_name <- paste0(file_name, ".", plot_format)
  file_path <- file.path(vis_path, file_name)
  
  if (save_plots) {
    ggsave(file_path, p, width = plot_width, height = plot_height, dpi = dpi)
  }
  
  if (show_plots) {
    print(p)
  }
  
  return(p)
}

plot_ibi <- function(df, groups, pure_IBI = FALSE, colors = NULL, add_line = TRUE, 
                     show_fewer_Rpeaks = TRUE, output_as_one = FALSE, 
                     vis_path, plot_width = 8, plot_height = 6, show_plots = FALSE, 
                     save_plots = TRUE, plot_format = "svg", show_title = TRUE, 
                     title_string = "", legend_title = "Condition", 
                     show_legend_title = TRUE, show_legend = TRUE, col_ = NULL, 
                     dpi = 150, say_ibi = FALSE) {
  
  # Check that at least one grouping variable is provided
  if (length(groups) < 1) {
    stop("Please provide at least one grouping variable in the 'groups' argument.")
  }
  
  # Determine which columns to pivot based on pure_IBI flag.
  prefix <- if (pure_IBI) "IBI" else "changeIBI"
  pivot_cols <- names(df)[startsWith(names(df), prefix) & !grepl("_2groups", names(df))]
  
  # Pivot data and always name the new column "ibi"
  data_long <- df %>%
    pivot_longer(
      cols = all_of(pivot_cols),
      names_to = "ibi",
      values_to = "value"
    )
  data_long$ibi <- factor(data_long$ibi, levels = unique(data_long$ibi))
  
  # Create summary statistics grouped by "ibi" and all variables in groups.
  stats_groups <- c("ibi", groups)
  data_stats <- data_long %>%
    group_by(across(all_of(stats_groups))) %>%
    summarise(
      mean = mean(value, na.rm = TRUE),
      lower_ci = mean - qt(0.975, df = n() - 1) * sd(value, na.rm = TRUE) / sqrt(n()),
      upper_ci = mean + qt(0.975, df = n() - 1) * sd(value, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  # Convert each grouping variable to factor and remove NA/empty levels.
  for (var in groups) {
    data_stats[[var]] <- as.factor(data_stats[[var]])
    data_stats <- data_stats %>% filter(!is.na(.data[[var]]) & .data[[var]] != "")
    data_stats[[var]] <- droplevels(data_stats[[var]])
  }
  
  # Define main group and extra groups.
  main_group <- groups[1]
  extra_groups <- if (length(groups) > 1) groups[-1] else NULL
  
  # Full mapping based on pure_IBI.
  if (pure_IBI) {
    full_labels <- c(-4, -3, -2, -1, 1, 2, 3, 4)
    full_vline <- 4.5
  } else {
    full_labels <- c(-3, -2, -1, 1, 2, 3, 4)
    full_vline <- 3.5
  }
  
  # If show_fewer_Rpeaks is TRUE, select a subset.
  if (show_fewer_Rpeaks) {
    if (pure_IBI) {
      subset_idx <- 4:7  # corresponds to -1, 1, 2, 3
    } else {
      subset_idx <- 3:6  # corresponds to -1, 1, 2, 3
    }
    full_levels <- levels(data_stats$ibi)
    subset_levels <- full_levels[subset_idx]
    data_stats <- data_stats %>% filter(ibi %in% subset_levels)
    data_stats$ibi <- factor(data_stats$ibi, levels = subset_levels)
    new_labels <- full_labels[subset_idx]
    vline_pos <- 1.5
  } else {
    new_labels <- full_labels
    vline_pos <- full_vline
  }
  
  # Determine final legend title.
  final_legend_title <- if (legend_title != "") legend_title else map_variable_name(main_group)
  if (!show_legend_title) final_legend_title <- NULL
  
  # Build base plot title.
  # For single plot (no extra groups): if title_string is empty, include main group.
  # For multiple plots: if title_string is empty, use just "RR interval" or "RR interval change" (without main group).
  base_title <- if (is.null(extra_groups)) {
    if (title_string != "") title_string else {
      if (pure_IBI) paste("RR interval for", map_variable_name(main_group))
      else paste("RR interval change for", map_variable_name(main_group))
    }
  } else {
    if (title_string != "") title_string else {
      if (pure_IBI) "RR interval for" else "RR interval change for"
    }
  }
  
  # Build default y-axis label.
  default_y_label <- if (pure_IBI) "RR interval (ms)" else "RR interval change (ms)"
  if (say_ibi) {
    default_y_label <- gsub("RR interval", "IBI", default_y_label, ignore.case = TRUE)
  }
  
  # Function to build extra grouping info from extra groups.
  make_extra_info <- function(subdat) {
    if (is.null(extra_groups)) return("")
    extra_info <- sapply(extra_groups, function(var) {
      paste(map_variable_name(var), as.character(unique(subdat[[var]]))[1])
    })
    paste(extra_info, collapse = if (length(extra_groups) > 1) ", " else " ")
  }
  
  # Helper plotting function.
  plot_fun <- function(dat, default_title) {
    # If extra groups exist, append their info.
    if (!is.null(extra_groups)) {
      extra_info <- make_extra_info(dat)
      default_title <- paste(default_title, extra_info)
    }
    
    # Apply title formatting:
    final_plot_title <- tolower(default_title)
    final_plot_title <- gsub(" pe ", " PE ", final_plot_title, ignore.case = TRUE)
    final_plot_title <- gsub("rr interval", "RR interval", final_plot_title, ignore.case = TRUE)
    final_plot_title <- paste0(toupper(substr(final_plot_title, 1, 1)), substr(final_plot_title, 2, nchar(final_plot_title)))
    if (say_ibi) {
      final_plot_title <- gsub("RR interval", "IBI", final_plot_title, ignore.case = TRUE)
    }
    
    p <- ggplot(dat, aes_string(x = "ibi", y = "mean", color = main_group))
    if (add_line) {
      p <- p + geom_vline(xintercept = vline_pos, linetype = "dashed", 
                          color = "lightgray", linewidth = 0.5)
    }
    p <- p + geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2, size = 1) +
      geom_line(aes_string(group = main_group), size = 1) +
      geom_point() +
      labs(
        title = if (show_title) final_plot_title else NULL,
        x = "R-peaks",
        y = default_y_label,
        color = final_legend_title
      ) +
      scale_x_discrete(labels = setNames(new_labels, levels(dat$ibi))) +
      theme_minimal() +
      theme(panel.grid = element_blank(),
            legend.position = if (show_legend) "right" else "none")
    
    if (!is.null(colors)) {
      lvls <- levels(dat[[main_group]])
      chosen_colors <- colors[seq_along(lvls)]
      p <- p + scale_color_manual(values = setNames(chosen_colors, lvls))
    }
    return(p)
  }
  
  # Container for plots.
  plot_list <- list()
  
  if (is.null(extra_groups)) {
    # Single-plot case.
    p <- plot_fun(data_stats, base_title)
    plot_list[[1]] <- p
    names(plot_list) <- main_group  # For file naming, use only main group name.
  } else {
    # Multiple plots: split data by extra groups.
    split_data <- split(data_stats, data_stats[, extra_groups, drop = FALSE])
    plot_list <- lapply(split_data, function(subdat) {
      plot_fun(subdat, base_title)
    })
    # For file naming, include main group (name only) plus each extra group's name and its level.
    plot_list_names <- lapply(split_data, function(subdat) {
      file_parts <- c(main_group)
      for (var in extra_groups) {
        file_parts <- c(file_parts, paste(var, as.character(unique(subdat[[var]]))[1], sep = "_"))
      }
      paste(file_parts, collapse = "_")
    })
    names(plot_list) <- unlist(plot_list_names)
  }
  
  # File prefix: if pure_IBI is TRUE, use "plot_pureibi_", otherwise "plot_ibi_"
  file_prefix <- if (pure_IBI) "plot_pureibi_" else "plot_ibi_"
  
  # Saving and/or showing plots.
  if (output_as_one && length(plot_list) > 1) {
    ncol_val <- if (!is.null(col_)) col_ else if (length(plot_list) == 2) 1 else 2
    combined_plot <- cowplot::plot_grid(plotlist = plot_list, ncol = ncol_val)
    # For combined file naming, use "combined_" plus all group names (without levels).
    combined_name <- paste(groups, collapse = "_")
    combined_file_name <- file.path(vis_path, paste0(file_prefix, "combined_", combined_name, ".", plot_format))
    if (save_plots) {
      ggsave(combined_file_name, combined_plot, width = plot_width, height = plot_height, dpi = dpi)
    }
    if (show_plots) print(combined_plot)
    return(combined_plot)
  } else {
    if (save_plots) {
      for (i in seq_along(plot_list)) {
        file_suffix <- names(plot_list)[i]
        file_name <- file.path(vis_path, paste0(file_prefix, file_suffix, ".", plot_format))
        ggsave(file_name, plot_list[[i]], width = plot_width, height = plot_height, dpi = dpi)
      }
    }
    if (show_plots) {
      for (p in plot_list) print(p)
    }
    return(plot_list)
  }
}

