library(ggplot2)
library(tidybayes)
library(dplyr)
library(rlang)
library(posterior)  # for as_draws_df()
library(purrr)
library(bayesplot)
library(RColorBrewer)

gather_coef_draws <- function(fit, coef_name) {
  # Remove "b_" prefix if already present to get the base coefficient name
  if (startsWith(coef_name, "b_")) {
    coef_name <- sub("^b_", "", coef_name)
  }
  
  # Get all parameter names from the model draws
  draws_names <- names(as_draws_df(fit))
  
  if (coef_name == "Intercept") {
    post_param <- "b_Intercept"       # Posterior intercept
    prior_param <- "prior_Intercept"    # Use prior_Intercept for intercept
  } else {
    # Posterior: always "b_" + coef_name
    post_param <- paste0("b_", coef_name)
    
    # Prior: check in order: "prior_" + coef_name, then "prior_b_" + coef_name, then general "prior_b"
    if (paste0("prior_", coef_name) %in% draws_names) {
      prior_param <- paste0("prior_", coef_name)
    } else if (paste0("prior_b_", coef_name) %in% draws_names) {
      prior_param <- paste0("prior_b_", coef_name)
    } else {
      prior_param <- "prior_b"
    }
  }
  
  post_draws <- fit %>%
    tidybayes::gather_draws(!!rlang::sym(post_param)) %>%
    dplyr::mutate(coefficient = coef_name,
                  type = "Posterior")
  
  prior_draws <- fit %>%
    tidybayes::gather_draws(!!rlang::sym(prior_param)) %>%
    dplyr::mutate(coefficient = coef_name,
                  type = "Prior")
  
  dplyr::bind_rows(post_draws, prior_draws)
}

# Revised single coefficient plotting function
plot_prior_posterior_rope <- function(fit, 
                                      coef_name, 
                                      rope = NULL, 
                                      x_limits = c(-0.21, 0.21),
                                      rope_color = NULL,
                                      xlab = "Possible parameter values (standardized coef.)",
                                      ylab = "Density",
                                      title = NULL,
                                      subtitle = NULL,
                                      legend_pos = c(1, 1),  # default: upper-right
                                      no_null = FALSE,            
                                      no_title = FALSE,
                                      short_title = FALSE,
                                      BF = NULL,
                                      summary_df = NULL,
                                      extract_rope_from_summary = FALSE,
                                      rope_overlap = FALSE) {
  
  # Determine title if not provided
  if (is.null(title)) {
    title_part1 <- map_variable_name(coef_name)
    title_part1 <- paste0(tolower(substr(title_part1, 1, 1)), substr(title_part1, 2, nchar(title_part1)))
    if (!short_title) {
      title_part2 <- map_variable_name(sub("^scale_", "", colnames(fit$data)[1]))
      title_part2 <- paste0(tolower(substr(title_part2, 1, 1)), substr(title_part2, 2, nchar(title_part2)))
      final_title <- paste("Effect of", title_part1, "on", title_part2)
    } else {
      final_title <- paste("Effect of", title_part1)
    }
  } else {
    final_title <- title
  }
  
  # Attempt to extract summary information if summary_df is provided.
  extraction_successful <- FALSE
  if (!is.null(summary_df)) {
    row_extract <- summary_df[summary_df$Model_parameter == coef_name, , drop = FALSE]
    if (nrow(row_extract) > 0) {
      extraction_successful <- TRUE
    }
  }
  
  # Override rope if extraction was successful and requested.
  if (extraction_successful && extract_rope_from_summary) {
    rope_extracted <- as.numeric(row_extract[, c("ROPE_low", "ROPE_high")])
    if (!all(is.na(rope_extracted))) {
      rope <- rope_extracted
    }
  }
  
  # Determine BF from summary_df if available.
  if (extraction_successful && is.null(BF)) {
    if (!is.null(rope)) {
      BF <- as.numeric(row_extract[,"rope_outside"])
    } else {
      BF <- as.numeric(row_extract[,"bf_against_null"])
    }
  }
  if (!is.null(BF)) BF <- round(BF, 3)
  
  # Build subtitle if not provided.
  if (is.null(subtitle) && extraction_successful) {
    if (!is.null(BF)) {
      if (rope_overlap) {
        rope_ov <- as.numeric(row_extract[,"ROPE_Percentage"])
        rope_ov <- round(rope_ov, 3)
        subtitle <- paste("BF", BF, ", ROPE overlap", rope_ov)
      } else {
        subtitle <- paste("BF", BF)
      }
    }
  }
  
  # Use the helper function to gather draws for the specified coefficient
  draws <- gather_coef_draws(fit, coef_name)
  
  # Plot the draws using ggplot2 with proper fill based on the 'type' column.
  p <- ggplot(draws, aes(x = .value, fill = type)) +
    geom_density(alpha = 0.7, linewidth = 0.4)+
    scale_fill_manual(values = c("Posterior" = "gray40", "Prior" = "gray90"),
                      name = "",
                      labels = c("Posterior", "Prior")) +
    xlim(x_limits) +
    labs(x = xlab, 
         y = ylab, 
         title = if (no_title) NULL else final_title,
         subtitle = subtitle) +
    theme(legend.position = legend_pos, legend.justification = c(1, 1))
  
  # Add ROPE display if desired.
  if (!no_null) {
    if (!is.null(rope)) {
      if (!is.null(rope_color)) {
        p <- p + geom_rect(aes(xmin = rope[1], xmax = rope[2], 
                               ymin = -Inf, ymax = Inf),
                           fill = rope_color, alpha = 0.01, inherit.aes = FALSE) +
          geom_col(data = data.frame(.value = 0), 
                   aes(x = .value, y = 0, fill = "ROPE"), 
                   show.legend = TRUE, alpha = 0) +
          scale_fill_manual(values = c("Posterior" = "gray40", 
                                       "Prior" = "gray90", 
                                       "ROPE" = rope_color), 
                            name = "", 
                            labels = c("Posterior", "Prior", "ROPE")) +
          guides(fill = guide_legend(override.aes = list(alpha = 1)))
      } else {
        p <- p + geom_vline(xintercept = rope, linetype = "dashed")
      }
    } else {
      p <- p + geom_vline(xintercept = 0, linetype = "dashed") 
    }
  }
  
  return(p)
}

plot_mcmc_chains_bayes <- function(fit,
                                   summary_df,
                                   display_intercept = TRUE,
                                   trace_title = NULL,
                                   trace_subtitle = NULL,
                                   palette = "PuOr",
                                   facet_ncol = 1,
                                   facet_nrow = NULL,
                                   legend_pos = "right") {
  # Filter out intercepts if requested.
  if (!display_intercept) {
    summary_df <- summary_df[!grepl("Intercept", summary_df$Model_parameter), ]
  }
  
  # Use the parameters from summary_df$Model_parameter.
  params <- summary_df$Model_parameter
  
  # Convert the fit to an array of draws.
  # Expected dimensions: iterations x chains x variable
  draws_array <- as.array(fit)
  
  # Subset to only those parameters present in summary_df.
  available_params <- intersect(params, dimnames(draws_array)$variable)
  if (length(available_params) == 0) {
    stop("None of the parameters in summary_df$Model_parameter were found in the fit.")
  }
  draws_sub <- draws_array[, , available_params, drop = FALSE]
  
  # Create a mapping from Model_parameter to friendly Parameter names.
  mapping_trace <- setNames(as.character(summary_df$Parameter), summary_df$Model_parameter)
  
  # Update the third-dimension names (the "variable" names) using the mapping.
  new_var_names <- mapping_trace[available_params]
  dimnames(draws_sub)$variable <- new_var_names
  
  # Determine the number of chains.
  n_chains <- dim(draws_sub)[2]
  
  # Process the palette.
  if (length(palette) == 1 && palette %in% rownames(brewer.pal.info)) {
    max_colors <- brewer.pal.info[palette, "maxcolors"]
    if (n_chains > max_colors) {
      stop("The chosen ColorBrewer palette cannot provide enough colors for all chains.")
    }
    pal <- brewer.pal(n = n_chains, name = palette)
  } else {
    if (length(palette) < n_chains) {
      stop("Not enough colors in the palette vector for the number of chains.")
    }
    pal <- palette[1:n_chains]
  }
  
  # Create the trace plot using bayesplot::mcmc_trace.
  p_trace <- bayesplot::mcmc_trace(draws_sub, 
                                   facet_args = list(ncol = facet_ncol, nrow = facet_nrow)) +
    labs(title = trace_title, subtitle = trace_subtitle,
         x = "Iteration", y = "Parameter value") +
    scale_color_manual(values = pal, guide = guide_legend(title = "Chain")) +
    theme(legend.position = legend_pos)
  
  return(p_trace)
}

plot_all_prior_posterior_rope <- function(fit, 
                                          summary_df,
                                          rope = NULL, 
                                          extract_rope_from_summary = FALSE,
                                          display_intercept = TRUE,
                                          no_null = FALSE,
                                          rope_overlap = FALSE,
                                          title = NULL,
                                          display_title = TRUE,
                                          subtitle = NULL,
                                          show_model_info = TRUE,
                                          model_info = NULL,
                                          x_limits = c(-0.21, 0.21),
                                          xlab = "Possible parameter values (standardised coef.)",
                                          ylab = "Density",
                                          legend_pos = c(1, 0.5),
                                          nrow = NULL,
                                          ncol = 1,
                                          flip_coord = FALSE,
                                          plot_trace = FALSE,
                                          palette = "PuOr") {
  
  # 1. Remove intercept rows if requested.
  if (!display_intercept) {
    summary_df <- summary_df[!grepl("Intercept", summary_df$hypothesis_parameter), ]
  }
  
  # 2. Choose BF column: use "rope_outside" if rope provided or extracting rope; otherwise, "bf_against_null".
  bf_col <- if (!is.null(rope) || extract_rope_from_summary) "rope_outside" else "bf_against_null"
  
  # 3. Format BF values: round to integer if >1000, otherwise to 3 decimals.
  bf_value <- as.numeric(summary_df[[bf_col]])
  bf_formatted <- ifelse(bf_value > 1000, round(bf_value, 0), round(bf_value, 3))
  
  # 4. Optionally append ROPE overlap info.
  overlap_str <- ""
  if (rope_overlap) {
    overlap_value <- round(as.numeric(summary_df$ROPE_Percentage), 3)
    overlap_str <- paste0(", RO ", overlap_value)
  }
  
  # 5. Create mapping for the prior-posterior plot.
  new_label <- paste0(as.character(summary_df$Parameter), ", BF ", bf_formatted, overlap_str)
  mapping <- setNames(new_label, summary_df$hypothesis_parameter)
  
  # 6. Extract coefficient names (using hypothesis_parameter for the main plot).
  coef_names <- summary_df$hypothesis_parameter
  
  # 7. Gather draws for each coefficient.
  all_draws <- map_dfr(coef_names, ~gather_coef_draws(fit, .x))
  
  # 8. Remap coefficient names in the draws to the friendly labels.
  all_draws$coefficient <- ifelse(all_draws$coefficient %in% names(mapping),
                                  mapping[all_draws$coefficient],
                                  all_draws$coefficient)
  all_draws$coefficient <- factor(all_draws$coefficient, 
                                  levels = mapping[unique(summary_df$hypothesis_parameter)])
  
  # 9. Handle plot title.
  if (display_title) {
    if (is.null(title)) {
      title <- paste("Dependent variable:", 
                     map_variable_name(sub("^scale_", "", colnames(fit$data)[1])))
    }
  } else {
    title <- NULL
  }
  
  # 10. Handle subtitle.
  sub_text <- NULL
  if (!is.null(subtitle)) {
    sub_text <- subtitle
  } else if (show_model_info && !is.null(model_info)) {
    sub_text <- paste0("Family: ", model_info$family, 
                       ", num. obs: ", model_info$n_obs, 
                       ", num. groups: ", model_info$n_groups)
  }
  
  # 11. Build the base faceted prior-posterior plot.
  p <- ggplot(all_draws, aes(x = .value, fill = type)) +
    geom_density(alpha = 0.7, linewidth = 0.4)+
    scale_fill_manual(values = c("Posterior" = "gray40", "Prior" = "gray90"),
                      name = "",
                      labels = c("Posterior", "Prior")) +
    xlim(x_limits) +
    labs(x = xlab, y = ylab, title = title, subtitle = sub_text) +
    theme(legend.position = legend_pos, legend.justification = c(1, 1)) +
    facet_wrap(~ coefficient, scales = "fixed", ncol = ncol, nrow = nrow)
  
  # 12. Add dashed lines for ROPE or null values (unless no_null is TRUE).
  if (!no_null) {
    if (!is.null(rope)) {
      p <- p + geom_vline(xintercept = rope, linetype = "dashed")
    } else if (extract_rope_from_summary) {
      rope_df <- summary_df
      rope_df$coefficient <- mapping[rope_df$hypothesis_parameter]
      rope_df$coefficient <- factor(rope_df$coefficient, levels = mapping[unique(summary_df$hypothesis_parameter)])
      p <- p +
        geom_vline(data = rope_df, mapping = aes(xintercept = ROPE_low), linetype = "dashed") +
        geom_vline(data = rope_df, mapping = aes(xintercept = ROPE_high), linetype = "dashed")
    } else {
      p <- p + geom_vline(xintercept = 0, linetype = "dashed")
    }
  }
  
  # 13. Optionally flip coordinates.
  if (flip_coord) {
    p <- p + coord_flip()
  }
  
  # 14. If plot_trace is TRUE, build the trace plot using bayesplot and join side-by-side.
  if (plot_trace) {
    # Determine trace plot title/subtitle: if the main plot shows them, pass "" to preserve spacing; otherwise, NULL.
    trace_title <- if (display_title && !is.null(title)) "" else NULL
    trace_subtitle <- if (!is.null(sub_text)) "" else NULL
    p_trace <- plot_mcmc_chains_bayes(fit, summary_df,
                                      display_intercept = display_intercept,
                                      trace_title = trace_title,
                                      trace_subtitle = trace_subtitle,
                                      palette = palette,
                                      facet_ncol = ncol,
                                      legend_pos = legend_pos)
    final_plot <- cowplot::plot_grid(p, p_trace, ncol = 2)
    return(final_plot)
  } else {
    return(p)
  }
}

compare_brms_model_plots <- function(fit_list, summary_list,
                                     rope = NULL,
                                     extract_rope_from_summary = FALSE,
                                     display_intercept = TRUE,
                                     no_null = FALSE,
                                     rope_overlap = FALSE,
                                     title = NULL,       
                                     display_title = TRUE,
                                     subtitle = NULL,    
                                     show_model_info = TRUE,
                                     model_info = NULL,
                                     x_limits = c(-0.21, 0.21),
                                     xlab = "Possible parameter values (standardized coef.)",
                                     ylab = "Density",
                                     legend_pos = c(1, 0.5),
                                     nrow = NULL,
                                     ncol = 1,
                                     flip_coord = FALSE,
                                     always_null_bf = TRUE) {
  ## 1. For each model, filter out intercept rows if requested.
  summary_list <- lapply(summary_list, function(df) {
    if (!display_intercept) {
      df <- df[!grepl("Intercept", df$hypothesis_parameter), ]
    }
    return(df)
  })
  
  ## 2. Compute the union of coefficients across the filtered summaries.
  #all_coefs <- unique(unlist(lapply(summary_list, function(df) df$hypothesis_parameter)))
  all_coefs <- Reduce(function(x, y) union(x, y), lapply(summary_list, function(df) df$hypothesis_parameter)) # NEW
  if (length(all_coefs) == 0) {
    warning("No coefficients to plot after filtering intercepts.")
    return(ggplot() + ggtitle("No coefficients to plot"))
  }
  
  ## 3. For each model, add dummy rows for any missing coefficients.
  add_dummy_coefs <- function(df, all_coefs) {
    missing_coefs <- setdiff(all_coefs, df$hypothesis_parameter)
    if (length(missing_coefs) > 0) {
      dummy_df <- data.frame(matrix(NA, nrow = length(missing_coefs), ncol = ncol(df)))
      colnames(dummy_df) <- colnames(df)
      dummy_df$hypothesis_parameter <- missing_coefs
      dummy_df$Parameter <- missing_coefs
      df <- rbind(df, dummy_df)
    }
    df$hypothesis_parameter <- factor(df$hypothesis_parameter, levels = all_coefs) # NEW
    df <- df[match(all_coefs, df$hypothesis_parameter), ]
    return(df)
  }
  summary_list_aligned <- lapply(summary_list, add_dummy_coefs, all_coefs = all_coefs)
  
  ## 4. Define a wrapper to gather coefficient draws.
  gather_coef_draws_wrapper <- function(fit, coef_name) {
    tryCatch({
      draws <- gather_coef_draws(fit, coef_name)
      draws
    }, error = function(e) {
      data.frame(.value = numeric(0), coefficient = character(0),
                 type = character(0), stringsAsFactors = FALSE)
    })
  }
  
  ## 5. Build a plot for each model.
  plot_list <- lapply(seq_along(fit_list), function(i) {
    fit <- fit_list[[i]]
    summary_df <- summary_list_aligned[[i]]
    
    ## Choose which Bayes factor column to use.
    bf_col <- if (always_null_bf) {
      "bf_against_null"
    } else {
      if (!is.null(rope) || extract_rope_from_summary) "rope_outside" else "bf_against_null"
    }
    bf_value <- as.numeric(summary_df[[bf_col]])
    
    ## If rope_overlap is TRUE, prepare overlap text.
    overlap_str <- ""
    if (rope_overlap) {
      overlap_str <- paste0(", RO ", round(as.numeric(summary_df$ROPE_Percentage), 3))
    }
    
    bf_formatted <- ifelse(!is.na(bf_value) & bf_value > 100,
                           round(bf_value, 0), round(bf_value, 3))
    
    ## Tag dummy rows (where BF is NA).
    is_dummy <- is.na(bf_value)
    
    ## Build facet labels – dummy rows get an empty label.
    label_vector <- ifelse(is_dummy,
                           "",
                           paste0(as.character(summary_df$Parameter),
                                  ", BF ", bf_formatted, overlap_str))
    mapping <- setNames(label_vector, summary_df$hypothesis_parameter)
    
    ## Gather draws for every coefficient (using raw names).
    all_draws <- do.call(rbind, lapply(all_coefs, function(coef) {
      draws <- gather_coef_draws_wrapper(fit, coef)
      if (nrow(draws) > 0) {
        draws$coefficient <- coef
      }
      draws
    }))
    
    ## Force empty facets: for any coefficient missing from draws, add a dummy row.
    missing_draws <- setdiff(all_coefs, unique(all_draws$coefficient))
    if (length(missing_draws) > 0) {
      dummy_rows <- data.frame(.value = NA, coefficient = missing_draws, type = NA, stringsAsFactors = FALSE)
      all_draws <- rbind(all_draws, dummy_rows)
    }
    all_draws$coefficient <- factor(all_draws$coefficient, levels = all_coefs)
    
    ## Build the base ggplot.
    p <- ggplot(all_draws, aes(x = .value, fill = type)) +
      geom_density(alpha = 0.7, linewidth = 0.4)+
      scale_fill_manual(values = c("Posterior" = "gray40", "Prior" = "gray90"),
                        name = "",
                        labels = c("Posterior", "Prior")) +
      xlim(x_limits) +
      labs(x = xlab, y = ylab) +
      theme(legend.position = legend_pos, legend.justification = c(1, 1)) +
      facet_wrap(~ coefficient, scales = "fixed", ncol = ncol, nrow = nrow, drop = FALSE,
                 labeller = as_labeller(mapping))
    
    ## Draw dashed vertical lines only on non-dummy facets.
    non_dummy_df <- summary_df[!is_dummy, ]
    non_dummy_df$coefficient <- non_dummy_df$hypothesis_parameter
    
    if (!no_null) {
      if (!is.null(rope)) {
        p <- p + geom_vline(data = non_dummy_df, aes(xintercept = rope),
                            linetype = "dashed", inherit.aes = FALSE)
      } else if (extract_rope_from_summary) {
        rope_df <- aggregate(cbind(ROPE_low, ROPE_high) ~ hypothesis_parameter,
                             data = non_dummy_df, FUN = function(x) x[1])
        rope_df$coefficient <- rope_df$hypothesis_parameter
        p <- p +
          geom_vline(data = rope_df, mapping = aes(xintercept = ROPE_low),
                     linetype = "dashed", inherit.aes = FALSE) +
          geom_vline(data = rope_df, mapping = aes(xintercept = ROPE_high),
                     linetype = "dashed", inherit.aes = FALSE)
      } else {
        p <- p + geom_vline(data = non_dummy_df, aes(xintercept = 0),
                            linetype = "dashed", inherit.aes = FALSE)
      }
    }
    
    if (flip_coord) {
      p <- p + coord_flip()
    }
    
    ## Add model-specific title/subtitle.
    if (display_title) {
      # If a list of titles is provided, use the i-th element; otherwise, fallback.
      if (!is.null(title)) {
        title_use <- if (is.list(title)) title[[i]] else title
      } else {
        title_use <- paste("Dependent variable:",
                           map_variable_name(sub("^scale_", "", colnames(fit$data)[1])))
      }
      p <- p + labs(title = title_use)
    }
    if (!is.null(subtitle)) {
      # If a list of subtitles is provided, use the i-th element; otherwise, fallback.
      subtitle_use <- if (is.list(subtitle)) subtitle[[i]] else subtitle
      p <- p + labs(subtitle = subtitle_use)
    } else if (show_model_info && !is.null(model_info)) {
      ## Use model_info[[i]] if model_info is a list.
      current_model_info <- if (is.list(model_info)) model_info[[i]] else model_info
      model_info_sub <- paste0("Family: ", current_model_info$family,
                               ", num. obs: ", as.character(current_model_info$n_obs),
                               ", num. groups: ", as.character(current_model_info$n_groups))
      p <- p + labs(subtitle = model_info_sub)
    }
    
    p
  })
  
  ## 6. Combine the individual model plots side-by-side.
  combined_plot <- cowplot::plot_grid(plotlist = plot_list, ncol = length(fit_list))
  return(combined_plot)
}

plot_coef_across_models <- function(coef_name, 
                                    fits, 
                                    summaries,
                                    rope = NULL, 
                                    extract_rope_from_summary = FALSE,
                                    labels = NULL,
                                    title = NULL,
                                    display_title = TRUE,
                                    x_limits = c(-0.21, 0.21),
                                    xlab = "Possible parameter values (standardized coef.)",
                                    ylab = "Density",
                                    legend_pos = c(1, 0.5),
                                    flip_coord = FALSE,
                                    rope_overlap = FALSE,
                                    no_null = FALSE) {
  # Check that the lists of fits and summaries have the same length.
  if (length(fits) != length(summaries)) {
    stop("The number of model fits and summaries must be equal.")
  }
  
  n_models <- length(fits)
  # If no labels provided, default to "M1", "M2", etc.
  if (is.null(labels)) {
    labels <- paste0("M", seq_len(n_models))
  }
  
  # Determine which BF column to use based on rope settings.
  bf_col <- if (!is.null(rope) || extract_rope_from_summary) "rope_outside" else "bf_against_null"
  
  # Initialize lists to store draws and summary info.
  all_draws_list <- vector("list", n_models)
  rope_info_list <- vector("list", n_models)
  
  for (i in seq_along(fits)) {
    fit <- fits[[i]]
    summary_df <- summaries[[i]]
    
    # Extract draws for the specified coefficient.
    # (Assumes the helper function `gather_coef_draws()` is available.)
    draw_df <- gather_coef_draws(fit, coef_name)
    
    # Extract the row in the summary that corresponds to the coefficient.
    # Matching is done using the "Model_parameter" column.
    coef_row <- summary_df[summary_df$Model_parameter == coef_name, ]
    
    if (nrow(coef_row) == 0) {
      warning(paste("Coefficient", coef_name, "not found in summary for model", labels[i]))
      computed_label <- labels[i]
    } else {
      bf_value <- as.numeric(coef_row[[bf_col]])
      bf_formatted <- ifelse(bf_value > 1000, round(bf_value, 0), round(bf_value, 3))
      
      if (rope_overlap && !is.na(coef_row$ROPE_Percentage)) {
        overlap_value <- round(as.numeric(coef_row$ROPE_Percentage), 3)
        computed_label <- paste0(labels[i], ", BF ", bf_formatted, ", ROPE ", overlap_value)
      } else {
        computed_label <- paste0(labels[i], ", BF ", bf_formatted)
      }
    }
    
    # Assign the computed label to the draw data.
    draw_df$model <- computed_label
    all_draws_list[[i]] <- draw_df
    
    # Save the summary info along with the computed label for vertical lines.
    if (nrow(coef_row) > 0) {
      coef_row$label <- computed_label
    } else {
      coef_row <- data.frame(label = computed_label, ROPE_low = NA, ROPE_high = NA, 
                             ROPE_Percentage = NA, stringsAsFactors = FALSE)
    }
    rope_info_list[[i]] <- coef_row
  }
  
  # Combine all draws and rope information into data frames.
  all_draws <- do.call(rbind, all_draws_list)
  computed_labels <- sapply(rope_info_list, function(x) x$label)
  all_draws$model <- factor(all_draws$model, levels = computed_labels)
  rope_df <- do.call(rbind, rope_info_list)
  rope_df$model <- factor(rope_df$label, levels = computed_labels)
  
  # Handle plot title.
  if (display_title && is.null(title)) {
    title <- map_variable_name(sub("^scale_", "", coef_name))
  }
  
  # Build the base faceted plot.
  p <- ggplot(all_draws, aes(x = .value, fill = type)) +
    geom_density(alpha = 0.7, linewidth = 0.4)+
    scale_fill_manual(values = c("Posterior" = "gray40", "Prior" = "gray90"),
                      name = "",
                      labels = c("Posterior", "Prior")) +
    xlim(x_limits) +
    labs(x = xlab, y = ylab, title = title) +
    theme(legend.position = legend_pos, legend.justification = c(1, 1))
  
  # Adjust facet layout based on flip_coord.
  if (flip_coord) {
    p <- p + facet_wrap(~ model, scales = "fixed", nrow = 1)
  } else {
    p <- p + facet_wrap(~ model, scales = "fixed", ncol = 1)
  }
  
  # Add dashed vertical lines.
  if (!no_null) {
    if (!is.null(rope)) {
      p <- p + geom_vline(xintercept = rope, linetype = "dashed")
    } else if (extract_rope_from_summary) {
      p <- p +
        geom_vline(data = rope_df, mapping = aes(xintercept = ROPE_low), linetype = "dashed") +
        geom_vline(data = rope_df, mapping = aes(xintercept = ROPE_high), linetype = "dashed")
    } else {
      p <- p + geom_vline(xintercept = 0, linetype = "dashed")
    }
  }
  
  # Option to flip coordinates.
  if (flip_coord) {
    p <- p + coord_flip()
  }
  
  return(p)
}

plot_brms_predictions <- function(fit, df, vcolumns,
                                  split_cont_groups = 3,
                                  colors = NULL,
                                  palette = "Dark2",
                                  fill = NULL,
                                  xlabs = NULL,
                                  ylabs = NULL,
                                  title = "Subject-level average values and expectation of posterior predictive",
                                  show_title = TRUE,
                                  ndraws = 100,
                                  convert_units = FALSE) {

  # If vcolumns is not a list, convert it into a list with one element.
  if (!is.list(vcolumns)) {
    vcolumns_list <- list(vcolumns)
  } else {
    vcolumns_list <- vcolumns
  }
  
  # Validate each element: each must be a character vector of length 1 or 2 and columns must exist in df.
  for (i in seq_along(vcolumns_list)) {
    if (!(length(vcolumns_list[[i]]) %in% c(1, 2))) {
      stop("Each element of vcolumns must be a character vector of length 1 or 2.")
    }
    missing_cols <- setdiff(vcolumns_list[[i]], names(df))
    if (length(missing_cols) > 0) {
      stop(paste("For plot", i, "the following columns are missing in df:", 
                 paste(missing_cols, collapse = ", ")))
    }
  }
  
  # Determine required columns: columns used in the model, all vcolumns (from all plots), plus "subID"
  model_cols <- names(fit$data)
  plot_cols  <- unique(unlist(vcolumns_list))
  required_cols <- unique(c(model_cols, plot_cols, "subID"))
  
  # If unit conversion is requested, compute conversion parameters from the originally provided df.
  # Assumes the dependent variable in the model is scaled and its unscaled version has a name
  # obtained by removing "scale_" from the response variable name.
  response_var <- colnames(fit$data)[1]
  if (convert_units) {
    base_name <- sub("^scale_", "", response_var)
    if (!(base_name %in% names(df))) {
      stop(paste("For unit conversion, column", base_name, "must be present in the original df."))
    }
    sd_val <- sd(df[[base_name]], na.rm = TRUE)
    mean_val <- mean(df[[base_name]], na.rm = TRUE)
  }
  
  # Subset the data frame to only required columns
  df_subset <- df[, required_cols, drop = FALSE]
  
  # Helper function: get bin labels for numeric variables
  get_labels <- function(n) {
    if(n == 2) {
      c("Low", "High")
    } else if(n == 3) {
      c("Low", "Medium", "High")
    } else {
      c(1:n) #paste("Group", 1:n)
    }
  }
  
  # Convert each plotting variable to a factor if not already
  for (v in plot_cols) {
    if (!is.factor(df_subset[[v]])) {
      if (is.numeric(df_subset[[v]])) {
        df_subset[[v]] <- cut(df_subset[[v]], 
                              breaks = split_cont_groups, 
                              labels = get_labels(split_cont_groups), 
                              include.lowest = TRUE)
      } else {
        df_subset[[v]] <- as.factor(df_subset[[v]])
      }
    }
  }
  
  # Get posterior predictive draws using the subsetted data and sample ndraws
  df_draws <- df_subset %>%
    add_epred_draws(fit) %>%
    sample_draws(ndraws)
  
  # Apply unit conversion to posterior draws if requested
  if (convert_units) {
    df_draws <- df_draws %>% dplyr::mutate(.epred = .epred * sd_val + mean_val)
  }

  if (is.null(ylabs)) {
    ylabs <- map_variable_name(sub("^scale_", "", response_var))
  }
  
  # Initialize list to store plots
  plot_list <- list()
  
  # Loop over each element in vcolumns_list and create a plot.
  for (i in seq_along(vcolumns_list)) {
    current_vcols <- vcolumns_list[[i]]
    primary_var <- current_vcols[1]
    
    # Set default labels for fill and x-axis based on primary variable if not provided.
    cur_fill <- if (is.null(fill)) map_variable_name(primary_var) else fill
    cur_xlab <- if (is.null(xlabs)) map_variable_name(primary_var) else xlabs
    
    # Aggregate the observed response at the subject level (grouping by subID plus the plotting variable(s))
    agg <- df_subset %>%
      group_by(subID, across(all_of(current_vcols))) %>%
      summarise(.value = mean(.data[[response_var]], na.rm = TRUE), .groups = "drop")
    
    # Apply unit conversion to aggregated values if requested
    if (convert_units) {
      agg <- agg %>% dplyr::mutate(.value = .value * sd_val + mean_val)
    }
    
    # Determine the factor levels and number of levels for the primary variable
    primary_levels <- levels(df_subset[[primary_var]])
    n_levels <- length(primary_levels)
    
    # Determine fill colors: use provided colors if enough, otherwise use the Brewer palette.
    if (!is.null(colors)) {
      if (length(colors) < n_levels) {
        stop(paste("Not enough colors provided for plot", i, 
                   ". Needed:", n_levels, "Provided:", length(colors)))
      }
      fill_colors <- colors[1:n_levels]
    } else {
      if (n_levels > RColorBrewer::brewer.pal.info[palette, "maxcolors"]) {
        stop(paste("The palette", palette, "cannot provide", n_levels, "colors for plot", i, 
                   ". Please supply a custom colors vector."))
      }
      fill_colors <- RColorBrewer::brewer.pal(n_levels, palette)
    }
    
    # Build the ggplot
    p <- ggplot(df_draws, aes(x = .data[[primary_var]], y = .epred, fill = .data[[primary_var]])) +
      stat_halfeye() +
      geom_point(data = agg,
                 mapping = aes(x = .data[[primary_var]], y = .value, fill = .data[[primary_var]]),
                 pch = 21, alpha = 0.5,
                 position = position_nudge(x = -0.1)) +
      scale_fill_manual(values = fill_colors) +
      labs(fill = cur_fill, x = cur_xlab, y = ylabs)
    
    # If two grouping variables are provided and the second is not NA, facet by the second.
    if (length(current_vcols) == 2 && !is.na(current_vcols[2])) {
      facet_formula <- as.formula(paste("~", current_vcols[2]))
      p <- p + facet_grid(facet_formula)
    }
    
    if (show_title) {
      p <- p + ggtitle(title)
    }
    
    plot_list[[i]] <- p
  }
  
  return(plot_list)
}

potential_description_for_plot_brms_predictions <- "Posterior Predictive Checks for [Outcome Variable]. The figure displays half-eye plots representing the posterior predictive distributions of [outcome variable] across bins of [predictor variable] (e.g., absolute prediction error levels) derived from our Bayesian model. Overlaid dots indicate the subject-level mean observed values, aggregated by subID. The close alignment between the model predictions (with their 95% credible intervals) and the empirical means suggests that the model adequately captures the data’s variability across conditions."

custom_conditional_effects <- function(fit, colors, xlab = NULL, ylab = NULL, legend_lab = NULL) {
  
  # Define a mapping function if not already available
  if (!exists("map_variable_name")) {
    map_variable_name <- function(x) { x }
  }
  
  # Get response variable (assumed to be the first column of fit$data)
  response_var <- colnames(fit$data)[1]
  # Default y-label: remove the "scale_" prefix and map the name
  default_ylab <- map_variable_name(sub("^scale_", "", response_var))
  
  # Get conditional effects object (without plotting)
  ce <- conditional_effects(fit, plot = FALSE)
  # Get a list of ggplot objects from the conditional effects.
  plots_list <- plot(ce, plot = FALSE)
  
  # Loop over each plot and modify aesthetics
  for (i in seq_along(plots_list)) {
    p <- plots_list[[i]]
    # Use the list name as the effect name; if not available, use an index.
    effect_name <- names(plots_list)[i]
    if (is.null(effect_name) || effect_name == "") {
      effect_name <- paste("Effect", i)
    }
    
    # Set default x-label: if not provided, use the effect name via map_variable_name()
    current_xlab <- if (is.null(xlab)) map_variable_name(effect_name) else xlab
    # Default legend title: if not provided, also use the effect name mapping
    current_legend_lab <- if (is.null(legend_lab)) map_variable_name(effect_name) else legend_lab
    # Default y-label: if not provided, use the unscaled response name
    current_ylab <- if (is.null(ylab)) default_ylab else ylab
    
    p <- p + labs(x = current_xlab, y = current_ylab, fill = current_legend_lab, color = current_legend_lab)
    
    # Update fill scale if the plot uses fill
    if (!is.null(p$mapping$fill)) {
      fill_var <- as_label(p$mapping$fill)
      n_groups <- length(unique(p$data[[fill_var]]))
      if (n_groups > length(colors)) {
        stop(paste("Not enough colors provided for fill aesthetic in effect", effect_name))
      }
      p <- p + scale_fill_manual(values = colors[1:n_groups])
    }
    
    # Update color scale if the plot uses color (or colour)
    if (!is.null(p$mapping$colour)) {
      color_var <- as_label(p$mapping$colour)
      n_groups <- length(unique(p$data[[color_var]]))
      if (n_groups > length(colors)) {
        stop(paste("Not enough colors provided for color aesthetic in effect", effect_name))
      }
      p <- p + scale_color_manual(values = colors[1:n_groups])
    } else if (!is.null(p$mapping$color)) {
      color_var <- as_label(p$mapping$color)
      n_groups <- length(unique(p$data[[color_var]]))
      if (n_groups > length(colors)) {
        stop(paste("Not enough colors provided for color aesthetic in effect", effect_name))
      }
      p <- p + scale_color_manual(values = colors[1:n_groups])
    }
    
    plots_list[[i]] <- p
  }
  
  # Ensure the list is named according to the effect names
  if (is.null(names(plots_list))) {
    names(plots_list) <- paste("Effect", seq_along(plots_list))
  }
  
  return(plots_list)
}

plot_empirical_binned <- function(df, varnames, bins = 400, color = 'darkgrey', ylim_vals = NULL) {
  # Check that varnames is a character vector of length 2:
  if (length(varnames) != 2) {
    stop("varnames must be a character vector of length 2: first for dependent (y), second for predictor (x).")
  }
  
  # Define dependent and predictor variable names
  dependent <- varnames[1]
  predictor <- varnames[2]
  
  # Use map_variable_name if available; otherwise, use identity function
  if (!exists("map_variable_name")) {
    map_variable_name <- function(x) x
  }
  
  # Get the range of the predictor variable
  predictor_range <- range(df[[predictor]], na.rm = TRUE)
  # Compute the bin width
  bin_width <- diff(predictor_range) / bins
  # Create labels for each bin: a sequence from the minimum to (max - bin_width)
  labels_seq <- seq(predictor_range[1], predictor_range[2] - bin_width, length.out = bins)
  
  # Name for the new binned variable (appended to predictor name)
  bin_var <- paste0(predictor, "_bins")
  
  # Create the binned variable using cut(); labels_seq will be a numeric sequence,
  # so the resulting factor levels can be converted back to numeric later.
  df[[bin_var]] <- cut(df[[predictor]], breaks = bins, labels = labels_seq)
  
  # Group the data by the binned predictor and calculate the mean of the dependent variable
  agg <- df %>%
    group_by(bin = .data[[bin_var]]) %>%
    summarise(mean_dep = mean(.data[[dependent]], na.rm = TRUE)) %>%
    ungroup()
  
  # Convert the binned factor to numeric so it can be plotted on a continuous x-axis
  agg$bin_numeric <- as.numeric(as.character(agg$bin))
  
  # Create the ggplot: points for the binned averages and a linear smooth with the given color.
  p <- ggplot(agg, aes(x = bin_numeric, y = mean_dep)) +
    geom_point(colour=colors[4],alpha=0.8, size=0.5) +
    geom_smooth(method = "lm", se = TRUE, color = color, fill = color) +
    labs(x = map_variable_name(predictor),
         y = map_variable_name(dependent))
  
  # Optionally set y-axis limits if provided
  if (!is.null(ylim_vals) && length(ylim_vals) == 2) {
    p <- p + ylim(ylim_vals[1], ylim_vals[2])
  }
  
  return(p)
}

plot_empirical_binned_with_conditional <- function(fit, df, varnames, bins = 400, color, 
                                                   ylim_vals = NULL, descale = TRUE,
                                                   xlab = NULL, ylab = NULL, plot_title = NULL,
                                                   dot_color = NULL, dot_alpha = 1) {
  # varnames: a character vector of length 2:
  #   first element: dependent variable, 
  #   second element: predictor variable.
  if (length(varnames) != 2) {
    stop("varnames must be a character vector of length 2: first for dependent (y) and second for predictor (x).")
  }
  
  # Save the user-supplied names
  user_dep <- varnames[1]
  user_pred <- varnames[2]
  
  # For the model (fit), the variables should be scaled.
  # If they don't start with "scale_", add it.
  fit_dep <- if (startsWith(user_dep, "scale_")) user_dep else paste0("scale_", user_dep)
  fit_pred <- if (startsWith(user_pred, "scale_")) user_pred else paste0("scale_", user_pred)
  
  # Check that the model contains these variables.
  if (!(fit_dep %in% names(fit$data))) {
    stop(paste("Dependent variable", fit_dep, "not found in fit$data."))
  }
  if (!(fit_pred %in% names(fit$data))) {
    stop(paste("Predictor variable", fit_pred, "not found in fit$data."))
  }
  
  # Helper: remove "scale_" from a name if present.
  get_original <- function(x) {
    if (startsWith(x, "scale_")) sub("^scale_", "", x) else x
  }
  
  # Determine the names to use for aggregation:
  # If descale is TRUE, we use the original (unscaled) names.
  if (descale) {
    agg_dep_var <- get_original(fit_dep)
    agg_pred_var <- get_original(fit_pred)
  } else {
    agg_dep_var <- fit_dep
    agg_pred_var <- fit_pred
  }
  
  # --- Axis Labels (with units) ---
  # Use map_variable_name on the original (descaled) name for the label.
  # For units, call map_variable_units on either the original name (if descale) or the scaled name.
  if (!exists("map_variable_name")) {
    map_variable_name <- function(x) x
  }
  if (!exists("map_variable_units")) {
    map_variable_units <- function(x) ""
  }
  
  # For the label text, always use the descaled name for the label.
  auto_xlab <- map_variable_name(get_original(fit_pred))
  auto_ylab <- map_variable_name(get_original(fit_dep))
  
  # Get the units string based on whether we're descaling.
  units_x <- if (descale) map_variable_units(get_original(fit_pred)) else map_variable_units(fit_pred)
  units_y <- if (descale) map_variable_units(get_original(fit_dep)) else map_variable_units(fit_dep)
  
  # If units string is non-empty, add it in parentheses.
  if (units_x != "") {
    auto_xlab <- paste0(auto_xlab, " (", units_x, ")")
  }
  if (units_y != "") {
    auto_ylab <- paste0(auto_ylab, " (", units_y, ")")
  }
  
  # Use provided axis labels if they are not NULL.
  final_xlab <- if (!is.null(xlab)) xlab else auto_xlab
  final_ylab <- if (!is.null(ylab)) ylab else auto_ylab
  
  # Use plot_title if provided; otherwise, leave it NULL (can be added later to the plot).
  
  # If descale is TRUE, compute conversion parameters from df using the original columns.
  if (descale) {
    if (!(agg_dep_var %in% names(df))) {
      stop(paste("Original dependent variable", agg_dep_var, "not found in df."))
    }
    if (!(agg_pred_var %in% names(df))) {
      stop(paste("Original predictor variable", agg_pred_var, "not found in df."))
    }
    mean_dep <- mean(df[[agg_dep_var]], na.rm = TRUE)
    sd_dep   <- sd(df[[agg_dep_var]], na.rm = TRUE)
    mean_pred <- mean(df[[agg_pred_var]], na.rm = TRUE)
    sd_pred   <- sd(df[[agg_pred_var]], na.rm = TRUE)
  }
  
  # ---- Binning the Predictor (for empirical averages) ----
  predictor_range <- range(df[[agg_pred_var]], na.rm = TRUE)
  bin_width <- diff(predictor_range) / bins
  labels_seq <- seq(predictor_range[1], predictor_range[2] - bin_width, length.out = bins)
  
  # Create a new column in df for the binned predictor.
  bin_var <- paste0(agg_pred_var, "_bins")
  df[[bin_var]] <- cut(df[[agg_pred_var]], breaks = bins, labels = labels_seq)
  
  # ---- Aggregate Dependent Variable by Binned Predictor ----
  library(dplyr)
  agg <- df %>%
    group_by(bin = .data[[bin_var]]) %>%
    summarise(mean_dep = mean(.data[[agg_dep_var]], na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(bin_numeric = as.numeric(as.character(bin)))
  
  # ---- Get Conditional Effects for the Predictor ----
  # Request the effect for the predictor variable used in the model (fit_pred).
  ce <- conditional_effects(fit, effects = fit_pred, plot = FALSE)
  ce_data <- ce[[1]]
  # ce_data should include a column with the predictor (named fit_pred) plus estimate__, lower__, and upper__.
  if (descale) {
    # Convert predictor values from scaled to original units.
    ce_data[[fit_pred]] <- ce_data[[fit_pred]] * sd_pred + mean_pred
    # Convert outcome predictions from scaled to original units.
    ce_data$estimate__ <- ce_data$estimate__ * sd_dep + mean_dep
    ce_data$lower__ <- ce_data$lower__ * sd_dep + mean_dep
    ce_data$upper__ <- ce_data$upper__ * sd_dep + mean_dep
  }
  
  # ---- Plotting ----
  library(ggplot2)
  # If dot_color is not provided, default to the input 'color'
  final_dot_color <- if (is.null(dot_color)) color else dot_color
  
  p <- ggplot() +
    # Empirical binned averages as points.
    geom_point(data = agg, aes(x = bin_numeric, y = mean_dep), 
               color = final_dot_color, alpha = dot_alpha) +
    # Conditional effect line. Use tidy evaluation for the x mapping.
    geom_line(data = ce_data, aes(x = !!sym(fit_pred), y = estimate__), 
              color = color, linewidth = 1) +
    # Uncertainty ribbon.
    geom_ribbon(data = ce_data, aes(x = !!sym(fit_pred), ymin = lower__, ymax = upper__), 
                fill = color, alpha = 0.3) +
    labs(x = final_xlab, y = final_ylab)
  
  # If a plot title is provided, add it.
  if (!is.null(plot_title)) {
    p <- p + ggtitle(plot_title)
  }
  
  if (!is.null(ylim_vals) && length(ylim_vals) == 2) {
    p <- p + ylim(ylim_vals[1], ylim_vals[2])
  }
  
  return(p)
}

plot_coef_orig_across_models <- function(coef_name,
                                         fits,
                                         summaries,
                                         rope = NULL,
                                         extract_rope_from_summary = FALSE,
                                         labels = NULL,
                                         title = NULL,
                                         display_title = TRUE,
                                         xlab = NULL,
                                         ylab = NULL,
                                         flip_coord = FALSE,
                                         rope_overlap = FALSE) {
  # Check that fits and summaries have the same length.
  if (length(fits) != length(summaries)) {
    stop("The number of model fits and summaries must be equal.")
  }
  
  n_models <- length(fits)
  
  # If no labels provided, default to "M1", "M2", etc.
  if (is.null(labels)) {
    labels <- paste0("M", seq_len(n_models))
  }
  
  # Determine which BF column to use based on rope settings.
  bf_col <- if (!is.null(rope) || extract_rope_from_summary) "rope_outside" else "bf_against_null"
  
  # Initialize a data frame to store coefficient info.
  df_coef <- data.frame(model = character(),
                        mode_orig = numeric(),
                        CI_low_orig = numeric(),
                        CI_high_orig = numeric(),
                        stringsAsFactors = FALSE)
  
  # We'll also collect the original_units from the summaries.
  units_list <- c()
  
  for (i in seq_along(fits)) {
    fit <- fits[[i]]
    summary_df <- summaries[[i]]
    
    # Extract the row corresponding to the coefficient.
    coef_row <- summary_df[summary_df$Model_parameter == coef_name, ]
    
    if (nrow(coef_row) == 0) {
      warning(paste("Coefficient", coef_name, "not found in summary for model", labels[i]))
      computed_label <- labels[i]
      mode_orig <- NA
      CI_low_orig <- NA
      CI_high_orig <- NA
    } else {
      bf_value <- as.numeric(coef_row[[bf_col]])
      bf_formatted <- ifelse(bf_value > 1000, round(bf_value, 0), round(bf_value, 3))
      
      if (rope_overlap && !is.na(coef_row$ROPE_Percentage)) {
        overlap_value <- round(as.numeric(coef_row$ROPE_Percentage), 3)
        computed_label <- paste0(labels[i], ", BF ", bf_formatted, ", ROPE ", overlap_value)
      } else {
        computed_label <- paste0(labels[i], ", BF ", bf_formatted)
      }
      
      mode_orig   <- coef_row$mode_orig
      CI_low_orig <- coef_row$CI_low_orig
      CI_high_orig<- coef_row$CI_high_orig
      
      # Save the original_units string from this summary row.
      units_list <- c(units_list, as.character(coef_row$original_units[1]))
    }
    
    df_coef <- rbind(df_coef,
                     data.frame(model = computed_label,
                                mode_orig = mode_orig,
                                CI_low_orig = CI_low_orig,
                                CI_high_orig = CI_high_orig,
                                stringsAsFactors = FALSE))
  }
  
  # Check that all original_units are the same.
  unique_units <- unique(units_list)
  if (length(unique_units) > 1) {
    stop("Original units differ across models. Cannot plot coefficients together in converted units.")
  }
  
  # Set the model column as a factor in the computed order.
  df_coef$model <- factor(df_coef$model, levels = df_coef$model)
  
  # Automatic title: if display_title is TRUE and title is NULL,
  # use map_variable_name on the coefficient name after removing any "scale_" prefix.
  if (display_title && is.null(title)) {
    if (!exists("map_variable_name")) {
      map_variable_name <- function(x) x
    }
    title <- map_variable_name(sub("^scale_", "", coef_name))
  }
  
  # For xlab: if xlab is NULL, we leave it blank.
  if (is.null(xlab)) {
    xlab <- ""
  }
  
  # For ylab, if not provided, use the original_units value.
  if (is.null(ylab)) {
    ylab <- unique_units[1]
  }
  
  # Build the basic plot.
  p <- ggplot(df_coef, aes(x = model, y = mode_orig)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = CI_low_orig, ymax = CI_high_orig), width = 0.2) +
    labs(x = xlab, y = ylab, title = title) 
  
  # If flip_coord is TRUE, flip the axes.
  if (flip_coord) {
    p <- p + coord_flip()
  }
  
  return(p)
}

visualize_prior_effects <- function(df, priors, vis_path, statresult_path, formula_vars) {
  # Initialize lists to store fits, summary data frames, and model info
  fits <- list()
  summary_dfs <- list()
  model_infos <- list()
  
  # Loop through each prior specification in the provided list
  for (i in seq_along(priors)) {
    cat("Running analysis for prior specification", i, "\n")
    
    # Run the prior predictive check with the specific prior
    brms_setup <- prior_pred_check(
      df = df, 
      formula_vars,
      vis_path = vis_path, 
      run_pp = FALSE, 
      random_slope = FALSE, 
      prior = priors[[i]]
    )
    
    # Fit the model with the specified setup (ensuring a fresh run each time)
    run_brms_model(
      df = df, 
      brms_setup = brms_setup, 
      statresult_path = statresult_path, 
      force_refit = TRUE
    )
    
    # Load the fitted model from file
    fit_file <- file.path(statresult_path, paste0(brms_setup$filename, ".rds"))
    fit <- readRDS(fit_file)
    fits[[i]] <- fit
    
    # Interpret the fit to extract summaries and model info
    info_new <- interpret_fit(
      fit, 
      skip_rope_bf = TRUE, 
      rope_upper = 0.005, 
      complete_verbal_description = FALSE
    )
    summary_dfs[[i]] <- info_new$summary_df
    model_infos[[i]] <- info_new$model_info
  }
  
  # Optionally return the fits and associated information for further inspection
  return(list(fits = fits, summary_dfs = summary_dfs, model_infos = model_infos))
}
