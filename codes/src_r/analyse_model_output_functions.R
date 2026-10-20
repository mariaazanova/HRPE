# Load required libraries
library(bayestestR)
library(brms)
library(effectsize)
library(dplyr)
library(tidybayes)

library(knitr)
library(kableExtra)
library(htmltools)

set_seed <- 423
set.seed(set_seed)

compute_rope <- function(df, var_name_ibi, var_name_beh, sd_ibi = NA, sd_beh = NA, use_ci = FALSE) {
  # If sd values are not provided, compute them from the data
  if(is.na(sd_ibi)) {
    sd_ibi <- sd(df[[var_name_ibi]], na.rm = T)
  }
  
  if(is.factor(df[[var_name_beh]])) {
    rope_upper <- 1 / sd_ibi
  } else {
    # For a continuous variable, compute sd of the behavioral variable if not provided
    if(is.na(sd_beh)) {
      sd_beh <- sd(df[[var_name_beh]], na.rm = T)
    }
  
  # Define the range of the behavioral variable:
  # either using the full range (max - min) or the 95% CI (quantile 0.975 - quantile 0.025)
  if(use_ci) {
    beh_range <- as.numeric(quantile(df[[var_name_beh]], 0.975, na.rm = T) - quantile(df[[var_name_beh]], 0.025, na.rm = T))
  } else {
    beh_range <- max(df[[var_name_beh]], na.rm = T) - min(df[[var_name_beh]], na.rm = T)
  }
  
  # Compute the upper bound for the scaled coefficient's ROPE
  rope_upper <- (1 / beh_range) * (sd_beh / sd_ibi)
  }
  
  # Round to 4 decimal places and return
  rope_upper <- round(rope_upper, 4)
  return(rope_upper)
}

compute_rope_model <- function(df, vars, sd_ibi = NA, sd_list = list(), use_ci = FALSE, assess_int = FALSE, include_intercept = FALSE) {
  # The first element is the IBI variable.
  var_ibi <- vars[1]
  if(is.na(sd_ibi)) {
    sd_ibi <- sd(df[[var_ibi]])
  }
  
  rope_list <- list()
  
  if(include_intercept) {
    rope_intercept <- 1 / sd_ibi
    rope_list[["Intercept"]] <- round(rope_intercept, 4)
  }
  
  # Helper: compute range of a variable.
  # For factors we set the effective range to 1.
  get_range <- function(x, use_ci) {
    if(is.factor(x)) {
      return(1)
    } else {
      if(use_ci) {
        return(as.numeric(quantile(x, 0.975) - quantile(x, 0.025)))
      } else {
        return(max(x) - min(x))
      }
    }
  }
  
  # First, extract the behavior variables (skipping "*")
  beh_vars <- vars[-1]
  beh_vars <- beh_vars[beh_vars != "*"]
  
  # Compute ROPE for main effects
  for(v in beh_vars) {
    rope_val <- compute_rope(df, var_name_ibi = var_ibi, var_name_beh = v, sd_ibi = sd_ibi, use_ci = use_ci)
    rope_list[[v]] <- rope_val
  }
  
  # If assess_int is TRUE, also compute ROPE for interactions
  if(assess_int) {
    star_positions <- which(vars == "*")
    for(pos in star_positions) {
      # Ensure there is a variable before and after the "*" token
      if(pos > 1 && pos < length(vars)) {
        var1 <- vars[pos - 1]
        var2 <- vars[pos + 1]
        # Skip if somehow either is "*" (should not happen)
        if(var1 == "*" || var2 == "*") next
        
        # Get range for each variable (using either full range or CI)
        range1 <- get_range(df[[var1]], use_ci)
        range2 <- get_range(df[[var2]], use_ci)
        
        # Get standard deviations, either from supplied sd_list or computed on the fly.
        sd1 <- if(!is.null(sd_list[[var1]])) sd_list[[var1]] else if(is.numeric(df[[var1]])) sd(df[[var1]]) else NA
        sd2 <- if(!is.null(sd_list[[var2]])) sd_list[[var2]] else if(is.numeric(df[[var2]])) sd(df[[var2]]) else NA
        
        # For factor variables, assume effective range and sd are 1.
        if(is.factor(df[[var1]])) { range1 <- 1; sd1 <- 1 }
        if(is.factor(df[[var2]])) { range2 <- 1; sd2 <- 1 }
        
        # Compute rope for interaction:
        # The interaction effect (in original IBI units) for a full-range change is:
        #   (range1 * range2) * beta_int_scaled * [sd(IBI)/(sd1*sd2)]
        # Requiring this to be less than 1 implies:
        #   |beta_int_scaled| < (1/(range1 * range2)) * ((sd1 * sd2)/sd(IBI))
        rope_int <- (1 / (range1 * range2)) * ((sd1 * sd2) / sd_ibi)
        name_int <- paste(var1, var2, sep = ":")
        rope_list[[name_int]] <- round(rope_int, 4)
      }
    }
  }
  
  return(rope_list)
}

compute_rope_model_all_int <- function(df, vars, sd_ibi = NA, sd_list = list(), use_ci = FALSE, assess_int = FALSE, include_intercept = FALSE) {
  # The first element is the IBI variable.
  var_ibi <- vars[1]
  if (is.na(sd_ibi)) {
    sd_ibi <- sd(df[[var_ibi]], na.rm = T)
  }
  
  rope_list <- list()
  
  if (include_intercept) {
    rope_intercept <- 1 / sd_ibi
    rope_list[["Intercept"]] <- round(rope_intercept, 4)
  }
  
  # Helper: compute range of a variable.
  # For factors we set the effective range to 1.
  get_range <- function(x, use_ci) {
    if (is.factor(x)) {
      return(1)
    } else {
      if (use_ci) {
        return(as.numeric(quantile(x, 0.975, na.rm = T) - quantile(x, 0.025, na.rm = T)))
      } else {
        return(max(x) - min(x))
      }
    }
  }
  
  # Extract behavior tokens (may include "*" tokens)
  behavior_tokens <- vars[-1]
  # Remove "*" tokens to compute main effects.
  beh_vars_clean <- behavior_tokens[behavior_tokens != "*"]
  
  # Compute ROPE for main effects.
  for (v in beh_vars_clean) {
    rope_val <- compute_rope(df, var_name_ibi = var_ibi, var_name_beh = v, sd_ibi = sd_ibi, use_ci = use_ci)
    rope_list[[v]] <- rope_val
  }
  
  # Only compute interactions if:
  # (a) assess_int is TRUE, AND 
  # (b) there is at least one "*" token in the behavior tokens (i.e., explicit chaining is provided),
  # (c) and there is more than one predictor.
  if (assess_int && any(behavior_tokens == "*") && length(beh_vars_clean) > 1) {
    # Loop over interaction orders (2-way, 3-way, etc.)
    for (k in 2:length(beh_vars_clean)) {
      combs <- combn(beh_vars_clean, k, simplify = FALSE)
      for (comb in combs) {
        # Get ranges and sds for each variable in the combination.
        ranges <- sapply(comb, function(v) get_range(df[[v]], use_ci))
        sds <- sapply(comb, function(v) {
          if (!is.null(sd_list[[v]])) {
            sd_list[[v]]
          } else if (is.numeric(df[[v]])) {
            sd(df[[v]])
          } else {
            NA
          }
        })
        # For factor variables, force range and sd to 1.
        for (i in seq_along(comb)) {
          v <- comb[i]
          if (is.factor(df[[v]])) {
            ranges[i] <- 1
            sds[i] <- 1
          }
        }
        # Compute the ROPE for the interaction:
        rope_int <- (1 / prod(ranges)) * (prod(sds) / sd_ibi)
        name_int <- paste(comb, collapse = ":")
        rope_list[[name_int]] <- round(rope_int, 4)
      }
    }
    
    # Group the predictors according to "*" tokens.
    groups <- list()
    current_group <- c()
    for (i in seq_along(behavior_tokens)) {
      token <- behavior_tokens[i]
      if (token == "*") {
        # "*" indicates the next variable should join the current group.
        next
      } else {
        if (i == 1) {
          current_group <- token
        } else {
          if (behavior_tokens[i - 1] == "*") {
            current_group <- c(current_group, token)
          } else {
            groups[[length(groups) + 1]] <- current_group
            current_group <- token
          }
        }
      }
    }
    if (length(current_group) > 0) {
      groups[[length(groups) + 1]] <- current_group
    }
    
    # Filter interactions: keep an interaction only if all its variables are in one group.
    rope_list <- rope_list[sapply(names(rope_list), function(term) {
      # Always keep main effects and the intercept.
      if (!grepl(":", term)) {
        return(TRUE)
      }
      term_vars <- unlist(strsplit(term, ":"))
      any(sapply(groups, function(g) all(term_vars %in% g)))
    })]
  }
  
  # --- SAFETY GATE: drop interactions with repeated variables (e.g., a:a or a:b:a) ---
  rope_list <- rope_list[sapply(names(rope_list), function(term) {
    # keep intercept + main effects
    if (!grepl(":", term, fixed = TRUE)) return(TRUE)
    
    vars_in_term <- strsplit(term, ":", fixed = TRUE)[[1]]
    # drop if any variable repeats within the interaction term
    length(unique(vars_in_term)) == length(vars_in_term)
  })]
  
  # --- FINAL SAFETY GATE: keep only interaction terms whose order matches the formula token order ---
  # Use first occurrence order from behavior tokens (ignore "*" and keep unique in appearance order)
  order_vars <- unique(behavior_tokens[behavior_tokens != "*"])
  order_idx  <- setNames(seq_along(order_vars), order_vars)
  
  rope_list <- rope_list[sapply(names(rope_list), function(term) {
    # Keep intercept + main effects
    if (!grepl(":", term, fixed = TRUE)) return(TRUE)
    
    term_vars <- strsplit(term, ":", fixed = TRUE)[[1]]
    
    # If an interaction contains a var not in order_idx (shouldn't happen), keep it (ultra-safe).
    if (any(!term_vars %in% names(order_idx))) return(TRUE)
    
    idx <- unname(order_idx[term_vars])
    
    # Keep only if order is strictly increasing (matches formula order)
    all(diff(idx) > 0)
  })]
  
  return(rope_list)
}

# bayestestR also takes a named list of ROPE but it has to have the same names (scale_ with cont and level with factors)
create_rope_ranges <- function(rope_list, named = FALSE) {
  rope_ranges <- lapply(rope_list, function(x) c(-x, x))
  if (named){
    return(rope_ranges)
  } else {
    unname(rope_ranges)
  }
}

map_rope_names_based_on_df <- function(df, rope_list, coef_names) {
  
  # Helper function to map a single (non-interaction) variable name
  map_single <- function(var) {
    candidate1 <- paste0("b_", var)
    candidate2 <- paste0("b_scale_", var)
    candidate3 <- if (var %in% names(df) && is.factor(df[[var]])) {
      paste0("b_", var, levels(df[[var]])[2])
    } else {
      NA
    }
    
    if (candidate1 %in% coef_names) {
      return(candidate1)
    } else if (candidate2 %in% coef_names) {
      return(candidate2)
    } else if (!is.na(candidate3) && candidate3 %in% coef_names) {
      return(candidate3)
    } else {
      # If no mapping found, return NA (or var)
      return(NA)
    }
  }
  
  # Process each name in the rope_list
  new_names <- sapply(names(rope_list), function(nm) {
    # Check if this is an interaction term (contains ":")
    if (grepl(":", nm)) {
      # Split the interaction term into its parts
      parts <- unlist(strsplit(nm, ":"))
      # Map each part individually and strip the initial "b_" for re-joining
      mapped_parts <- sapply(parts, function(part) {
        single_map <- map_single(part)
        if (!is.na(single_map)) {
          sub("^b_", "", single_map)
        } else {
          part
        }
      })
      # Recombine with colon and add a single "b_" at the front
      final_candidate <- paste0("b_", paste(mapped_parts, collapse = ":"))
      if (final_candidate %in% coef_names) {
        return(final_candidate)
      } else {
        return(nm)
      }
    } else {
      # Non-interaction term
      single_map <- map_single(nm)
      if (!is.na(single_map)) {
        return(single_map)
      } else {
        return(nm)
      }
    }
  }, USE.NAMES = FALSE)
  
  names(rope_list) <- new_names
  return(rope_list)
}

compute_coef_info <- function(fit, coef_name, rope_upper, complete_verbal_description = TRUE) {

  if (startsWith(coef_name, "b_")) {
    bayes_coef <- coef_name
    hypothesis_coef <- sub("^b_", "", coef_name)
  } else {
    bayes_coef <- paste0("b_", coef_name)
    hypothesis_coef <- coef_name
  }
  
  ## 1. Comprehensive summary via describe_posterior()
  dp <- describe_posterior(
    fit,
    parameters = bayes_coef,
    centrality = "mode",
    ci_method = "hdi",
    test = "all",
    rope_range = c(-rope_upper, rope_upper)
  )
  
  ## 2. HDI at multiple levels (50%, 75%, 89%, 95%)
  hdi_res <- hdi(fit, parameters = bayes_coef, ci = c(0.5, 0.75, 0.89, 0.95))
  
  ## 3. Support intervals (SI) for BF thresholds 1, 3, and 10.
  si_res <- si(fit, BF = c(1, 3, 10), parameters = bayes_coef)
  
  ## 4. p_direction and corresponding p-value, and get other required info
  # pd_val <- p_direction(fit, parameter = bayes_coef)
  get_coef_mask <- as.character(dp$Parameter) == bayes_coef
  mode <- dp$Mode[get_coef_mask]
  hdi_ci <- c(dp$CI_low[get_coef_mask], dp$CI_high[get_coef_mask])
  pd_val <- dp$pd[get_coef_mask]
  p_val_pd  <- pd_to_p(pd_val)
  rope_overlap <- dp$ROPE_Percentage[get_coef_mask]
  p_val_rope <- dp$p_ROPE[get_coef_mask]
  rope_decision <- dp$ROPE_Equivalence[get_coef_mask]
  
  ## 5. Bayes Factors
  bf_against_null <- exp(dp$log_BF[get_coef_mask])
  bf_for_null  <-1 / bf_against_null 
  bf_rope_df <- bayesfactor_parameters(fit, null = c(-rope_upper, rope_upper), parameters = bayes_coef)
  bf_outside <- exp(bf_rope_df$log_BF[get_coef_mask])
  bf_inside <- 1 / bf_outside
  
  # Interpret the BF values using effectsize::interpret_bf()
  interp_alt    <- if(!is.na(bf_against_null)) effectsize::interpret_bf(bf_against_null) else NA
  interp_null     <- if(!is.na(bf_for_null)) effectsize::interpret_bf(bf_for_null) else NA
  interp_inside  <- if(!is.na(bf_inside)) effectsize::interpret_bf(bf_inside) else NA
  interp_outside <- if(!is.na(bf_outside)) effectsize::interpret_bf(bf_outside) else NA
  
  bayes_factors <- list(
    bf_against_null = bf_against_null,
    interpretation_for_alternative = interp_alt,
    bf_for_null = bf_for_null,
    interpretation_for_null = interp_null,
    rope_outside = bf_outside,
    interpretation_outside = interp_outside,
    rope_inside = bf_inside,
    interpretation_inside = interp_inside
  )
  
  # Assemble all outputs into one list
  result <- list(
    describe_posterior = dp,
    hdi = hdi_res,
    support_intervals = si_res,
    mode = mode,
    hdi_ci = hdi_ci, 
    pd_val = pd_val,
    p_val_pd  = p_val_pd, 
    rope_percent_overlap = rope_overlap, 
    p_val_rope = p_val_rope, 
    rope_decision = rope_decision,
    bayes_factors = bayes_factors
  )
  
  if (complete_verbal_description) {
    sig_text <- if(p_val_pd < 0.05) ", suggesting significant direction of effect" else ""
    verbal_description <- paste0(
      "The resulting coefficient has a mode of ", round(mode, 3),
      " and a 95% highest density interval of [", round(hdi_ci[1], 3), ", ", round(hdi_ci[2], 3), "]. ",
      "The probability of direction is ", round(pd_val, 4)*100,
      "% with a corresponding p-value of ", round(p_val_pd, 4), sig_text, ". ",
      "The Bayes Factor of ", round(bayes_factors$bf_against_null, 3),
      " shows ", bayes_factors$interpretation_for_alternative, " the alternative hypothesis, while the Bayes Factor for the point null is ", round(bayes_factors$bf_for_null, 3),
      " (", bayes_factors$interpretation_for_null, " the null). ",
      "The region of practical equivalence [-", rope_upper, ", ", rope_upper, "] overlaps with the posterior by ", 
      round(rope_overlap, 4)*100, "% with a p-value of ", round(p_val_rope, 4), ". ",
      "The Bayes Factor for the effect being outside the ROPE is ", round(bayes_factors$rope_outside, 3),
      " (", bayes_factors$interpretation_outside, " the alternative), and for the effect being negligible is ", 
      round(bayes_factors$rope_inside, 3), " (", bayes_factors$interpretation_inside, " the null)."
    )
    result$verbal_description <- verbal_description
  }
  
  return(result)
}

compute_coef_info_many <- function(fit, coef_name = NULL, rope_upper, complete_verbal_description = TRUE) {
    
    # If no parameters are provided, extract all variables that start with "b_"
    if (is.null(coef_name)) {
      all_vars <- get_variables(fit)
      coef_name <- all_vars[grepl("^b_", all_vars)]
    }
    coef_name <- as.character(coef_name)
    
    # Define two versions of parameter names:
    # - bayes_coef: with "b_" (for posterior summaries)
    # - hypothesis_coef: without "b_" (matching bayesfactor_parameters() output)
    bayes_coef <- sapply(coef_name, function(x) ifelse(startsWith(x, "b_"), x, paste0("b_", x)))
    hypothesis_coef <- sapply(coef_name, function(x) ifelse(startsWith(x, "b_"), sub("^b_", "", x), x))
    
    # 1. Comprehensive summary via describe_posterior()
    dp <- describe_posterior(
      fit,
      parameters = bayes_coef,
      centrality = "mode",
      ci_method = "hdi",
      test = "all",
      rope_range = c(-rope_upper, rope_upper)
    )
    # Add a helper column: parameter name without the "b_" prefix
    #dp <- dp %>% mutate(hypothesis_parameter = sub("^b_", "", Parameter))
    
    # 2. Compute HDIs at multiple levels (50%, 75%, 89%, 95%)
    hdi_res <- hdi(fit, parameters = bayes_coef, ci = c(0.5, 0.75, 0.89, 0.95))
    
    # 3. Compute support intervals (SI) for BF thresholds 1, 3, and 10.
    si_res <- si(fit, BF = c(1, 3, 10), parameters = bayes_coef)
    
    # 4. Build a summary data frame from dp (one row per parameter)
    summary_df <- dp %>% mutate(hypothesis_parameter = sub("^b_", "", Parameter)) %>% transmute(
      Parameter,
      Mode,
      CI_low,
      CI_high,
      pd,
      p_value_pd = pd_to_p(pd),
      ROPE_Percentage,
      p_ROPE,
      ROPE_Equivalence,
      hypothesis_parameter,
      bf_against_null = exp(log_BF),
      bf_for_null = 1/exp(log_BF)
    )
    
    # 6. Bayes Factors for the ROPE test (null = c(-rope_upper, rope_upper))
    bf_rope_df <- bayesfactor_parameters(fit, null = c(-rope_upper, rope_upper), effects = "all", parameters = bayes_coef)
    bf_rope_df <- bf_rope_df %>% 
      mutate(Parameter = sub("^b_", "", Parameter),
             BF = exp(log_BF))
    summary_df <- summary_df %>%
      left_join(bf_rope_df %>% select(Parameter, BF), 
                by = c("hypothesis_parameter" = "Parameter")) %>%
      rename(rope_outside = BF) %>%
      mutate(rope_inside = ifelse(!is.na(rope_outside) & rope_outside != 0, 1 / rope_outside, NA))
    
    # 7. Interpret the Bayes Factors using effectsize::interpret_bf()
    summary_df <- summary_df %>%
      mutate(
        interpretation_for_alternative = ifelse(!is.na(bf_against_null), effectsize::interpret_bf(bf_against_null), NA),
        interpretation_for_null = ifelse(!is.na(bf_for_null), effectsize::interpret_bf(bf_for_null), NA),
        interpretation_outside = ifelse(!is.na(rope_outside), effectsize::interpret_bf(rope_outside), NA),
        interpretation_inside = ifelse(!is.na(rope_inside), effectsize::interpret_bf(rope_inside), NA)
      )
    
    # 8. Build a verbal description for each parameter
    verbal_texts <- summary_df %>% rowwise() %>% mutate(
      sig_text = ifelse(p_value_pd < 0.05, ", suggesting significant direction of effect", ""),
      text = paste0(
        "The resulting coefficient for ", hypothesis_parameter, " has a mode of ", round(Mode, 3),
        " and a 95% highest density interval of [", round(CI_low, 3), ", ", round(CI_high, 3), "]. ",
        "The probability of direction is ", round(pd, 4)*100, "% with a corresponding p-value of ", round(p_value_pd, 4), sig_text, ". ",
        "The Bayes Factor of ", round(bf_against_null, 3), " shows ", interpretation_for_alternative, " the alternative hypothesis, while the Bayes Factor for the point null is ", round(bf_for_null, 3),
        " (", interpretation_for_null, " the null). ",
        "The region of practical equivalence [-", rope_upper, ", ", rope_upper, "] overlaps with the posterior by ", 
        round(ROPE_Percentage, 4)*100, "% with a p-value of ", round(p_ROPE, 4), ". ",
        "The Bayes Factor for the effect being outside the ROPE is ", round(rope_outside, 3),
        " (", interpretation_outside, " the alternative), and for the effect being negligible is ", 
        round(rope_inside, 3), " (", interpretation_inside, " the null)."
      )
    ) %>% pull(text)
    
    verbal_description <- paste(verbal_texts)
    
    # Assemble all outputs into one final result list
    result <- list(
      describe_posterior = dp,
      hdi = hdi_res,
      support_intervals = si_res,
      summary_df = summary_df,
      verbal_description = if (complete_verbal_description) verbal_description else NULL
    )
    
    return(result)
}

map_rope_names_based_on_fit <- function(rope_list, coef_names) {
  
  # Helper function to map a single (non-interaction) variable name using approximate matching.
  map_single <- function(var) {
    candidate1 <- paste0("b_", var)
    candidate2 <- paste0("b_scale_", var)
    
    if (candidate1 %in% coef_names) {
      return(candidate1)
    } else if (candidate2 %in% coef_names) {
      return(candidate2)
    } else {
      # If exact matches not found, search for any coefficient that contains 'var'
      matches <- grep(var, coef_names, value = TRUE)
      if (length(matches) > 0) {
        return(matches[1])
      } else {
        return(NA)
      }
    }
  }
  
  # Process each name in rope_list
  new_names <- sapply(names(rope_list), function(nm) {
    # Check if this is an interaction term (contains ":")
    if (grepl(":", nm)) {
      # Split the interaction term into its parts
      parts <- unlist(strsplit(nm, ":"))
      # Map each part individually and remove the "b_" prefix for re-joining
      mapped_parts <- sapply(parts, function(part) {
        single_map <- map_single(part)
        if (!is.na(single_map)) {
          sub("^b_", "", single_map)
        } else {
          part
        }
      })
      # Recombine with colon and add a single "b_" at the front
      final_candidate <- paste0("b_", paste(mapped_parts, collapse = ":"))
      if (final_candidate %in% coef_names) {
        return(final_candidate)
      } else {
        # If final candidate is not found exactly, try an approximate match
        matches <- grep(paste(mapped_parts, collapse = ":"), coef_names, value = TRUE)
        if (length(matches) > 0) {
          return(matches[1])
        } else {
          return(nm)
        }
      }
    } else {
      # Non-interaction term: use map_single
      single_map <- map_single(nm)
      if (!is.na(single_map)) {
        return(single_map)
      } else {
        return(nm)
      }
    }
  }, USE.NAMES = FALSE)
  
  names(rope_list) <- new_names
  return(rope_list)
}

interpret_fit <- function(fit, rope_upper=NULL, complete_verbal_description = TRUE, skip_rope_bf = FALSE) {
  
  all_vars <- get_variables(fit)
  coef_names <- all_vars[grepl("^b_", all_vars)]
  
  # Compute rope_range:
  # - If rope_upper is a single number, use the same range for all parameters.
  # - If rope_upper is a list, compute individual ROPE ranges.
  mapping_df <- NULL
  if (is.numeric(rope_upper) && length(rope_upper) == 1) {
    rope_range <- c(-rope_upper, rope_upper)
  } else if (is.list(rope_upper)) {
    # Save the original rope names.
    original_rope_names <- names(rope_upper)
    # Map each element using the map_rope_names function.
    mapped_rope_list <- map_rope_names_based_on_fit(rope_upper, coef_names = coef_names)
    # Create a mapping data frame: mapped names (names of mapped_rope_list) and the original names.
    mapping_df <- data.frame(
      mapped = names(mapped_rope_list),
      initial_names = original_rope_names,
      stringsAsFactors = FALSE
    )
    # Create individual ROPE ranges (named list) using the mapped rope list.
    rope_range <- create_rope_ranges(mapped_rope_list, named = TRUE)
  } else {
    rope_range <- NULL
  }
  
  # 1. Comprehensive summary via describe_posterior()
  dp <- describe_posterior(
    fit,
    centrality = "mode",
    ci_method = "hdi",
    test = "all",
    rope_range = rope_range
  )
  
  # 2. Compute HDIs at multiple levels (50%, 75%, 89%, 95%)
  hdi_res <- bayestestR::hdi(fit, ci = c(0.5, 0.75, 0.89, 0.95))
  
  # 3. Compute support intervals (SI) for BF thresholds 1, 3, and 10.
  si_res <- si(fit, BF = c(1, 3, 10))
  
  # 4. Build a summary data frame from dp 
  summary_df <- dp %>%
    mutate(hypothesis_parameter = sub("^b_", "", Parameter)) %>%
    transmute(
      Parameter,
      Mode,
      CI_low,
      CI_high,
      pd,
      p_value_pd = pd_to_p(pd),
      ROPE_Percentage,
      p_ROPE,
      ROPE_Equivalence,
      ROPE_low,
      ROPE_high,
      hypothesis_parameter,
      bf_against_null = exp(log_BF),
      bf_for_null = 1 / exp(log_BF)
    )
  
  # 5. Bayes Factors for the ROPE test (using rope_range as null)
  if (is.list(rope_range)) {
    if (skip_rope_bf) {
      # Create a data frame with NA for BF for each parameter in summary_df.
      bf_rope_df <- summary_df %>%
        select(hypothesis_parameter) %>%
        rename(Parameter = hypothesis_parameter) %>%
        mutate(BF = NA)
    } else {
    bf_rope_list <- lapply(names(rope_range), function(par) {
      current_range <- rope_range[[par]]
      bf_temp <- bayesfactor_parameters(fit, null = current_range, effects = "all", parameters = par)
      # Filter to ensure only the exact matching parameter is retained.
      bf_temp <- bf_temp[grepl(paste0("^", par, "$"), bf_temp$Parameter), ]
      bf_temp
    })
    bf_rope_df <- do.call(rbind, bf_rope_list)
    bf_rope_df <- bf_rope_df %>% 
      mutate(Parameter = sub("^b_", "", Parameter),
             BF = exp(log_BF))
    }
  } else {
    if (is.null(rope_range)){
      bf_rope_df <- summary_df %>%
        select(hypothesis_parameter) %>%
        rename(Parameter = hypothesis_parameter) %>%
        mutate(BF = NA)
    } else {
      bf_rope_df <- bayesfactor_parameters(fit, null = rope_range, effects = "all")
      bf_rope_df <- bf_rope_df %>% 
        mutate(Parameter = sub("^b_", "", Parameter),
               BF = exp(log_BF))
    }
  }
  
  summary_df <- summary_df %>%
    left_join(bf_rope_df %>% select(Parameter, BF), 
              by = c("hypothesis_parameter" = "Parameter")) %>%
    rename(rope_outside = BF) %>%
    mutate(rope_inside = ifelse(!is.na(rope_outside) & rope_outside != 0, 1 / rope_outside, NA))
  
  # 6. Interpret the Bayes Factors using effectsize::interpret_bf()
  summary_df <- summary_df %>%
    mutate(
      interpretation_for_alternative = ifelse(!is.na(bf_against_null),
                                              effectsize::interpret_bf(bf_against_null),
                                              NA),
      interpretation_for_null = ifelse(!is.na(bf_for_null),
                                       effectsize::interpret_bf(bf_for_null),
                                       NA),
      interpretation_outside = ifelse(!is.na(rope_outside),
                                      effectsize::interpret_bf(rope_outside),
                                      NA),
      interpretation_inside = ifelse(!is.na(rope_inside),
                                     effectsize::interpret_bf(rope_inside),
                                     NA)
    )
  
  if (!is.null(mapping_df)) {
    summary_df <- summary_df %>%
      left_join(mapping_df, by = c("Parameter" = "mapped"))
  }
  
  # 7. Add diagnostic info from summary(fit)$fixed.
  fix_sum <- as.data.frame(summary(fit)$fixed)
  fix_sum$hypothesis_parameter <- rownames(fix_sum)
  fix_sum <- fix_sum %>% select(hypothesis_parameter, Rhat, Bulk_ESS, Tail_ESS)
  summary_df <- summary_df %>%
    left_join(fix_sum, by = "hypothesis_parameter")
  
  # 8. Add prior information from prior_summary(fit).
  # First, fill down the 'prior' column so that missing values are replaced by the last non-empty value.
  prior_df <- prior_summary(fit) %>%
    mutate(prior = ifelse(prior == "", NA, prior)) %>%
    tidyr::fill(prior, .direction = "down")
  # Extract intercept prior and join others based on the coef column.
  intercept_prior <- prior_df %>% filter(class == "Intercept") %>% pull(prior)
  other_priors <- prior_df %>% filter(class == "b")
  summary_df <- summary_df %>%
    left_join(other_priors %>% select(coef, prior), by = c("hypothesis_parameter" = "coef")) %>%
    mutate(prior_info = ifelse(hypothesis_parameter == "Intercept",
                               ifelse(length(intercept_prior) > 0, intercept_prior, NA),
                               prior)) %>%
    select(-prior)
  
  ## Check for random slopes and add random effects summary columns
  random_slopes_exist <- FALSE
  # Check using the second row’s coefficient from summary_df (if available) to see if random slopes exist
  if (nrow(summary_df) >= 2) {
    test_coef <- summary_df$hypothesis_parameter[2]
    random_slopes_exist <- tryCatch({
      # Attempt to extract the estimated slopes using the provided indexing structure
      tmp <- ranef(fit)$subID[, , test_coef]
      TRUE
    }, error = function(e) FALSE)
  }
  if (random_slopes_exist) {
    # Extract random effects structures
    vc <- VarCorr(fit)
    ranef_group <- ranef(fit)$subID
    # Compute the additional metrics row by row using tryCatch to handle any issues per row individually
    summary_df <- summary_df %>%
      rowwise() %>%
      mutate(
        slope_variance = tryCatch({
          if (hypothesis_parameter %in% rownames(vc$subID$sd)) {
            # Extract the standard deviation and square it
            sd_val <- vc$subID$sd[hypothesis_parameter, "Estimate"]
            sd_val^2
          } else {
            NA_real_
          }
        }, error = function(e) NA_real_),
        min_slope = tryCatch({
          # Attempt to extract the estimated slopes using the three-dimensional indexing.
          if (hypothesis_parameter %in% dimnames(ranef_group)[[3]]) {
            slope_estimates <- ranef_group[, , hypothesis_parameter][, "Estimate"]
            min(slope_estimates, na.rm = TRUE)
          } else {
            NA_real_
          }
        }, error = function(e) NA_real_),
        max_slope = tryCatch({
          if (hypothesis_parameter %in% dimnames(ranef_group)[[3]]) {
            slope_estimates <- ranef_group[, , hypothesis_parameter][, "Estimate"]
            max(slope_estimates, na.rm = TRUE)
          } else {
            NA_real_
          }
        }, error = function(e) NA_real_)
      ) %>%
      ungroup()
  }
  
  # 9. Build a verbal description for each parameter.
  summary_df <- summary_df %>%
    rename(Model_parameter = Parameter) %>%
    mutate(Parameter = sapply(Model_parameter, map_variable_name))
  verbal_texts <- summary_df %>%
    rowwise() %>%
    mutate(
      sig_text = ifelse(p_value_pd < 0.05, ", suggesting significant direction of effect", ""),
      base_text = paste0(
        "The resulting coefficient for ", Parameter, " has a mode of ", round(Mode, 3),
        " and a 95% highest density interval of [", round(CI_low, 3), ", ", round(CI_high, 3), "]. ",
        "The probability of direction is ", round(pd * 100, 4), "% with a corresponding p-value of ", round(p_value_pd, 4), sig_text, ". ",
        "The Bayes Factor of ", round(bf_against_null, 3), " shows ", interpretation_for_alternative, " the alternative hypothesis, while the Bayes Factor for the point null is ", round(bf_for_null, 3),
        " (", interpretation_for_null, " the null). ",
        "The region of practical equivalence [", round(ROPE_low, 3), ", ", round(ROPE_high, 3), "] overlaps with the posterior by ", round(ROPE_Percentage * 100, 4), "% with a p-value of ", round(p_ROPE, 4), ". "
      ),
      rope_text = if (!skip_rope_bf) {
        paste0(
          "The Bayes Factor for the effect being outside the ROPE is ", round(rope_outside, 3),
          " (", interpretation_outside, " the alternative), and for the effect being negligible is ", round(rope_inside, 3),
          " (", interpretation_inside, " the null)."
        )
      } else {
        ""
      },
      text = paste0(base_text, rope_text)
    ) %>%
    ungroup() %>%
    select(hypothesis_parameter, text)
  
  # Create a named list where each coefficient's verbal description is accessible by its name.
  verbal_description <- setNames(as.list(verbal_texts$text), verbal_texts$hypothesis_parameter)
  
  # compute model fit metrics to add them to model info and to list of results
  loo_fit <- loo(fit)
  waic_fit <- waic(fit)
  bR2 <- bayes_R2(fit)
  
  # measures of correlation among posteriors
  posterior_samples <- as.matrix(fit, pars = "^b_")
  if (ncol(posterior_samples) > 1) {
    cor_matrix <- cor(posterior_samples)
    abs_cor <- abs(cor_matrix)
    diag(abs_cor) <- NA
    max_abs_corr <- max(abs_cor, na.rm = TRUE)
    eigen_values <- eigen(cor_matrix)$values
    condition_index <- sqrt(max(eigen_values) / min(eigen_values))
  } else {
    max_abs_corr <- NA
    condition_index <- NA
  }
  
  mod_sum <- summary(fit)
  n_groups <- if ("ngrps" %in% names(mod_sum)) mod_sum$ngrps else NA
  model_info <- list(
    family = fit$family$family,
    n_obs = nrow(fit$data),
    n_groups = n_groups, 
    looic = loo_fit$estimates["looic","Estimate"],
    p_loo = loo_fit$estimates["p_loo","Estimate"],
    waic = waic_fit$estimates["waic","Estimate"],
    p_waic = waic_fit$estimates["p_waic","Estimate"],
    bR2 = bR2["R2","Estimate"],
    max_abs_corr = max_abs_corr,
    condition_index = condition_index
  )
  
  # Assemble all outputs into one final result list
  result <- list(
    describe_posterior = dp,
    hdi = hdi_res,
    support_intervals = si_res,
    summary_df = summary_df,
    verbal_description = if (complete_verbal_description) verbal_description else NULL,
    model_info = model_info,
    loo_fit = loo_fit,
    waic_fit = waic_fit,
    bR2 = bR2
  )
  
  return(result)
}

convert_coef_to_orig_units <- function(coef, var_names, df, sd_dep = NA, sd_list = list()) {
  # The first element is the dependent variable.
  dep_var <- var_names[1]
  if (is.na(sd_dep)) {
    sd_dep <- sd(df[[dep_var]])
  }
  
  # Extract independent variables, ignoring any "*" tokens.
  indep_vars <- var_names[-1]
  indep_vars <- indep_vars[indep_vars != "*"]
  indep_vars <- indep_vars[indep_vars != ":"]
  
  # Compute standard deviations for the independent variables.
  sds <- sapply(indep_vars, function(v) {
    if (!is.null(sd_list[[v]])) {
      sd_list[[v]]
    } else if (is.numeric(df[[v]])) {
      sd(df[[v]])
    } else {
      # For factors we assume no scaling, effectively sd = 1.
      1
    }
  })
  
  # Convert the coefficient back to original units.
  # For interactions, multiply by the inverse of the product of the independent sds.
  orig_coef <- coef * (sd_dep / prod(sds))
  
  return(orig_coef)
}

convert_summary_to_orig_units <- function(summary_df, fit, df) {
  # Determine the base variable: first column name from fit$data, with "scale_" removed if present.
  base_var <- sub("^scale_", "", colnames(fit$data)[1])
  # Determine the original units for the base variable.
  orig_units <- map_variable_units(base_var)
  
  summary_df <- summary_df %>%
    rowwise() %>%
    mutate(
      mode_orig = if (Parameter == "Intercept") {
        NA
      } else {
        convert_coef_to_orig_units(
          Mode,
          c(
            base_var,
            if (grepl(":", initial_names)) {
              strsplit(initial_names, ":")[[1]]
            } else {
              initial_names
            }
          ),
          df
        )
      },
      CI_low_orig = if (Parameter == "Intercept") {
        NA
      } else {
        convert_coef_to_orig_units(
          CI_low,
          c(
            base_var,
            if (grepl(":", initial_names)) {
              strsplit(initial_names, ":")[[1]]
            } else {
              initial_names
            }
          ),
          df
        )
      },
      CI_high_orig = if (Parameter == "Intercept") {
        NA
      } else {
        convert_coef_to_orig_units(
          CI_high,
          c(
            base_var,
            if (grepl(":", initial_names)) {
              strsplit(initial_names, ":")[[1]]
            } else {
              initial_names
            }
          ),
          df
        )
      },
      original_units = orig_units
    ) %>%
    ungroup()
  
  return(summary_df)
}

create_combined_summary_table <- function(summary_df, model_info,
                                          cols = c("Parameter", "Mode", "CI_low", "CI_high", "pd", 
                                                   "ROPE_Percentage", "bf_against_null", "Rhat", "Bulk_ESS", "Tail_ESS"),
                                          params = NULL,
                                          readable_names = NULL) {
  if (!is.null(params)) {
    potential_cols <- c("Parameter", "Model_parameter", "hypothesis_parameter", "initial_names")
    existing_cols <- intersect(potential_cols, colnames(summary_df))
    # Keep rows if any of the existing parameter columns has an exact match to any of the provided params.
    summary_df <- summary_df %>% 
      filter(if_any(all_of(existing_cols), ~ . %in% params))
  }
  
  # Subset summary_df to the desired columns.
  if (!is.null(cols)) {
    short_table <- summary_df %>% select(all_of(cols)) # or any_of() if we don't want errors for missing cols
  }
  
  # Round numeric values
  short_table <- short_table %>% 
    mutate(across(where(is.numeric), ~ as.character(ifelse(abs(.) > 1000, round(.), round(., 3))))) %>%
    mutate(across(where(~ !is.numeric(.)), as.character))
  
  # If readable_names are provided and the length matches, rename columns.
  if (!is.null(readable_names) && length(readable_names) == length(cols)) {
    colnames(short_table) <- readable_names
    display_cols <- readable_names
  } else {
    display_cols <- cols
  }
  
  # Determine which column will hold the model info values; typically, the second column.
  value_col <- if (length(display_cols) >= 2) display_cols[2] else display_cols[1]
  
  # Construct a list for model info.
  model_info_list <- list(Family = model_info$family, n_obs = model_info$n_obs, n_groups = model_info$n_groups)
  n_model_info <- length(model_info_list)
  n_cols <- length(display_cols)
  
  # Create an empty data frame for model info rows, filled with empty strings.
  model_info_df <- as.data.frame(matrix("", nrow = n_model_info, ncol = n_cols), stringsAsFactors = FALSE)
  colnames(model_info_df) <- display_cols
  
  # Wrap model info labels in <b> tags using htmltools::HTML.
  model_info_labels <- sapply(names(model_info_list), function(x) as.character(HTML(paste0("<b>", x, "</b>"))))
  model_info_df[[display_cols[1]]] <- model_info_labels
  
  # Fill the value column with the corresponding model info values.
  model_info_df[[value_col]] <- sapply(model_info_list, function(x) {
    if (is.numeric(x)) formatC(x, format = "fg", digits = 3) else as.character(x)
  })
  
  # Combine the parameter estimates and model info rows.
  combined_df <- bind_rows(short_table, model_info_df)
  
  # Replace any remaining NA with an empty string.
  combined_df[is.na(combined_df)] <- ""
  
  # Identify the first row of the model info rows.
  first_model_info_row <- nrow(short_table) + 1
  
  # Create the final table.
  combined_table <- combined_df %>% 
    kable(caption = "Combined Parameter and Model Info Summary", align = "c", escape = FALSE, format = "html") %>%
    kable_styling(bootstrap_options = c("condensed"), full_width = FALSE) %>%
    row_spec(first_model_info_row, extra_css = "border-top: 2px solid black;")
  
  return(combined_table)
}

library(flextable)
library(magrittr)
library(officer)

create_combined_summary_table_flex <- function(summary_df, model_info,
                                               cols = c("Parameter", "Mode", "CI_low", "CI_high", "pd", 
                                                        "ROPE_Percentage", "bf_against_null", "Rhat", "Bulk_ESS", "Tail_ESS"),
                                               params = NULL,
                                               readable_names = NULL) {
  if (!is.null(params)) {
    potential_cols <- c("Parameter", "Model_parameter", "hypothesis_parameter", "initial_names")
    existing_cols <- intersect(potential_cols, colnames(summary_df))
    # Keep rows if any of the existing parameter columns has an exact match to any of the provided params.
    summary_df <- summary_df %>% 
      filter(if_any(all_of(existing_cols), ~ . %in% params))
  }
  
  # Subset summary_df to the desired columns.
  if (!is.null(cols)) {
    short_table <- summary_df %>% select(all_of(cols)) # or any_of() if we don't want errors for missing cols
  }
  
  # Round numeric values
  short_table <- short_table %>% 
    mutate(across(where(is.numeric), ~ as.character(ifelse(abs(.) > 1000, round(.), round(., 3))))) %>%
    mutate(across(where(~ !is.numeric(.)), as.character))
  
  # Map readable names if provided.
  if (!is.null(readable_names)) {
    num_rename <- min(length(readable_names), ncol(short_table))
    colnames(short_table)[1:num_rename] <- readable_names[1:num_rename]
  }
  display_cols <- colnames(short_table)
  
  # Determine which column will hold the model info values (typically the second column).
  value_col <- if (length(display_cols) >= 2) display_cols[2] else display_cols[1]
  
  # Build model info rows.
  model_info_list <- list(Family = model_info$family, n_obs = model_info$n_obs, n_groups = model_info$n_groups)
  n_model_info <- length(model_info_list)
  n_cols <- length(display_cols)
  
  # Create a data frame for model info rows, filled with empty strings.
  model_info_df <- as.data.frame(matrix("", nrow = n_model_info, ncol = n_cols), stringsAsFactors = FALSE)
  colnames(model_info_df) <- display_cols
  
  # In the first column, set the model info labels.
  model_info_df[[ display_cols[1] ]] <- names(model_info_list)
  
  # In the value column, set the corresponding model info values (rounded if numeric).
  model_info_df[[ value_col ]] <- sapply(model_info_list, function(x) {
    if (is.numeric(x)) formatC(x, format = "fg", digits = 3) else as.character(x)
  })
  
  # Combine the parameter summary and model info rows.
  combined_df <- bind_rows(short_table, model_info_df)
  combined_df[is.na(combined_df)] <- ""
  
  # Determine the first row index of model info rows.
  first_model_info_row <- nrow(short_table) + 1
  
  # Create the flextable.
  ft <- flextable(combined_df)
  
  # Use a built-in theme for a clean, white background.
  ft <- theme_box(ft)
  
  # Remove any alternate row shading by forcing a white background on all cells.
  ft <- bg(ft, bg = "white", part = "all")
  
  # Bold only the first column in the model info rows.
  ft <- bold(ft, i = first_model_info_row:nrow(combined_df), j = display_cols[1], bold = TRUE)
  
  # Add a horizontal line above the model info rows.
  # (Set border on the row just above the model info block.)
  ft <- hline(ft, i = first_model_info_row - 1, border = fp_border(color = "black", width = 2))
  
  # Autofit for a cleaner look.
  ft <- autofit(ft)
  
  return(ft)
}

library(dplyr)
library(flextable)
library(magrittr)
library(officer)  

# Helper function to process a single model's summary_df
process_summary_df <- function(summary_df, 
                               cols = c("Parameter", "Mode", "CI_low", "CI_high", "pd", 
                                        "ROPE_Percentage", "bf_against_null", "Rhat", "Bulk_ESS", "Tail_ESS"),
                               params = NULL,
                               readable_names = c("Parameter","Mode", "lower 95% CI", "higher 95% CI", "Prob.Dir.", "ROPE overlap", "BF", "RHat","Bulk ESS","Tail ESS"),
                               convert_CI_range = TRUE) {
  
  if (!is.null(params)) {
    potential_cols <- c("Parameter", "Model_parameter", "hypothesis_parameter", "initial_names")
    existing_cols <- intersect(potential_cols, colnames(summary_df))
    # Keep rows if any of the existing parameter columns has an exact match to any of the provided params.
    summary_df <- summary_df %>% 
      filter(if_any(all_of(existing_cols), ~ . %in% params))
  }
  
  # Subset summary_df to the desired columns.
  if (!is.null(cols)) {
  short_table <- summary_df %>% select(all_of(cols)) # or any_of() if we don't want errors for missing cols
  }
  
  # Round numeric values
  short_table <- short_table %>% 
    mutate(across(where(is.numeric), ~ as.character(ifelse(abs(.) > 100, round(.), round(., 3))))) %>%
    mutate(across(where(~ !is.numeric(.)), as.character))
  
  if (convert_CI_range && all(c("CI_low", "CI_high") %in% names(short_table))) {
    # Create the new column using rounded values.
    short_table <- short_table %>% 
      mutate(`95% CI` = paste0("[", CI_low, ",", CI_high, "]"))
    
    # Remove the original CI_low and CI_high columns.
    short_table <- short_table %>% select(-any_of(c("CI_low", "CI_high")))
    
    # Reconstruct the column order: remove CI_low and CI_high from the original 'cols'
    # and insert "95% CI" at the position of the first occurrence among them.
    if ("CI_low" %in% cols && "CI_high" %in% cols) {
      pos <- min(match("CI_low", cols), match("CI_high", cols))
      new_cols <- cols[!cols %in% c("CI_low", "CI_high")]
      new_cols <- append(new_cols, "95% CI", after = pos - 1)
      # Ensure the new order includes only columns present in short_table.
      new_cols <- new_cols[new_cols %in% names(short_table)]
      short_table <- short_table[, new_cols, drop = FALSE]
      
      # Adjust the readable_names similarly.
      if (!is.null(readable_names)) {
        idx_low <- match("CI_low", cols)
        idx_high <- match("CI_high", cols)
        # Remove both CI-related entries.
        readable_names <- readable_names[-c(idx_low, idx_high)]
        # Insert "95% CI" at the original position.
        readable_names <- append(readable_names, "95% CI", after = pos - 1)
      }
    }
  }
  
  # Map readable names if provided.
  if (!is.null(readable_names)) {
    num_rename <- min(length(readable_names), ncol(short_table))
    colnames(short_table)[1:num_rename] <- readable_names[1:num_rename]
  }
  
  return(short_table)
}

combine_model_summaries_flex <- function(list_summaries, list_model_info, 
                                         param_col = "Parameter",
                                         cols = c("Parameter", "Mode", "CI_low", "CI_high", "pd", 
                                                  "ROPE_Percentage", "bf_against_null", "Rhat", "Bulk_ESS", "Tail_ESS"),
                                         params = NULL,
                                         readable_names = c("Parameter","Mode", "lower 95% CI", "higher 95% CI", "Prob.Dir.", "ROPE overlap", "BF", "RHat","Bulk ESS","Tail ESS"),
                                         model_names = NULL,
                                         dependent_vars = NULL,
                                         convert_CI_range = TRUE,
                                         add_model_names = TRUE) {
  # Check: dependent_vars must be provided and match the number of models.
  if (is.null(dependent_vars)) {
    stop("Please supply a dependent_vars vector with one element per model.")
  }
  if (length(dependent_vars) != length(list_summaries)) {
    stop("Length of dependent_vars must equal the number of models.")
  }
  
  # If model_names not provided, default to "Model 1", "Model 2", etc.
  if (is.null(model_names)) {
    model_names <- paste("Model", seq_along(list_summaries))
  }
  
  # Process each model's summary_df.
  list_processed <- lapply(list_summaries, process_summary_df, cols = cols, params = params, readable_names = readable_names, convert_CI_range = convert_CI_range)
  
  # For each processed summary, combine all columns except the parameter column into one multi-line string.
  list_combined <- lapply(list_processed, function(df) {
    available_measures <- setdiff(colnames(df), param_col)
    if (length(available_measures) > 0) {
      df$combined <- apply(df[, available_measures, drop = FALSE], 1, function(x) {
        paste(paste0(available_measures, ": ", x), collapse = "\n")
      })
    } else {
      df$combined <- ""
    }
    df %>% select(!!sym(param_col), combined)
  })
  
  # Merge all processed summaries by the parameter column (full outer join).
  master_df <- Reduce(function(x, y) merge(x, y, by = param_col, all = TRUE, sort = FALSE), list_combined)
  master_df[is.na(master_df)] <- ""
  # Rename merged columns (other than the parameter) with the unique model_names.
  if (ncol(master_df) > 1) {
    colnames(master_df)[-1] <- model_names
  }
  
  # Create flextable from master_df.
  ft <- flextable(master_df)
  ft <- autofit(ft)
  #ft <- set_table_properties(ft, width = 0.8, layout = "autofit")
  #ft <- fontsize(ft, size = 8, part = "all")
  
  # Remove default header.
  ft <- delete_part(ft, part = "header")
  
  # --- Add custom header row ---
  # We want one header row where:
  # - The first cell is "Parameter"
  # - The remaining cells display a two-line label:
  #    "Dependent variable:" (line 1) and the supplied dependent_vars string (line 2)
  header_row <- c("Parameter", sapply(dependent_vars, function(x) paste("Dependent variable:", x, sep = "\n")))
  ft <- add_header_row(ft, values = header_row, colwidths = rep(1, length(header_row)))
  if (add_model_names){
    header_row2 <- c("", model_names)
    ft <- add_header_row(ft, values = header_row2, colwidths = rep(1, length(header_row2)))
  }
  
  # Add a bold horizontal line (thick bottom border) below the header row.
  ft <- border(ft, part = "header", border.bottom = fp_border(color = "black", width = 2))
  
  # --- Add footer row with model-specific global info ---
  # Build a string for each model (e.g., "Family: X\nn_obs: Y\nn_groups: Z")
  model_info_strings <- sapply(list_model_info, function(mi) {
    paste0("Family: ", mi$family, "\n",
           "Num. obs: ", mi$n_obs, "\n",
           "Num. groups: ", mi$n_groups, "\n",
           "LOOIC: ", round(mi$looic, 2), "\n",
           "pLOOIC: ", round(mi$p_loo, 2), "\n",
           "WAIC: ", round(mi$waic, 2), "\n",
           "pWAIC: ", round(mi$p_waic, 2), "\n",
           "Bayes R2: ", round(mi$bR2, 3))
  })
  footer_values <- c("", model_info_strings)
  ft <- add_footer_row(ft, values = footer_values, colwidths = rep(1, length(footer_values)))
  ft <- align(ft, part = "footer", align = "left")
  ft <- border(ft, part = "footer", border.top = fp_border(color = "black", width = 2))
  
  # Add a bold horizontal line above the footer in the body.
  #ft <- hline(ft, i = nrow(master_df), part = "body", border = fp_border(color = "black", width = 2))
  
  return(ft)
}

group_model_formulas <- function(model_formulas, ignore_last = FALSE) {
  
  # Build tail keys for each formula
  tail_keys <- sapply(model_formulas, function(formula) {
    # If ignore_last is TRUE, skip the first and last element, else skip only the first element
    if (ignore_last) {
      tail_elements <- formula[-c(1, length(formula))]
    } else {
      tail_elements <- formula[-1]
    }
    paste(tail_elements, collapse = " ")
  })
  
  # Get the unique tail keys
  unique_keys <- unique(tail_keys)
  
  groups <- list()
  
  for (key in unique_keys) {
    # Get indices of formulas with this tail
    indices <- which(tail_keys == key)
    
    # The common tail is taken from the first matching formula
    common_tail <- model_formulas[[indices[1]]][-1]
    
    # Collect the first elements (the changeIBIX_X parts)
    first_elements <- sapply(model_formulas[indices], function(formula) formula[1])
    
    # Optionally, you might want to keep the full formulas as well:
    group_formulas <- model_formulas[indices]
    
    groups[[key]] <- list(common = common_tail,
                          firsts = first_elements,
                          formulas = group_formulas)
  }
  
  return(groups)
}

formula_list_difference <- function(list1, list2) {
  # Helper function: Convert a formula vector to a single string
  formula_to_string <- function(formula_vec) {
    paste(formula_vec, collapse = " ")
  }
  
  # Convert both lists to vectors of strings
  list1_str <- sapply(list1, formula_to_string)
  list2_str <- sapply(list2, formula_to_string)
  
  # Identify formulas that are unique to list1
  unique1 <- list1[!(list1_str %in% list2_str)]
  # Identify formulas that are unique to list2
  unique2 <- list2[!(list2_str %in% list1_str)]
  
  # Combine the unique formulas from both lists
  difference <- c(unique1, unique2)
  return(difference)
}

name_fits <- function(fits_list, model_formulas) {
  # Check if fits_list contains both fits and infos
  if (is.list(fits_list) && !is.null(fits_list$fits)) {
    fits <- fits_list$fits
    infos <- fits_list$infos
  } else {
    fits <- fits_list
    infos <- NULL
  }
  named_fits <- list()
  named_infos <- list()
  
  for (i in seq_along(model_formulas)) {
    formula_str <- paste(model_formulas[[i]], collapse = " ")
    named_fits[[formula_str]] <- fits[[i]]
    if (!is.null(infos)) {
      named_infos[[formula_str]] <- infos[[i]]
    }
  }
  if (!is.null(infos)) {
    return(list(fits = named_fits, infos = named_infos))  # CHANGE: return a list containing both
  } else {
    return(named_fits)
  }
}

process_grouped_fits <- function(grouped_models, named_fits, df, rope_upper = NULL, disable_plots = TRUE) {
  group_outputs <- list()
  
  nms <- names(grouped_models)
  if (is.null(nms) || any(nms == "")) {
    nms <- ifelse(is.na(nms) | nms == "",
                  as.character(seq_along(grouped_models)),
                  nms)
    names(grouped_models) <- nms
  }
  
  # Iterate over each group (the names come from the grouping key, e.g., the common tail)
  for (group_key in names(grouped_models)) {
    group <- grouped_models[[group_key]]
    # group$formulas is a list of formulas in this group
    formulas <- group$formulas
    
    # Prepare lists to collect outputs for this group
    group_fits <- list()
    summary_list <- list()
    model_info_list <- list()
    
    # Process each formula in the group
    for (j in seq_along(formulas)) {
      formula <- formulas[[j]]
      formula_str <- paste(formula, collapse = " ")
      
      # Retrieve the corresponding fit from the named fits list
      if (is.list(named_fits) && !is.null(named_fits$fits)) {
        fit <- named_fits$fits[[formula_str]]
        info <- named_fits$infos[[formula_str]]
      } else {
        fit <- named_fits[[formula_str]]
        info <- NULL
      }
      
      if (is.null(fit)) {
        warning(paste("Fit not found for formula:", formula_str))
        next
      }
      
      group_fits[[j]] <- fit
      
      # If info is already available (pre-computed), use it
      if (is.null(info)) {
        if (is.factor(df[[formula[1]]])) {
          info <- interpret_fit(fit, 0.005, skip_rope_bf = TRUE)
        } else if (!is.null(rope_upper)){
          info <- interpret_fit(fit, rope_upper)
          info$summary_df <- convert_summary_to_orig_units(summary_df = info$summary_df, fit = fit, df)
        } else {
          rope <- compute_rope_model_all_int(df, formula, assess_int = TRUE, include_intercept = TRUE)
          info <- interpret_fit(fit, rope)
          info$summary_df <- convert_summary_to_orig_units(summary_df = info$summary_df, fit = fit, df)
        }
      }
      
      summary_list[[j]] <- info$summary_df
      model_info_list[[j]] <- info$model_info
      
      while (sink.number(type = "output") > 0) {
        sink(NULL, type = "output")
      }
      
    }
    
    # Only if we have at least one valid fit do we continue:
    if (length(group_fits) > 0) {
      # For dependent_vars, extract from the first fit's data (adjust if needed)
      dep_vars <- lapply(group_fits, function(fit) {
        map_variable_name(sub("^scale_", "", colnames(fit$data)[1]))
      })
      
      # Combine summaries into a flex table
      summary_table <- combine_model_summaries_flex(
        list_summaries = summary_list,
        list_model_info = model_info_list,
        readable_names = c("Parameter", "Mode", "lower 95% CI", "higher 95% CI", "BF", "ROPE Bound", "% in ROPE", "BF in ROPE", "RHat"), # "% in ROPE"  # maybe also rope bounds then? # , "Bulk ESS", "Tail ESS"
        dependent_vars = dep_vars,
        convert_CI_range = TRUE,
        add_model_names = FALSE, 
        cols = c("Parameter", "Mode", "CI_low", "CI_high", "bf_against_null", "ROPE_high", "ROPE_Percentage", "rope_inside", "Rhat") # "ROPE_Percentage" # maybe also rope bounds then? # , "Bulk_ESS", "Tail_ESS"
      )
      
      # Generate a comparison plot
      if (disable_plots) {
        comparison_plot <- NULL
      } else {
        comparison_plot <- compare_brms_model_plots(
          group_fits,                  # list of fits in this group
          summary_list,                # corresponding list of summary data frames
          show_model_info = TRUE,
          model_info = model_info_list, # list of model_info objects
          extract_rope_from_summary = TRUE,
          display_title = TRUE,
          legend_pos = 'bottom'
        )
      }
    } else {
      summary_table <- NULL
      comparison_plot <- NULL
    }
    
    # Save the outputs for this group in a named list entry
    group_outputs[[group_key]] <- list(
      summary_table = summary_table,
      comparison_plot = comparison_plot,
      fits = group_fits,
      summaries = summary_list,
      model_info = model_info_list
    )
  }
  
  return(group_outputs)
}

run_model_pipeline <- function(model_formulas, df, statresult_path, modelsummary_path = NULL, 
                               prior_pred_args = list(run_pp = FALSE, statresult_path=statresult_path), ignore_last = FALSE,
                               force_refit = FALSE, force_info = FALSE, iter = 4000, rope_upper = NULL, disable_plots = TRUE) {  
  model_formulas <- lapply(model_formulas, function(vec) unlist(lapply(vec, split_star)))
  model_formulas <- model_formulas[!duplicated(sapply(model_formulas, paste, collapse = " "))]
  run_brms_models(model_formulas, df, statresult_path, modelsummary_path, prior_pred_args, force_refit = force_refit, force_info = force_info, iter = iter, rope_upper = rope_upper)  
  fits_list <- load_fits(model_formulas, df, statresult_path, modelsummary_path, prior_pred_args) 
  named_fits <- name_fits(fits_list, model_formulas)
  fits_list <- NULL
  grouped_models <- group_model_formulas(model_formulas, ignore_last)
  group_results <- process_grouped_fits(grouped_models, named_fits, df, rope_upper = rope_upper, disable_plots = disable_plots)
  return(group_results)
}

split_df_by_cont_var <- function(df, var_name) {
  if (!var_name %in% names(df)) {
    stop("Variable not found in data frame")
  }
  # Calculate the median value of the specified variable
  med <- median(df[[var_name]], na.rm = TRUE)
  
  # Create a grouping variable with names based on var_name and group ("low" or "high")
  group_names <- ifelse(df[[var_name]] <= med, paste0(var_name, "_low"), paste0(var_name, "_high"))
  
  # Split the data frame based on the grouping variable
  split_df <- split(df, group_names)
  
  return(split_df)
}

split_df_by_group_var <- function(df, var_name, assign_env = .GlobalEnv) {
  # 1) sanity check
  if (!var_name %in% names(df)) {
    stop("Variable '", var_name, "' not found in data frame")
  }
  # 2) extract the vector of levels (dropping any NAs)
  grp <- df[[var_name]]
  if (!is.factor(grp)) grp <- as.character(grp)
  lvls <- unique(na.omit(grp))
  # 3) split into a named list
  split_list <- split(df, grp)
  # 4) rename each element to "VarName_Level"
  names(split_list) <- paste0(var_name, "_", names(split_list))

  return(split_list)
}

create_word_from_groups <- function(group_results, chapter_path, filename) {
  # ensure output directory exists
  if (!dir.exists(chapter_path)) {
    dir.create(chapter_path, recursive = TRUE)
  }
  
  # start a fresh Word document
  doc <- read_docx()
  
  # for each group, add a heading + its flextable
  for (group_name in names(group_results)) {
    title_group_name <- map_variable_name(group_name)
    
    doc <- body_add_par(doc, title_group_name, style = "heading 1")
    
    ft <- group_results[[group_name]]$summary_table
    ft <- autofit(ft)
    ft <- set_table_properties(ft, layout = "autofit", width = 1)
    ft <- fontsize(ft, size = 9, part = "all")
    
    doc <- body_add_flextable(doc, value = ft)
    doc <- body_add_par(doc, "", style = "Normal")
  }
  
  # compose full output path (add .docx if missing)
  if (!grepl("\\.docx$", filename, ignore.case = TRUE)) {
    filename <- paste0(filename, ".docx")
  }
  out_file <- file.path(chapter_path, filename)
  out_file <- handle_long_filename(out_file)
  
  # write file to disk
  print(doc, target = out_file)
  message("Written: ", out_file)
}
