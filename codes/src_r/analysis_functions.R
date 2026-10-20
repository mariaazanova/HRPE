Sys.setenv(
  SDKROOT = "/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk",
  CXXFLAGS = "-isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk",
  CFLAGS   = "-isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk",
  CPPFLAGS = "-isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk",
  LDFLAGS  = "-isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk"
)

library(brms)
options(brms.backend = "cmdstanr")

set_seed <- 423
set.seed(set_seed)

generate_model_formulas <- function(main = TRUE,
                                    pure_IBI = FALSE,
                                    pure_IBI_prev = FALSE,
                                    int_prev_ibi = FALSE,
                                    phase_ms = FALSE,
                                    other_r = FALSE,
                                    only_other_r = FALSE,
                                    no_phase_int = FALSE,
                                    only_phase_int = FALSE,
                                    add_controls = FALSE,
                                    only_controls = FALSE,
                                    star_controls = FALSE,
                                    controls = NULL) {
  # Initialize the list to store formulas
  model_formulas <- list()
  
  # Check that at least one primary option is TRUE
  if (!main && !pure_IBI) {
    warning("At least one of the options 'main' or 'pure_IBI' should be set to TRUE")
    return(model_formulas)
  }
  
  # ---- Define First Options ----
  first_options <- list(c("expectation"), c("outcome"), c("PE"), c("absPE"), c("expectation", "*", "outcome"), c("absPE", "*", "outcome"))
  
  # ---- Define Second Options ----
  second_options_with_star <- list(c("*", "phaseFeedb"), c("*", "sysDiaFeedb"))
  second_options_without_star <- list(c("phaseFeedb"), c("sysDiaFeedb"))
  if (phase_ms) {
    # Add the additional msDistToFeedb option to both with and without star options
    second_options_with_star <- c(second_options_with_star, list(c("*", "msDistToFeedb")))
    second_options_without_star <- c(second_options_without_star, list(c("msDistToFeedb")))
  }
  
  # Container for formulas before adding controls
  formulas_no_control <- list()
  
  # ---- Build Main Formulas (using "changeIBI") ----
  if (main) {
    if (other_r) {
      if (only_other_r) {
        # Only the additional ones (3 versions)
        main_versions <- c("2_1", "3_2", "8_7")
      } else {
        # Full extended list (7 versions)
        main_versions <- c("2_1", "3_2", "4_3", "5_4", "6_5", "7_6", "8_7")
      }
    } else {
      # Base case (4 versions)
      main_versions <- c("4_3", "5_4", "6_5", "7_6")
    }
    
    for (v in main_versions) {
      base <- paste0("changeIBI", v)
      for (f in first_options) {
        if (!no_phase_int) {
          for (s in second_options_with_star) {
            formulas_no_control[[length(formulas_no_control) + 1]] <- c(base, f, s)
          }
        } 
        if (!only_phase_int) {
          for (s in second_options_without_star) {
            formulas_no_control[[length(formulas_no_control) + 1]] <- c(base, f, s)
          }
        }
      }
    }
  }
  
  # ---- Build Pure IBI Formulas (using "IBI") ----
  if (pure_IBI) {
    if (other_r) {
      if (only_other_r) {
        # Only the additional ones (3 versions)
        pure_versions <- c("2", "3", "8")
      } else {
        # Full extended list (7 versions)
        pure_versions <- c("2", "3", "4", "5", "6", "7", "8")
      }
    } else {
      # Base case (4 versions)
      pure_versions <- c("4", "5", "6", "7")
    }
    
    for (v in pure_versions) {
      base <- paste0("IBI", v)
      for (f in first_options) {
        if (!no_phase_int) {
          for (s in second_options_with_star) {
            formulas_no_control[[length(formulas_no_control) + 1]] <- c(base, f, s)
          }
        }
        for (s in second_options_without_star) {
          formulas_no_control[[length(formulas_no_control) + 1]] <- c(base, f, s)
        }
      }
    }
  }
  
  # ---- Add Controls ----
  formulas_with_controls <- list()
  if (add_controls && !is.null(controls) && length(controls) > 0) {
    for (formula_vec in formulas_no_control) {
      for (control in controls) {
        if (star_controls) {
          new_formula <- c(formula_vec, "*", control)
        } else {
          new_formula <- c(formula_vec, control)
        }
        formulas_with_controls[[length(formulas_with_controls) + 1]] <- new_formula
      }
    }
  }
  
  # ---- Final Output Selection ----
  if (only_controls) {
    model_formulas <- formulas_with_controls
  } else {
    model_formulas <- c(formulas_no_control, formulas_with_controls)
  }
  
  # ---- Append Pure IBI Previous (pure_IBI_prev) Last ----
  # This applies only to pure IBI formulas (those starting with "IBI" but not "changeIBI")
  if (pure_IBI_prev) {
    model_formulas <- lapply(model_formulas, function(formula_vec) {
      if (startsWith(formula_vec[1], "IBI") && !startsWith(formula_vec[1], "changeIBI")) {
        version_num <- as.numeric(sub("IBI", "", formula_vec[1]))
        if (!is.na(version_num) && version_num > 1) {
          prev_version <- as.character(version_num - 1)
          if (int_prev_ibi) {
            # Append "*" then the previous IBI variable at the very end
            formula_vec <- c(formula_vec, "*", paste0("IBI", prev_version))
          } else {
            formula_vec <- c(formula_vec, paste0("IBI", prev_version))
          }
        }
      }
      formula_vec
    })
  }
  
  model_formulas <- lapply(model_formulas, function(vec) unlist(lapply(vec, split_star_space)))
  return(model_formulas)
}

generate_model_formulas_alternative <- function(phase_ms = FALSE,
                                                no_phase_int = FALSE,
                                                only_phase_int = FALSE,
                                                add_controls = FALSE,
                                                only_controls = FALSE,
                                                star_controls = FALSE,
                                                controls = NULL, 
                                                add_stim = FALSE) {
  # Use the full set of provided first-variable names
  first_vars <- c("afterFeedbAveIBI", "afterFeedbRMSSD", "afterFeedbSDNN",
                  "beforeFeedbAveIBI", "beforeFeedbRMSSD", "beforeFeedbSDNN")
  if (add_stim){
    first_vars <- c("beforeStimAveIBI", "beforeStimRMSSD", "beforeStimSDNN",
                    "afterStimAveIBI", "afterStimRMSSD", "afterStimSDNN",
                    "afterFeedbAveIBI", "afterFeedbRMSSD", "afterFeedbSDNN",
                    "beforeFeedbAveIBI", "beforeFeedbRMSSD", "beforeFeedbSDNN")
  }
  
  # ---- Define First Options ----
  first_options <- list(c("expectation"), c("outcome"), c("PE"), c("absPE"), c("expectation", "*", "outcome"), c("absPE", "*", "outcome"))
  
  # ---- Define Second Options ----
  second_options_with_star <- list(c("*", "phaseFeedb"), c("*", "sysDiaFeedb"))
  second_options_without_star <- list(c("phaseFeedb"), c("sysDiaFeedb"))
  if (phase_ms) {
    # Add the additional msDistToFeedb option to both with and without star options
    second_options_with_star <- c(second_options_with_star, list(c("*", "msDistToFeedb")))
    second_options_without_star <- c(second_options_without_star, list(c("msDistToFeedb")))
  }
  
  # Initialize container for formulas
  model_formulas <- list()
  
  # ---- Build Formulas using all first-variable names ----
  for (base in first_vars) {
    for (f in first_options) {
      if (!no_phase_int) {
        for (s in second_options_with_star) {
          model_formulas[[length(model_formulas) + 1]] <- c(base, f, s)
        }
      }
      if (!only_phase_int) {
        for (s in second_options_without_star) {
          model_formulas[[length(model_formulas) + 1]] <- c(base, f, s)
        }
      }
    }
  }
  
  # ---- Add Controls if specified ----
  formulas_with_controls <- list()
  if (add_controls && !is.null(controls) && length(controls) > 0) {
    for (formula_vec in model_formulas) {
      for (control in controls) {
        if (star_controls) {
          new_formula <- c(formula_vec, "*", control)
        } else {
          new_formula <- c(formula_vec, control)
        }
        formulas_with_controls[[length(formulas_with_controls) + 1]] <- new_formula
      }
    }
  }
  
  # ---- Final Output Selection ----
  if (only_controls) {
    final_formulas <- formulas_with_controls
  } else {
    final_formulas <- c(model_formulas, formulas_with_controls)
  }
  
  final_formulas <- lapply(final_formulas, function(vec) unlist(lapply(vec, split_star_space)))
  return(final_formulas)
}

generate_model_formulas_performance <- function(phase_factor = FALSE,
                                                phase_ms = FALSE,
                                                no_phase_int = FALSE,
                                                only_phase_int = FALSE,
                                                add_controls = FALSE,
                                                only_controls = FALSE,
                                                star_controls = FALSE,
                                                controls = NULL,
                                                no_beh_control = TRUE,
                                                beh_control_int = FALSE,
                                                only_correct_nextTrial = TRUE) {
  # Define dependent variables (first element in the formula)
  if (only_correct_nextTrial){
    dep_vars <- c("trueCorrect_nextTrial")
  } else {
    dep_vars <- c("trueCorrect_nextTrial", 
                "switchStay_nextTrialAction", 
                "switchStay_sameCueAction", 
                "switchStay_nextTrialStrategy")
  }
  
  # Define the main versions (base case) for changeIBI formulas
  main_versions <- c("4_3", "5_4", "6_5", "7_6")
  
  # ---- Define First Options ----
  first_options <- list(c("expectation", "*", "outcome"), c("absPE", "*", "outcome"))
  if (no_beh_control) {
    first_options <- c(first_options, list(c()))
  }
  
  # Add beh_control_int option if true and no_beh_control is false.
  if (beh_control_int) {
    first_options <- lapply(first_options, function(x) {
      if (length(x) > 0) c(x, "*") else x
    })
  }
  
  # ---- Define Second Options ----
  second_options_with_star <- list(c("*", "phaseFeedb"))
  second_options_without_star <- list(c("phaseFeedb"))
  if (phase_factor) {
    # Add the additional msDistToFeedb option to both with-star and without-star lists
    second_options_with_star <- c(second_options_with_star, list(c("*", "sysDiaFeedb")))
    second_options_without_star <- c(second_options_without_star, list(c("sysDiaFeedb")))
  }
  if (phase_ms) {
    # Add the additional msDistToFeedb option to both with-star and without-star lists
    second_options_with_star <- c(second_options_with_star, list(c("*", "msDistToFeedb")))
    second_options_without_star <- c(second_options_without_star, list(c("msDistToFeedb")))
  }
  
  # Container for base formulas
  formulas <- list()
  
  # ---- Build Base Formulas ----
  for (dv in dep_vars) {
    for (v in main_versions) {
      base <- paste0("changeIBI", v)
      for (f in first_options) {
        if (!no_phase_int) {
          for (s in second_options_with_star) {
            formulas[[length(formulas) + 1]] <- c(dv, f, base, s)
          }
        }
        if (!only_phase_int){
          for (s in second_options_without_star) {
            formulas[[length(formulas) + 1]] <- c(dv, f, base, s)
          }
        }
      }
    }
  }
  
  # ---- Add Controls if Specified ----
  formulas_with_controls <- list()
  if (add_controls && !is.null(controls) && length(controls) > 0) {
    for (formula_vec in formulas) {
      for (control in controls) {
        if (star_controls) {
          new_formula <- c(formula_vec, "*", control)
        } else {
          new_formula <- c(formula_vec, control)
        }
        formulas_with_controls[[length(formulas_with_controls) + 1]] <- new_formula
      }
    }
  }
  
  # ---- Final Output Selection ----
  if (only_controls) {
    final_formulas <- formulas_with_controls
  } else {
    final_formulas <- c(formulas, formulas_with_controls)
  }
  
  final_formulas <- lapply(final_formulas, function(vec) unlist(lapply(vec, split_star_space)))
  return(final_formulas)
}

replace_double_underscore_in_colnames <- function(df) {
  colnames(df) <- gsub("__", "_", colnames(df))
  return(df)
}

generate_model_formulas_neural <- function(phase_factor = TRUE,
                                           phase_ms = FALSE,
                                           no_phase_int = FALSE,
                                           only_phase_int = FALSE,
                                           no_phase_at_all = FALSE,
                                           add_controls = FALSE,
                                           only_controls = FALSE,
                                           star_controls = FALSE,
                                           controls = NULL,
                                           no_beh_control = TRUE,
                                           beh_control_int = FALSE, 
                                           simple_beh_control = TRUE,
                                           add_tfr = TRUE,
                                           only_tfr = FALSE,
                                           withBase = TRUE, 
                                           only_pz = TRUE,
                                           baseline_dep_var = TRUE,
                                           control_for_baseline = FALSE,
                                           ecg_control = FALSE,
                                           next_trial= FALSE,
                                           reverse_ibi_and_p3 = FALSE) {
  # ---- Define Dependent Variables ----
  # Base dependent variables
  if (only_pz){
    base_dep_vars <- c("Pz__p3b_Win", "Pz__p3b_Lose") #, "Pz__p3b"
  } else {
    base_dep_vars <- c("FCz__p3a", "Pz__p3b_Win", "Pz__p3b_Lose") #, "Pz__p3b" "FCz__frn",
  }
  
  # If add_tfr is TRUE, define additional TFR options
  if (add_tfr) {
    tfr_dep_vars <- c("Pz__high_alpha__baseline__TFR", "Pz__low_alpha__baseline__TFR")
  } else {
    tfr_dep_vars <- character(0)
  }
  
  # Decide on the final list based on the only_tfr flag
  if (only_tfr) {
    dep_vars <- tfr_dep_vars
  } else {
    dep_vars <- c(base_dep_vars, tfr_dep_vars)
  }
  
  # If withBase is TRUE, add additional options for each dep_var
  if (withBase) {
    dep_vars <- c(dep_vars, paste0(dep_vars, "__withBase"))
  }
  
  if (baseline_dep_var){
    if (only_pz){
      dep_vars <- c(dep_vars, "Pz__baseline")
    } else {
      dep_vars <- c(dep_vars, "Pz__baseline", "FCz__baseline")
    }
  }
  
  # ---- Define Main Versions ----
  main_versions <- c("4_3", "5_4", "6_5", "7_6")
  
  # ---- Define First Options ----
  if (simple_beh_control){
    first_options <- list(c("expectation", "outcome"))
  } else {
    first_options <- list(c("expectation", "*", "outcome"), c("absPE", "*", "outcome"))
  }
  if (no_beh_control) {
    first_options <- c(first_options, list(c()))
  }
  
  # Append "*" to non-empty first_options entries if beh_control_int is TRUE
  if (beh_control_int) {
    first_options <- lapply(first_options, function(x) {
      if (length(x) > 0) c(x, "*") else x
    })
  }
  
  # ---- Define Second Options ----
  second_options_with_star <- list(c("*", "phaseFeedb"))
  second_options_without_star <- list(c("phaseFeedb"))
  if (phase_factor) {
    second_options_with_star <- c(second_options_with_star, list(c("*", "sysDiaFeedb")))
    second_options_without_star <- c(second_options_without_star, list(c("sysDiaFeedb")))
  }
  if (phase_ms) {
    second_options_with_star <- c(second_options_with_star, list(c("*", "msDistToFeedb")))
    second_options_without_star <- c(second_options_without_star, list(c("msDistToFeedb")))
  }
  if (no_phase_at_all){
    second_options_with_star <- list(c())
    second_options_without_star <- list(c())
  }
  
  # Container for base formulas
  formulas <- list()
  
  # ---- Build Base Formulas ----
  for (dv in dep_vars) {
    for (v in main_versions) {
      base <- paste0("changeIBI", v)
      for (f in first_options) {
        if (!no_phase_int) {
          for (s in second_options_with_star) {
            formulas[[length(formulas) + 1]] <- c(dv, f, base, s)
          }
        }
        if (!only_phase_int) {
          for (s in second_options_without_star) {
            formulas[[length(formulas) + 1]] <- c(dv, f, base, s)
          }
        }
      }
    }
  }
  
  # ---- Add Controls if Specified ----
  formulas_with_controls <- list()
  if (add_controls && !is.null(controls) && length(controls) > 0) {
    for (formula_vec in formulas) {
      for (control in controls) {
        if (star_controls) {
          new_formula <- c(formula_vec, "*", control)
        } else {
          new_formula <- c(formula_vec, control)
        }
        formulas_with_controls[[length(formulas_with_controls) + 1]] <- new_formula
      }
    }
  }
  
  # ---- Final Output Selection ----
  if (only_controls) {
    final_formulas <- formulas_with_controls
  } else {
    final_formulas <- c(formulas, formulas_with_controls)
  }
  
  # ---- Add Baseline Control ----
  if (control_for_baseline) {
    final_formulas <- lapply(final_formulas, function(vec) {
      # skip if dep var already a baseline var
      if (vec[1] %in% c("Pz__baseline", "FCz__baseline")) {
        return(vec)
      }
      if (next_trial) {
        baseline_var <- ifelse(startsWith(vec[1], "Pz"), "Pz__baseline_next_time", "FCz__baseline_next_time")
      } else {
        baseline_var <- ifelse(startsWith(vec[1], "Pz"), "Pz__baseline", "FCz__baseline")
      }
      c(vec, baseline_var)
    })
  }
  
  # ---- Add ECG Control ----
  if (ecg_control) {
    final_formulas <- lapply(final_formulas, function(vec) {
      parts <- unlist(strsplit(vec[1], "__"))
      parts <- parts[!parts %in% c("TFR", "low_alpha", "high_alpha", "withBase")]
      second_part <- parts[2]
      if (next_trial) {
        ekg_var <- paste0("EKG1__", second_part, '_next_time')
      } else {
        ekg_var <- paste0("EKG1__", second_part)
      }
      c(vec, ekg_var)
    })
  }
  
  # ---- Add next trial ----
  if (next_trial) {
    final_formulas <- lapply(final_formulas, function(vec) {
      vec[1] <- paste0(vec[1], "_next_time")
      vec
    })
  }
  
  # ---- Reverse IBI and DV ----
  if (reverse_ibi_and_p3) {
    final_formulas <- lapply(final_formulas, function(vec) {
      base_idx <- which(grepl("^changeIBI", vec))
      if (length(base_idx) == 1 && base_idx != 1) {
        vec[c(1, base_idx)] <- vec[c(base_idx, 1)]
      }
      vec
    })
  }
  
  # ---- Clean Up and Return ----
  final_formulas <- lapply(final_formulas, function(vec) {
    sapply(vec, function(x) gsub("__", "_", x), USE.NAMES = FALSE)
  })
  final_formulas <- lapply(final_formulas, function(vec) unlist(lapply(vec, split_star_space)))
  return(final_formulas)
}

create_formula <- function(df, var_list, random_slope = FALSE, simple_model = FALSE, 
                           out_filename = FALSE, additional_filename_string = "") {
  # Check that at least one variable is provided.
  if (length(var_list) < 1) {
    stop("var_list must contain at least one variable.")
  }
  
  # Create a mapping for each non-operator token:
  # if a variable is continuous (i.e. not a factor or a binary 0/1 variable) then 
  # prepend "scale_" to indicate that we use a scaled version.
  var_mapping <- list()
  for (tok in var_list) {
    if (tok == "*") next  # Skip operator tokens
    if (!tok %in% names(df)) {
      stop(paste("Variable", tok, "not found in dataframe."))
    }
    if (!is.null(var_mapping[[tok]])) next  # Already processed
    if (length(var_list) == 1) {
      var_mapping[[tok]] <- tok
    } else if (is.factor(df[[tok]])) {
      var_mapping[[tok]] <- tok
    } else if (all(unique(df[[tok]]) %in% c(0, 1))) {
      warning(paste("Variable", tok, "appears to be binary. Converting to factor."))
      df[[tok]] <- as.factor(df[[tok]])
      var_mapping[[tok]] <- tok
    } else {
      # For continuous variables, expect a scaled version in df.
      scaled_var <- paste0("scale_", tok)
      if (!(scaled_var %in% names(df))) {
        warning(paste("Scaled variable", scaled_var, "not found in dataframe."))
      }
      var_mapping[[tok]] <- scaled_var
    }
  }
  
  # Build the formula string.
  # Assume the first element is the outcome variable.
  outcome_orig <- var_list[1]
  outcome <- var_mapping[[outcome_orig]]
  
  # Process the rest of the tokens to build the predictor part.
  predictor_tokens <- var_list[-1]
  terms <- c()
  current_term <- ""
  prev_was_star <- FALSE
  for (tok in predictor_tokens) {
    if (tok == "*") {
      prev_was_star <- TRUE
      next
    }
    if (prev_was_star) {
      # Append interaction to the current term.
      current_term <- paste0(current_term, "*", var_mapping[[tok]])
      prev_was_star <- FALSE
    } else {
      # Finish the previous term if any, and start a new one.
      if (current_term != "") {
        terms <- c(terms, current_term)
      }
      current_term <- var_mapping[[tok]]
    }
  }
  if (current_term != "") {
    terms <- c(terms, current_term)
  }
  
  # If no predictor is provided, default to 1 (one-sample t-test analogue).
  predictors_str <- if (length(terms) == 0) "1" else paste(terms, collapse = " + ")
  
  # Build the fixed-effects part.
  formula_str <- paste0(outcome, " ~ ", predictors_str)
  
  # Append random effects structure if not a simple model.
  if (!simple_model) {
    if (random_slope) {
      # Identify the unique predictors (ignoring duplicates and "*" tokens).
      unique_predictors <- unique(sapply(predictor_tokens[predictor_tokens != "*"],
                                         function(x) var_mapping[[x]]))
      if (length(unique_predictors) > 0 && !(length(unique_predictors) == 1 && unique_predictors == "1")) {
        random_str <- paste0(" (1 + ", paste(unique_predictors, collapse = " + "), " | subID)")
      } else {
        random_str <- " (1 | subID)"
      }
    } else {
      random_str <- " (1 | subID)"
    }
    formula_str <- paste0(formula_str, " +", random_str)
  }
  
  # If out_filename is TRUE, create a filename based on the original variable names.
  if (out_filename) {
    # Replace "*" with "int" for readability while preserving token order.
    filename_tokens <- sapply(var_list, function(x) if(x == "*") "int" else x)
    filename <- paste(filename_tokens, collapse = "_")
    
    # Append suffix: use "_simple" if simple_model is TRUE, or "_random_slope" if random_slope is TRUE.
    if (simple_model) {
      filename <- paste0(filename, "_simple")
    } else if (random_slope) {
      filename <- paste0(filename, "_random_slope")
    }
    
    # Append additional filename string if provided (with an underscore separator if non-empty).
    if (additional_filename_string != "") {
      filename <- paste0(filename, "_", additional_filename_string)
    }
    
    return(filename)
  } else {
    return(formula_str)
  }
}

prior_pred_check <- function(df, 
                             array_of_vars, 
                             result_path = ".",
                             family = NULL, 
                             prior = NULL, 
                             run_pp = TRUE,
                             default_prior = c(set_prior("normal(0,0.1)", class = "b")),
                             random_slope = FALSE, 
                             simple_model = FALSE,
                             vis_path = ".", 
                             plot_width = 8, 
                             plot_height = 6,
                             show_plots = FALSE, 
                             save_plots = TRUE,
                             plot_format = "svg", 
                             dpi = 150,
                             normal_intercept = FALSE, 
                             pick_phase = NULL,
                             additional_filename_string = ""
                             ) {
  
  if (!is.null(pick_phase) && pick_phase %in% c("Systole", "Diastole")) {
    df <- df[(df$sysDiaFeedb == pick_phase) & (!is.na(df$sysDiaFeedb)), ]
  }
  
  # 1. Build the formula using create_formula with out_filename = FALSE
  formula <- create_formula(df, 
                            array_of_vars, 
                            random_slope = random_slope, 
                            simple_model = simple_model, 
                            out_filename = FALSE)
  
  # Determine the output variable (assumed to be the first element)
  output_var <- array_of_vars[1]
  
  # If either family or prior is missing, warn and derive defaults.
  if (is.null(family) || is.null(prior)) {
    warning("Family and/or prior not provided. Using suggested family and default prior settings.")
    
    # If family is missing, use the suggestion function on the output variable.
    if (is.null(family)) {
      family <- suggest_brms_family(df[[output_var]])
    }
    
    # If prior is missing, choose based on output variable name.
    if (is.null(prior)) {
    if (length(array_of_vars) == 1) {
      if (normal_intercept) {
        prior <- c(set_prior(paste0("normal(0,", round(sd(df[[output_var]]) * 1.5, 2), ")"), class = "Intercept"))
      } else {
        prior <- NULL
      }
    } else {
      if (grepl("IBI", output_var, ignore.case = TRUE)) {
        prior <- c(set_prior("normal(0,0.05)", class = "b"))
      } else if (grepl("Pz|FCz", output_var, ignore.case = TRUE)) {
        prior <- c(set_prior("normal(0,0.1)", class = "b"))
      } else {
        warning("Output variable could not be mapped to a specific prior. Using default prior.")
        prior <- default_prior
      }
      }
    }
  }
  
  if (!is.null(pick_phase) && pick_phase %in% c("Systole", "Diastole") && (additional_filename_string != "")) {
    additional_filename_string <- paste0(additional_filename_string,"_df",pick_phase) # used to be: paste0("df",pick_phase), should be: paste0(additional_filename_string,"_",pick_phase)
  } else if (!is.null(pick_phase) && pick_phase %in% c("Systole", "Diastole")) {
    additional_filename_string <- paste0("df", pick_phase) # used to be: paste0("df", pick_phase), should be: pick_phase
  }
  
  # Use create_formula to generate a file prefix for saving the figure.
  file_prefix <- create_formula(df, 
                                array_of_vars, 
                                random_slope = random_slope, 
                                simple_model = simple_model, 
                                out_filename = TRUE,
                                additional_filename_string = additional_filename_string)
  
  file_prefix_path <- file.path(result_path, file_prefix)
  file_prefix_path <- handle_long_filename(file_prefix_path)
  file_prefix <- basename(file_prefix_path)
  
  # 2. If run_pp is TRUE, fit the model and produce prior predictive check plots.
  if (run_pp) {
    fit <- brm(formula = formula,
               data = df,
               family = family,
               warmup = 1000,
               sample_prior = "only",
               prior = prior,
               iter = 2000,
               chains = 4,
               cores = 2,
               seed = set_seed)
    
    # Generate various pp_check plots.
    p4 <- pp_check(fit, ndraws = 100) +
      xlab(map_variable_name(output_var)) + ggtitle('Prior predictive distribution')
    p1 <- pp_check(fit, type = 'stat', stat = 'min', ndraws = 100) +
      xlab(map_variable_name(output_var)) + ggtitle('Prior predictive distribution of minimum values')
    p2 <- pp_check(fit, type = 'stat', stat = 'mean', ndraws = 100) +
      xlab(map_variable_name(output_var)) + ggtitle('Prior predictive distribution of means')
    p3 <- pp_check(fit, type = 'stat', stat = 'max', ndraws = 100) +
      xlab(map_variable_name(output_var)) + ggtitle('Prior predictive distribution of maximum values')
    
    # Combine the individual plots.
    p <- combine_plots(list(p4, p1, p2, p3))
    
    # Save and/or display the combined plot.
    if (save_plots) {
      file_name <- file.path(vis_path, paste0(file_prefix, "_priorpc.", plot_format))
      ggsave(file_name, p, width = plot_width, height = plot_height, dpi = dpi)
    }
    if (show_plots) {
      print(p)
    }
  }
  
  # Return a list containing the brms setup (formula, family, and prior)
  return(list(formula = formula, family = family, prior = prior, filename = file_prefix))
}

run_brms_model <- function(df_to_model, brms_setup, statresult_path, force_refit = FALSE, iter = 4000,
                           run_pp = FALSE,
                           save_plots = TRUE,
                           show_plots = FALSE,
                           vis_path = statresult_path,
                           plot_format = "svg",
                           plot_width = 8,
                           plot_height = 6,
                           dpi = 150) {
  
  # Check if filename indicates phase filtering is needed
  pick_phase <- NULL
  if (grepl("_dfSystole$", brms_setup$filename)) {
    pick_phase <- "Systole"
  } else if (grepl("_dfDiastole$", brms_setup$filename)) {
    pick_phase <- "Diastole"
  }
  
  # but would it work with _simple or random_slope models?
  # later replace with correct option
  #if (grepl("_Systole$", brms_setup$filename)) {
  #  pick_phase <- "Systole"
  #} else if (grepl("_Diastole$", brms_setup$filename)) {
  #  pick_phase <- "Diastole"
  #}
  
  
  # Filter df based on phase if needed
  if (!is.null(pick_phase)) {
    df_to_model <- df_to_model[(df_to_model$sysDiaFeedb == pick_phase) & (!is.na(df_to_model$sysDiaFeedb)), ]
    warning(paste("Data was filtered for", pick_phase, "prior to model fitting."))
  }
  
  # Define file path for model storage
  model_file <- file.path(statresult_path, paste0(brms_setup$filename, ".rds"))
  
  # Set file_refit option based on force_refit
  refit_option <- if (force_refit) "always" else "never"
  
  fit <- brm(
    formula = brms_setup$formula,
    data = df_to_model,
    family = brms_setup$family,
    warmup = 2000,
    sample_prior = TRUE,
    prior = brms_setup$prior,
    iter = iter, # ideally should be 12 000 which yields 40 000 posterior samples
    chains = 4,
    cores = 2,
    seed = set_seed,
    control = list(adapt_delta = 0.95),
    file = model_file,
    save_pars = save_pars(all = TRUE),
    file_refit = refit_option
  )
  
  # If run_pp is TRUE, run posterior predictive checks and optionally save/show the figure
  ## here we need to improve xlab naming based on the map_var_name(first term in the formula)
  if (run_pp) {
    # Generate the individual pp_check plots
    p4 <- pp_check(fit, ndraws = 100) +
      xlab('changeHR') + ggtitle('Posterior predictive distribution')
    p1 <- pp_check(fit, type = 'stat', stat = 'min', ndraws = 100) +
      xlab('changeHR') + ggtitle('Posterior predictive distribution of minimum values')
    p2 <- pp_check(fit, type = 'stat', stat = 'mean', ndraws = 100) +
      xlab('changeHR') + ggtitle('Posterior predictive distribution of means')
    p3 <- pp_check(fit, type = 'stat', stat = 'max', ndraws = 100) +
      xlab('changeHR') + ggtitle('Posterior predictive distribution of maximum values')
    
    # Combine the individual plots into a single figure
    pp_plot <- combine_plots(list(p4, p1, p2, p3))
    
    # Save the combined plot if requested
    if (save_plots) {
      file_name <- file.path(vis_path, paste0(brms_setup$filename, "_postpc.", plot_format))
      ggsave(file_name, pp_plot, width = plot_width, height = plot_height, dpi = dpi)
    }
    # Display the plot if requested
    if (show_plots) {
      print(pp_plot)
    }
  }
  
  return(fit)
}

split_star <- function(s) {
  pos <- regexpr("\\*", s)  # find the position of "*"
  if (pos > 0) {
    # Create parts: before "*", the "*" itself, and after "*"
    parts <- c(substr(s, 1, pos - 1), "*", substr(s, pos + 1, nchar(s)))
    # Remove any empty strings from the result
    parts[parts != ""]
  } else {
    s  # return the string unchanged if "*" is not found
  }
}

split_space <- function(s) {
  pos <- regexpr(" ", s)  # find the position of the first space
  if (pos > 0) {
    # Create parts: before the space and after the space
    parts <- c(substr(s, 1, pos - 1), substr(s, pos + 1, nchar(s)))
    # Remove any empty strings from the result
    parts[parts != ""]
  } else {
    s  # Return the string unchanged if no space is found
  }
}

split_star_space <- function(s) {
  if (grepl("\\*", s)) {
    split_star(s)
  } else if (grepl(" ", s)) {
    split_space(s)
  } else {
    s  # Return unchanged if neither "*" nor " " is found
  }
}

remove_expression <- function(formulas, expr, remove_int = TRUE) {
  lapply(formulas, function(x) {
    # Create a logical vector for removal flags, initially all FALSE.
    remove <- rep(FALSE, length(x))
    # Loop through each token in the expression.
    for (i in seq_along(x)) {
      if (x[i] == expr) {
        remove[i] <- TRUE  # Mark the expression for removal.
        # If remove_int is TRUE and there's an asterisk preceding the expression, mark that too.
        if (remove_int && i > 1 && x[i-1] == "*") {
          remove[i-1] <- TRUE
        }
      }
    }
    # Return the expression vector with the marked elements removed.
    x[!remove]
  })
}

compute_subject_averages <- function(df, test_col, factor_col = NULL) {
  # Ensure that the factor columns are treated as factors.
  if (!is.null(factor_col) && length(factor_col) > 0) {
    df <- df %>% mutate(across(all_of(factor_col), as.factor))
  }
  
  ### 1. Overall Averages ###
  overall <- df %>%
    group_by(subID) %>%
    summarise(across(all_of(test_col), ~mean(.x, na.rm = TRUE)), .groups = "drop") #%>%
  #rename_with(~paste0(.x, "_overall"), all_of(test_col))
  
  result <- overall
  
  if (!is.null(factor_col) && length(factor_col) > 0) {
    ### 2. Main Effects ###
    for (f in factor_col) {
      group_avg <- df %>%
        group_by(subID, !!sym(f)) %>%
        summarise(across(all_of(test_col), ~mean(.x, na.rm = TRUE)), .groups = "drop")
      
      # Pivot so that each factor level becomes its own column.
      # Column names will be: <testvar>_<f>_<level>
      group_avg_wide <- group_avg %>%
        pivot_wider(
          names_from = !!sym(f),
          values_from = all_of(test_col),
          names_sep = paste0("_", f, "_"),
          values_fill = NA
        )
      
      # Merge into the result by subID (each subject appears only once).
      result <- left_join(result, group_avg_wide, by = "subID")
    }
    
    ### 3. Interaction Averages ###
    # Compute averages for each unique combination of factor levels.
    interaction_df <- df %>%
      group_by(subID, across(all_of(factor_col))) %>%
      summarise(across(all_of(test_col), ~mean(.x, na.rm = TRUE)), .groups = "drop")
    
    # Unite factor columns into one combined interaction label.
    # Remove the original factor columns so they don't duplicate.
    interaction_df <- interaction_df %>%
      unite("interaction", all_of(factor_col), sep = "_", remove = TRUE)
    
    # Pivot so that each unique interaction becomes its own set of columns.
    interaction_wide <- interaction_df %>%
      pivot_wider(
        names_from = interaction,
        values_from = all_of(test_col),
        names_sep = "_",
        values_fill = NA
      )
    
    # Merge the interaction averages.
    result <- left_join(result, interaction_wide, by = "subID")
    
  }
  
  return(result)
}

run_brms_models <- function(model_formulas, df, statresult_path, modelsummary_path = NULL,
                            prior_pred_args = list(run_pp = FALSE, random_slope = FALSE, result_path=statresult_path),
                            parallel = FALSE, force_refit = FALSE, force_info = FALSE, iter = 4000, rope_upper = NULL) {
  
  model_formulas <- lapply(model_formulas, function(vec) unlist(lapply(vec, split_star)))
  
  # Define an inner function to process a single model
  process_model <- function(i, formula) {
    formula_str <- paste(formula, collapse = " ")
    log_info("Processing model {i}: {formula_str}")
    
    # Combine default prior_pred_args with the formula and data
    brms_setup <- tryCatch({
      do.call(prior_pred_check, c(list(df = df, array_of_vars = formula), prior_pred_args))
    }, error = function(e) {
      log_error("Error in prior_pred_check for model {i}: {e$message}")
      return(NULL)
    })
    
    if (is.null(brms_setup)) return(NULL)
    
    # Determine the file name for the model (.rds file)
    model_file <- file.path(statresult_path, paste0(brms_setup$filename, ".rds"))
    info_file <- if (!is.null(modelsummary_path)) file.path(modelsummary_path, paste0(brms_setup$filename, ".rds")) else NULL 
    
    # no encoding and decoding for long filenames here?
    
    # Skip model fitting if the file already exists
    if (file.exists(model_file) && !force_refit) {
      if ((!is.null(info_file) && !file.exists(info_file))||(force_info)) {
        log_info("Model file {model_file} exists but info file {info_file} missing or force refit is required. Loading fit to compute info.")
        fit <- tryCatch({
          readRDS(model_file)
        }, error = function(e) {
          log_error("Error reading model file for model {i}: {e$message}")
          return(NULL)
        })
        if (is.null(fit)) return(NULL)
        # Compute the info and save it
        if (is.factor(df[[formula[1]]])) {
          info <- interpret_fit(fit, rope_upper = NULL, skip_rope_bf = TRUE)
        } else if (!is.null(rope_upper)){
          info <- interpret_fit(fit, rope_upper)
          #info$summary_df <- convert_summary_to_orig_units(summary_df = info$summary_df, fit = fit, df)
        } else if (length(formula)==1) {
          info <- interpret_fit(fit, rope_upper = NULL, skip_rope_bf = TRUE)
        } else {
          rope <- compute_rope_model_all_int(df, formula, assess_int = TRUE, include_intercept = TRUE)
          info <- interpret_fit(fit, rope)
          info$summary_df <- convert_summary_to_orig_units(summary_df = info$summary_df, fit = fit, df)
        }
        saveRDS(info, info_file)
        log_info("Successfully computed and saved info for model {i}: {info_file}")
        rm(fit)
        gc(verbose = FALSE)
      } else {
        log_info("Model file {model_file} already exists and no refit required. Skipping model {i}.")
      }
      return(NULL)
    }
    
    # Run the brms model
    fit <- tryCatch({
      run_brms_model(df_to_model = df, brms_setup = brms_setup, statresult_path = statresult_path, force_refit = force_refit, iter = iter)
    }, error = function(e) {
      log_error("Error in run_brms_model for model {i}: {e$message}")
      return(NULL)
    })
    
    # After obtaining fit, compute and save the interpretation info if modelsummary_path is provided
    if (!is.null(modelsummary_path)) {
      info_file <- file.path(modelsummary_path, paste0(brms_setup$filename, ".rds"))
      if (file.exists(info_file)) {
        log_info("Info file {info_file} already exists. Skipping interpretation for model {i}.")
      } else {
        if (is.factor(df[[formula[1]]])) {
          info <- interpret_fit(fit, rope_upper = NULL, skip_rope_bf = TRUE)
        } else if (!is.null(rope_upper)){
          info <- interpret_fit(fit, rope_upper)
          #info$summary_df <- convert_summary_to_orig_units(summary_df = info$summary_df, fit = fit, df)
        } else if (length(formula)==1) {
          info <- interpret_fit(fit, rope_upper = NULL, skip_rope_bf = TRUE)
        } else {
          rope <- compute_rope_model_all_int(df, formula, assess_int = TRUE, include_intercept = TRUE)
          info <- interpret_fit(fit, rope)
          info$summary_df <- convert_summary_to_orig_units(summary_df = info$summary_df, fit = fit, df)
        }
        saveRDS(info, info_file)
        log_info("Successfully computed and saved info for model {i}: {info_file}")
      }
    }
    
    log_info("Finished processing model {i}: {formula_str}")
    
    # Clean up fit from memory and force garbage collection
    rm(fit)
    gc(verbose = FALSE)
    return(NULL)
  }
  
  # Process models either in parallel or sequentially
  if (parallel) {
    library(future)
    library(future.apply)
    plan(multisession, workers = 2)
    future_lapply(seq_along(model_formulas), function(i) {
      process_model(i, model_formulas[[i]])
    })
  } else {
    for (i in seq_along(model_formulas)) {
      process_model(i, model_formulas[[i]])
    }
  }
}

load_fits <- function(model_formulas, df, statresult_path, modelsummary_path = NULL,
                      prior_pred_args = list(run_pp = FALSE, random_slope = FALSE, result_path=statresult_path)) {
  fits <- list()
  infos <- list()
  
  for (i in seq_along(model_formulas)) {
    formula <- model_formulas[[i]]
    formula_str <- paste(formula, collapse = " ")
    log_info("Loading model {i}: {formula_str}")
    
    # Use prior_pred_check to obtain brms_setup (which includes the filename)
    brms_setup <- tryCatch({
      do.call(prior_pred_check, c(list(df = df, array_of_vars = formula), prior_pred_args))
    }, error = function(e) {
      log_error("Error in prior_pred_check for model {i}: {e$message}")
      return(NULL)
    })
    
    if (is.null(brms_setup)) {
      fits[[i]] <- NULL
      infos[[i]] <- NULL
      next
    }
    
    # Construct the model file path based on the brms_setup filename
    model_file <- file.path(statresult_path, paste0(brms_setup$filename, ".rds"))
    if (!file.exists(model_file)) {
      log_error("File not found for model {i}: {model_file}")
      fits[[i]] <- NULL
      infos[[i]] <- NULL
      next
    }
    
    # Try to load the model fit from the file
    fit <- tryCatch({
      readRDS(model_file)
    }, error = function(e) {
      log_error("Error reading model file for model {i}: {e$message}")
      return(NULL)
    })
    
    fits[[i]] <- fit
    log_info("Successfully loaded model {i}: {model_file}")
    # If modelsummary_path is provided, try to load the pre-computed info file
    if (!is.null(modelsummary_path)) {
      info_file <- file.path(modelsummary_path, paste0(brms_setup$filename, ".rds"))
      if (file.exists(info_file)) {
        info <- tryCatch({
          readRDS(info_file)
        }, error = function(e) {
          log_error("Error reading info file for model {i}: {e$message}")
          return(NULL)
        })
        infos[[i]] <- info
      } else {
        infos[[i]] <- NULL
        log_info("No info file for model {i}: {info_file}")
      }
    }
  }
  
  if (!is.null(modelsummary_path) && !is_empty(infos)) {
    return(list(fits = fits, infos = infos)) 
  } else {
    return(fits)
  }
}

