# Run in RStudio
# Input:  data/data.csv
# Output: data/StatResults/*.rds and data/ModelSummaries/*.rds

# Setup ####
library(logger)
library(loo)

force_refit <- FALSE
force_info <- FALSE

code_path <- dirname(rstudioapi::getSourceEditorContext()$path)
data_path <- file.path(dirname(dirname(code_path)), "data")
data_filename <- file.path(data_path, "data.csv")
log_path <- file.path(data_path, "derivatives", "logs")
statresult_path <- file.path(data_path, "StatResults")
modelsummary_path <- file.path(data_path, "ModelSummaries")
chapter_path <- file.path(data_path, "Chapters")

dir.create(log_path, recursive = TRUE, showWarnings = FALSE)
dir.create(statresult_path, recursive = TRUE, showWarnings = FALSE)
dir.create(modelsummary_path, recursive = TRUE, showWarnings = FALSE)
dir.create(chapter_path, recursive = TRUE, showWarnings = FALSE)
log_appender(appender_file(file.path(log_path, "prerun.txt")))

source(file.path(code_path, "preprocess.R"))
source(file.path(code_path, "plotting_helpers.R"))
source(file.path(code_path, "analyse_column_distributions.R"))
source(file.path(code_path, "analysis_functions.R"))
source(file.path(code_path, "analyse_model_output_functions.R"))
source(file.path(code_path, "plot_model_results_functions.R"))
source(file.path(code_path, "power_simulation_functions.R"))
source(file.path(code_path, "neural_plotting_functions.R"))
source(file.path(code_path, "rpeak_plotting_functions.R"))
source(file.path(code_path, "average_plotting_functions.R"))

# Formula helper
add_control <- function(formulas, var_name, int = FALSE) {
  if (!is.list(formulas) || length(formulas) == 0) return(formulas)

  add_tokens <- if (length(var_name) == 1L) {
    c(var_name)
  } else {
    as.character(var_name)
  }

  join_token <- if (isTRUE(int)) "*" else NULL

  lapply(formulas, function(tokens) {
    c(tokens, join_token, add_tokens)
  })
}

# Behavioural/cardiac data ####
df <- read.csv(data_filename)
result <- preprocess_df(df, rescale_after_outliers = TRUE)
df <- result$df_filtered
rm(result)

# Four heartbeats x four behavioural predictors, continuous cardiac phase.
main_formulas <- list(
  c("changeIBI4_3", "expectation", "*", "phaseFeedb"),
  c("changeIBI5_4", "expectation", "*", "phaseFeedb"),
  c("changeIBI6_5", "expectation", "*", "phaseFeedb"),
  c("changeIBI7_6", "expectation", "*", "phaseFeedb"),
  c("changeIBI4_3", "outcome", "*", "phaseFeedb"),
  c("changeIBI5_4", "outcome", "*", "phaseFeedb"),
  c("changeIBI6_5", "outcome", "*", "phaseFeedb"),
  c("changeIBI7_6", "outcome", "*", "phaseFeedb"),
  c("changeIBI4_3", "PE", "*", "phaseFeedb"),
  c("changeIBI5_4", "PE", "*", "phaseFeedb"),
  c("changeIBI6_5", "PE", "*", "phaseFeedb"),
  c("changeIBI7_6", "PE", "*", "phaseFeedb"),
  c("changeIBI4_3", "absPE", "*", "phaseFeedb"),
  c("changeIBI5_4", "absPE", "*", "phaseFeedb"),
  c("changeIBI6_5", "absPE", "*", "phaseFeedb"),
  c("changeIBI7_6", "absPE", "*", "phaseFeedb")
)

group_results <- run_model_pipeline(
  main_formulas,
  df,
  statresult_path,
  modelsummary_path,
  prior_pred_args = list(run_pp = FALSE, random_slope = FALSE),
  force_refit = force_refit,
  force_info = force_info
)
create_word_from_groups(group_results, chapter_path, "formulas_main")
rm(group_results)

# R2: systole-only tests ####
main_systole_formulas <- list(
  c("changeIBI6_5", "outcome"),
  c("changeIBI6_5", "PE"),
  c("changeIBI6_5", "absPE")
)

group_results <- run_model_pipeline(
  main_systole_formulas,
  df,
  statresult_path,
  modelsummary_path,
  prior_pred_args = list(run_pp = FALSE, pick_phase = "Systole", random_slope = FALSE),
  force_refit = force_refit,
  force_info = force_info
)
create_word_from_groups(group_results, chapter_path, "main_systole_formulas")
rm(group_results, df)

# Neural/P3b data: identical cleaning and renaming to original ####
df <- read.csv(data_filename)
result <- preprocess_df(df, rescale_after_outliers = TRUE, clean_neural = TRUE)
df <- result$df_filtered
rm(result)
df <- replace_double_underscore_in_colnames(df)

main_systole_formulas_neural <- list(
  c("changeIBI5_4", "Pz_p3b_Win"),
  c("changeIBI5_4", "Pz_p3b_Lose")
)

main_diastole_formulas_neural <- list(
  c("changeIBI6_5", "Pz_p3b_Win"),
  c("changeIBI6_5", "Pz_p3b_Lose")
)

group_results <- run_model_pipeline(
  main_systole_formulas_neural,
  df,
  statresult_path,
  modelsummary_path,
  prior_pred_args = list(run_pp = FALSE, pick_phase = "Systole", random_slope = FALSE),
  force_refit = force_refit,
  force_info = force_info
)
create_word_from_groups(group_results, chapter_path, "formulas_erps_main_systole")
rm(group_results)

group_results <- run_model_pipeline(
  main_diastole_formulas_neural,
  df,
  statresult_path,
  modelsummary_path,
  prior_pred_args = list(run_pp = FALSE, pick_phase = "Diastole", random_slope = FALSE),
  force_refit = force_refit,
  force_info = force_info
)
create_word_from_groups(group_results, chapter_path, "formulas_erps_main_diastole")
rm(group_results)

message("Main cardiac and P3b models processed.")

# Control analyses ####

# Cardiac controls ####
df <- read.csv(data_filename)
result <- preprocess_df(df, rescale_after_outliers = TRUE)
df <- result$df_filtered
rm(result)

# Main effects of interest
overall_formulas <- list(
  c("changeIBI4_3", "expectation", "*", "phaseFeedb"),
  c("changeIBI7_6", "expectation", "*", "phaseFeedb"),
  c("changeIBI6_5", "outcome", "*", "phaseFeedb"),
  c("changeIBI6_5", "PE", "*", "phaseFeedb"),
  c("changeIBI6_5", "absPE", "*", "phaseFeedb")
)

# Approach versus avoid
control_formulas_approach <- add_control(
  overall_formulas, var_name = c("choiceApproach", "*", "phaseFeedb"), int = FALSE
)
group_results <- run_model_pipeline(
  control_formulas_approach, df, statresult_path, modelsummary_path,
  prior_pred_args = list(run_pp = FALSE), force_info = force_info, force_refit = force_refit
)
create_word_from_groups(group_results, chapter_path, "control_formulas_approach")

# Trial number within block
control_formulas_numtrial <- add_control(
  overall_formulas, var_name = c("trialInBlock", "*", "phaseFeedb"), int = FALSE
)
group_results <- run_model_pipeline(
  control_formulas_numtrial, df, statresult_path, modelsummary_path,
  prior_pred_args = list(run_pp = FALSE), force_info = force_info, force_refit = force_refit
)
create_word_from_groups(group_results, chapter_path, "control_formulas_numtrial")

# Resting HR and HRV
control_formulas_restingHR <- add_control(
  overall_formulas, var_name = c("RestAveIBI", "*", "phaseFeedb"), int = FALSE
)
control_formulas_restingHRV <- add_control(
  overall_formulas, var_name = c("RestRMSSD", "*", "phaseFeedb"), int = FALSE
)
control_formulas_resting <- c(control_formulas_restingHR, control_formulas_restingHRV)
group_results <- run_model_pipeline(
  control_formulas_resting, df, statresult_path, modelsummary_path,
  prior_pred_args = list(run_pp = FALSE), force_info = force_info, force_refit = force_refit
)
create_word_from_groups(group_results, chapter_path, "control_formulas_resting")

# Average IBI and RMSSD within trial, before/after feedback
control_formulas_beforeHR <- add_control(
  overall_formulas, var_name = c("beforeFeedbAveIBI", "*", "phaseFeedb"), int = FALSE
)
control_formulas_beforeHRV <- add_control(
  overall_formulas, var_name = c("beforeFeedbRMSSD", "*", "phaseFeedb"), int = FALSE
)
control_formulas_afterHR <- add_control(
  overall_formulas, var_name = c("afterFeedbAveIBI", "*", "phaseFeedb"), int = FALSE
)
control_formulas_afterHRV <- add_control(
  overall_formulas, var_name = c("afterFeedbRMSSD", "*", "phaseFeedb"), int = FALSE
)
control_formulas_withinTrial <- c(
  control_formulas_beforeHR, control_formulas_beforeHRV,
  control_formulas_afterHR, control_formulas_afterHRV
)
group_results <- run_model_pipeline(
  control_formulas_withinTrial, df, statresult_path, modelsummary_path,
  prior_pred_args = list(run_pp = FALSE), force_info = force_info, force_refit = force_refit
)
create_word_from_groups(group_results, chapter_path, "control_formulas_withinTrial")

# Prior sensitivity
prior1 <- c(set_prior("normal(0,0.1)", class = "b"))
prior2 <- c(set_prior("normal(0,0.03)", class = "b"))
prior3 <- c(set_prior("normal(0,1)", class = "b"))

group_results <- run_model_pipeline(
  overall_formulas, df, statresult_path, modelsummary_path,
  prior_pred_args = list(run_pp = FALSE, prior = prior1,
                         additional_filename_string = "prior1_01"),
  force_info = force_info, force_refit = force_refit
)
create_word_from_groups(group_results, chapter_path, "control_prior_sensitivity1_01")
group_results <- run_model_pipeline(
  overall_formulas, df, statresult_path, modelsummary_path,
  prior_pred_args = list(run_pp = FALSE, prior = prior2,
                         additional_filename_string = "prior2_003"),
  force_info = force_info, force_refit = force_refit
)
create_word_from_groups(group_results, chapter_path, "control_prior_sensitivity2_003")
group_results <- run_model_pipeline(
  overall_formulas, df, statresult_path, modelsummary_path,
  prior_pred_args = list(run_pp = FALSE, prior = prior3,
                         additional_filename_string = "prior3_1"),
  force_info = force_info, force_refit = force_refit
)
create_word_from_groups(group_results, chapter_path, "control_prior_sensitivity3_1")

# Reaction-time control
RT_df <- df[!is.na(df$RT), 1:410]
control_formulas_RT <- add_control(
  overall_formulas, var_name = c("RT", "*", "phaseFeedb"), int = FALSE
)
group_results <- run_model_pipeline(
  control_formulas_RT, RT_df, statresult_path, modelsummary_path,
  prior_pred_args = list(run_pp = FALSE),
  force_info = force_info, force_refit = force_refit
)
create_word_from_groups(group_results, chapter_path, "control_formulas_RT")

# Separate approach- and avoid-only fits
RT_df <- df[df$choiceApproach == "Approach", 1:410]
group_results <- run_model_pipeline(
  overall_formulas, RT_df, statresult_path, modelsummary_path,
  prior_pred_args = list(run_pp = FALSE, additional_filename_string = "onlyApproach"),
  force_info = force_info, force_refit = force_refit
)
create_word_from_groups(group_results, chapter_path, "control_formulas_Approach_only")
RT_df <- df[df$choiceApproach == "Avoid", 1:410]
group_results <- run_model_pipeline(
  overall_formulas, RT_df, statresult_path, modelsummary_path,
  prior_pred_args = list(run_pp = FALSE, additional_filename_string = "onlyAvoid"),
  force_info = force_info, force_refit = force_refit
)
create_word_from_groups(group_results, chapter_path, "control_formulas_Avoid_only")

# R2 during diastole, sensitivity: the equivalent systole models already run above
group_results <- run_model_pipeline(
  main_systole_formulas, df, statresult_path, modelsummary_path,
  prior_pred_args = list(run_pp = FALSE, pick_phase = "Diastole"),
  force_info = force_info, force_refit = force_refit
)
create_word_from_groups(group_results, chapter_path, "control_formulas_main_Diastole")

# Strict HRV/heart-rate outlier removal
df_strict <- read.csv(data_filename)
result <- preprocess_df(
  df_strict, rescale_after_outliers = TRUE, very_strict_cardiac_clean = TRUE
)
df_strict <- result$df_filtered
rm(result)
group_results <- run_model_pipeline(
  overall_formulas, df_strict, statresult_path, modelsummary_path,
  prior_pred_args = list(run_pp = FALSE, additional_filename_string = "very_strict"),
  force_info = force_info, force_refit = force_refit
)
create_word_from_groups(group_results, chapter_path, "control_formulas_main_very_strict")
rm(df_strict, RT_df, group_results)

# Uncorrected (pure) IBI models perform worse ####
# Set FALSE to skip these extensive secondary fits.
run_pure_IBI_controls <- TRUE
if (run_pure_IBI_controls) {
  formulas_pure_with_prev <- generate_model_formulas(
    main = FALSE, other_r = TRUE, pure_IBI = TRUE,
    pure_IBI_prev = TRUE, int_prev_ibi = TRUE
  )
  formulas_pure_without_prev <- generate_model_formulas(
    main = FALSE, other_r = TRUE, pure_IBI = TRUE,
    pure_IBI_prev = FALSE, int_prev_ibi = FALSE
  )

  has_subseq <- function(tokens, pattern) {
    n <- length(tokens); m <- length(pattern)
    if (n < m) return(FALSE)
    for (i in 1:(n - m + 1)) {
      if (all(tokens[i:(i + m - 1)] == pattern)) return(TRUE)
    }
    FALSE
  }

  should_remove <- function(tokens) {
    bad1 <- "sysDiaFeedb" %in% tokens
    phase_idx <- which(tokens == "phaseFeedb")
    bad2 <- length(phase_idx) > 0 && any(
      phase_idx == 1L | tokens[phase_idx - 1L] != "*"
    )
    bad3 <- has_subseq(tokens, c("expectation", "*", "outcome")) ||
      has_subseq(tokens, c("absPE", "*", "outcome"))
    bad1 || bad2 || bad3
  }

  formulas_pure_with_prev <- formulas_pure_with_prev[
    !vapply(formulas_pure_with_prev, should_remove, logical(1))
  ]
  formulas_pure_without_prev <- formulas_pure_without_prev[
    !vapply(formulas_pure_without_prev, should_remove, logical(1))
  ]
  formulas_pure <- c(formulas_pure_with_prev, formulas_pure_without_prev)

  group_results <- run_model_pipeline(
    formulas_pure, df, statresult_path, modelsummary_path,
    prior_pred_args = list(run_pp = FALSE),
    ignore_last = TRUE, force_info = force_info, force_refit = force_refit
  )
  create_word_from_groups(group_results, chapter_path, "control_formulas_pure")
  rm(group_results)
}
rm(df)

# Neural/P3b controls ####

# Apply the original neural cleaning before every control type.
df <- read.csv(data_filename)
result <- preprocess_df(df, rescale_after_outliers = TRUE, clean_neural = TRUE)
df <- result$df_filtered
rm(result)
df <- replace_double_underscore_in_colnames(df)

# The main neural formulas have already been defined and run above:
# main_systole_formulas_neural and main_diastole_formulas_neural.
# Here as a sensitivity check we run the same two neural outcomes in each phase, including opposite-phase checks.
sysdia_formulas <- c(main_systole_formulas_neural, main_diastole_formulas_neural)
for (phase in levels(df$sysDiaFeedb)) {
  group_results <- run_model_pipeline(
    sysdia_formulas, df, statresult_path, modelsummary_path,
    prior_pred_args = list(run_pp = FALSE, pick_phase = phase),
    force_info = force_info, force_refit = force_refit
  )
  create_word_from_groups(group_results, chapter_path, paste0("control_neuro_formulas_", phase))
}

# Prior sensitivity
for (phase in c("Systole", "Diastole")) {
  formulas <- if (phase == "Systole") {
    main_systole_formulas_neural
  } else {
    main_diastole_formulas_neural
  }

  group_results <- run_model_pipeline(
    formulas, df, statresult_path, modelsummary_path,
    prior_pred_args = list(run_pp = FALSE, prior = prior1,
                           additional_filename_string = "prior1_01", pick_phase = phase),
    force_info = force_info, force_refit = force_refit
  )
  create_word_from_groups(group_results, chapter_path, paste0("control_neuro_prior_sensitivity1_01_", tolower(phase)))
  group_results <- run_model_pipeline(
    formulas, df, statresult_path, modelsummary_path,
    prior_pred_args = list(run_pp = FALSE, prior = prior2,
                           additional_filename_string = "prior2_003", pick_phase = phase),
    force_info = force_info, force_refit = force_refit
  )
  create_word_from_groups(group_results, chapter_path, paste0("control_neuro_prior_sensitivity2_003_", tolower(phase)))
  group_results <- run_model_pipeline(
    formulas, df, statresult_path, modelsummary_path,
    prior_pred_args = list(run_pp = FALSE, prior = prior3,
                           additional_filename_string = "prior3_1", pick_phase = phase),
    force_info = force_info, force_refit = force_refit
  )
  create_word_from_groups(group_results, chapter_path, paste0("control_neuro_prior_sensitivity3_1_", tolower(phase)))
}

# Neural controls 
var_names <- c(
  "outcome", "PE", "expectation", "RestAveIBI", "RestRMSSD",
  "trialInBlock", "choiceApproach", "Pz_baseline",
  "beforeFeedbAveIBI", "beforeFeedbRMSSD",
  "afterFeedbAveIBI", "afterFeedbRMSSD"
)
for (var_name in var_names) {
  word_suffix <- if (var_name %in% c("beforeFeedbAveIBI", "beforeFeedbRMSSD",
                                     "afterFeedbAveIBI", "afterFeedbRMSSD")) "" else "_noint"
  control_formulas_systole <- add_control(
    main_systole_formulas_neural, var_name = c(var_name), int = FALSE
  )
  group_results <- run_model_pipeline(
    control_formulas_systole, df, statresult_path, modelsummary_path,
    prior_pred_args = list(run_pp = FALSE, pick_phase = "Systole"),
    force_info = force_info, force_refit = force_refit
  )
  create_word_from_groups(group_results, chapter_path, paste0("control_neuro_formulas_", var_name, word_suffix, "_Systole"))
  control_formulas_diastole <- add_control(
    main_diastole_formulas_neural, var_name = c(var_name), int = FALSE
  )
  group_results <- run_model_pipeline(
    control_formulas_diastole, df, statresult_path, modelsummary_path,
    prior_pred_args = list(run_pp = FALSE, pick_phase = "Diastole"),
    force_info = force_info, force_refit = force_refit
  )
  create_word_from_groups(group_results, chapter_path, paste0("control_neuro_formulas_", var_name, word_suffix, "_Diastole"))
}

# ECG amplitude in the neural P3b time window
control_systole_formulas <- list(
  c("changeIBI5_4", "Pz_p3b_Win", "EKG1_p3b_Win"),
  c("changeIBI5_4", "Pz_p3b_Lose", "EKG1_p3b_Lose")
)
control_diastole_formulas <- list(
  c("changeIBI6_5", "Pz_p3b_Win", "EKG1_p3b_Win"),
  c("changeIBI6_5", "Pz_p3b_Lose", "EKG1_p3b_Lose")
)
group_results <- run_model_pipeline(
  control_systole_formulas, df, statresult_path, modelsummary_path,
  prior_pred_args = list(run_pp = FALSE, pick_phase = "Systole"),
  force_info = force_info, force_refit = force_refit
)
create_word_from_groups(group_results, chapter_path, paste0("control_neuro_formulas_", "EKG", "_Systole"))
group_results <- run_model_pipeline(
  control_diastole_formulas, df, statresult_path, modelsummary_path,
  prior_pred_args = list(run_pp = FALSE, pick_phase = "Diastole"),
  force_info = force_info, force_refit = force_refit
)
create_word_from_groups(group_results, chapter_path, paste0("control_neuro_formulas_", "EKG", "_Diastole"))

# Neural RT control
RT_df <- df[!is.na(df$RT), 1:410]
control_formulas_systole <- add_control(
  main_systole_formulas_neural, var_name = c("RT"), int = FALSE
)
group_results <- run_model_pipeline(
  control_formulas_systole, RT_df, statresult_path, modelsummary_path,
  prior_pred_args = list(run_pp = FALSE, pick_phase = "Systole"),
  force_info = force_info, force_refit = force_refit
)
create_word_from_groups(group_results, chapter_path, "control_neuro_formulas_RT_Systole")
control_formulas_diastole <- add_control(
  main_diastole_formulas_neural, var_name = c("RT"), int = FALSE
)
group_results <- run_model_pipeline(
  control_formulas_diastole, RT_df, statresult_path, modelsummary_path,
  prior_pred_args = list(run_pp = FALSE, pick_phase = "Diastole"),
  force_info = force_info, force_refit = force_refit
)
create_word_from_groups(group_results, chapter_path, "control_neuro_formulas_RT_Diastole")

# Approach-only and avoid-only neural models in systole and diastole
RT_df <- df[df$choiceApproach == "Approach", 1:410]
group_results <- run_model_pipeline(
  main_systole_formulas_neural, RT_df, statresult_path, modelsummary_path,
  prior_pred_args = list(run_pp = FALSE, additional_filename_string = "onlyApproach",
                         pick_phase = "Systole"),
  force_info = force_info, force_refit = force_refit
)
create_word_from_groups(group_results, chapter_path, "control_neuro_Approach_only_Systole")
group_results <- run_model_pipeline(
  main_diastole_formulas_neural, RT_df, statresult_path, modelsummary_path,
  prior_pred_args = list(run_pp = FALSE, additional_filename_string = "onlyApproach",
                         pick_phase = "Diastole"),
  force_info = force_info, force_refit = force_refit
)
create_word_from_groups(group_results, chapter_path, "control_neuro_Approach_only_Diastole")
RT_df <- df[df$choiceApproach == "Avoid", 1:410]
group_results <- run_model_pipeline(
  main_systole_formulas_neural, RT_df, statresult_path, modelsummary_path,
  prior_pred_args = list(run_pp = FALSE, additional_filename_string = "onlyAvoid",
                         pick_phase = "Systole"),
  force_info = force_info, force_refit = force_refit
)
create_word_from_groups(group_results, chapter_path, "control_neuro_Avoid_only_Systole")
group_results <- run_model_pipeline(
  main_diastole_formulas_neural, RT_df, statresult_path, modelsummary_path,
  prior_pred_args = list(run_pp = FALSE, additional_filename_string = "onlyAvoid",
                         pick_phase = "Diastole"),
  force_info = force_info, force_refit = force_refit
)
create_word_from_groups(group_results, chapter_path, "control_neuro_Avoid_only_Diastole")

# Strict HRV/outlier cleaning of neural data
df_strict <- read.csv(data_filename)
result <- preprocess_df(
  df_strict, rescale_after_outliers = TRUE,
  clean_neural = TRUE, very_strict_cardiac_clean = TRUE
)
df_strict <- result$df_filtered
rm(result)
df_strict <- replace_double_underscore_in_colnames(df_strict)
group_results <- run_model_pipeline(
  main_systole_formulas_neural, df_strict, statresult_path, modelsummary_path,
  prior_pred_args = list(run_pp = FALSE, additional_filename_string = "very_strict",
                         pick_phase = "Systole"),
  force_info = force_info, force_refit = force_refit
)
create_word_from_groups(group_results, chapter_path, "control_neuro_very_strict_systole")
group_results <- run_model_pipeline(
  main_diastole_formulas_neural, df_strict, statresult_path, modelsummary_path,
  prior_pred_args = list(run_pp = FALSE, additional_filename_string = "very_strict",
                         pick_phase = "Diastole"),
  force_info = force_info, force_refit = force_refit
)
create_word_from_groups(group_results, chapter_path, "control_neuro_very_strict_diastole")
rm(group_results, RT_df, df_strict, df)

message("Cardiac and neural control models processed.")

# Main-model random-intercept vs random-slope comparisons ####

run_random_slope_comparisons <- TRUE
if (run_random_slope_comparisons) {
  
  df_cardiac <- read.csv(data_filename)
  prep <- preprocess_df(df_cardiac, rescale_after_outliers = TRUE)
  df_cardiac <- prep$df_filtered
  rm(prep)

  slope_overall_formulas <- overall_formulas
  slope_systole_formulas <- main_systole_formulas

  group_results <- run_model_pipeline(
    slope_overall_formulas, df_cardiac, statresult_path, modelsummary_path,
    prior_pred_args = list(run_pp = FALSE, random_slope = TRUE),
    force_refit = force_refit, force_info = force_info
  )
  create_word_from_groups(group_results, chapter_path, "main_overall_formulas_try_random_slopes")
  rm(group_results)
  group_results <- run_model_pipeline(
    slope_systole_formulas, df_cardiac, statresult_path, modelsummary_path,
    prior_pred_args = list(run_pp = FALSE, pick_phase = "Systole", random_slope = TRUE),
    force_refit = force_refit, force_info = force_info
  )
  create_word_from_groups(group_results, chapter_path, "main_systole_formulas_try_random_slopes")
  rm(group_results)

  df_neural <- read.csv(data_filename)
  prep <- preprocess_df(df_neural, rescale_after_outliers = TRUE, clean_neural = TRUE)
  df_neural <- replace_double_underscore_in_colnames(prep$df_filtered)
  rm(prep)

  group_results <- run_model_pipeline(
    main_systole_formulas_neural, df_neural, statresult_path, modelsummary_path,
    prior_pred_args = list(run_pp = FALSE, pick_phase = "Systole", random_slope = TRUE),
    force_refit = force_refit, force_info = force_info
  )
  create_word_from_groups(group_results, chapter_path, "main_systole_formulas_neural_random_slopes")
  rm(group_results)
  group_results <- run_model_pipeline(
    main_diastole_formulas_neural, df_neural, statresult_path, modelsummary_path,
    prior_pred_args = list(run_pp = FALSE, pick_phase = "Diastole", random_slope = TRUE),
    force_refit = force_refit, force_info = force_info
  )
  create_word_from_groups(group_results, chapter_path, "main_diastole_formulas_neural_random_slopes")
  rm(group_results)

  # Read the saved leave-one-out results
  # Intercept-only fits were already run above and are not refitted 
  get_loo_from_rds <- function(df, formula_vec, random_slope, pick_phase = NULL) {
    setup <- prior_pred_check(
      df, formula_vec, run_pp = FALSE,
      random_slope = random_slope, pick_phase = pick_phase
    )
    info_file <- file.path(modelsummary_path, paste0(setup$filename, ".rds"))
    if (!file.exists(info_file)) stop("Missing model summary: ", info_file)
    info <- readRDS(info_file)
    if (is.null(info$loo_fit)) stop("Missing loo_fit in: ", info_file)
    info$loo_fit
  }

  # 3 * SE threshold
  compare_slope_pair <- function(df, formula_vec, label, pick_phase = NULL) {
    loo_with <- get_loo_from_rds(df, formula_vec, TRUE, pick_phase)
    loo_without <- get_loo_from_rds(df, formula_vec, FALSE, pick_phase)

    looic_with <- loo_with$estimates["looic", "Estimate"]
    looic_without <- loo_without$estimates["looic", "Estimate"]
    p_loo_with <- loo_with$estimates["p_loo", "Estimate"]
    p_loo_without <- loo_without$estimates["p_loo", "Estimate"]
    delta_looic <- looic_without - looic_with

    elpd_pw_with <- loo_with$pointwise[, "elpd_loo"]
    elpd_pw_without <- loo_without$pointwise[, "elpd_loo"]
    stopifnot(length(elpd_pw_with) == length(elpd_pw_without))
    pw_diff <- elpd_pw_with - elpd_pw_without
    elpd_diff <- sum(pw_diff)
    se_diff <- sqrt(length(pw_diff)) * stats::sd(pw_diff)
    k <- 3
    conclusion <- if (abs(elpd_diff) < k * se_diff) {
      "No meaningful difference"
    } else if (elpd_diff > 0) {
      "Random slopes model preferred"
    } else {
      "Intercept-only model preferred"
    }

    data.frame(
      Model_Pair = label,
      LOOIC_with_slopes = round(looic_with, 2),
      LOOIC_without_slopes = round(looic_without, 2),
      p_LOO_with = round(p_loo_with, 2),
      p_LOO_without = round(p_loo_without, 2),
      Delta_LOOIC_without_minus_with = round(delta_looic, 2),
      ELPD_diff_with_minus_without = round(elpd_diff, 2),
      SE_diff = round(se_diff, 2),
      Conclusion = conclusion,
      check.names = FALSE
    )
  }

  # This includes optional exploratory main-model variants, but no controls
  cardiac_overall_names <- c(
    "Overall_R-1_expectation", "Overall_R3_expectation",
    "Overall_R2_outcome", "Overall_R2_PE", "Overall_R2_absPE"
  )
  cardiac_sys_names <- c(
    "Systole_R2_outcome", "Systole_R2_PE", "Systole_R2_absPE"
  )
  neural_sys_names <- c("Systole_R1_P3b_Win", "Systole_R1_P3b_Lose")
  neural_dia_names <- c("Diastole_R2_P3b_Win", "Diastole_R2_P3b_Lose")

  compare_group <- function(df, formulas, labels, pick_phase = NULL) {
    do.call(rbind, lapply(seq_along(formulas), function(i) {
      compare_slope_pair(df, formulas[[i]], labels[[i]], pick_phase)
    }))
  }

  comparison_table <- rbind(
    compare_group(df_cardiac, slope_overall_formulas, cardiac_overall_names),
    compare_group(df_cardiac, slope_systole_formulas, cardiac_sys_names, "Systole"),
    compare_group(df_neural, main_systole_formulas_neural, neural_sys_names, "Systole"),
    compare_group(df_neural, main_diastole_formulas_neural, neural_dia_names, "Diastole")
  )
  comparison_file <- file.path(modelsummary_path, "main_random_slope_LOO_comparisons.csv")
  write.csv(comparison_table, comparison_file, row.names = FALSE)
  message("Saved main-model random-slope LOO comparisons: ", comparison_file)

  # Word comparison table
  comparison_table_word <- comparison_table
  names(comparison_table_word)[names(comparison_table_word) == "Delta_LOOIC_without_minus_with"] <-
    "ΔLOOIC (without - with)"
  ft <- flextable::flextable(comparison_table_word) %>%
    flextable::autofit() %>%
    flextable::set_table_properties(layout = "autofit", width = 1) %>%
    flextable::fontsize(size = 9, part = "all")
  doc <- officer::read_docx()
  doc <- officer::body_add_par(doc, "Model Comparisons (Random Intercepts vs. Slopes)", style = "heading 1")
  doc <- flextable::body_add_flextable(doc, value = ft)
  comparison_word_file <- file.path(chapter_path, "model_comparisons.docx")
  print(doc, target = comparison_word_file)
  message("Written: ", comparison_word_file)

  rm(df_cardiac, df_neural, comparison_table)
}
