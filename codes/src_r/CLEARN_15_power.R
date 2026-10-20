# Posterior predictive evidential power for the main models
# Run CLEARN_13_run.R before this script to ensure fits exist.
# Default: 34 participants
# New results go into data/ModelSummaries/CLEARN_power/

# Setup ####
library(logger)

code_path <- dirname(rstudioapi::getSourceEditorContext()$path)
data_path <- file.path(dirname(dirname(code_path)), "data")
data_filename <- file.path(data_path, "data.csv")
statresult_path <- file.path(data_path, "StatResults")
modelsummary_path <- file.path(data_path, "ModelSummaries")
power_path <- file.path(modelsummary_path, "CLEARN_power")
dir.create(power_path, recursive = TRUE, showWarnings = FALSE)

# Original helper files, unchanged and in their original load order.
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

# Simulation settings ####
n_subs <- c(34, 100)  
n_sims <- 200         
gate_scope <- "per_term"

# Reusable wrapper
run_power <- function(df, formula, n_s, n_trials, pick_phase = NULL) {
  brms_setup <- prior_pred_check(
    df, formula, run_pp = FALSE, random_slope = FALSE, pick_phase = pick_phase
  )
  model_file <- file.path(statresult_path, paste0(brms_setup$filename, ".rds"))
  output_file <- file.path(power_path, paste0("pwr_", n_s, "_", brms_setup$filename, ".rds"))

  if (file.exists(output_file)) {
    message("Existing CLEARN power result left unchanged: ", output_file)
    return(invisible(NULL))
  }
  if (!file.exists(model_file)) stop("Missing fitted model: ", model_file)
  fit <- readRDS(model_file)

  pilot_df <- df
  if (!is.null(pick_phase)) {
    pilot_df <- df[(df$sysDiaFeedb == pick_phase) & (!is.na(df$sysDiaFeedb)), ]
  }

  power_res <- simulate_power_BF_pp(
    brms_fit = fit,
    pilot_df = pilot_df,
    n_sim = n_sims,
    n_sub = n_s,
    n_trials = n_trials,
    BF_threshold = 10,
    refit_iter = 2000,
    chains = 2,
    cores = 2,
    hypothesis_direction = "equal",
    param_source = "single_draw",
    base_seed = 1234,
    diagnostics_gate = TRUE,
    adapt_delta = 0.95,
    max_treedepth = 12,
    contrasts_spec = NULL,
    include_intercept = FALSE,
    rhat_threshold = 1.01,
    gating_scope = gate_scope,
    eps = 1e-12
  )

  saveRDS(power_res, output_file)
  message("Saved: ", output_file)
  print(power_res$evidence_effect)  # BF > 10, after global diagnostic gating
  print(colMeans(as.matrix(power_res$BF_matrix_final) > 3, na.rm = TRUE))
  invisible(power_res)
}

# Behavioural/cardiac models ####
df <- read.csv(data_filename)
result <- preprocess_df(df, rescale_after_outliers = TRUE)
df <- result$df_filtered
rm(result)

# Exact original five overall models
main_overall_formulas <- list(
  c("changeIBI4_3", "expectation", "*", "phaseFeedb"),
  c("changeIBI7_6", "expectation", "*", "phaseFeedb"),
  c("changeIBI6_5", "outcome", "*", "phaseFeedb"),
  c("changeIBI6_5", "PE", "*", "phaseFeedb"),
  c("changeIBI6_5", "absPE", "*", "phaseFeedb")
)

# Exact original three systole-only models
main_systole_formulas <- list(
  c("changeIBI6_5", "outcome"),
  c("changeIBI6_5", "PE"),
  c("changeIBI6_5", "absPE")
)

for (n_s in n_subs) {
  for (formula in main_overall_formulas) {
    run_power(df, formula, n_s, n_trials = 300)
  }
  for (formula in main_systole_formulas) {
    run_power(df, formula, n_s, n_trials = 100, pick_phase = "Systole")
  }
}
rm(df)

# Neural/P3b models ####
df <- read.csv(data_filename)
result <- preprocess_df(df, rescale_after_outliers = TRUE, clean_neural = TRUE)
df <- result$df_filtered
rm(result)
df <- replace_double_underscore_in_colnames(df)

# Models for each cardiac phase.
main_systole_formulas_neural <- list(
  c("changeIBI5_4", "Pz_p3b_Win"),
  c("changeIBI5_4", "Pz_p3b_Lose")
)
main_diastole_formulas_neural <- list(
  c("changeIBI6_5", "Pz_p3b_Win"),
  c("changeIBI6_5", "Pz_p3b_Lose")
)

for (n_s in n_subs) {
  for (formula in main_systole_formulas_neural) {
    run_power(df, formula, n_s, n_trials = 100, pick_phase = "Systole")
  }
  for (formula in main_diastole_formulas_neural) {
    run_power(df, formula, n_s, n_trials = 100, pick_phase = "Diastole")
  }
}

message("Power simulations processed. Outputs: ", power_path)
