# This script LOADS existing model summaries instead of running models again.
# Run CLEARN_13_run.R first to ensure model summaries exist.
# Output: data/ModelSummaries/CLEARN_FDR/

# Setup ####
library(logger)

code_path <- dirname(rstudioapi::getSourceEditorContext()$path)
data_path <- file.path(dirname(dirname(code_path)), "data")
data_filename <- file.path(data_path, "data.csv")
statresult_path <- file.path(data_path, "StatResults")
modelsummary_path <- file.path(data_path, "ModelSummaries")
fdr_path <- file.path(modelsummary_path, "CLEARN_FDR")
dir.create(fdr_path, recursive = TRUE, showWarnings = FALSE)

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

# Behavioural data preprocessing ####
df <- read.csv(data_filename)
result <- preprocess_df(df, rescale_after_outliers = TRUE)
df <- result$df_filtered
rm(result)
distances_df <- compute_distances(df)
df <- cbind(df, distances_df)
rm(distances_df)

# FDR function ####
bayes_fdr_group <- function(
    group_block,
    coef = NULL,
    coef_col = c("hypothesis_parameter", "Model_parameter"),
    q = 0.05,
    lfdr_col = "ROPE_Percentage",
    exclude_intercept = TRUE
) {
  coef_col <- match.arg(coef_col)

  sums <- group_block$summaries
  df <- do.call(rbind, lapply(seq_along(sums), function(i) {
    x <- sums[[i]]
    x$model_index <- i
    x
  }))

  if (exclude_intercept) {
    df <- df[df$hypothesis_parameter != "Intercept", , drop = FALSE]
  }

  if (!is.null(coef)) {
    if (coef_col == "hypothesis_parameter") {
      df <- df[df$hypothesis_parameter %in% coef, , drop = FALSE]
    } else {
      df <- df[df$Model_parameter %in% coef, , drop = FALSE]
    }
  }

  # The local false discovery measure is ROPE_Percentage, NOT p_ROPE.
  df$lfdr <- df[[lfdr_col]]
  ord <- order(df$lfdr, na.last = TRUE)
  df <- df[ord, , drop = FALSE]

  df$rank <- seq_len(nrow(df))
  df$cum_fdr <- cumsum(df$lfdr) / df$rank

  ok <- which(df$cum_fdr <= q)
  k <- if (length(ok) == 0) 0 else max(ok)

  df$discover <- FALSE
  if (k > 0) df$discover[seq_len(k)] <- TRUE

  list(
    q = q,
    n_tested = nrow(df),
    n_discoveries = sum(df$discover),
    lfdr_cutoff = if (k > 0) max(df$lfdr[seq_len(k)], na.rm = TRUE) else NA_real_,
    results = df
  )
}

# Load the 16 main model summaries in their original order ####
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

load_model_summary <- function(formula) {
  brms_setup <- prior_pred_check(df, formula, run_pp = FALSE)
  summary_file <- file.path(modelsummary_path, paste0(brms_setup$filename, ".rds"))
  if (!file.exists(summary_file)) {
    stop("Missing summary from CLEARN_13: ", summary_file)
  }
  info <- readRDS(summary_file)
  if (is.null(info$summary_df)) stop("summary_df is absent in: ", summary_file)
  info$summary_df
}

# Four matched models per behavioural predictor
formulas_by_predictor <- split(main_formulas, rep(c("expectation", "outcome", "PE", "absPE"), each = 4))
summary_groups <- lapply(formulas_by_predictor, function(formulas) {
  list(summaries = lapply(formulas, load_model_summary))
})

# Exact seven tests from try_fdr_correction_example.R ####
fdr_results <- list(
  expectation = bayes_fdr_group(summary_groups$expectation,
    coef = "scale_expectation", coef_col = "hypothesis_parameter", q = 0.05),
  outcome = bayes_fdr_group(summary_groups$outcome,
    coef = "outcomeLose", coef_col = "hypothesis_parameter", q = 0.05),
  outcome_phase = bayes_fdr_group(summary_groups$outcome,
    coef = "outcomeLose:scale_phaseFeedb", coef_col = "hypothesis_parameter", q = 0.05),
  PE = bayes_fdr_group(summary_groups$PE,
    coef = "scale_PE", coef_col = "hypothesis_parameter", q = 0.05),
  PE_phase = bayes_fdr_group(summary_groups$PE,
    coef = "scale_PE:scale_phaseFeedb", coef_col = "hypothesis_parameter", q = 0.05),
  absPE = bayes_fdr_group(summary_groups$absPE,
    coef = "scale_absPE", coef_col = "hypothesis_parameter", q = 0.05),
  absPE_phase = bayes_fdr_group(summary_groups$absPE,
    coef = "scale_absPE:scale_phaseFeedb", coef_col = "hypothesis_parameter", q = 0.05)
)

# Save FDR reports ####
for (name in names(fdr_results)) {
  result_fdr <- fdr_results[[name]]
  if (result_fdr$n_tested != 4) {
    stop("Expected exactly four heartbeat-wise tests for ", name,
         "; found ", result_fdr$n_tested,
         ". Check summary coefficient names before interpreting the result.")
  }
  # Model index is within its four-formula predictor group (not the 16-item list)
  model_type <- sub("_phase$", "", name)
  results_df <- result_fdr$results
  predictor_formulas <- formulas_by_predictor[[model_type]]
  results_df$formula <- vapply(results_df$model_index, function(i) {
    paste(predictor_formulas[[i]], collapse = " ")
  }, character(1))
  output_file <- file.path(fdr_path, paste0("fdr_", name, ".csv"))
  if (file.exists(output_file)) stop("Output already exists (not overwritten): ", output_file)
  write.csv(results_df, output_file, row.names = FALSE)
  cat("\n", name, ": ", result_fdr$n_discoveries, "/", result_fdr$n_tested,
      " discoveries\n", sep = "")
  print(results_df[, c("formula", "hypothesis_parameter", "lfdr", "cum_fdr", "discover")])
}

message("FDR results saved under: ", fdr_path)
