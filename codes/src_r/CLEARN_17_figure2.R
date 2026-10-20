# Figure 2 — expectation and signed prediction-error associations
# Run this file in RStudio from codes/src_r/
# Input: data/data.csv
# Output: vis/origart_fig2.png

# Paths and helpers ####
code_path <- dirname(rstudioapi::getSourceEditorContext()$path)
project_path <- dirname(dirname(code_path))
data_filename <- file.path(project_path, "data", "data.csv")
vis_path <- file.path(project_path, "vis")
dir.create(vis_path, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(data_filename)) stop("Missing data: ", data_filename)

source(file.path(code_path, "preprocess.R"))
source(file.path(code_path, "plotting_helpers.R"))
source(file.path(code_path, "plot_model_results_functions.R"))
source(file.path(code_path, "analysis_functions.R"))

colors <- get_colors()

# Load ####
df <- read.csv(data_filename)
result <- preprocess_df(df, rescale_after_outliers = TRUE)
df <- result$df_filtered
rm(result)
df <- replace_double_underscore_in_colnames(df)

# A–D: model-derived expectation ####
# Only the ten-bin grouping is needed for Figure 2 (within participant).
df <- df %>%
  group_by(subID) %>%
  mutate(expectation_10groups = factor(ntile(expectation, 10))) %>%
  dplyr::ungroup()

pA <- plot_empirical_binned(
  df, c("changeIBI4_3", "expectation"), bins = 500
) + labs(y = "Change of IBI -1 (ms)") + ylim(-50, 50)

pB <- plot_empirical_binned(
  df, c("changeIBI7_6", "expectation"), bins = 500
) + labs(y = "Change of IBI 3 (ms)") + ylim(-50, 50)

pC <- plot_subject_summary(
  df, "changeIBI4_3", "expectation_10groups", type = "bar",
  dotplot_binwidth = 1, dotplot_size = 1, dot_move = -0.3,
  fill_color = colors[1]
) + ylim(-25, 15) +
  labs(y = "Change of IBI -1 (ms)", x = "Expectation bins") +
  plot_theme_default

pD <- plot_subject_summary(
  df, "changeIBI7_6", "expectation_10groups", type = "bar",
  dotplot_binwidth = 1, dotplot_size = 1, dot_move = -0.3,
  fill_color = colors[1]
) + ylim(-25, 15) +
  labs(y = "Change of IBI 3 (ms)", x = "Expectation bins") +
  plot_theme_default

# E–H: signed PE at IBI 2, early versus late cardiac cycle ####
df <- df %>%
  group_by(subID) %>%
  mutate(PE_10groups = factor(ntile(PE, 10))) %>%
  dplyr::ungroup()

sysdf <- df[!is.na(df$sysDiaFeedb) & df$sysDiaFeedb == "Systole", , drop = FALSE]
diadf <- df[!is.na(df$sysDiaFeedb) & df$sysDiaFeedb == "Diastole", , drop = FALSE]

pE <- plot_empirical_binned(
  sysdf, c("changeIBI6_5", "PE"), bins = 500
) + labs(y = "Change of IBI 2 (ms)", subtitle = "Early cardiac cycle") +
  ylim(-50, 50)

pF <- plot_empirical_binned(
  diadf, c("changeIBI6_5", "PE"), bins = 500
) + labs(y = "Change of IBI 2 (ms)", subtitle = "Late cardiac cycle") +
  ylim(-50, 50)

pG <- plot_subject_summary(
  sysdf, "changeIBI6_5", "PE_10groups", type = "bar",
  dotplot_binwidth = 1, dotplot_size = 1, dot_move = -0.3,
  fill_color = colors[1]
) + ylim(-27, 17) +
  labs(y = "Change of IBI 2 (ms)", x = "PE bins",
       subtitle = "Early cardiac cycle") +
  plot_theme_default

pH <- plot_subject_summary(
  diadf, "changeIBI6_5", "PE_10groups", type = "bar",
  dotplot_binwidth = 1, dotplot_size = 1, dot_move = -0.3,
  fill_color = colors[1]
) + ylim(-27, 17) +
  labs(y = "Change of IBI 2 (ms)", x = "PE bins",
       subtitle = "Late cardiac cycle") +
  plot_theme_default

# Assemble A–H ####
first_row <- plot_grid(
  ggdraw(pA), ggdraw(pB), ggdraw(pC), ggdraw(pD),
  ncol = 4, labels = c("A", "B", "C", "D"), align = "h"
)
second_row <- plot_grid(
  ggdraw(pE), ggdraw(pF), ggdraw(pG), ggdraw(pH),
  ncol = 4, labels = c("E", "F", "G", "H"), align = "h"
)
fig <- plot_grid(first_row, second_row, ncol = 1, rel_heights = c(1, 1))
save_figure(fig, aspect_ratio = 0.75, vis_path,
            filename = "origart_fig2", ext = "png")

message("Saved ", file.path(vis_path, "origart_fig2.png"))
