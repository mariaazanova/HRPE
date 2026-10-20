# Figure 3 — pseudotrial-corrected Pz responses and difference topographies
# Run this file in RStudio from codes/src_r/
# Input: data/GA/ (grand-average *_Channels_CI.csv and *_topo.csv; filename_mapping.csv if files have shortened names)
# Output: vis/origart_fig3.png
# Figures are read from previously computed grand averages (in src_python)

library(logger)

# Paths and helpers ####
code_path <- dirname(rstudioapi::getSourceEditorContext()$path)
project_path <- dirname(dirname(code_path))
GA_path <- file.path(project_path, "data", "GA")
vis_path <- file.path(project_path, "vis")
dir.create(vis_path, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(GA_path)) stop("Missing neural grand-average folder: ", GA_path)

source(file.path(code_path, "plotting_helpers.R"))
source(file.path(code_path, "neural_plotting_functions.R"))
colors <- get_colors()

# A–D: four pseudotrial-corrected grand averages at Pz ####
versions <- c("5_4", "6_5")
phases <- c("Systole", "Diastole")
responses <- list()

make_dif_file <- function(phase, v) {
  sprintf(
    paste0(
      "dif_sysDiaFeedb_%s_changeIBI%s_2groups_high_corrected_",
      "sysDiaFeedb_%s_changeIBI%s_2groups_high__",
      "sysDiaFeedb_%s_changeIBI%s_2groups_low_corrected_",
      "sysDiaFeedb_%s_changeIBI%s_2groups_low"
    ),
    phase, v, phase, v, phase, v, phase, v
  )
}

# Require the CI time courses used in the figure
mapping <- load_filename_mapping(file.path(GA_path, "filename_mapping.csv"))$orig2short
check_ga_ci <- function(condition) {
  virtual <- paste0(condition, "_Channels_CI.csv")
  actual <- if (virtual %in% names(mapping)) mapping[[virtual]] else virtual
  path <- file.path(GA_path, actual)
  if (!file.exists(path)) stop("Missing pseudotrial-corrected CI time course: ", path)
}

for (v in versions) {
  for (phase in phases) {
    dif_file <- make_dif_file(phase, v)
    processed <- process_dif_file(dif_file, GA_path, prefixes)
    cleaned <- clean_condition(
      processed$condition1, processed$condition2,
      processed$label_condition1, processed$label_condition2,
      condition_name = TRUE, remove_condition_name = FALSE,
      Rpeak_in_name = TRUE
    )
    title <- decode_file_name(cleaned$cleaned_condition, Rpeak_in_name = FALSE)
    keywords <- c(processed$condition1, processed$condition2)
    invisible(lapply(keywords, check_ga_ci))

    responses[[paste(phase, v, sep = "_")]] <- plot_evoked_responses(
      keywords, GA_path, vis_path,
      labels = c("High", "Low"),
      legend_title = cleaned$condition_name,
      show_title = TRUE, title_string = title,
      plot_ci = TRUE, chosen_channels = c("Pz"),
      show_plots = FALSE, save_plots = FALSE,
      plot_format = "pdf", palette = "PuOr", colors = colors,
      remove_grid = TRUE, add_0 = TRUE,
      add_shaded_areas = FALSE, pick_roi = "all",
      latencies_file = NULL,
      plot_width = 8, plot_height = 6,
      legend_position = "right", tmin = -0.3, tmax = 0.7
    )
  }
}

pA <- responses$Systole_5_4 +
  labs(title = "", subtitle = "Pz, IBI 1 Change early in cycle", colour = "")
pB <- responses$Diastole_5_4 +
  labs(title = "", subtitle = "Pz, IBI 1 Change late in cycle", colour = "")
pC <- responses$Systole_6_5 +
  labs(title = "", subtitle = "Pz, IBI 2 Change early in cycle", colour = "")
pD <- responses$Diastole_6_5 +
  labs(title = "", subtitle = "Pz, IBI 2 Change late in cycle", colour = "")

# E/F: difference topographies in the phase-specific P300 time windows ####
pE <- plot_topographies(
  make_dif_file("Systole", "5_4"), GA_path, vis_path,
  palette = "PuOr", show_plots = FALSE, save_plots = FALSE,
  plot_format = "pdf", chosen_roi = "p3b_Win",
  fill_title = NULL, title = NULL,
  plot_width = 6, plot_height = 6, adjust_colorscale = FALSE
)
pF <- plot_topographies(
  make_dif_file("Diastole", "6_5"), GA_path, vis_path,
  palette = "PuOr", show_plots = FALSE, save_plots = FALSE,
  plot_format = "pdf", chosen_roi = "p3b_Lose",
  fill_title = NULL, title = NULL,
  plot_width = 6, plot_height = 6, adjust_colorscale = FALSE
)

# Adjust point-geom size (electrode positions)
reduce_sensor_point_size <- function(p) {
  point_layers <- which(vapply(
    p$layers,
    function(layer) inherits(layer$geom, "GeomPoint"),
    logical(1)
  ))
  if (length(point_layers) == 0) {
    stop("The topography plot has no point layer for sensor markers.")
  }
  p$layers[[point_layers[1]]]$aes_params$size <- 0.6
  p
}
pE <- reduce_sensor_point_size(pE) +
  labs(title = "", subtitle = "IBI 1 in Systole", colour = "")
pF <- reduce_sensor_point_size(pF) +
  labs(title = "", subtitle = "IBI 2 in Diastole", colour = "")

pE <- pE +
  guides(fill = guide_colourbar(
    barheight = grid::unit(2.5, "cm"),
    barwidth = grid::unit(0.25, "cm")
  )) +
  labs(fill = expression(mu * V))
pF <- pF +
  guides(fill = guide_colourbar(
    barheight = grid::unit(2.5, "cm"),
    barwidth = grid::unit(0.25, "cm")
  )) +
  labs(fill = expression(mu * V))

# Assemble A–F using the 2-column layout ####
first_row_first_column <- plot_grid(
  ggdraw(pA), ggdraw(pB), ggdraw(pC), ggdraw(pD),
  ncol = 2, labels = c("A", "B", "C", "D"),
  rel_widths = c(1, 1, 1, 1), align = "h"
)
first_row_second_column <- plot_grid(
  ggdraw(pE), ggdraw(pF), ncol = 1,
  labels = c("E", "F"), rel_widths = c(1, 1, 1, 1), align = "h"
)
first_row <- plot_grid(
  first_row_first_column, first_row_second_column,
  ncol = 2, rel_widths = c(1, 0.35), align = "h"
)
fig <- plot_grid(first_row, ncol = 1)
save_figure(fig, aspect_ratio = 0.75, vis_path,
            filename = "origart_fig3", ext = "png")
message("Saved ", file.path(vis_path, "origart_fig3.png"))
