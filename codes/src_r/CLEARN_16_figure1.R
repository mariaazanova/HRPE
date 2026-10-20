# Figure 1 — experimental paradigm and beat-to-beat cardiac dynamics
# Run this file in RStudio from codes/src_r/
# Inputs: data/data.csv; vis/design_figures/Slide1.png, Slide11.png, Slide13.png
# Output: vis/origart_fig1.png

library(magick)

# Paths and existing functions ####
code_path <- dirname(rstudioapi::getSourceEditorContext()$path)
project_path <- dirname(dirname(code_path))
data_filename <- file.path(project_path, "data", "data.csv")
vis_path <- file.path(project_path, "vis")
png_dir <- file.path(vis_path, "design_figures")

dir.create(vis_path, recursive = TRUE, showWarnings = FALSE)
source(file.path(code_path, "preprocess.R"))
source(file.path(code_path, "plotting_helpers.R"))
source(file.path(code_path, "plot_model_results_functions.R"))
source(file.path(code_path, "rpeak_plotting_functions.R"))

required_files <- c(
  data_filename,
  file.path(png_dir, c("Slide1.png", "Slide11.png", "Slide13.png"))
)
if (any(!file.exists(required_files))) {
  stop("Missing Figure 1 input(s):\n", paste(required_files[!file.exists(required_files)], collapse = "\n"))
}

colors <- get_colors()

# Data for the empirical plots ####

df <- read.csv(data_filename)
result <- preprocess_df(df, rescale_after_outliers = TRUE)
df <- result$df_filtered
rm(result)

df_plot <- df %>%
  mutate(
    sysDiaFeedb = recode(
      as.character(sysDiaFeedb),
      "Systole" = "Early",
      "Diastole" = "Late"
    ),
    sysDiaFeedb = factor(sysDiaFeedb, levels = c("Early", "Late")),
    phaseFeedb_01 = phaseFeedb / 360
  )

# Three existing design panels ####
read_design_panel <- function(slide, border) {
  img <- image_read(file.path(png_dir, paste0("Slide", slide, ".png")))
  img <- image_trim(img)
  img <- image_border(img, color = "white", geometry = border)
  ggdraw() + draw_image(img)
}

design_A <- read_design_panel(1, "30x30")
design_B <- read_design_panel(11, "150x150")
design_C <- read_design_panel(13, "30x30")

# F/G: Average IBI and consecutive IBI change by early/late phase ####

phase_group <- c("sysDiaFeedb")
ibi_args <- list(
  df = df_plot, groups = phase_group, colors = colors,
  add_line = TRUE, show_fewer_Rpeaks = FALSE,
  output_as_one = TRUE, vis_path = vis_path,
  plot_width = 8, plot_height = 6,
  show_plots = FALSE, save_plots = FALSE, plot_format = "pdf",
  show_title = FALSE, title_string = "",
  legend_title = "Cardiac cycle", show_legend_title = TRUE,
  col_ = NULL, dpi = 150, say_ibi = TRUE
)

ibi_args$pure_IBI <- TRUE
ibi_args$show_legend <- FALSE
pF <- do.call(plot_ibi, ibi_args)$sysDiaFeedb + plot_theme_default

ibi_args$pure_IBI <- FALSE
ibi_args$show_legend <- TRUE
pG <- do.call(plot_ibi, ibi_args)$sysDiaFeedb + plot_theme_default

# D/E: Empirical phase effects on the changes in IBI 1 and IBI 2 ####
pD <- plot_empirical_binned(
  df_plot, c("changeIBI5_4", "phaseFeedb_01"), bins = 500
) +
  labs(x = "Cardiac cycle phase", y = "Change in IBI 1 (ms)") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))

pE <- plot_empirical_binned(
  df_plot, c("changeIBI6_5", "phaseFeedb_01"), bins = 500
) +
  labs(x = "Cardiac cycle phase", y = "Change in IBI 2 (ms)") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))

# Assemble A–G with panel ordering and relative sizes ####
first_row <- plot_grid(
  design_A, design_B,
  ncol = 2, labels = c("A", "B"), align = "h",
  rel_widths = c(1, 0.6)
)
second_row_first_column <- plot_grid(
  design_C, ncol = 1, labels = c("C"), align = "h"
)
second_row_second_column <- plot_grid(
  ggdraw(pD), ggdraw(pE), ncol = 2,
  labels = c("D", "E"), align = "h"
)
second_row <- plot_grid(
  second_row_first_column, second_row_second_column,
  ncol = 2, rel_widths = c(0.5, 1), align = "h"
)
third_row <- plot_grid(
  ggdraw(pF), ggdraw(pG), ncol = 2,
  labels = c("F", "G"), rel_widths = c(1, 1), align = "h"
)
fig <- plot_grid(
  first_row, second_row, third_row,
  ncol = 1, rel_heights = c(0.8, 1, 1)
)
save_figure(fig, aspect_ratio = 1, vis_path,
            filename = "origart_fig1", ext = "png")

message("Saved ", file.path(vis_path, "origart_fig1.png"))
