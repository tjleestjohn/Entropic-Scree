# ==============================================================================
# APPENDIX B: ABLATION SPECTRA VISUALIZATION
# ==============================================================================

# Ensure required plotting libraries are loaded
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
library(ggplot2)
library(patchwork)

cat("\n[1/3] Extracting baseline eigenvalues (17 bins)...\n")
# 1. Grab the baseline from your original primary results object
eig_17 <- results$eigenvalues

# 2. Re-run the engine quickly for 71 and 152 bins just to grab their eigenspectra 
cat("[2/3] Fetching eigenspectrum for 71 bins...\n")
capture.output({
  res_71 <- calculate_entropic_scree(data = observed_data
                                     , num_bins = 71 
                                     , interactive_mode = FALSE 
                                     , purge_constants = FALSE 
                                     , check_collinearity = FALSE)
})
eig_71 <- res_71$eigenvalues

cat("[3/3] Fetching eigenspectrum for 152 bins...\n")
capture.output({
  res_152 <- calculate_entropic_scree(data = observed_data
                                      , num_bins = 152
                                      , interactive_mode = FALSE
                                      , purge_constants = FALSE
                                      , check_collinearity = FALSE)
})
eig_152 <- res_152$eigenvalues

# ==============================================================================
# PLOTTING ENGINE
# ==============================================================================

# Helper function to guarantee perfectly uniform styling across all 3 panels
create_panel <- function(eig_vec, title_text, hide_y_label = FALSE) {
  # Zoom window: Ranks 1 to 75 to isolate the boundary and the artifact
  df <- data.frame(Rank = 1:75, Eigenvalue = eig_vec[1:75])
  
  p <- ggplot(df, aes(x = Rank, y = Eigenvalue)) +
    geom_line(color = "dodgerblue", linewidth = 1) +
    geom_point(color = "dodgerblue", size = 1.5) +
    # True generative boundary (Invariant)
    geom_vline(xintercept = 10, color = "#D55E00", linetype = "dashed", linewidth = 1.2) +
    # Visual guide for the noise artifact at 55
    geom_vline(xintercept = 55, color = "gray60", linetype = "dotted", linewidth = 1) +
    scale_y_continuous(trans = 'log10') +
    labs(
      title = title_text, 
      x = "Eigenvalue Index", 
      y = if (hide_y_label) "" else "Log(Eigenvalue)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      panel.grid.minor = element_blank(),
      axis.title.y = element_text(margin = margin(r = 10))
    )
  
  return(p)
}

# Generate the 3 panels (bquote allows for clean mathematical exponents in the titles)
p1 <- create_panel(eig_17, bquote("Conservative" ~ (N^{1/3}) ~ "[17 Bins]"), hide_y_label = FALSE)
p2 <- create_panel(eig_71, bquote("Dense" ~ (N^{1/2}) ~ "[71 Bins]"), hide_y_label = TRUE)
p3 <- create_panel(eig_152, "Freedman-Diaconis [152 Bins]", hide_y_label = TRUE)

# Combine using patchwork
final_plot <- p1 + p2 + p3 + 
  plot_layout(ncol = 3) + 
  plot_annotation(
    title = "Eigenspectra Across Discretization Regimes",
    subtitle = "Generative signal (K=10) remains invariant while the internal noise artifact (K=55) dissolves",
    theme = theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5, color = "gray30", margin = margin(b = 15))
    )
  )

# Render to the viewer
print(final_plot)