# ==============================================================================
# APPENDIX B: ABLATION SPECTRA VISUALIZATION
# ==============================================================================

# Ensure required plotting libraries are loaded
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
library(ggplot2)
library(patchwork)

# ------------------------------------------------------------------------------
# 1. DYNAMICALLY CALCULATE BIN TARGETS BASED ON CURRENT 'N'
# ------------------------------------------------------------------------------
n_rows <- nrow(observed_data)

bin_default <- ceiling(n_rows^(1/3))
bin_dense   <- ceiling(sqrt(n_rows))

# Dynamically estimate Freedman-Diaconis across a sample of continuous variables
# (Filtering for variables with >10 unique values as a proxy for continuous)
cont_cols <- names(observed_data)[sapply(observed_data, function(x) data.table::uniqueN(x) > 10)]

if(length(cont_cols) > 0) {
  set.seed(42)
  sample_cols <- sample(cont_cols, min(100, length(cont_cols)))
  # Use base R's nclass.FD to get the Freedman-Diaconis bins per variable
  fd_counts <- sapply(observed_data[, ..sample_cols], nclass.FD)
  bin_extreme <- round(mean(fd_counts))
} else {
  bin_extreme <- ceiling(n_rows^(2/3)) # Safe fallback if no continuous vars exist
}

cat("\n=================================================================\n")
cat(" DYNAMIC ABLATION TARGETS\n")
cat("=================================================================\n")
cat(sprintf(" -> Conservative (N^(1/3)) : %d bins\n", bin_default))
cat(sprintf(" -> Dense (N^(1/2))        : %d bins\n", bin_dense))
cat(sprintf(" -> Extreme (FD Avg)       : %d bins\n", bin_extreme))
cat("=================================================================\n")

# ------------------------------------------------------------------------------
# 2. RUN EXTRACTIONS
# ------------------------------------------------------------------------------
cat(sprintf("\n[1/3] Extracting baseline eigenvalues (%d bins)...\n", bin_default))
# Grab the baseline from original primary results object
eig_default <- results$eigenvalues

cat(sprintf("[2/3] Fetching eigenspectrum for %d bins...\n", bin_dense))
capture.output({
  res_dense <- calculate_entropic_scree(data = observed_data
                                        , num_bins = bin_dense 
                                        , interactive_mode = FALSE 
                                        , purge_constants = FALSE 
                                        , check_collinearity = FALSE)
})
eig_dense <- res_dense$eigenvalues

cat(sprintf("[3/3] Fetching eigenspectrum for %d bins...\n", bin_extreme))
capture.output({
  res_extreme <- calculate_entropic_scree(data = observed_data
                                          , num_bins = bin_extreme
                                          , interactive_mode = FALSE
                                          , purge_constants = FALSE
                                          , check_collinearity = FALSE)
})
eig_extreme <- res_extreme$eigenvalues

# ==============================================================================
# 3. OUTPUT VALUES FOR TABLE 2
# ==============================================================================
cat("\n=================================================================\n")
cat(" TABLE 2 UPDATE VALUES (R_eff)\n")
cat("=================================================================\n")
cat(sprintf(" -> %d Bins (Conservative) : %.3f\n", bin_default, results$R_eff))
cat(sprintf(" -> %d Bins (Dense)        : %.3f\n", bin_dense, res_dense$R_eff))
cat(sprintf(" -> %d Bins (Extreme)      : %.3f\n", bin_extreme, res_extreme$R_eff))
cat("=================================================================\n")

# ==============================================================================
# 4. PLOTTING ENGINE
# ==============================================================================

# Helper function to guarantee perfectly uniform styling across all 3 panels
create_panel <- function(eig_vec, title_text, hide_y_label = FALSE) {
  # Zoom window: Ranks 1 to 75 to isolate the boundary
  df <- data.frame(Rank = 1:75, Eigenvalue = eig_vec[1:75])
  
  # Calculate 1.025x of the second eigenvalue to use as the upper limit
  y_upper_limit <- eig_vec[2] * 1.025
  
  p <- ggplot(df, aes(x = Rank, y = Eigenvalue)) +
    geom_line(color = "dodgerblue", linewidth = 1) +
    geom_point(color = "dodgerblue", size = 1.5) +
    # True generative boundary (Dynamically points to K_TRUE from global environment)
    geom_vline(xintercept = K_TRUE, color = "#D55E00", linetype = "dashed", linewidth = 1.2) +
    scale_y_continuous(trans = 'log10') +
    # Use coord_cartesian to visually zoom without discarding the first data point
    coord_cartesian(ylim = c(NA, y_upper_limit)) +
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

# Generate the 3 panels dynamically (bquote and sprintf insert the exact dynamic bin counts)
p1 <- create_panel(eig_default, bquote(atop(bold("Conservative Default" ~ (N^{1/3})), bold(.(sprintf("[%d Bins]", bin_default))))), hide_y_label = FALSE)
p2 <- create_panel(eig_dense, bquote(atop(bold("Dense" ~ (N^{1/2})), bold(.(sprintf("[%d Bins]", bin_dense))))), hide_y_label = TRUE)
p3 <- create_panel(eig_extreme, bquote(atop(bold("Freedman-Diaconis" ~ phantom(N^{1/2})), bold(.(sprintf("[%d Bins]", bin_extreme))))), hide_y_label = TRUE)

# Combine using patchwork
final_plot <- p1 + p2 + p3 + 
  plot_layout(ncol = 3) + 
  plot_annotation(
    title = "Eigenspectra Across Discretization Regimes",
    subtitle = sprintf("Generative signal (K=%d) remains invariant despite inflation of the continuous noise tail", K_TRUE),
    theme = theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5, color = "gray30", margin = margin(b = 15))
    )
  )

# Render to the viewer
print(final_plot)
