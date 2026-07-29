# ==============================================================================
# ENTROPIC SCREE: Simulation & Demonstration Script
#
# Author: Terrence J. Lee-St. John, PhD
# Organization: Enli (www.enli.com.au)
# 
# Description: Generates a high-dimensional, mixed-type, noisy synthetic 
# dataset to demonstrate the systematic degradation of standard PCA, and utilizes 
# the Entropic Scree to extract the Latent Generative Rank (r).
# ==============================================================================

options(max.print = 999999)

rm(list = ls())
gc(verbose = FALSE)

# List of all required packages (added pkgbuild for robust Rtools checking)
required_packages <- c("Rcpp", "data.table", "infotheo", "ggplot2", "patchwork", "MASS", "stringr", "pkgbuild")

# Find out which ones are missing
missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

# Install the missing ones
if(length(missing_packages) > 0) {
  cat("Installing missing dependencies:", paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages)
} else {
  cat("All package dependencies are already installed!\n")
}

# ==============================================================================
# RTOOLS C++ COMPILER CHECK (WINDOWS ONLY)
# ==============================================================================
if (.Platform$OS.type == "windows") {
  # suppressWarnings hides the native pkgbuild warning so our clean custom block prints alone
  if (!suppressWarnings(pkgbuild::has_rtools())) {
    stop(
      "\n=================================================================\n",
      " [!] MISSING OR INCOMPATIBLE C++ COMPILER (Rtools)\n",
      "=================================================================\n",
      " Rtools is required to build the C++ backend on Windows.\n",
      " It is either missing, not on your PATH, or your Rtools version \n",
      " does not match your R version.\n\n",
      sprintf(" Your current R version is: %s\n", getRversion()),
      " You must install the version of Rtools that matches this R version.\n\n",
      " Please download and install the correct Rtools here:\n",
      " https://cran.r-project.org/bin/windows/Rtools/\n\n",
      " Note: After installing, you MUST restart your R session before \n",
      " running this script again.\n",
      "=================================================================\n",
      call. = FALSE
    )
  } else {
    cat("Compatible Rtools C++ compiler found. Ready to build backend.\n")
  }
}

# ==============================================================================
# 0. C++ OPENMP MUTUAL INFORMATION ENGINE (Rcpp)
# ==============================================================================
Rcpp::sourceCpp(code = '
#include <Rcpp.h>
#include <omp.h>
#include <cmath>
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix fast_parallel_MI(IntegerMatrix mat, int num_bins, int cores) {
    int n = mat.nrow();
    int p = mat.ncol();
    NumericMatrix MI(p, p);
    
    std::vector<std::vector<int>> margins(p, std::vector<int>(num_bins, 0));
    std::vector<double> H(p, 0.0);

    for(int j = 0; j < p; ++j) {
        for(int i = 0; i < n; ++i) {
            int val = mat(i, j) - 1; 
            if(val >= 0 && val < num_bins) margins[j][val]++;
        }
        double entropy = 0.0;
        for(int b = 0; b < num_bins; ++b) {
            if(margins[j][b] > 0) {
                double prob = (double)margins[j][b] / n;
                entropy -= prob * log(prob);
            }
        }
        H[j] = entropy;
        MI(j, j) = entropy;
    }

    #pragma omp parallel for num_threads(cores) schedule(dynamic)
    for(int i = 0; i < p; ++i) {
        for(int j = i + 1; j < p; ++j) {
            std::vector<int> joint(num_bins * num_bins, 0);
            
            for(int row = 0; row < n; ++row) {
                int val1 = mat(row, i) - 1;
                int val2 = mat(row, j) - 1;
                if(val1 >= 0 && val1 < num_bins && val2 >= 0 && val2 < num_bins) {
                    joint[val1 * num_bins + val2]++;
                }
            }
            
            double joint_entropy = 0.0;
            for(int k = 0; k < num_bins * num_bins; ++k) {
                if(joint[k] > 0) {
                    double prob = (double)joint[k] / n;
                    joint_entropy -= prob * log(prob);
                }
            }
            
            double mi_val = H[i] + H[j] - joint_entropy;
            if (mi_val < 0) mi_val = 0; 
            
            MI(i, j) = mi_val;
            MI(j, i) = mi_val;
        }
    }
    return MI;
}
')

# ==============================================================================
# 1. ENTROPIC SCREE FUNCTION
# ==============================================================================
calculate_entropic_scree <- function(data
                                     , low_entropy_thresh = 0.05
                                     , num_bins = NULL
                                     , bin_multiplier = 1.0
                                     , num_cores = parallel::detectCores() - 2
                                     , interactive_mode = TRUE
                                     , purge_constants = TRUE
                                     , check_collinearity = TRUE) {
  
  start_time <- Sys.time()
  dt <- data.table::as.data.table(data)
  
  # ----------------------------------------------------------------------------
  # [0/9] INITIAL DIMENSION CHECK
  # ----------------------------------------------------------------------------
  if (ncol(dt) < 2) {
    stop("Execution Halted: The input dataset must contain at least 2 columns to calculate mutual information.")
  }
  
  # ----------------------------------------------------------------------------
  # [1/9] PURGE CONSTANTS & DUPLICATES
  # ----------------------------------------------------------------------------
  if (purge_constants) {
    cat("[1/9] Purging constants and identical duplicates...\n")
    const_cols <- names(dt)[sapply(dt, function(x) data.table::uniqueN(x, na.rm = TRUE) <= 1)]
    if (length(const_cols) > 0) dt[, (const_cols) := NULL]
    
    dup_cols <- duplicated(as.list(dt))
    if (any(dup_cols)) dt <- dt[, !dup_cols, with = FALSE]
  } else {
    cat("[1/10] Skipping constant and duplicate purge (user requested)...\n")
  }
  
  # ----------------------------------------------------------------------------
  # [2/9] MULTIVARIATE COLLINEARITY CHECK
  # ----------------------------------------------------------------------------
  if (check_collinearity) {
    cat("[2/10] Checking for perfect multivariate linear combinations...\n")
    num_cols <- names(dt)[sapply(dt, is.numeric)]
    if (length(num_cols) > 1 && nrow(dt) > length(num_cols)) {
      
      p_cols <- length(num_cols)
      target_rows <- min(nrow(dt), p_cols + 500)
      
      if ((as.numeric(target_rows) * as.numeric(p_cols)) < 2147000000) {
        set.seed(42) 
        sample_idx <- sample(seq_len(nrow(dt)), target_rows)
        num_mat <- as.matrix(na.omit(dt[sample_idx, ..num_cols]))
        
        if (nrow(num_mat) > 0) {
          qr_mat <- cbind(Intercept = 1, num_mat)
          qr_decomp <- qr(qr_mat, tol = 1e-7)
          
          if (qr_decomp$rank < ncol(qr_mat)) {
            drop_indices <- qr_decomp$pivot[(qr_decomp$rank + 1):ncol(qr_mat)]
            lin_combos <- setdiff(colnames(qr_mat)[drop_indices], "Intercept")
            if (length(lin_combos) > 0) {
              cat(sprintf("      -> Purged %d perfectly collinear variables to protect downstream SIA.\n", length(lin_combos)))
              dt[, (lin_combos) := NULL]
            }
          }
        }
      } else {
        cat("      -> Matrix exceeds absolute LINPACK bounds even when sub-sampled. Bypassing...\n")
      }
    }
  } else {
    cat("[2/10] Skipping collinearity check (user requested)...\n")
  }
  
  cat("[3/10] Discretizing continuous data and dense-ranking categoricals...\n")
  if (is.null(num_bins)) {
    num_bins <- max(2, ceiling(bin_multiplier * nrow(dt)^(1/3))) 
  }
  
  cols <- names(dt)
  dt[, (cols) := lapply(.SD, function(x) {
    if (is.numeric(x) && data.table::uniqueN(x) > num_bins) {
      infotheo::discretize(x, disc = "equalfreq", nbins = num_bins)[[1]]
    } else {
      data.table::frank(x, ties.method = "dense")
    }
  })]
  
  cat("[4/10] Purging non-linear monotonic duplicates and locking types...\n")
  dup_cols_post <- duplicated(as.list(dt))
  if (any(dup_cols_post)) dt <- dt[, !dup_cols_post, with = FALSE]
  
  valid_cols <- names(dt)
  dt[, (valid_cols) := lapply(.SD, function(x) as.integer(as.factor(x)))]
  
  cat("[5/10] Calculating marginal entropies and purging near-constants...\n")
  H_vec <- sapply(dt, infotheo::entropy)
  valid_vars <- names(H_vec)[H_vec >= low_entropy_thresh]
  if (length(H_vec) > length(valid_vars)) {
    dt <- dt[, ..valid_vars]
    H_vec <- H_vec[valid_vars]
  }
  
  p <- ncol(dt)
  if (p < 2) stop("Execution Halted: Less than 2 valid variables remain.")
  
  cat(sprintf("[6/10] Computing %d x %d Mutual Information Matrix (C++ OpenMP)...\n", p, p))
  
  mat_data <- as.matrix(dt)
  bin_sample_sizes <- apply(mat_data, 2, tabulate)
  rm(dt) 
  gc(verbose = FALSE)
  
  safe_target_cores <- max(1, num_cores)
  MI_mat <- fast_parallel_MI(mat_data, num_bins = num_bins, cores = safe_target_cores)
  rownames(MI_mat) <- colnames(MI_mat) <- valid_vars
  
  cat("[7/10] Applying Joint Entropy (Jaccard) Normalization...\n")
  sum_H_mat <- outer(H_vec, H_vec, FUN = "+")
  joint_H_mat <- sum_H_mat - MI_mat
  joint_H_mat[joint_H_mat < 1e-9] <- 1e-9
  
  NMI_mat <- MI_mat / joint_H_mat
  diag(NMI_mat) <- 1.0
  
  cat("[8/10] Applying Double-Centering...\n")
  m_valid <- length(valid_vars)
  
  # Vectorized Double-Centering
  row_means <- rowMeans(NMI_mat)
  grand_mean <- mean(row_means)
  
  NMI_mat_c <- NMI_mat - outer(row_means, row_means, FUN = "+") + grand_mean
  
  # Calculate trace of the centered matrix
  Tr_Mc <- sum(diag(NMI_mat_c))
  mean_trace <- Tr_Mc / m_valid
  
  cat("[9/10] Extracting Entropic Latent Factors (Eigen Decomposition)...\n")
  eigen_res <- eigen(NMI_mat_c, symmetric = TRUE)
  raw_eig_vals <- eigen_res$values
  
  # Calculate SCDR before clipping
  m_plus <- sum(raw_eig_vals[raw_eig_vals > 0])
  m_minus <- sum(abs(raw_eig_vals[raw_eig_vals < 0]))
  SCDR <- (m_minus / m_plus) * 100
  
  eig_vals <- pmax(raw_eig_vals, 1e-9)
  
  # Constructive Spectral Mass (sum of positive clipped eigenvalues)
  m_plus <- sum(eig_vals)
  
  cat("[10/10] Calculating R_eff and Estimating Elbow...\n")
  sig_vals <- eig_vals[eig_vals > 0]
  if (length(sig_vals) > 0) {
    p_vals <- sig_vals / sum(sig_vals)
    H_spec <- -sum(p_vals * log(p_vals)) 
    R_eff <- exp(H_spec)
  } else {
    R_eff <- 1
  }
  
  # ==========================================================================
  # DIAGNOSTIC ONLY: MACRO GAP (NOISE Cliff BOUNDARY)
  # ==========================================================================
  n_total <- length(eig_vals)
  valid_k <- sum(eig_vals > mean_trace) 
  
  macro_max_noise_gap <- NA_real_
  macro_actual_gap <- NA_real_
  macro_gap_ratio <- NA_real_
  top_of_bulk_idx <- NA_integer_
  elbow_ratio <- NA_real_
  elbow_pct_drop <- NA_real_
  
  valid_search_space <- eig_vals[eig_vals > 1e-8]
  
  if (length(valid_search_space) > 10) {
    all_gaps_diag <- abs(diff(valid_search_space))
    n_active <- length(valid_search_space)
    noise_start_idx <- min(valid_k + max(3, floor(n_active * 0.05)), n_active - 5)
    noise_tail_idx <- noise_start_idx:(n_active - 1)
    
    if(length(noise_tail_idx) > 0) {
      noise_gaps <- all_gaps_diag[noise_tail_idx]
      max_noise_gap <- max(noise_gaps)
      
      macro_multiplier <- 20
      gap_threshold <- max(1e-6, max_noise_gap * macro_multiplier) 
      
      macroscopic_gap_indices <- which(all_gaps_diag > gap_threshold)
      if (length(macroscopic_gap_indices) > 0) {
        top_of_bulk_idx <- max(macroscopic_gap_indices) + 1
        macro_max_noise_gap <- max_noise_gap
        macro_actual_gap <- all_gaps_diag[top_of_bulk_idx - 1]
        if (max_noise_gap > 1e-9) macro_gap_ratio <- macro_actual_gap / max_noise_gap
      }
    }
  }
  
  # ==========================================================================
  # PRIMARY ENGINE: MAXIMUM SECONDARY EIGENVALUE RATIO (LOG-GAP)
  # ==========================================================================
  
  valid_search_space <- eig_vals[eig_vals > 1e-8]
  n_valid_search <- length(valid_search_space)
  
  if (n_valid_search >= 3) {
    # Transform to Log-Space to evaluate the relative Ratio (percentage drop)
    all_gaps <- abs(diff(log(valid_search_space)))
    
    # --- EDGE CASE: Ultra-low rank (Macro gap boundary severely truncates search space) ---
    if (!is.na(top_of_bulk_idx) && top_of_bulk_idx < 4) {
      ordered_gaps <- order(all_gaps, decreasing = TRUE)
      second_largest_gap_idx <- ordered_gaps[2]
      
      # >= ensures any gap at or after the start of the noise is bypassed
      if (second_largest_gap_idx >= top_of_bulk_idx) {
        K_elbow <- max(1, top_of_bulk_idx - 1)
        elbow_method <- "Macro Gap Boundary - 1"
      } else {
        K_elbow <- second_largest_gap_idx
        elbow_method <- "Second Largest Eigenvalue Ratio"
      }
    } else {
      # --- STANDARD ENGINE (Maximum Secondary Eigenvalue Ratio) ---
      # Bound the search space by the top of the bulk if it exists
      if (!is.na(top_of_bulk_idx) && top_of_bulk_idx > 2) {
        search_limit <- top_of_bulk_idx - 1
        secondary_gaps <- all_gaps[2:search_limit]
      } else {
        secondary_gaps <- all_gaps[2:length(all_gaps)] 
      }
      K_elbow <- which.max(secondary_gaps) + 1
      
      # Final safety clamp
      if (!is.na(top_of_bulk_idx) && K_elbow >= top_of_bulk_idx) {
        K_elbow <- max(1, top_of_bulk_idx - 1)
      }
      elbow_method <- "Maximum Secondary Eigenvalue Ratio"
    }
  } else {
    K_elbow <- max(1, valid_k)
    elbow_method <- "Kaiser Criterion (> Mean Trace)"
  }
  
  if (K_elbow < length(eig_vals)) {
    val_current <- eig_vals[K_elbow]
    val_next <- eig_vals[K_elbow + 1]
    if (val_next > 1e-12) {
      elbow_ratio <- val_current / val_next
      elbow_pct_drop <- (1 - (val_next / val_current)) * 100
    }
  }
  
  # ============================================================================
  # WAVE 1: INITIAL OUTPUT & METRICS
  # ============================================================================
  pct_prob_volume <- (R_eff / m_valid) * 100
  pct_redundant_signal <- (1 - (R_eff / m_valid)) * 100
  redundant_signal_volume <- m_valid - R_eff
  
  n_eigen_gt_mean <- valid_k
  n_eigen_le_mean <- n_total - valid_k
  
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    # ZOOMED VIEW
    zoom_start <- max(1, K_elbow - 5)
    zoom_end <- min(length(eig_vals), K_elbow + 15)
    plot_df_zoom <- data.frame(Rank = zoom_start:zoom_end, Eigenvalue = eig_vals[zoom_start:zoom_end])
    
    p_scree_zoom <- ggplot2::ggplot(plot_df_zoom, ggplot2::aes(x = Rank, y = Eigenvalue)) +
      ggplot2::geom_line(color = "dodgerblue", linewidth = 1) +
      ggplot2::geom_point(color = "dodgerblue", size = 2) +
      ggplot2::geom_vline(xintercept = K_elbow, color = "#D55E00", linetype = "dashed", linewidth = 1.2) +
      ggplot2::scale_y_continuous(trans = 'log10') +
      ggplot2::scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(min(x), max(x)))))) +
      ggplot2::labs(title = "Zoomed View", x = "Eigenvalue Index", y = "Log(Eigenvalue)") +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 12))
    
    # MACRO VIEW (Dynamically extended to show top of noise bulk)
    macro_base <- max(50, K_elbow * 10)
    if (!is.na(top_of_bulk_idx)) {
      macro_base <- max(macro_base, top_of_bulk_idx + 25) # 25 index visual cushion
    }
    macro_end <- min(length(eig_vals), macro_base)
    
    plot_df_macro <- data.frame(Rank = 1:macro_end, Eigenvalue = eig_vals[1:macro_end])
    macro_y_max <- if(length(eig_vals) >= 2) eig_vals[2] * 1.1 else max(eig_vals)
    
    p_scree_macro <- ggplot2::ggplot(plot_df_macro, ggplot2::aes(x = Rank, y = Eigenvalue)) +
      ggplot2::geom_line(color = "dodgerblue", linewidth = 1) +
      ggplot2::geom_point(color = "dodgerblue", size = 2) +
      ggplot2::geom_vline(xintercept = K_elbow, color = "#D55E00", linetype = "dashed", linewidth = 1.2) +
      ggplot2::scale_y_continuous(trans = 'log10') +
      ggplot2::coord_cartesian(ylim = c(NA, macro_y_max)) + 
      ggplot2::labs(title = "Macro View", x = "Eigenvalue Index", y = "Log(Eigenvalue)") +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 12))
    
    if (requireNamespace("patchwork", quietly = TRUE)) {
      combined_plot <- (p_scree_macro + p_scree_zoom) +
        patchwork::plot_annotation(
          title = "Entropic Scree Results",
          subtitle = sprintf("Automated Elbow = %d", K_elbow),
          theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 16, hjust = 0.5),
                                 plot.subtitle = ggplot2::element_text(size = 14, hjust = 0.5))
        )
      print(combined_plot)
    } else {
      print(p_scree_macro)
      print(p_scree_zoom)
    }
  }
  
  cat("\n=================================================================\n")
  cat(" STRUCTURAL COMPOSITION\n")
  cat("=================================================================\n")
  cat(sprintf(" -> %-50s : %d\n", "Valid Variables (m)", m_valid))
  cat(sprintf(" -> %-50s : %.3f\n", "Centered Trace (Tr_Mc)", Tr_Mc))
  cat(sprintf(" -> %-50s : %.3f%%\n", "Synergistic Curvature Deficit Ratio (SCDR)", SCDR))
  cat("-----------------------------------------------------------------\n")
  cat(sprintf(" -> %-50s : %.3f\n", "Total Unique Probabilistic Volume (R_eff)", R_eff))
  cat(sprintf(" -> %-50s : %.3f%%\n", "%", pct_prob_volume))
  cat("      (Unique Signal + Structural Uncertainty + \u22A5 Measurement Error)\n")
  cat("-----------------------------------------------------------------\n")
  cat(sprintf(" -> %-50s : %.3f\n", "Redundant Signal Volume (Tr_Mc - R_eff)", redundant_signal_volume))
  cat(sprintf(" -> %-50s : %.3f%%\n", "%", pct_redundant_signal))
  cat("-----------------------------------------------------------------\n")
  cat(sprintf(" -> %-50s : %d\n", "Eigenvalues > Mean Trace", n_eigen_gt_mean))
  cat(sprintf(" -> %-50s : %d\n", "Eigenvalues <= Mean Trace", n_eigen_le_mean))
  cat("=================================================================\n\n")
  
  cat("=================================================================\n")
  cat(" AUTOMATED ELBOW DETECTION (HEURISTIC)\n")
  cat("=================================================================\n")
  cat(sprintf(" -> %-43s : %d\n", "Automated Extracted Elbow Rank (K_elbow)", K_elbow))
  cat(sprintf(" -> %-43s : %s\n", "Extraction Method Tripped", elbow_method))
  if (!is.na(elbow_ratio)) {
    cat(sprintf(" -> %-43s : %.2fx (%.1f%% Drop)\n", "Elbow Magnitude (Topological Variance Drop)", elbow_ratio, elbow_pct_drop))
  }
  if (!is.na(macro_gap_ratio)) {
    cat("-----------------------------------------------------------------\n")
    cat(" [Diagnostic: Macro Gap (Noise Cliff)]\n")
    cat(sprintf(" -> %-43s : %d\n", "Identified Top of Noise Bulk (Index)", top_of_bulk_idx)) 
    cat(sprintf(" -> %-43s : %.6f\n", "Macro Gap Baseline (Max Noise Gap)", macro_max_noise_gap))
    cat(sprintf(" -> %-43s : %.2fx Baseline\n", "Actual Macro Gap Magnitude", macro_gap_ratio))
  }
  cat("=================================================================\n")
  
  # ============================================================================
  # WAVE 2: INTERACTIVE USER OVERRIDE
  # ============================================================================
  K_final <- K_elbow
  
  if (interactive_mode) {
    cat("\n[WARNING]: The automated elbow extractor relies on statistical heuristics and\n")
    cat("may not perfectly align with the true structural elbow of your specific dataset.\n")
    cat("Please visually examine the generated scree plot.\n\n")
    
    first_prompt <- TRUE
    while (TRUE) {
      if (first_prompt) {
        prompt_msg <- sprintf("Do you want to keep the Extracted Elbow Rank of %d? (Type 'Y' to keep, or enter custom rank): ", K_final)
      } else {
        prompt_msg <- sprintf("Do you want to keep the updated rank of %d? (Type 'Y' to finalize, or enter a new custom rank): ", K_final)
      }
      
      ans <- trimws(readline(prompt = prompt_msg))
      if (tolower(ans) %in% c("y", "yes")) {
        cat("\n[+] Finalizing rank selection.\n")
        break
      } else {
        parsed_k <- suppressWarnings(as.integer(ans))
        if (!is.na(parsed_k) && parsed_k > 0 && parsed_k <= m_valid) {
          K_final <- parsed_k
          cat(sprintf("\n[+] Rank manually updated to %d.\n", K_final))
          
          # ====================================================================
          # INTERACTIVE GRAPH PREVIEW & METRICS
          # ====================================================================
          if (requireNamespace("ggplot2", quietly = TRUE)) {
            # 1. ZOOMED VIEW
            zoom_start_upd <- max(1, min(K_elbow, K_final) - 5)
            zoom_end_upd <- min(length(eig_vals), max(K_elbow, K_final) + 15)
            
            plot_df_zoom_upd <- data.frame(
              Rank = zoom_start_upd:zoom_end_upd, 
              Eigenvalue = eig_vals[zoom_start_upd:zoom_end_upd]
            )
            
            p_scree_zoom_upd <- ggplot2::ggplot(plot_df_zoom_upd, ggplot2::aes(x = Rank, y = Eigenvalue)) +
              ggplot2::geom_line(color = "dodgerblue", linewidth = 1) +
              ggplot2::geom_point(color = "dodgerblue", size = 2) +
              ggplot2::geom_vline(xintercept = K_elbow, color = "gray60", linetype = "dashed", linewidth = 1) +
              ggplot2::geom_vline(xintercept = K_final, color = "forestgreen", linetype = "solid", linewidth = 1.2) +
              ggplot2::scale_y_continuous(trans = 'log10') +
              ggplot2::scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(min(x), max(x)))))) +
              ggplot2::labs(
                title = "Zoomed View",
                x = "Eigenvalue Index (m)", 
                y = "Log(Eigenvalue)"
              ) +
              ggplot2::theme_minimal(base_size = 14) +
              ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 12))
            
            # 2. MACRO VIEW (Dynamically extended)
            macro_base_upd <- max(50, K_final * 10, K_elbow * 10)
            if (!is.na(top_of_bulk_idx)) {
              macro_base_upd <- max(macro_base_upd, top_of_bulk_idx + 25) # 25 index visual cushion
            }
            macro_end_upd <- min(length(eig_vals), macro_base_upd)
            
            plot_df_macro_upd <- data.frame(
              Rank = 1:macro_end_upd,
              Eigenvalue = eig_vals[1:macro_end_upd]
            )
            
            macro_y_max <- if(length(eig_vals) >= 2) eig_vals[2] * 1.1 else max(eig_vals)
            
            p_scree_macro_upd <- ggplot2::ggplot(plot_df_macro_upd, ggplot2::aes(x = Rank, y = Eigenvalue)) +
              ggplot2::geom_line(color = "dodgerblue", linewidth = 1) +
              ggplot2::geom_point(color = "dodgerblue", size = 2) +
              ggplot2::geom_vline(xintercept = K_elbow, color = "gray60", linetype = "dashed", linewidth = 1) +
              ggplot2::geom_vline(xintercept = K_final, color = "forestgreen", linetype = "solid", linewidth = 1.2) +
              ggplot2::scale_y_continuous(trans = 'log10') +
              ggplot2::coord_cartesian(ylim = c(NA, macro_y_max)) + 
              ggplot2::labs(
                title = "Macro View",
                x = "Eigenvalue Index (m)", 
                y = "Log(Eigenvalue)"
              ) +
              ggplot2::theme_minimal(base_size = 14) +
              ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 12))
            
            # Render Side-by-Side with global titles
            if (requireNamespace("patchwork", quietly = TRUE)) {
              combined_plot_upd <- (p_scree_macro_upd + p_scree_zoom_upd) +
                patchwork::plot_annotation(
                  title = "Entropic Scree Results",
                  subtitle = sprintf("User Confirmed Elbow = %d (Auto: %d)", K_final, K_elbow),
                  theme = ggplot2::theme(
                    plot.title = ggplot2::element_text(face = "bold", size = 16, hjust = 0.5),
                    plot.subtitle = ggplot2::element_text(size = 14, hjust = 0.5)
                  )
                )
              print(combined_plot_upd)
            } else {
              print(p_scree_macro_upd)
              print(p_scree_zoom_upd)
            }
          }
          
          # --- PREVIEW GRAVITY CALCULATIONS ---
          signal_variance_prev <- sum(eig_vals[1:K_final])
          signal_weight_prev <- signal_variance_prev / m_plus
          unique_signal_volume_prev <- R_eff * signal_weight_prev
          total_signal_volume_prev <- unique_signal_volume_prev + redundant_signal_volume
          
          AIG_prev <- total_signal_volume_prev / K_final
          core_eigenvals_prev <- eig_vals[1:K_final]
          p_core_prev <- core_eigenvals_prev / sum(core_eigenvals_prev)
          FSIG_prev <- p_core_prev * total_signal_volume_prev
          
          cat("\n=================================================================\n")
          cat(sprintf(" (PREVIEW) ELBOW LATENT METRICS FOR K_elbow = %d\n", K_final))
          cat("=================================================================\n")
          cat(sprintf(" -> Preview Rank (K_elbow)      : %d\n", K_final))
          cat(sprintf(" -> Avg Info Gravity (AIG)      : %.3f\n", AIG_prev))
          cat(" -> Factor-Specific Informational Gravity (FSIG):\n")
          print(round(FSIG_prev, 3))
          cat("=================================================================\n\n")
          
          first_prompt <- FALSE
        } else {
          cat("[-] Invalid input. Please enter 'Y' to finalize, or a valid positive integer.\n\n")
        }
      }
    }
  }
  
  # --- FINAL GRAVITY CALCULATIONS ---
  signal_variance <- sum(eig_vals[1:K_final])
  signal_weight <- signal_variance / m_plus
  unique_signal_volume <- R_eff * signal_weight
  total_signal_volume <- unique_signal_volume + redundant_signal_volume
  
  AIG <- total_signal_volume / K_final
  core_eigenvals <- eig_vals[1:K_final]
  p_core <- core_eigenvals / sum(core_eigenvals)
  FSIG_final <- p_core * total_signal_volume
  
  cat("\n=================================================================\n")
  cat(" (FINAL) ELBOW LATENT METRICS (based user-confirmed K_elbow)\n")
  cat("=================================================================\n")
  cat(sprintf(" -> Final Retained Rank (K_elbow)     : %d\n", K_final))
  cat(sprintf(" -> Avg Info Gravity (AIG)            : %.3f\n", AIG))
  cat(" -> Factor-Specific Informational Gravity (FSIG):\n")
  print(round(FSIG_final, 3))
  cat("=================================================================\n\n")
  
  # ============================================================================
  # WAVE 3: FINAL TRIPARTITE STRUCTURAL COMPOSITION
  # ============================================================================
  noise_variance <- if (K_final < m_valid) sum(eig_vals[(K_final + 1):m_valid]) else 0
  noise_weight <- noise_variance / m_plus
  idiosyncratic_noise_volume <- R_eff * noise_weight
  
  cat("=================================================================\n")
  cat(sprintf(" (FINAL) TRIPARTITE STRUCTURAL COMPOSITION (based user-confirmed K_elbow = %d)\n", K_final))
  cat("=================================================================\n")
  cat(sprintf(" -> %-50s : %d\n", "Valid Variables (m)", m_valid))
  cat(sprintf(" -> %-50s : %.2f\n", "Centered Trace (Tr_Mc)", Tr_Mc))
  cat("-----------------------------------------------------------------\n")
  cat(sprintf(" -> %-50s : %.3f\n", "Total Signal Volume (Unique + Redundant)", total_signal_volume))
  cat(sprintf(" -> %-40s : %.3f\n", "   (Unique Signal Volume)", unique_signal_volume))
  cat(sprintf(" -> %-40s : %.3f\n", "   (Redundant Signal Volume)", redundant_signal_volume))
  cat("-----------------------------\n")
  cat(sprintf(" -> %-50s : %.3f\n", "Idiosyncratic Noise Volume", idiosyncratic_noise_volume))
  cat("      (Structural Uncertainty + \u22A5 Measurement Error)\n")
  cat("-----------------------------------------------------------------\n")
  cat(sprintf(" -> %-50s : %.3f%%\n", "% Total Signal", (total_signal_volume / m_valid) * 100))
  cat(sprintf(" -> %-40s : %.3f%%\n", "   (% Unique Signal)", (unique_signal_volume / m_valid) * 100))
  cat(sprintf(" -> %-40s : %.3f%%\n", "   (% Redundant Signal)", (redundant_signal_volume / m_valid) * 100))
  cat("-----------------------------\n")
  cat(sprintf(" -> %-50s : %.3f%%\n", "% Idiosyncratic Noise", (idiosyncratic_noise_volume / m_valid) * 100))
  cat("      (Structural Uncertainty + \u22A5 Measurement Error)\n")
  cat("=================================================================\n\n")
  
  # ============================================================================
  # EXTENDED Factor-Specific Informational Gravity (FSIG) CALCULATIONS
  # ============================================================================
  
  # --- EXTENDED MODEL A: MACRO BULK BOUNDARY ---
  top_bulk_safe <- if(!is.na(top_of_bulk_idx)) top_of_bulk_idx else valid_k
  extended_bulk_k <- max(K_final, top_bulk_safe - 1)
  extended_eigenvals_bulk <- eig_vals[1:extended_bulk_k]
  p_extended_bulk <- extended_eigenvals_bulk / sum(extended_eigenvals_bulk)
  sig_var_bulk <- sum(extended_eigenvals_bulk)
  total_sig_vol_bulk <- (R_eff * (sig_var_bulk / m_plus)) + redundant_signal_volume
  FSIG_extended_bulk <- p_extended_bulk * total_sig_vol_bulk
  
  # --- EXTENDED MODEL B: KAISER RULE BOUNDARY ---
  extended_kaiser_k <- max(K_final, valid_k)
  extended_eigenvals_kaiser <- eig_vals[1:extended_kaiser_k]
  p_extended_kaiser <- extended_eigenvals_kaiser / sum(extended_eigenvals_kaiser)
  sig_var_kaiser <- sum(extended_eigenvals_kaiser)
  total_sig_vol_kaiser <- (R_eff * (sig_var_kaiser / m_plus)) + redundant_signal_volume
  FSIG_extended_kaiser <- p_extended_kaiser * total_sig_vol_kaiser
  
  return(list(
    eigenvalues = eig_vals,
    similarity_matrix = NMI_mat_c,
    retained_features = valid_vars,
    bin_distributions = bin_sample_sizes,
    R_eff = R_eff,
    K_auto_extracted = K_elbow,
    extraction_method = elbow_method,
    K_final = K_final,
    top_of_bulk = top_bulk_safe,
    total_signal_volume = total_signal_volume,
    unique_signal_volume = unique_signal_volume,
    redundant_signal_volume = redundant_signal_volume, 
    idiosyncratic_noise_volume = idiosyncratic_noise_volume,
    AIG = AIG,
    FSIG_final = FSIG_final,
    FSIG_extended_bulk = FSIG_extended_bulk,
    FSIG_extended_kaiser = FSIG_extended_kaiser,
    eigenvectors = eigen_res$vectors
  ))
}

################################################################################
################################################################################
# SIMULATION TO TEST ENTROPIC SCREE FUNCTION 
################################################################################
################################################################################

# ==============================================================================
# 2. DATA GENERATION & MEASUREMENT ERROR LOGIC
# ==============================================================================

generate_random_corr_matrix <- function(k) {
  mat <- matrix(rnorm(k * k), nrow = k)
  cov_mat <- crossprod(mat)
  return(cov2cor(cov_mat))
}

# STEP 2A: Generate the pure, uncorrupted Ground Truth Proxies
generate_true_mixed_proxies <- function(s1_continuous, m_proxies, max_interaction_order = 3, max_polynomial_order = 3, int_scaling = 1, continuous_ratio = 1.0) {
  k <- ncol(s1_continuous)
  n <- nrow(s1_continuous)
  
  # Allow dynamic control over the continuous vs binary split
  m_cont <- floor(m_proxies * continuous_ratio)
  m_bin <- m_proxies - m_cont
  
  # 1. Expand the Continuous Latent Space
  s1_df <- as.data.frame(s1_continuous)
  colnames(s1_df) <- paste0("X", 1:k)
  
  # ENGINE 1: Cross-interactions (Capped at k)
  eff_cross_order <- min(k, max_interaction_order)
  formula_str <- if(eff_cross_order == 1) as.formula("~ .") else as.formula(paste0("~ .^", eff_cross_order))
  
  cat(sprintf("      -> Expanding continuous roots up to interaction order %d and polynomial order %d...\n", eff_cross_order, max_polynomial_order))
  design_mat <- model.matrix(formula_str, data = s1_df)[, -1, drop = FALSE]
  
  # ENGINE 2: Pure Polynomial Powers (Uncapped)
  if (max_polynomial_order > 1) {
    power_list <- list()
    for (i in 1:k) {
      for (p in 2:max_polynomial_order) {
        power_col <- s1_df[[i]]^p
        col_name <- paste(rep(paste0("X", i), p), collapse = ":")
        power_list[[col_name]] <- power_col
      }
    }
    if (length(power_list) > 0) {
      design_mat <- cbind(design_mat, as.matrix(as.data.frame(power_list)))
    }
  }
  
  # --- STANDARDIZE THE ENTIRE DESIGN MATRIX ---
  # Every term (main effects and all interactions/powers) now has Variance = 1.0
  design_mat <- scale(design_mat)
  # -------------------------------------------------
  
  term_names <- colnames(design_mat)
  interaction_orders <- stringr::str_count(term_names, ":") + 1
  n_terms <- ncol(design_mat)
  
  # --- THE DIALS ---
  # To use a factorial penalty, swap the comments on the next two lines. This makes the underlying generation more linear:
  # term_sds <- sqrt(1 / factorial(interaction_orders))
  term_sds <- rep(1.0, n_terms) 
  # int_scaling is passed as a function argument
  # -----------------
  
  # 2. Build the Weight Matrix
  coeffs <- matrix(0, nrow = n_terms, ncol = m_proxies)
  is_main <- interaction_orders == 1
  is_int <- interaction_orders > 1
  n_main <- sum(is_main)
  n_int <- sum(is_int)
  
  for(j in 1:m_proxies) {
    # Main effects
    if (n_main > 0) coeffs[is_main, j] <- rnorm(n_main, mean = 0, sd = term_sds[is_main])
    # Interactions and Polynomials
    if (n_int > 0)  coeffs[is_int, j]  <- rnorm(n_int, mean = 0, sd = term_sds[is_int] * int_scaling)
  }
  mask <- matrix(rbinom(n_terms * m_proxies, 1, 0.05), nrow = n_terms, ncol = m_proxies)
  coeffs <- coeffs * mask
  
  # 3. Generate the Raw Structural Signal
  raw_signal <- design_mat %*% coeffs
  
  # 4A. The True Continuous Signal
  cat(sprintf("      -> Generating %d True Continuous Proxies...\n", m_cont))
  true_cont <- raw_signal[, 1:m_cont, drop = FALSE]
  
  # 4B. The True Binary Signal (SAFELY BYPASSED IF 0)
  if (m_bin > 0) {
    cat(sprintf("      -> Generating %d True Binary Proxies...\n", m_bin))
    signal_bin <- raw_signal[, (m_cont + 1):m_proxies, drop = FALSE]
    
    apply_copula_mapping <- function(scores) {
      if (sd(scores) < 1e-9) return(rep(0.5, length(scores)))
      z_scores <- as.vector(scale(scores))
      probs <- pnorm(z_scores) 
      if (runif(1) > 0.5) probs <- 1 - probs 
      return(probs)
    }
    
    prob_mat <- apply(signal_bin, 2, apply_copula_mapping)
    true_bin <- matrix(rbinom(length(prob_mat), 1, prob_mat), nrow = n, ncol = m_bin)
    
    # Combine mixed types
    true_proxies <- cbind(true_cont, true_bin)
    is_continuous <- c(rep(TRUE, m_cont), rep(FALSE, m_bin))
  } else {
    # Pure continuous universe
    true_proxies <- true_cont
    is_continuous <- rep(TRUE, m_cont)
  }
  
  # 5. Randomly shuffle the columns
  mix_idx <- sample(m_proxies)
  true_proxies <- true_proxies[, mix_idx, drop = FALSE]
  is_continuous <- is_continuous[mix_idx]
  
  # Calculate K_rlzd before returning
  active_terms <- sum(rowSums(abs(coeffs)) > 0)
  
  return(list(
    data_matrix = true_proxies,
    is_continuous = is_continuous,
    active_terms = active_terms
  ))
}

# STEP 2B: Apply independent Measurement Error to the True Data
apply_measurement_error <- function(true_universe, snr_continuous = 2.0, binary_error_rate = 0.15) {
  cat(sprintf("      -> Applying Measurement Error (Continuous SNR = %.2f, Binary Bit-Flip Rate = %.3f)...\n", snr_continuous, binary_error_rate))
  
  obs_mat <- true_universe$data_matrix
  is_cont <- true_universe$is_continuous
  n <- nrow(obs_mat)
  m <- ncol(obs_mat)
  
  for (j in 1:m) {
    if (is_cont[j]) {
      # Add Gaussian Noise mapped to target SNR
      true_var <- var(obs_mat[, j])
      if (true_var < 1e-9) true_var <- 1e-9 
      
      noise <- rnorm(n, mean = 0, sd = 1)
      noise_var <- var(noise)
      
      scaling_factor <- sqrt(true_var / (noise_var * snr_continuous))
      obs_mat[, j] <- obs_mat[, j] + (noise * scaling_factor)
      
    } else {
      # Add Bit-Flip Measurement Error (Misreading the true state)
      # Flips 1s to 0s, and 0s to 1s with probability = binary_error_rate
      flip_mask <- rbinom(n, 1, binary_error_rate)
      obs_mat[, j] <- abs(obs_mat[, j] - flip_mask) # Equivalent to XOR
    }
  }
  
  obs_dt <- data.table::as.data.table(obs_mat)
  data.table::setnames(obs_dt, paste0("V", 1:m))
  return(obs_dt)
}

# ==============================================================================
# 3. SIMULATE DATA
# ==============================================================================
set.seed(19862026)

K_TRUE <- 10
N_ROWS <- 5000
M_PROXIES <- 10000
CONTINUOUS_RATIO <- .80  # SET TO 1.0 FOR PURE CONTINUOUS, 0.0 FOR PURE BINARY, OR ANYWHERE IN BETWEEN

# --- DIALS ---
MAX_INTERACTION_ORDER <- round(K_TRUE/2)    # Controls cross-interactions (e.g., X1:X2). Maxes out at K_TRUE.
MAX_POLYNOMIAL_ORDER <- 4    # Controls pure powers (e.g., X1^2, X1^3). Uncapped.
# ----------------------

# ERROR KNOBS
CONTINUOUS_SNR <- 2  # Lower is dirtier (e.g., 0.5 is garbage, 10 is clean)
BINARY_ERROR_RATE <- 0.15   # 0.0 is perfect sensor, 0.50 is pure static coin-flip

cat(sprintf("Generating Mixed Synthetic Universe: %s rows, %d Proxies, %d Latent Drivers...\n", format(N_ROWS, big.mark=","), M_PROXIES, K_TRUE))

# 1. Generate S(1) Latent Space
medium_corr_matrix <- generate_random_corr_matrix(K_TRUE)
identity_matrix <- diag(K_TRUE)
Z_latent_continuous <- MASS::mvrnorm(n = N_ROWS
                                     , mu = rep(0, K_TRUE)
                                     , Sigma = identity_matrix # Change to medium_corr_matrix for correlated states
)

# 2. Calculate Theoretical Ceiling
s1_df <- as.data.frame(Z_latent_continuous)
colnames(s1_df) <- paste0("X", 1:K_TRUE)

# ENGINE 1: Cross-interactions (Capped at K_TRUE)
eff_cross_order <- min(K_TRUE, MAX_INTERACTION_ORDER)
formula_str <- as.formula(paste0("~ .^", eff_cross_order))
design_mat <- model.matrix(formula_str, data = s1_df)[, -1, drop = FALSE]

# ENGINE 2: Pure Polynomial Powers (Uncapped)
if (MAX_POLYNOMIAL_ORDER > 1) {
  power_list <- list()
  for (i in 1:K_TRUE) {
    for (p in 2:MAX_POLYNOMIAL_ORDER) {
      power_col <- s1_df[[i]]^p
      col_name <- paste(rep(paste0("X", i), p), collapse = ":")
      power_list[[col_name]] <- power_col
    }
  }
  if (length(power_list) > 0) {
    design_mat <- cbind(design_mat, as.matrix(as.data.frame(power_list)))
  }
}

# --- STANDARDIZE THE THEORETICAL MATRIX TO MATCH GENERATOR ---
design_mat <- scale(design_mat)
# ------------------------------------------------------------------

term_names <- colnames(design_mat)
interaction_orders <- stringr::str_count(term_names, ":") + 1

# --- THE DIALS ---
term_sds <- rep(1.0, length(term_names)) # Use sqrt(1 / factorial(interaction_orders)) instead to shrink sds for interactions to move the systemt towards linearity.
int_scaling <- 1.0 
# -----------------

term_weights <- ifelse(interaction_orders == 1, term_sds, term_sds * int_scaling)
muzzled_design_mat <- t(t(design_mat) * term_weights)

expanded_cov <- cov(muzzled_design_mat)
expanded_eigen <- pmax(eigen(expanded_cov, symmetric = TRUE)$values, 0)
p_expanded <- expanded_eigen[expanded_eigen > 1e-9] / sum(expanded_eigen[expanded_eigen > 1e-9])
true_continuous_ceiling <- exp(-sum(p_expanded * log(p_expanded)))

cat(sprintf("\n[***] R_conf (Effective Latent Configurational Rank): %.2f\n\n", true_continuous_ceiling))

# 3. Generate the Data (Uncoupled Truth and Error)
true_universe <- generate_true_mixed_proxies(
  Z_latent_continuous, 
  M_PROXIES, 
  max_interaction_order = MAX_INTERACTION_ORDER, 
  max_polynomial_order = MAX_POLYNOMIAL_ORDER,
  int_scaling = int_scaling,
  continuous_ratio = CONTINUOUS_RATIO
)

# Print the K_rlzd extracted safely from the function
cat(sprintf("[***] K_rlzd (Realized Latent Configurational Rank): %d\n", true_universe$active_terms))

observed_data <- apply_measurement_error(true_universe, snr_continuous = CONTINUOUS_SNR, binary_error_rate = BINARY_ERROR_RATE)

cat("\nDataset is ready. Starting pipeline...\n\n")

# ==============================================================================
# 4. RUN ENTROPIC SCREE
# ==============================================================================
results <- calculate_entropic_scree(observed_data
                                    , purge_constants = FALSE
                                    , check_collinearity = FALSE
                                    , interactive_mode = FALSE)

# ==============================================================================
# 5. STANDARD, SPEARMAN, & KERNEL PCA EXTRACTION (FOR COMPARISON)
# ==============================================================================
m_total <- ncol(observed_data)
n_rows <- nrow(observed_data)

# --- 5A. STANDARD PCA ---
cat("\nExtracting Standard PCA for comparison...\n")
start_pca <- Sys.time()
# Use prcomp with scaling to mirror a Pearson Correlation Matrix extraction
pca_res <- prcomp(observed_data, center = TRUE, scale. = TRUE)
pca_eigenvalues <- pca_res$sdev^2
pca_time <- round(as.numeric(difftime(Sys.time(), start_pca, units = "secs")), 2)
cat(sprintf("Standard PCA completed in %.2f seconds.\n", pca_time))
pca_kaiser <- sum(pca_eigenvalues > 1.0)

# --- 5B. SPEARMAN RANK PCA ---
cat("\nExtracting Spearman Rank PCA for comparison...\n")
start_spearman <- Sys.time()
# Apply dense rank to mimic Spearman correlation structure
ranked_data <- apply(observed_data, 2, data.table::frank)
spearman_res <- prcomp(ranked_data, center = TRUE, scale. = TRUE)
spearman_eigenvalues <- spearman_res$sdev^2
spearman_time <- round(as.numeric(difftime(Sys.time(), start_spearman, units = "secs")), 2)
cat(sprintf("Spearman Rank PCA completed in %.2f seconds.\n", spearman_time))
rm(ranked_data); gc(verbose = FALSE) # Clear memory
spearman_kaiser <- sum(spearman_eigenvalues > 1.0)

# --- 5C. KERNEL PCA (RBF / GAUSSIAN) ---
cat("\nExtracting Kernel PCA (RBF) for comparison...\n")
start_kpca <- Sys.time()
# Fast squared distance and Kernel matrix computation to protect RAM
X_scaled <- scale(as.matrix(observed_data))
X_norm <- rowSums(X_scaled^2)
D2 <- outer(X_norm, X_norm, "+") - 2 * tcrossprod(X_scaled)
gamma <- 1 / m_total
K_mat <- exp(-gamma * D2)
rm(X_scaled, D2); gc(verbose = FALSE)

# Double center the Kernel matrix
one_n <- matrix(1/n_rows, n_rows, n_rows)
K_c <- K_mat - one_n %*% K_mat - K_mat %*% one_n + one_n %*% K_mat %*% one_n
rm(K_mat, one_n); gc(verbose = FALSE)

# Extract Eigenvalues
kpca_eigenvalues <- eigen(K_c, symmetric = TRUE, only.values = TRUE)$values
kpca_eigenvalues <- pmax(kpca_eigenvalues, 0) # Clip numeric artifacts
kpca_time <- round(as.numeric(difftime(Sys.time(), start_kpca, units = "secs")), 2)
cat(sprintf("Kernel PCA (RBF) completed in %.2f seconds.\n", kpca_time))
rm(K_c); gc(verbose = FALSE)

# --- 5D. BUILD COMPARISON DATAFRAME ---
# Pad N-length vectors with zeros to match m_total
pca_padded <- c(pca_eigenvalues, rep(0, m_total - length(pca_eigenvalues)))
spearman_padded <- c(spearman_eigenvalues, rep(0, m_total - length(spearman_eigenvalues)))
kpca_padded <- c(kpca_eigenvalues, rep(0, m_total - length(kpca_eigenvalues)))

df_compare <- data.frame(
  Rank = rep(1:m_total, 4),
  Eigenvalue = c(pca_padded, spearman_padded, kpca_padded, results$eigenvalues),
  Method = factor(rep(c("Standard PCA", "Spearman PCA", "Kernel PCA", "Entropic Scree"), each = m_total),
                  levels = c("Standard PCA", "Spearman PCA", "Kernel PCA", "Entropic Scree"))
)

# ==============================================================================
# 6. SIDE-BY-SIDE VISUAL PROOF (4-PANEL GRID)
# ==============================================================================
if (requireNamespace("ggplot2", quietly = TRUE)) {
  
  # Calculate caps for visual clarity
  pca_y_max <- if(length(pca_eigenvalues) >= 2) pca_eigenvalues[2] * 1.1 else max(pca_eigenvalues)
  spear_y_max <- if(length(spearman_eigenvalues) >= 2) spearman_eigenvalues[2] * 1.1 else max(spearman_eigenvalues)
  kpca_y_max <- if(length(kpca_eigenvalues) >= 2) kpca_eigenvalues[2] * 1.1 else max(kpca_eigenvalues)
  ent_y_max <- if(length(results$eigenvalues) >= 2) results$eigenvalues[2] * 1.1 else max(results$eigenvalues)
  
  ent_y_min <- min(results$eigenvalues[results$eigenvalues > 1e-8])
  pca_y_min <- 0 
  kpca_y_min <- 0
  
  # Dynamic log breaks generator
  log10_breaks <- function(x) {
    10^seq(floor(log10(min(x))), ceiling(log10(max(x))))
  }
  
  # Extract K_rlzd and R_conf
  K_rlzd <- true_universe$active_terms
  R_conf <- true_continuous_ceiling 
  ref_size <- 3.0 # Slightly scaled down for 4-panel density
  
  # Helper function to generate staggered Y heights for annotations
  get_heights <- function(y_min, y_max) {
    list(
      h90 = y_min + (y_max - y_min) * 0.90,
      h75 = y_min + (y_max - y_min) * 0.75,
      h60 = y_min + (y_max - y_min) * 0.60,
      h52 = y_min + (y_max - y_min) * 0.525, # Midpoint between K_rlzd (h75) and Rank Ceiling (h30)
      h45 = y_min + (y_max - y_min) * 0.45,
      h30 = y_min + (y_max - y_min) * 0.30
    )
  }
  
  hpca <- get_heights(pca_y_min, pca_y_max)
  hspear <- get_heights(pca_y_min, spear_y_max)
  hkpca <- get_heights(kpca_y_min, kpca_y_max)
  hent <- get_heights(ent_y_min, ent_y_max)
  
  # Split the data
  pca_data <- df_compare[df_compare$Method == "Standard PCA", ]
  spear_data <- df_compare[df_compare$Method == "Spearman PCA", ]
  kpca_data <- df_compare[df_compare$Method == "Kernel PCA", ]
  ent_data <- df_compare[df_compare$Method == "Entropic Scree" & df_compare$Rank < M_PROXIES, ]
  
  # ---------------------------------------------------------
  # PANEL 1: STANDARD PCA
  # ---------------------------------------------------------
  p_pca <- ggplot2::ggplot(pca_data, ggplot2::aes(x = Rank, y = Eigenvalue)) +
    ggplot2::geom_line(data = pca_data[pca_data$Rank < n_rows, ], color = "firebrick", linewidth = 1) +
    ggplot2::geom_line(data = pca_data[pca_data$Rank >= (n_rows - 1), ], color = "gray60", linetype = "dashed", linewidth = 1) +
    ggplot2::geom_vline(xintercept = K_TRUE, color = "forestgreen", linetype = "solid", linewidth = 1.2) +
    ggplot2::annotate("text", x = K_TRUE, y = hpca$h90, label = sprintf(" r \n (%d)", K_TRUE), hjust = -0.1, color = "forestgreen", fontface = "bold", size = ref_size) +
    ggplot2::geom_vline(xintercept = K_rlzd, color = "purple", linetype = "dotdash", linewidth = 1) +
    ggplot2::annotate("text", x = K_rlzd, y = hpca$h75, label = sprintf(" K_rlzd \n (%d)", K_rlzd), hjust = 1.05, color = "purple", fontface = "italic", size = ref_size) +
    ggplot2::geom_vline(xintercept = pca_kaiser, color = "darkorange", linetype = "dashed", linewidth = 1) +
    ggplot2::annotate("text", x = pca_kaiser, y = hpca$h52, label = sprintf("Kaiser, E>1 \n(%d) ", pca_kaiser), hjust = 1.05, color = "darkorange", fontface = "italic", size = ref_size) +
    ggplot2::geom_vline(xintercept = (n_rows - 1), color = "black", linetype = "dotted", linewidth = 1) +
    ggplot2::annotate("text", x = (n_rows - 1), y = hpca$h30, label = "Rank \nCeiling (N-1)", hjust = 1.05, size = ref_size) +
    ggplot2::scale_x_log10(breaks = log10_breaks) +
    ggplot2::coord_cartesian(ylim = c(pca_y_min, pca_y_max)) +
    ggplot2::labs(title = "Standard PCA", x = "Eigenvalue Index [Log Scale]", y = "Eigenvalue") +
    ggplot2::theme_minimal(base_size = 12) + ggplot2::theme(plot.title = ggplot2::element_text(face="bold"))
  
  # ---------------------------------------------------------
  # PANEL 2: SPEARMAN RANK PCA
  # ---------------------------------------------------------
  p_spear <- ggplot2::ggplot(spear_data, ggplot2::aes(x = Rank, y = Eigenvalue)) +
    ggplot2::geom_line(data = spear_data[spear_data$Rank < n_rows, ], color = "darkgoldenrod", linewidth = 1) +
    ggplot2::geom_line(data = spear_data[spear_data$Rank >= (n_rows - 1), ], color = "gray60", linetype = "dashed", linewidth = 1) +
    ggplot2::geom_vline(xintercept = K_TRUE, color = "forestgreen", linetype = "solid", linewidth = 1.2) +
    ggplot2::annotate("text", x = K_TRUE, y = hspear$h90, label = sprintf(" r \n (%d)", K_TRUE), hjust = -0.1, color = "forestgreen", fontface = "bold", size = ref_size) +
    ggplot2::geom_vline(xintercept = K_rlzd, color = "purple", linetype = "dotdash", linewidth = 1) +
    ggplot2::annotate("text", x = K_rlzd, y = hspear$h75, label = sprintf(" K_rlzd \n (%d)", K_rlzd), hjust = 1.05, color = "purple", fontface = "italic", size = ref_size) +
    ggplot2::geom_vline(xintercept = spearman_kaiser, color = "darkorange", linetype = "dashed", linewidth = 1) +
    ggplot2::annotate("text", x = spearman_kaiser, y = hspear$h52, label = sprintf("Kaiser, E>1 \n(%d) ", spearman_kaiser), hjust = 1.05, color = "darkorange", fontface = "italic", size = ref_size) +
    ggplot2::geom_vline(xintercept = (n_rows - 1), color = "black", linetype = "dotted", linewidth = 1) +
    ggplot2::annotate("text", x = (n_rows - 1), y = hspear$h30, label = "Rank \nCeiling (N-1)", hjust = 1.05, size = ref_size) +
    ggplot2::scale_x_log10(breaks = log10_breaks) +
    ggplot2::coord_cartesian(ylim = c(pca_y_min, spear_y_max)) +
    ggplot2::labs(title = "Spearman Rank PCA", x = "Eigenvalue Index [Log Scale]", y = "Eigenvalue") +
    ggplot2::theme_minimal(base_size = 12) + ggplot2::theme(plot.title = ggplot2::element_text(face="bold"))
  
  # ---------------------------------------------------------
  # PANEL 3: KERNEL PCA (RBF)
  # ---------------------------------------------------------
  # Calculate KPCA Kaiser (Eigenvalues > Mean Trace) for visual consistency
  kpca_kaiser <- sum(kpca_eigenvalues > mean(kpca_eigenvalues))
  
  p_kpca <- ggplot2::ggplot(kpca_data, ggplot2::aes(x = Rank, y = Eigenvalue)) +
    ggplot2::geom_line(data = kpca_data[kpca_data$Rank < n_rows, ], color = "mediumorchid4", linewidth = 1) +
    ggplot2::geom_line(data = kpca_data[kpca_data$Rank >= (n_rows - 1), ], color = "gray60", linetype = "dashed", linewidth = 1) +
    ggplot2::geom_vline(xintercept = K_TRUE, color = "forestgreen", linetype = "solid", linewidth = 1.2) +
    ggplot2::annotate("text", x = K_TRUE, y = hkpca$h90, label = sprintf(" r \n (%d)", K_TRUE), hjust = -0.1, color = "forestgreen", fontface = "bold", size = ref_size) +
    ggplot2::geom_vline(xintercept = K_rlzd, color = "purple", linetype = "dotdash", linewidth = 1) +
    ggplot2::annotate("text", x = K_rlzd, y = hkpca$h75, label = sprintf(" K_rlzd \n (%d)", K_rlzd), hjust = 1.05, color = "purple", fontface = "italic", size = ref_size) +
    ggplot2::geom_vline(xintercept = kpca_kaiser, color = "darkorange", linetype = "dashed", linewidth = 1) +
    ggplot2::annotate("text", x = kpca_kaiser, y = hkpca$h52, label = sprintf("Kaiser, E>Mean \n(%d) ", kpca_kaiser), hjust = 1.05, color = "darkorange", fontface = "italic", size = ref_size) +
    ggplot2::geom_vline(xintercept = (n_rows - 1), color = "black", linetype = "dotted", linewidth = 1) +
    ggplot2::annotate("text", x = (n_rows - 1), y = hkpca$h30, label = "Rank \nCeiling (N-1)", hjust = 1.05, size = ref_size) +
    ggplot2::scale_x_log10(breaks = log10_breaks) +
    ggplot2::coord_cartesian(ylim = c(kpca_y_min, kpca_y_max)) +
    ggplot2::labs(title = "Kernel PCA (RBF)", x = "Eigenvalue Index [Log Scale]", y = "Eigenvalue") +
    ggplot2::theme_minimal(base_size = 12) + ggplot2::theme(plot.title = ggplot2::element_text(face="bold"))
  
  # ---------------------------------------------------------
  # PANEL 4: ENTROPIC SCREE
  # ---------------------------------------------------------
  # Calculate Entropic Kaiser (Eigenvalues > Mean Trace)
  ent_kaiser <- sum(results$eigenvalues > mean(results$eigenvalues))
  
  p_ent <- ggplot2::ggplot(ent_data, ggplot2::aes(x = Rank, y = Eigenvalue)) +
    ggplot2::geom_line(color = "dodgerblue", linewidth = 1) +
    ggplot2::geom_vline(xintercept = K_TRUE, color = "forestgreen", linetype = "solid", linewidth = 1.2) +
    ggplot2::annotate("text", x = K_TRUE, y = hent$h90, label = sprintf(" r \n (%d)", K_TRUE), hjust = -0.1, color = "forestgreen", fontface = "bold", size = ref_size) +
    ggplot2::geom_vline(xintercept = K_rlzd, color = "purple", linetype = "dotdash", linewidth = 1) +
    ggplot2::annotate("text", x = K_rlzd, y = hent$h75, label = sprintf(" K_rlzd \n (%d)", K_rlzd), hjust = 1.05, color = "purple", fontface = "italic", size = ref_size) +
    ggplot2::geom_vline(xintercept = ent_kaiser, color = "darkorange", linetype = "dashed", linewidth = 1) +
    ggplot2::annotate("text", x = ent_kaiser, y = hent$h52, label = sprintf("Kaiser, E>Mean \n(%d) ", ent_kaiser), hjust = 1.05, color = "darkorange", fontface = "italic", size = ref_size) +
    ggplot2::geom_vline(xintercept = (M_PROXIES - 1), color = "black", linetype = "dotted", linewidth = 1) +
    ggplot2::annotate("text", x = (M_PROXIES - 1), y = hent$h30, label = "Rank \nCeiling (m-1)", hjust = 1.05, size = ref_size) +
    ggplot2::scale_x_log10(breaks = log10_breaks) +
    ggplot2::coord_cartesian(ylim = c(ent_y_min, ent_y_max)) +
    ggplot2::labs(title = "Entropic Scree", x = "Eigenvalue Index [Log Scale]", y = "Eigenvalue") +
    ggplot2::theme_minimal(base_size = 12) + ggplot2::theme(plot.title = ggplot2::element_text(face="bold"))
  
  # ---------------------------------------------------------
  # RENDER 4-PANEL GRID
  # ---------------------------------------------------------
  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined <- (p_pca + p_spear) / (p_kpca + p_ent)
    print(combined)
  } else {
    cat("\n[!] To view the side-by-side 4-panel plot, please install the 'patchwork' package.\n")
    print(p_pca)
    print(p_spear)
    print(p_kpca)
    print(p_ent)
  }
}
