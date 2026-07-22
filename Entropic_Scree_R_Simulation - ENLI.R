# ==============================================================================
# ENTROPIC SCREE: Simulation & Demonstration Script
#
# Author: Terrence J. Lee-St. John, PhD
# Organization: Enli (www.enli.com.au)
# 
# Description: Generates a high-dimensional, mixed-type, noisy synthetic 
# dataset to demonstrate the structural collapse of standard PCA, and utilizes 
# the Entropic Scree to flawlessly extract the true Latent Generative Rank.
# ==============================================================================

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
  # [1/9] PURGE CONSTANTS & DUPLICATES
  # ----------------------------------------------------------------------------
  if (purge_constants) {
    cat("[1/9] Purging constants and identical duplicates...\n")
    const_cols <- names(dt)[sapply(dt, function(x) data.table::uniqueN(x, na.rm = TRUE) <= 1)]
    if (length(const_cols) > 0) dt[, (const_cols) := NULL]
    
    dup_cols <- duplicated(as.list(dt))
    if (any(dup_cols)) dt <- dt[, !dup_cols, with = FALSE]
  } else {
    cat("[1/9] Skipping constant and duplicate purge (user requested)...\n")
  }
  
  # ----------------------------------------------------------------------------
  # [2/9] MULTIVARIATE COLLINEARITY CHECK
  # ----------------------------------------------------------------------------
  if (check_collinearity) {
    cat("[2/9] Checking for perfect multivariate linear combinations...\n")
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
    cat("[2/9] Skipping collinearity check (user requested)...\n")
  }
  
  cat("[3/9] Discretizing continuous data and dense-ranking categoricals...\n")
  if (is.null(num_bins)) {
    num_bins <- max(2, floor(bin_multiplier * nrow(dt)^(1/3))) 
  }
  
  cols <- names(dt)
  dt[, (cols) := lapply(.SD, function(x) {
    if (is.numeric(x) && data.table::uniqueN(x) > num_bins) {
      infotheo::discretize(x, disc = "equalfreq", nbins = num_bins)[[1]]
    } else {
      data.table::frank(x, ties.method = "dense")
    }
  })]
  
  cat("[4/9] Purging non-linear monotonic duplicates and locking types...\n")
  dup_cols_post <- duplicated(as.list(dt))
  if (any(dup_cols_post)) dt <- dt[, !dup_cols_post, with = FALSE]
  
  valid_cols <- names(dt)
  dt[, (valid_cols) := lapply(.SD, function(x) as.integer(as.factor(x)))]
  
  cat("[5/9] Calculating marginal entropies and purging near-constants...\n")
  H_vec <- sapply(dt, infotheo::entropy)
  valid_vars <- names(H_vec)[H_vec >= low_entropy_thresh]
  if (length(H_vec) > length(valid_vars)) {
    dt <- dt[, ..valid_vars]
    H_vec <- H_vec[valid_vars]
  }
  
  p <- ncol(dt)
  if (p < 2) stop("Execution Halted: Less than 2 valid variables remain.")
  
  cat(sprintf("[6/9] Computing %d x %d Mutual Information Matrix (C++ OpenMP)...\n", p, p))
  
  mat_data <- as.matrix(dt)
  bin_sample_sizes <- apply(mat_data, 2, tabulate)
  rm(dt) 
  gc(verbose = FALSE)
  
  safe_target_cores <- max(1, num_cores)
  MI_mat <- fast_parallel_MI(mat_data, num_bins = num_bins, cores = safe_target_cores)
  rownames(MI_mat) <- colnames(MI_mat) <- valid_vars
  
  cat("[7/9] Applying Joint Entropy (Jaccard) Normalization...\n")
  sum_H_mat <- outer(H_vec, H_vec, FUN = "+")
  joint_H_mat <- sum_H_mat - MI_mat
  joint_H_mat[joint_H_mat < 1e-9] <- 1e-9
  
  NMI_mat <- MI_mat / joint_H_mat
  diag(NMI_mat) <- 1.0
  
  cat("[8/9] Extracting Entropic Latent Factors (Eigen Decomposition)...\n")
  eigen_res <- eigen(NMI_mat, symmetric = TRUE)
  eig_vals <- pmax(eigen_res$values, 1e-9)
  
  cat("[9/9] Calculating R_eff and Backward Scan Spectral Elbow...\n")
  sig_vals <- eig_vals[eig_vals > 0]
  if (length(sig_vals) > 0) {
    p_vals <- sig_vals / sum(sig_vals)
    H_spec <- -sum(p_vals * log(p_vals)) 
    R_eff <- exp(H_spec)
  } else {
    R_eff <- 1
  }
  
  # ==========================================================================
  # DIAGNOSTIC ONLY: MACRO GAP (NOISE CLIFF BOUNDARY)
  # ==========================================================================
  n_total <- length(eig_vals)
  valid_k <- sum(eig_vals > 1.0) 
  
  macro_max_noise_gap <- NA_real_
  macro_actual_gap <- NA_real_
  macro_gap_ratio <- NA_real_
  top_of_bulk_idx <- NA_integer_
  
  valid_search_space <- eig_vals[eig_vals > 1e-8]
  
  if (length(valid_search_space) > 10) {
    all_gaps_diag <- abs(diff(valid_search_space))
    n_active <- length(valid_search_space)
    noise_start_idx <- min(valid_k + max(3, floor(n_active * 0.05)), n_active - 5)
    noise_tail_idx <- noise_start_idx:(n_active - 1)
    
    if(length(noise_tail_idx) > 0) {
      noise_gaps <- all_gaps_diag[noise_tail_idx]
      max_noise_gap <- max(noise_gaps)
      
      macro_multiplier <- 10
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
  # LOGIC FOR MAXIMUM SECONDARY SPECTRAL GAP & TRIPLE TAP
  # ==========================================================================
  
  # --- BASE FALLBACK: MAXIMUM SECONDARY SPECTRAL GAP ---
  if (n_total >= 3) {
    all_gaps <- abs(diff(eig_vals))
    # Exclude the first gap (between lambda_1 and lambda_2) due to Perron-Frobenius
    secondary_gaps <- all_gaps[2:length(all_gaps)] 
    # Index is shifted by 1 relative to the true gap
    fallback_k <- which.max(secondary_gaps) + 1
    elbow_method <- "Maximum Secondary Spectral Gap"
  } else {
    fallback_k <- max(1, valid_k)
    elbow_method <- "Kaiser Criterion (> 1.0)"
  }
  
  K_elbow <- fallback_k
  tripped_sigma <- NA_real_
  tripped_breakout <- NA_real_
  
  # --- PRIMARY: TRIPLE-TAP SCANNER (OVERRIDE) ---
  if (valid_k >= 4 && n_total > 5) {
    window_size <- max(5, floor(n_total * 0.05)) 
    log_vals <- log(eig_vals)
    min_sigma <- 1e-4
    sigma_multiplier <- 10
    
    start_k <- n_total - 2
    
    for (k in seq(start_k, 3, by = -1)) {
      window_end <- min(n_total, k + window_size)
      tail_idxs <- k:window_end
      if (length(tail_idxs) < 3) next 
      
      tail_log_vals <- log_vals[tail_idxs]
      fit <- lm(tail_log_vals ~ tail_idxs)
      sigma <- summary(fit)$sigma
      if(is.nan(sigma) || is.na(sigma) || sigma < min_sigma) sigma <- min_sigma
      
      cand_idx <- k - 1
      pred_cand <- predict(fit, newdata = data.frame(tail_idxs = cand_idx))
      actual_cand <- log_vals[cand_idx]
      
      if (actual_cand > (pred_cand + (sigma_multiplier * sigma))) {
        confirmed <- TRUE
        
        # Verify it wasn't a fluke by checking the next two points up the curve
        if (cand_idx > 1) {
          pred_1 <- predict(fit, newdata = data.frame(tail_idxs = cand_idx - 1))
          if (log_vals[cand_idx - 1] <= (pred_1 + (sigma_multiplier * sigma))) confirmed <- FALSE
        }
        if (confirmed && cand_idx > 2) {
          pred_2 <- predict(fit, newdata = data.frame(tail_idxs = cand_idx - 2))
          if (log_vals[cand_idx - 2] <= (pred_2 + (sigma_multiplier * sigma))) confirmed <- FALSE
        }
        
        if (confirmed) {
          K_elbow <- cand_idx
          tripped_sigma <- sigma    
          tripped_breakout <- (actual_cand - pred_cand) / sigma
          elbow_method <- "Triple-Tap Scanner"
          break
        }
      }
    }
  }
  
  # ============================================================================
  # WAVE 1: INITIAL OUTPUT & METRICS
  # ============================================================================
  m_valid <- length(valid_vars)
  pct_prob_volume <- (R_eff / m_valid) * 100
  pct_redundant_signal <- (1 - (R_eff / m_valid)) * 100
  redundant_signal_volume <- m_valid - R_eff
  
  n_eigen_gt_1 <- valid_k
  n_eigen_le_1 <- n_total - valid_k
  
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
      ggplot2::labs(title = "Zoomed View", x = "Eigenvalue Index (m)", y = "Log(Eigenvalue)") +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 12))
    
    # MACRO VIEW
    macro_end <- min(length(eig_vals), max(50, K_elbow * 3))
    plot_df_macro <- data.frame(Rank = 1:macro_end, Eigenvalue = eig_vals[1:macro_end])
    macro_y_max <- if(length(eig_vals) >= 2) eig_vals[2] * 1.1 else max(eig_vals)
    
    p_scree_macro <- ggplot2::ggplot(plot_df_macro, ggplot2::aes(x = Rank, y = Eigenvalue)) +
      ggplot2::geom_line(color = "dodgerblue", linewidth = 1) +
      ggplot2::geom_point(color = "dodgerblue", size = 2) +
      ggplot2::geom_vline(xintercept = K_elbow, color = "#D55E00", linetype = "dashed", linewidth = 1.2) +
      ggplot2::scale_y_continuous(trans = 'log10') +
      ggplot2::coord_cartesian(ylim = c(NA, macro_y_max)) + 
      ggplot2::labs(title = "Macro View", x = "Eigenvalue Index (m)", y = "Log(Eigenvalue)") +
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
  cat("-----------------------------------------------------------------\n")
  cat(sprintf(" -> %-50s : %.2f\n", "Total Unique Probabilistic Volume (R_eff)", R_eff))
  cat(sprintf(" -> %-50s : %.1f%%\n", "%", pct_prob_volume))
  cat("      (Unique Signal + Structural Uncertainty + \u22A5 Measurement Error)\n")
  cat("-----------------------------------------------------------------\n")
  cat(sprintf(" -> %-50s : %.2f\n", "Redundant Signal Volume (m - R_eff)", redundant_signal_volume))
  cat(sprintf(" -> %-50s : %.1f%%\n", "%", pct_redundant_signal))
  cat("-----------------------------------------------------------------\n")
  cat(sprintf(" -> %-50s : %d\n", "Eigenvalues > 1.0", n_eigen_gt_1))
  cat(sprintf(" -> %-50s : %d\n", "Eigenvalues <= 1.0", n_eigen_le_1))
  cat("=================================================================\n\n")
  
  cat("=================================================================\n")
  cat(" AUTOMATED ELBOW DETECTION (HEURISTIC)\n")
  cat("=================================================================\n")
  cat(sprintf(" -> %-43s : %d\n", "Automated Extracted Elbow Rank (K_elbow)", K_elbow))
  cat(sprintf(" -> %-43s : %s\n", "Extraction Method Tripped", elbow_method))
  if (!is.na(tripped_breakout)) {
    cat("-----------------------------------------------------------------\n")
    cat(" [Triple-Tap Scanner Details]\n")
    cat(sprintf(" -> %-43s : %.6f\n", "Scanner Regression Baseline Sigma", tripped_sigma))
    cat(sprintf(" -> %-43s : %.2f-Sigma\n", "Actual Breakout Magnitude", tripped_breakout))
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
          # RESTORED: INTERACTIVE GRAPH PREVIEW & METRICS
          # ====================================================================
          if (requireNamespace("ggplot2", quietly = TRUE)) {
            # 1. ZOOMED VIEW (Updated)
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
            
            # 2. MACRO VIEW (Updated)
            macro_end_upd <- min(length(eig_vals), max(50, K_final * 3, K_elbow * 3))
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
          signal_weight_prev <- signal_variance_prev / m_valid
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
          cat(sprintf(" -> Avg Info Gravity (AIG)      : %.2f\n", AIG_prev))
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
  signal_weight <- signal_variance / m_valid
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
  cat(sprintf(" -> Avg Info Gravity (AIG)            : %.2f\n", AIG))
  cat(" -> Factor-Specific Informational Gravity (FSIG):\n")
  print(round(FSIG_final, 3))
  cat("=================================================================\n\n")
  
  # ============================================================================
  # WAVE 3: FINAL TRIPARTITE STRUCTURAL COMPOSITION
  # ============================================================================
  noise_variance <- sum(eig_vals[(K_final + 1):m_valid])
  noise_weight <- noise_variance / m_valid
  idiosyncratic_noise_volume <- R_eff * noise_weight
  
  cat("=================================================================\n")
  cat(sprintf(" (FINAL) TRIPARTITE STRUCTURAL COMPOSITION (based user-confirmed K_elbow = %d)\n", K_final))
  cat("=================================================================\n")
  cat(sprintf(" -> %-50s : %d\n", "Valid Variables (m)", m_valid))
  cat("-----------------------------------------------------------------\n")
  cat(sprintf(" -> %-50s : %.2f\n", "Total Signal Volume (Unique + Redundant)", total_signal_volume))
  cat(sprintf(" -> %-40s : %.3f\n", "   (Unique Signal Volume)", unique_signal_volume))
  cat(sprintf(" -> %-40s : %.3f\n", "   (Redundant Signal Volume)", redundant_signal_volume))
  cat("-----------------------------\n")
  cat(sprintf(" -> %-50s : %.2f\n", "Idiosyncratic Noise Volume", idiosyncratic_noise_volume))
  cat("      (Structural Uncertainty + \u22A5 Measurement Error)\n")
  cat("-----------------------------------------------------------------\n")
  cat(sprintf(" -> %-50s : %.2f%%\n", "% Total Signal", (total_signal_volume / m_valid) * 100))
  cat(sprintf(" -> %-40s : %.3f%%\n", "   (% Unique Signal)", (unique_signal_volume / m_valid) * 100))
  cat(sprintf(" -> %-40s : %.3f%%\n", "   (% Redundant Signal)", (redundant_signal_volume / m_valid) * 100))
  cat("-----------------------------\n")
  cat(sprintf(" -> %-50s : %.2f%%\n", "% Idiosyncratic Noise", (idiosyncratic_noise_volume / m_valid) * 100))
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
  total_sig_vol_bulk <- (R_eff * (sig_var_bulk / m_valid)) + redundant_signal_volume
  FSIG_extended_bulk <- p_extended_bulk * total_sig_vol_bulk
  
  # --- EXTENDED MODEL B: KAISER RULE BOUNDARY ---
  extended_kaiser_k <- max(K_final, valid_k)
  extended_eigenvals_kaiser <- eig_vals[1:extended_kaiser_k]
  p_extended_kaiser <- extended_eigenvals_kaiser / sum(extended_eigenvals_kaiser)
  sig_var_kaiser <- sum(extended_eigenvals_kaiser)
  total_sig_vol_kaiser <- (R_eff * (sig_var_kaiser / m_valid)) + redundant_signal_volume
  FSIG_extended_kaiser <- p_extended_kaiser * total_sig_vol_kaiser
  
  return(list(
    eigenvalues = eig_vals,
    similarity_matrix = NMI_mat,
    retained_features = valid_vars,
    bin_distributions = bin_sample_sizes,
    R_eff = R_eff,
    K_auto_extracted = K_elbow,
    extraction_method = elbow_method,
    tripped_sigma = tripped_sigma,
    tripped_breakout = tripped_breakout,
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
  
  # NEW: Allow dynamic control over the continuous vs binary split
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
  # To restore the factorial penalty later, swap the comments on the next two lines:
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
  mask <- matrix(rbinom(n_terms * m_proxies, 1, 0.25), nrow = n_terms, ncol = m_proxies)
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
  
  # Calculate Algebraic K_rlzd before returning
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

cat(sprintf("Generating Mixed G2G Universe: %s rows, %d Proxies, %d Latent Drivers...\n", format(N_ROWS, big.mark=","), M_PROXIES, K_TRUE))

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
# term_sds <- sqrt(1 / factorial(interaction_orders))
term_sds <- rep(1.0, length(term_names))
int_scaling <- 1.0 
# -----------------

term_weights <- ifelse(interaction_orders == 1, term_sds, term_sds * int_scaling)
muzzled_design_mat <- t(t(design_mat) * term_weights)

expanded_cov <- cov(muzzled_design_mat)
expanded_eigen <- pmax(eigen(expanded_cov, symmetric = TRUE)$values, 0)
p_expanded <- expanded_eigen[expanded_eigen > 1e-9] / sum(expanded_eigen[expanded_eigen > 1e-9])
true_continuous_ceiling <- exp(-sum(p_expanded * log(p_expanded)))

cat(sprintf("\n[***] R_alg (Effective Latent Configuration Rank): %.2f\n\n", true_continuous_ceiling))

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
                                    )

# ==============================================================================
# 5. STANDARD PCA EXTRACTION (FOR COMPARISON)
# ==============================================================================
cat("\nExtracting Standard PCA for comparison...\n")
start_pca <- Sys.time()

# We use prcomp with scaling to mirror a Pearson Correlation Matrix extraction
pca_res <- prcomp(observed_data, center = TRUE, scale. = TRUE)

# Calculate standard PCA eigenvalues
pca_eigenvalues <- pca_res$sdev^2

pca_time <- round(as.numeric(difftime(Sys.time(), start_pca, units = "secs")), 2)
cat(sprintf("Standard PCA completed in %.2f seconds.\n", pca_time))

# Create a comparison data frame
# PCA will only return N eigenvalues, padding the rest with exactly 0
m_total <- ncol(observed_data)
pca_padded <- c(pca_eigenvalues, rep(0, m_total - length(pca_eigenvalues)))

df_compare <- data.frame(
  Rank = rep(1:m_total, 2),
  Eigenvalue = c(pca_padded, results$eigenvalues),
  Method = factor(rep(c("Standard PCA", "Entropic Scree"), each = m_total),
                  levels = c("Standard PCA", "Entropic Scree"))
)

# Standard Kaiser Rule for comparison
pca_kaiser <- sum(pca_eigenvalues > 1.0)

# ==============================================================================
# 6. SIDE-BY-SIDE VISUAL PROOF
# ==============================================================================
if (requireNamespace("ggplot2", quietly = TRUE)) {
  
  # Calculate caps for visual clarity
  pca_y_max <- if(length(pca_eigenvalues) >= 2) pca_eigenvalues[2] * 1.1 else max(pca_eigenvalues)
  ent_y_max <- if(length(results$eigenvalues) >= 2) results$eigenvalues[2] * 1.1 else max(results$eigenvalues)
  
  # Strict lower bound for Entropic Scree (ignoring the 1e-9 artificial clamp)
  ent_y_min <- min(results$eigenvalues[results$eigenvalues > 1e-8])
  pca_y_min <- 0 # Standard PCA floor
  
  # Dynamic log breaks generator function (1, 10, 100, 1000...)
  log10_breaks <- function(x) {
    10^seq(floor(log10(min(x))), ceiling(log10(max(x))))
  }
  
  # Extract K_rlzd and R_alg
  K_rlzd <- true_universe$active_terms
  R_alg <- true_continuous_ceiling 
  
  # Global annotation font size
  ref_size <- 3.5
  
  # --- Linear-Space Positions for Annotations ---
  # Staggering every label across unique vertical heights guarantees 
  # zero visual overlap, regardless of how close the lines are on the x-axis.
  pca_y_90 <- pca_y_min + (pca_y_max - pca_y_min) * 0.90
  pca_y_75 <- pca_y_min + (pca_y_max - pca_y_min) * 0.75
  pca_y_60 <- pca_y_min + (pca_y_max - pca_y_min) * 0.60
  pca_y_45 <- pca_y_min + (pca_y_max - pca_y_min) * 0.45
  pca_y_30 <- pca_y_min + (pca_y_max - pca_y_min) * 0.30
  
  ent_y_90 <- ent_y_min + (ent_y_max - ent_y_min) * 0.90
  ent_y_75 <- ent_y_min + (ent_y_max - ent_y_min) * 0.75
  ent_y_60 <- ent_y_min + (ent_y_max - ent_y_min) * 0.60
  
  # Split the PCA data to visually distinguish the Null Space
  pca_data <- df_compare[df_compare$Method == "Standard PCA", ]
  
  # Plot PCA (Linear Scale Y, Log Scale X)
  p_pca <- ggplot2::ggplot(pca_data, ggplot2::aes(x = Rank, y = Eigenvalue)) +
    
    # 1. Calculated Variance (Solid Red)
    ggplot2::geom_line(data = pca_data[pca_data$Rank < N_ROWS, ], color = "firebrick", linewidth = 1) +
    
    # 2. Algebraic Null Space (Dashed Gray - exactly 0 due to N-1 ceiling)
    ggplot2::geom_line(data = pca_data[pca_data$Rank >= (N_ROWS - 1), ], color = "gray60", linetype = "dashed", linewidth = 1) +
    
    # Add r (Latent Generative Rank)
    ggplot2::geom_vline(xintercept = K_TRUE, color = "forestgreen", linetype = "solid", linewidth = 1.2) +
    ggplot2::annotate("text", x = K_TRUE, y = pca_y_90, label = sprintf("r (%d)", K_TRUE), hjust = -0.1, color = "forestgreen", fontface = "bold", size = ref_size) +
    
    # Add R_alg (Effective Latent Configurational Rank)
    ggplot2::geom_vline(xintercept = R_alg, color = "magenta", linetype = "longdash", linewidth = 1) +
    ggplot2::annotate("text", x = R_alg, y = pca_y_75, label = sprintf("R_alg (%.1f)", R_alg), hjust = 1.05, color = "magenta", fontface = "italic", size = ref_size) +
    
    # Add K_rlzd (Realized Latent Configurational Rank)
    ggplot2::geom_vline(xintercept = K_rlzd, color = "purple", linetype = "dotdash", linewidth = 1) +
    ggplot2::annotate("text", x = K_rlzd, y = pca_y_60, label = sprintf("K_rlzd (%d)", K_rlzd), hjust = -0.05, color = "purple", fontface = "italic", size = ref_size) +
    
    # Add PCA Kaiser Rule (E > 1.0)
    ggplot2::geom_vline(xintercept = pca_kaiser, color = "darkorange", linetype = "dashed", linewidth = 1) +
    ggplot2::annotate("text", x = pca_kaiser, y = pca_y_45, label = sprintf("Kaiser Cutoff (E>1): %d", pca_kaiser), hjust = 1.05, color = "darkorange", fontface = "italic", size = ref_size) +
    
    # Add PCA N-1 Ceiling (Dynamic)
    ggplot2::geom_vline(xintercept = (N_ROWS - 1), color = "black", linetype = "dotted", linewidth = 1) +
    ggplot2::annotate("text", x = (N_ROWS - 1), y = pca_y_30, label = "PCA Rank Ceiling (N-1)", hjust = 1.05, size = ref_size) +
    
    # Linear Y-axis, Log X-axis
    ggplot2::scale_x_log10(breaks = log10_breaks) +
    ggplot2::coord_cartesian(ylim = c(pca_y_min, pca_y_max)) +
    ggplot2::labs(
      title = "Standard PCA",
      x = "Eigenvalue Index (m) [Log Scale]", y = "Eigenvalue (from Correlation Matrix)"
    ) +
    ggplot2::theme_minimal(base_size = 14)
  
  # Plot Entropic Scree (Linear Scale Y, Log Scale X)
  p_ent <- ggplot2::ggplot(df_compare[df_compare$Method == "Entropic Scree", ], 
                           ggplot2::aes(x = Rank, y = Eigenvalue)) +
    ggplot2::geom_line(color = "dodgerblue", linewidth = 1) +
    
    # Add r (Latent Generative Rank)
    ggplot2::geom_vline(xintercept = K_TRUE, color = "forestgreen", linetype = "solid", linewidth = 1.2) +
    ggplot2::annotate("text", x = K_TRUE, y = ent_y_90, label = sprintf("r (%d)", K_TRUE), hjust = -0.1, color = "forestgreen", fontface = "bold", size = ref_size) +
    
    # Add R_alg (Effective Latent Configurational Rank)
    ggplot2::geom_vline(xintercept = R_alg, color = "magenta", linetype = "longdash", linewidth = 1) +
    ggplot2::annotate("text", x = R_alg, y = ent_y_75, label = sprintf("R_alg (%.1f)", R_alg), hjust = 1.05, color = "magenta", fontface = "italic", size = ref_size) +
    
    # Add K_rlzd (Realized Latent Configurational Rank)
    ggplot2::geom_vline(xintercept = K_rlzd, color = "purple", linetype = "dotdash", linewidth = 1) +
    ggplot2::annotate("text", x = K_rlzd, y = ent_y_60, label = sprintf("K_rlzd (%d)", K_rlzd), hjust = -0.05, color = "purple", fontface = "italic", size = ref_size) +
    
    # Linear Y-axis, Log X-axis
    ggplot2::scale_x_log10(breaks = log10_breaks) +
    
    # Bounded to crop out the 1e-9 clamp tail for a clean linear floor
    ggplot2::coord_cartesian(ylim = c(ent_y_min, ent_y_max)) +
    ggplot2::labs(
      title = "Entropic Scree",
      x = "Eigenvalue Index (m) [Log Scale]", y = "Eigenvalue (from NMI Matrix)"
    ) +
    ggplot2::theme_minimal(base_size = 14)
  
  if (requireNamespace("patchwork", quietly = TRUE)) {
    print(p_pca + p_ent)
  } else {
    print(p_pca)
    print(p_ent)
  }
}