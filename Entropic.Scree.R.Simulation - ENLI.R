options(max.print = 9999999)
rm(list = ls())
gc(verbose = FALSE)

# ==============================================================================
# ENTROPIC SCREE Function (v1.0.0)
#
# Author: Terrence J. Lee-St. John
# Organization: Enli (www.enli.com.au)
#
# Description: An information-theoretic diagnostic technique for estimating the
# intrinsic dimensionality of tabular datasets. Evaluates shared probability mass
# via a transformed mutual information metric. Aims to extract the Intrinsic
# Generative Rank (r) and structural topology.
# ==============================================================================
#
# Copyright 2026 Terrence J. Lee-St. John (Enli)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ==============================================================================

# List of all required packages (added pkgbuild for robust Rtools checking)
required_packages <- c("Rcpp", "data.table", "infotheo", "ggplot2", "patchwork", "pkgbuild")

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
      "\n===================================================================================\n",
      " [!] MISSING OR INCOMPATIBLE C++ COMPILER (Rtools)\n",
      "===================================================================================\n",
      " Rtools is required to build the C++ backend on Windows.\n",
      " It is either missing, not on your PATH, or your Rtools version \n",
      " does not match your R version.\n\n",
      sprintf(" Your current R version is: %s\n", getRversion()),
      " You must install the version of Rtools that matches this R version.\n\n",
      " Please download and install the correct Rtools here:\n",
      " https://cran.r-project.org/bin/windows/Rtools/\n\n",
      " Note: After installing, you MUST restart your R session before \n",
      " running this script again.\n",
      "===================================================================================\n",
      call. = FALSE
    )
  } else {
    cat("Compatible Rtools C++ compiler found. Ready to build backend.\n")
  }
}

# ==============================================================================
# C++ COMPILATION (OpenMP Degradation)
# ==============================================================================
cat("[+] Attempting to compile parallelized C++ backend (OpenMP)...\n")

cpp_parallel <- '
#include <Rcpp.h>
#include <omp.h>
#include <cmath>
#include <vector>
#include <algorithm>
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
    
    #pragma omp parallel num_threads(cores)
    {
        std::vector<int> joint(num_bins * num_bins, 0);
        
        #pragma omp for schedule(dynamic)
        for(int i = 0; i < p; ++i) {
            for(int j = i + 1; j < p; ++j) {
                std::fill(joint.begin(), joint.end(), 0);
                
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
    } 
    return MI;
}
'

cpp_single <- '
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>

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
    
    std::vector<int> joint(num_bins * num_bins, 0);
    
    for(int i = 0; i < p; ++i) {
        for(int j = i + 1; j < p; ++j) {
            std::fill(joint.begin(), joint.end(), 0);
            
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
'

# Attempt to compile the parallel version first. If it fails (missing OpenMP), catch the error and run the single-threaded fallback.
compile_success <- tryCatch({
  # suppressWarnings hides the red compilation text from cluttering the user's console if it fails
  suppressWarnings(suppressMessages(Rcpp::sourceCpp(code = cpp_parallel)))
  TRUE
}, error = function(e) {
  FALSE
})

if (compile_success) {
  cat("    -> Success! OpenMP detected. Parallel processing enabled.\n")
} else {
  cat("    -> OpenMP not detected (common on Mac). Compiling single-threaded fallback...\n")
  cat("       (Note: For parallel processing on Mac, 'libomp' must be manually installed).\n")
  Rcpp::sourceCpp(code = cpp_single)
}

# ==============================================================================
# 1. ENTROPIC SCREE FUNCTION
# ==============================================================================
Entropic.Scree <- function(data
                           , low_entropy_thresh = 0.05
                           , num_bins = NULL
                           , bin_multiplier = 1.0
                           , num_cores = parallel::detectCores() - 2
                           , interactive_mode = TRUE
                           , purge_constants = TRUE
                           , check_collinearity = TRUE
                           , triple_tap_window = 20
                           , extract_eigenvectors = FALSE
                           , extract_bipolar_modules = FALSE
                           , bipolar_top_n = 10
                           , return_processed_data = FALSE) {
  
  # Enforce data.table requirement
  if (!data.table::is.data.table(data)) {
    stop("Error: 'data' must be a data.table. Please convert your dataset using data.table::as.data.table() before running Entropic.Scree.")
  }
  
  # Enforce eigenvector extraction if bipolar modules are requested
  if (extract_bipolar_modules) {
    extract_eigenvectors <- TRUE
  }
  
  start_time <- Sys.time()
  dt <- data.table::copy(data)
  
  # Enforce a minimum window size of 3 for linear regression (df >= 1)
  triple_tap_window <- max(3, as.integer(triple_tap_window))
  
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
    const_cols <- names(dt)[sapply(dt, function(x)
      data.table::uniqueN(x, na.rm = TRUE) <= 1)]
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
        num_mat <- as.matrix(na.omit(dt[sample_idx, num_cols, with = FALSE]))
        
        if (nrow(num_mat) > 0) {
          qr_mat <- cbind(Intercept = 1, num_mat)
          qr_decomp <- qr(qr_mat, tol = 1e-7)
          
          if (qr_decomp$rank < ncol(qr_mat)) {
            drop_indices <- qr_decomp$pivot[(qr_decomp$rank + 1):ncol(qr_mat)]
            lin_combos <- setdiff(colnames(qr_mat)[drop_indices], "Intercept")
            if (length(lin_combos) > 0) {
              cat(sprintf("      -> Purged %d perfectly collinear variables.\n", length(lin_combos)))
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
    dt <- dt[, valid_vars, with = FALSE]
    H_vec <- H_vec[valid_vars]
  }
  
  p <- ncol(dt)
  if (p < 2) stop("Execution Halted: Less than 2 valid variables remain.")
  
  cat(sprintf("[6/10] Computing %d x %d Mutual Information Matrix (C++ OpenMP)...\n", p, p))
  
  if (return_processed_data) {
    processed_data_out <- data.table::copy(dt)
  } else {
    processed_data_out <- NULL
  }
  
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
  if (extract_eigenvectors) {
    eigen_res <- eigen(NMI_mat_c, symmetric = TRUE)
    raw_eig_vals <- eigen_res$values
    vectors_out <- eigen_res$vectors
  } else {
    eigen_res <- eigen(NMI_mat_c, symmetric = TRUE, only.values = TRUE)
    raw_eig_vals <- eigen_res$values
    vectors_out <- NULL
  }
  
  # Calculate SCDR before clipping
  m_plus <- sum(raw_eig_vals[raw_eig_vals > 0])
  m_minus <- sum(abs(raw_eig_vals[raw_eig_vals < 0]))
  SCDR <- (m_minus / m_plus) * 100
  
  eig_vals <- pmax(raw_eig_vals, 1e-9)
  
  # Count true positive dimensions before the geometric zero to bound plots
  max_plot_idx <- sum(raw_eig_vals > 1e-12)
  
  # Constructive Spectral Mass (sum of positive clipped eigenvalues)
  m_plus <- sum(eig_vals)
  
  cat("[10/10] Calculating Total Unique Probabilistic Volume (R_eff) and Estimating Phase Transitions...\n")
  sig_vals <- eig_vals[eig_vals > 0]
  if (length(sig_vals) > 0) {
    p_vals <- sig_vals / sum(sig_vals)
    H_spec <- -sum(p_vals * log(p_vals))
    R_eff <- exp(H_spec)
  } else {
    R_eff <- 1
  }
  
  # ==========================================================================
  # STEP 1: MACRO GAP BOUNDARY (IDIOSYNCRATIC INFORMATIONAL VARIANCE CLIFF DETECTION)
  # ==========================================================================
  n_total <- length(eig_vals)
  valid_k <- sum(eig_vals > mean_trace)
  
  macro_max_bulk_gap <- NA_real_
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
      
      macro_multiplier <- 1.5
      gap_threshold <- max(1e-6, max_noise_gap * macro_multiplier)
      
      macroscopic_gap_indices <- which(all_gaps_diag > gap_threshold)
      if (length(macroscopic_gap_indices) > 0) {
        top_of_bulk_idx <- max(macroscopic_gap_indices) + 1
        macro_max_bulk_gap <- max_noise_gap
        macro_actual_gap <- all_gaps_diag[top_of_bulk_idx - 1]
        if (max_noise_gap > 1e-9) macro_gap_ratio <- macro_actual_gap / max_noise_gap
      }
    }
  }
  
  search_start_idx <- if (!is.na(top_of_bulk_idx)) max(2, top_of_bulk_idx - 1) else valid_k
  
  # ==========================================================================
  # STEP 2: ENGINE A - MAXIMUM SECONDARY EIGENVALUE RATIO (LOG-GAP)
  # ==========================================================================
  # This engine isolates the Observed Generative Rank (K_roots)
  K_log_gap <- NA_integer_
  log_gap_pct_drop <- NA_real_
  log_gap_ratio <- NA_real_
  log_gap_fallback <- FALSE
  n_valid_search <- length(valid_search_space)
  
  if (n_valid_search >= 3) {
    all_log_gaps <- abs(diff(log(valid_search_space)))
    
    if (!is.na(top_of_bulk_idx) && top_of_bulk_idx < 4) {
      # Ultra-low rank edge case
      ordered_gaps <- order(all_log_gaps, decreasing = TRUE)
      second_largest_gap_idx <- ordered_gaps[2]
      if (second_largest_gap_idx >= top_of_bulk_idx) {
        K_log_gap <- max(1, top_of_bulk_idx - 1)
        log_gap_fallback <- TRUE
      } else {
        K_log_gap <- second_largest_gap_idx
      }
    } else {
      # Standard Evaluation
      search_limit <- if (!is.na(top_of_bulk_idx)) max(2, top_of_bulk_idx - 1) else valid_k
      if (search_limit >= 2) {
        secondary_gaps <- all_log_gaps[2:search_limit]
        K_log_gap <- which.max(secondary_gaps) + 1
      } else {
        K_log_gap <- 1
        log_gap_fallback <- TRUE
      }
    }
  } else {
    K_log_gap <- max(1, valid_k)
  }
  
  if (K_log_gap < length(eig_vals)) {
    val_current <- eig_vals[K_log_gap]
    val_next <- eig_vals[K_log_gap + 1]
    if (val_next > 1e-12) {
      log_gap_ratio <- val_current / val_next
      log_gap_pct_drop <- (1 - (val_next / val_current)) * 100
    }
  }
  
  # ==========================================================================
  # STEP 3: ENGINE B - TRIPLE-TAP (DYNAMIC LINEAR WITH MACRO-STITCH)
  # ==========================================================================
  K_triple_tap <- NA_integer_
  triple_tap_multiplier <- NA_real_
  triple_tap_actual_sigma <- NA_real_
  triple_tap_expected_val <- NA_real_
  triple_tap_actual_val <- NA_real_
  triple_tap_local_sd <- NA_real_
  prob_target <- NA_real_
  stitch_applied <- FALSE
  triple_tap_fallback <- FALSE
  
  if (search_start_idx >= 3) {
    # 1. Apply the Topological Stitch directly to the eigenvalues
    stitched_eig_vals <- eig_vals
    if (!is.na(top_of_bulk_idx) && !is.na(macro_actual_gap)) {
      stitch_constant <- 0.5 * macro_actual_gap
      # Add constant to all eigenvalues at or inside the noise bulk
      stitched_eig_vals[top_of_bulk_idx:length(stitched_eig_vals)] <- stitched_eig_vals[top_of_bulk_idx:length(stitched_eig_vals)] + stitch_constant
      stitch_applied <- TRUE
    }
    
    # 2. Dynamically scale sigma using t-dist to target a family-wise false positive rate of 1/100
    # Bonferroni correction: alpha = 0.01 / number of tests (search_start_idx)
    prob_target <- 1 - (1 / (100 * search_start_idx))
    
    # Pre-compute design matrix for the dynamic reference window - LINEAR ONLY
    x_ref <- 0:(triple_tap_window - 1)
    X_mat <- cbind(1, x_ref)
    
    # Scan backwards from the top of the bounded space down to 1
    for (i in seq(search_start_idx, 1, by = -1)) {
      target_val <- stitched_eig_vals[i]
      
      # Reference window: next N eigenvalues to calculate linear expectation
      ref_start <- i + 1
      ref_end <- min(length(stitched_eig_vals), i + triple_tap_window)
      
      if (ref_start > length(stitched_eig_vals)) next # Safety check at the absolute tail
      
      ref_vals <- stitched_eig_vals[ref_start:ref_end]
      
      if (length(ref_vals) == triple_tap_window) {
        # Fit 1st degree polynomial (2 params). df = window - 2
        df_local <- triple_tap_window - 2
        dynamic_t <- qt(prob_target, df = df_local)
        current_multiplier <- max(2.0, dynamic_t)
        
        fit <- lm.fit(x = X_mat, y = ref_vals)
        
        # Intercept is the prediction at the first reference point (x = 0).
        # We predict the target eigenvalue one step backward (x = -1).
        expected_val <- fit$coefficients[1] - fit$coefficients[2]
        
        # Protect against linear overshoot pulling the expected value negative
        expected_val <- max(1e-12, expected_val)
        
        # Calculate local standard deviation of residuals (RSS / df)
        rss <- sum(fit$residuals^2)
        local_sd <- max(1e-12, sqrt(rss / df_local))
        
        if (target_val > (expected_val + current_multiplier * local_sd)) {
          K_triple_tap <- i
          triple_tap_multiplier <- current_multiplier
          triple_tap_expected_val <- expected_val
          triple_tap_actual_val <- target_val
          triple_tap_actual_sigma <- (target_val - expected_val) / local_sd
          triple_tap_local_sd <- local_sd
          break
        }
      } else if (length(ref_vals) >= 2) {
        # Fallback for the absolute tail where full window points are not available. df = n - 1
        df_local <- length(ref_vals) - 1
        dynamic_t <- qt(prob_target, df = df_local)
        current_multiplier <- max(2.0, dynamic_t)
        
        local_mean <- mean(ref_vals)
        local_sd <- max(1e-12, sd(ref_vals))
        
        if (target_val > (local_mean + current_multiplier * local_sd)) {
          K_triple_tap <- i
          triple_tap_multiplier <- current_multiplier
          triple_tap_expected_val <- local_mean
          triple_tap_actual_val <- target_val
          triple_tap_actual_sigma <- (target_val - local_mean) / local_sd
          triple_tap_local_sd <- local_sd
          break
        }
      }
    }
  }
  
  # Fallback if Triple-Tap never triggers (or if search_start_idx < 3)
  if (is.na(K_triple_tap)) {
    K_triple_tap <- K_log_gap
    triple_tap_fallback <- TRUE
  }
  
  # Safely capture the multiplier for the printout if it never broke but we had space
  if (is.na(triple_tap_multiplier) && search_start_idx >= 3) {
    triple_tap_multiplier <- max(2.0, qt(1 - (1 / (100 * search_start_idx)), df = triple_tap_window - 2))
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
    
    # Define dual-lines for GGplot dynamically based on convergence
    if (K_log_gap == K_triple_tap) {
      lines_geom <- ggplot2::geom_vline(xintercept = K_log_gap, color = "forestgreen", linetype = "dashed", linewidth = 1.2)
      sub_title_text <- sprintf("Converged: Observed Generative Rank & Extended Signal Tail = %d", K_log_gap)
    } else {
      lines_geom <- list(
        ggplot2::geom_vline(xintercept = K_log_gap, color = "#D55E00", linetype = "dashed", linewidth = 1.2),
        ggplot2::geom_vline(xintercept = K_triple_tap, color = "purple", linetype = "twodash", linewidth = 1.2)
      )
      sub_title_text <- sprintf("Observed Generative Rank (Orange) = %d | Extended Signal Tail (Purple) = %d", K_log_gap, K_triple_tap)
    }
    
    # ZOOMED VIEW
    max_k <- max(K_log_gap, K_triple_tap)
    zoom_start <- max(1, min(K_log_gap, K_triple_tap) - 5)
    zoom_end <- min(max_plot_idx, max_k + 15)
    plot_df_zoom <- data.frame(Rank = zoom_start:zoom_end, Eigenvalue = eig_vals[zoom_start:zoom_end])
    
    p_scree_zoom <- ggplot2::ggplot(plot_df_zoom, ggplot2::aes(x = Rank, y = Eigenvalue)) +
      ggplot2::geom_line(color = "dodgerblue", linewidth = 1) +
      ggplot2::geom_point(color = "dodgerblue", size = 2) +
      lines_geom +
      ggplot2::scale_y_continuous(trans = 'log10') +
      ggplot2::scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(min(x), max(x)))))) +
      ggplot2::labs(title = "Zoomed View", x = "Eigenvalue Index", y = "Log(Eigenvalue)") +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 12))
    
    # MACRO VIEW (Dynamically extended)
    macro_base <- max(50, max_k * 10)
    if (!is.na(top_of_bulk_idx)) {
      macro_base <- max(macro_base, top_of_bulk_idx + 25)
    }
    macro_end <- min(max_plot_idx, macro_base)
    
    plot_df_macro <- data.frame(Rank = 1:macro_end, Eigenvalue = eig_vals[1:macro_end])
    macro_y_max <- if(length(eig_vals) >= 2) eig_vals[2] * 1.1 else max(eig_vals)
    
    p_scree_macro <- ggplot2::ggplot(plot_df_macro, ggplot2::aes(x = Rank, y = Eigenvalue)) +
      ggplot2::geom_line(color = "dodgerblue", linewidth = 1) +
      ggplot2::geom_point(color = "dodgerblue", size = 2) +
      lines_geom +
      ggplot2::scale_y_continuous(trans = 'log10') +
      ggplot2::coord_cartesian(ylim = c(NA, macro_y_max)) +
      ggplot2::labs(title = "Macro View", x = "Eigenvalue Index", y = "Log(Eigenvalue)") +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 12))
    
    if (requireNamespace("patchwork", quietly = TRUE)) {
      combined_plot <- (p_scree_macro + p_scree_zoom) +
        patchwork::plot_annotation(
          title = "Entropic Scree Results",
          subtitle = sub_title_text,
          theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 16, hjust = 0.5),
                                 plot.subtitle = ggplot2::element_text(size = 14, hjust = 0.5))
        )
      print(combined_plot)
    } else {
      print(p_scree_macro)
      print(p_scree_zoom)
    }
  }
  
  cat("\n===================================================================================\n")
  cat(" STRUCTURAL COMPOSITION\n")
  cat("===================================================================================\n")
  cat(sprintf(" -> %-50s : %d\n", "Valid Variables (m)", m_valid))
  cat(sprintf(" -> %-50s : %.3f\n", "Centered Trace (Tr_Mc)", Tr_Mc))
  cat(sprintf(" -> %-50s : %.3f%%\n", "Synergistic Curvature Deficit Ratio (SCDR)", SCDR))
  cat("-----------------------------------------------------------------\n")
  cat(sprintf(" -> %-50s : %.3f\n", "Total Unique Probabilistic Volume (R_eff)", R_eff))
  cat(sprintf(" -> %-50s : %.3f%%\n", "%", pct_prob_volume))
  cat("      (Unique Signal Volume + Structural Uncertainty\n")
  cat("      + Independent Measurement Error + Unshared Signal Geometry)\n")
  cat("-----------------------------------------------------------------\n")
  cat(sprintf(" -> %-50s : %.3f\n", "Redundant Signal Volume (m - R_eff)", redundant_signal_volume))
  cat(sprintf(" -> %-50s : %.3f%%\n", "%", pct_redundant_signal))
  cat("-----------------------------------------------------------------\n")
  cat(sprintf(" -> %-50s : %d\n", "Eigenvalues > Mean Trace", n_eigen_gt_mean))
  cat(sprintf(" -> %-50s : %d\n", "Eigenvalues <= Mean Trace", n_eigen_le_mean))
  cat("===================================================================================\n\n")
  
  cat("===================================================================================\n")
  cat(" AUTOMATED ELBOW DETECTION (DUAL-DIAGNOSTIC ENSEMBLE)\n")
  cat("===================================================================================\n")
  cat(" [Diagnostic: Continuous Bulk Boundary (Macro Cliff)]\n")
  if (!is.na(top_of_bulk_idx)) {
    cat(sprintf(" -> %-43s : Index %d\n", "Identified Top of Idiosyncratic", top_of_bulk_idx))
    cat("    Informational Variance Bulk\n")
    cat(sprintf(" -> %-43s : %.6f\n", "Macro Gap Baseline (Max Bulk Gap)", macro_max_bulk_gap))
    cat(sprintf(" -> %-43s : %.2fx Baseline\n", "Actual Macro Gap Magnitude", macro_gap_ratio))
  } else {
    cat(sprintf(" -> %-43s : %s\n", "Identified Top of Idiosyncratic Informational Variance Bulk", "Failed (Defaulted to Kaiser)"))
  }
  cat("-----------------------------------------------------------------\n")
  cat(" [Engine A: Log-Gap]\n")
  cat(sprintf(" -> %-43s : %d%s\n", "Observed Generative Rank (K_roots)", K_log_gap, ifelse(log_gap_fallback, " [Forced by Macro Boundary]", "")))
  if (!is.na(log_gap_ratio)) {
    cat(sprintf(" -> %-43s : %.2fx (%.1f%% Drop)\n", "Log-Gap Magnitude", log_gap_ratio, log_gap_pct_drop))
  }
  cat("-----------------------------------------------------------------\n")
  cat(sprintf(" [Engine B: Triple-Tap (%d-Point Linear)]\n", triple_tap_window))
  cat(sprintf(" -> %-43s : %d%s\n", "Extended Signal Tail Rank (K_extended)", K_triple_tap, ifelse(triple_tap_fallback, " [Failed to trigger; defaulted to Log-Gap]", "")))
  cat(sprintf(" -> %-43s : %s\n", "Macro-Stitch Applied", ifelse(stitch_applied, "Yes (+0.5x Macro Gap to Bulk)", "No")))
  if (!is.na(triple_tap_expected_val)) {
    cat(sprintf(" -> %-43s : %.5f\n", "Expected Local Eigenvalue", triple_tap_expected_val))
    cat(sprintf(" -> %-43s : %.5f\n", "Actual Local Eigenvalue", triple_tap_actual_val))
    cat(sprintf(" -> %-43s : %.6f\n", "Calculated Local Baseline Sigma (SD)", triple_tap_local_sd))
  }
  if (!is.na(triple_tap_multiplier)) {
    cat(sprintf(" -> %-43s : %.2f (1-alpha = %.6f)\n", "Required Local t-Multiplier", triple_tap_multiplier, prob_target))
  }
  if (!is.na(triple_tap_actual_sigma)) {
    cat(sprintf(" -> %-43s : %.2f Sigma Breakout\n", "Actual vs Expected", triple_tap_actual_sigma))
  }
  cat("===================================================================================\n")
  
  # ============================================================================
  # WAVE 2: INTERACTIVE USER OVERRIDE (OR DEFAULTING)
  # ============================================================================
  # Default Assignments for non-interactive execution
  K_roots <- K_log_gap
  K_extended <- max(K_log_gap, K_triple_tap)
  
  if (interactive_mode) {
    cat("\n[WARNING]: The automated extractors rely on statistical heuristics and\n")
    cat("may not perfectly align with the true structural elbow of your specific dataset.\n")
    cat("Please visually examine the generated scree plot.\n\n")
    
    while (TRUE) {
      # 1. Prompt for Observed Generative Rank
      ans_r <- trimws(readline(prompt = sprintf("Enter the Observed Generative Rank (K_roots) [Press Enter to keep %d]: ", K_roots)))
      if (ans_r != "") {
        parsed_r <- suppressWarnings(as.integer(ans_r))
        if (!is.na(parsed_r) && parsed_r > 0 && parsed_r <= m_valid) K_roots <- parsed_r
      }
      
      # 2. Prompt for Extended Signal Tail Rank
      ans_ext <- trimws(readline(prompt = sprintf("Enter the rank of the Extended Signal Tail [Press Enter to keep %d]: ", K_extended)))
      if (ans_ext != "") {
        parsed_ext <- suppressWarnings(as.integer(ans_ext))
        if (!is.na(parsed_ext) && parsed_ext > 0 && parsed_ext <= m_valid) K_extended <- parsed_ext
      }
      
      if (K_roots > K_extended) {
        cat("[-] Warning: The Extended Signal Tail Rank should be >= the Observed Generative Rank.\n")
      }
      
      # ====================================================================
      # INTERACTIVE GRAPH PREVIEW & METRICS
      # ====================================================================
      if (requireNamespace("ggplot2", quietly = TRUE)) {
        
        # Setup previous engine lines (faded)
        if (K_log_gap == K_triple_tap) {
          prev_lines <- ggplot2::geom_vline(xintercept = K_log_gap, color = "gray60", linetype = "dashed", linewidth = 1)
        } else {
          prev_lines <- list(
            ggplot2::geom_vline(xintercept = K_log_gap, color = "gray60", linetype = "dashed", linewidth = 1),
            ggplot2::geom_vline(xintercept = K_triple_tap, color = "gray70", linetype = "twodash", linewidth = 1)
          )
        }
        
        # Draw user lines: Solid for roots, dashed for extended tail
        root_line <- ggplot2::geom_vline(xintercept = K_roots, color = "forestgreen", linetype = "solid", linewidth = 1.2)
        if (K_roots != K_extended) {
          ext_line <- ggplot2::geom_vline(xintercept = K_extended, color = "forestgreen", linetype = "dashed", linewidth = 1.2)
        } else {
          ext_line <- NULL
        }
        
        final_lines <- list(root_line, ext_line)
        
        # 1. ZOOMED VIEW
        zoom_start_upd <- max(1, min(K_log_gap, K_triple_tap, K_roots) - 5)
        zoom_end_upd <- min(max_plot_idx, max(K_log_gap, K_triple_tap, K_extended) + 15)
        
        plot_df_zoom_upd <- data.frame(Rank = zoom_start_upd:zoom_end_upd, Eigenvalue = eig_vals[zoom_start_upd:zoom_end_upd])
        
        p_scree_zoom_upd <- ggplot2::ggplot(plot_df_zoom_upd, ggplot2::aes(x = Rank, y = Eigenvalue)) +
          ggplot2::geom_line(color = "dodgerblue", linewidth = 1) +
          ggplot2::geom_point(color = "dodgerblue", size = 2) +
          prev_lines + final_lines +
          ggplot2::scale_y_continuous(trans = 'log10') +
          ggplot2::scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(min(x), max(x)))))) +
          ggplot2::labs(title = "Zoomed View", x = "Eigenvalue Index (m)", y = "Log(Eigenvalue)") +
          ggplot2::theme_minimal(base_size = 14) +
          ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 12))
        
        # 2. MACRO VIEW
        macro_base_upd <- max(50, K_extended * 10, K_log_gap * 10)
        if (!is.na(top_of_bulk_idx)) {
          macro_base_upd <- max(macro_base_upd, top_of_bulk_idx + 25)
        }
        macro_end_upd <- min(max_plot_idx, macro_base_upd)
        
        plot_df_macro_upd <- data.frame(Rank = 1:macro_end_upd, Eigenvalue = eig_vals[1:macro_end_upd])
        macro_y_max <- if(length(eig_vals) >= 2) eig_vals[2] * 1.1 else max(eig_vals)
        
        p_scree_macro_upd <- ggplot2::ggplot(plot_df_macro_upd, ggplot2::aes(x = Rank, y = Eigenvalue)) +
          ggplot2::geom_line(color = "dodgerblue", linewidth = 1) +
          ggplot2::geom_point(color = "dodgerblue", size = 2) +
          prev_lines + final_lines +
          ggplot2::scale_y_continuous(trans = 'log10') +
          ggplot2::coord_cartesian(ylim = c(NA, macro_y_max)) +
          ggplot2::labs(title = "Macro View", x = "Eigenvalue Index (m)", y = "Log(Eigenvalue)") +
          ggplot2::theme_minimal(base_size = 14) +
          ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 12))
        
        # Render Side-by-Side
        if (requireNamespace("patchwork", quietly = TRUE)) {
          combined_plot_upd <- (p_scree_macro_upd + p_scree_zoom_upd) +
            patchwork::plot_annotation(
              title = "Entropic Scree Results",
              subtitle = sprintf("User Confirmed: Observed Generative Rank (Solid) = %d | Extended Signal Tail (Dashed) = %d", K_roots, K_extended),
              theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 16, hjust = 0.5),
                                     plot.subtitle = ggplot2::element_text(size = 14, hjust = 0.5))
            )
          print(combined_plot_upd)
        } else {
          print(p_scree_macro_upd)
          print(p_scree_zoom_upd)
        }
      }
      
      # --- PREVIEW GRAVITY CALCULATIONS ---
      # 1. Total Volume is aggregated up to K_extended
      signal_variance_prev <- sum(eig_vals[1:K_extended])
      signal_weight_prev <- signal_variance_prev / m_plus
      unique_signal_volume_prev <- R_eff * signal_weight_prev
      total_signal_volume_prev <- unique_signal_volume_prev + redundant_signal_volume
      
      # 2. Rebundle Gravity purely into K_roots
      AIG_prev <- total_signal_volume_prev / K_roots
      core_eigenvals_prev <- eig_vals[1:K_roots]
      p_core_prev <- core_eigenvals_prev / sum(core_eigenvals_prev)
      FSIG_prev <- p_core_prev * total_signal_volume_prev
      
      cat("\n===================================================================================\n")
      cat(sprintf(" (PREVIEW) ELBOW METRICS (K_roots = %d, K_extended = %d) [K_extended Volume Rebundled to K_roots]\n", K_roots, K_extended))
      cat("===================================================================================\n")
      STP_prev <- FSIG_prev / FSIG_prev[1]
      cat(sprintf(" -> Average Informational Gravity (AIG): %.3f Variable Equivalents\n", AIG_prev))
      cat(" -> Factor-Specific Informational Gravity (FSIG):\n")
      print(round(FSIG_prev, 3))
      cat(" -> Structural Topology Profile (Relative to FSIG_1):\n")
      print(round(STP_prev, 3))
      cat("===================================================================================\n\n")
      
      ans_confirm <- trimws(readline(prompt = "Do you want to finalize these ranks? (Y/N): "))
      if (tolower(ans_confirm) %in% c("y", "yes")) {
        cat("\n[+] Finalizing rank selections.\n")
        break
      }
    }
  }
  
  # --- FINAL GRAVITY CALCULATIONS ---
  # 1. Capture Total Shared Signal Volume strictly up to the Extended Signal Tail
  signal_variance <- sum(eig_vals[1:K_extended])
  signal_weight <- signal_variance / m_plus
  unique_signal_volume <- R_eff * signal_weight
  total_signal_volume <- unique_signal_volume + redundant_signal_volume
  
  # 2. Rebundle AIG and FSIG strictly into the Observed Generative Rank (K_roots)
  AIG <- total_signal_volume / K_roots
  core_eigenvals <- eig_vals[1:K_roots]
  p_core <- core_eigenvals / sum(core_eigenvals)
  FSIG_final <- p_core * total_signal_volume
  
  elbow_label <- ifelse(interactive_mode, "User-Confirmed", "Automated Default")
  
  cat("\n===================================================================================\n")
  cat(sprintf(" (FINAL) ENTROPIC SCREE METRICS (based on %s)\n", elbow_label))
  cat("===================================================================================\n")
  STP_final <- FSIG_final / FSIG_final[1]
  cat(" [Rank Estimates]\n")
  cat(sprintf(" -> Observed Generative Rank (K_roots)     : %d\n", K_roots))
  cat(sprintf(" -> Extended Signal Tail Rank (K_extended) : %d\n", K_extended))
  cat("-----------------------------------------------------------------\n")
  cat(" [Informational Gravity Metrics Estimates] (K_extended Volume Rebundled to K_roots)\n")
  cat(sprintf(" -> Average Informational Gravity (AIG): %.3f Variable Equivalents\n", AIG))
  cat(" -> Factor-Specific Informational Gravity (FSIG):\n")
  print(round(FSIG_final, 3))
  cat(" -> Structural Topology Profile (Relative to FSIG_1):\n")
  print(round(STP_final, 3))
  cat("===================================================================================\n\n")
  
  # ============================================================================
  # WAVE 3: FINAL TRIPARTITE STRUCTURAL COMPOSITION
  # ============================================================================
  idiosyncratic_variance <- if (K_extended < m_valid) sum(eig_vals[(K_extended + 1):m_valid]) else 0
  idiosyncratic_weight <- idiosyncratic_variance / m_plus
  idiosyncratic_volume <- R_eff * idiosyncratic_weight
  
  cat("===================================================================================\n")
  cat(sprintf(" (FINAL) TRIPARTITE STRUCTURAL COMPOSITION (Rebundling Based on K_extended = %d)\n", K_extended))
  cat("===================================================================================\n")
  cat(sprintf(" -> %-50s : %d\n", "Valid Variables (m)", m_valid))
  cat(sprintf(" -> %-50s : %.2f\n", "Centered Trace (Tr_Mc)", Tr_Mc))
  cat("-----------------------------------------------------------------\n")
  cat(sprintf(" -> %-50s : %.3f\n", "Total Shared Signal Volume (Unique + Redundant)", total_signal_volume))
  cat(sprintf(" -> %-40s : %.3f\n", "   (Unique Signal Volume)", unique_signal_volume))
  cat(sprintf(" -> %-40s : %.3f\n", "   (Redundant Signal Volume)", redundant_signal_volume))
  cat("-----------------------------\n")
  cat(sprintf(" -> %-50s : %.3f\n", "Idiosyncratic Informational Variance", idiosyncratic_volume))
  cat("      (Structural Uncertainty + Independent Measurement Error\n")
  cat("      + Unshared Signal Geometry)\n")
  cat("-----------------------------------------------------------------\n")
  cat(sprintf(" -> %-50s : %.3f%%\n", "% Total Shared Signal", (total_signal_volume / m_valid) * 100))
  cat(sprintf(" -> %-40s : %.3f%%\n", "   (% Unique Signal)", (unique_signal_volume / m_valid) * 100))
  cat(sprintf(" -> %-40s : %.3f%%\n", "   (% Redundant Signal)", (redundant_signal_volume / m_valid) * 100))
  cat("-----------------------------\n")
  cat(sprintf(" -> %-50s : %.3f%%\n", "% Idiosyncratic Informational Variance", (idiosyncratic_volume / m_valid) * 100))
  cat("      (Structural Uncertainty + Independent Measurement Error\n")
  cat("      + Unshared Signal Geometry)\n")
  cat("===================================================================================\n\n")
  
  # ============================================================================
  # WAVE 4: METHODOLOGICAL REFERENCE
  # ============================================================================
  cat("===================================================================================\n")
  cat(" ENTROPIC SCREE (v1.0.0) - METHODOLOGICAL REFERENCE & LICENSE\n")
  cat("===================================================================================\n")
  cat(" -> Framework developed by Terrence J. Lee-St. John (Enli)\n")
  cat(" -> Released under the Apache License 2.0 (Open Source)\n")
  cat(" -> For full methods and metric definitions, see:\n")
  cat("    The Entropic Scree (2026) - https://doi.org/10.5281/zenodo.22028087\n")
  cat("===================================================================================\n\n")
  
  # ============================================================================
  # EXTENDED Factor-Specific Informational Gravity (FSIG) CALCULATIONS
  # ============================================================================
  
  # --- EXTENDED MODEL A: MACRO BULK BOUNDARY ---
  top_bulk_safe <- if(!is.na(top_of_bulk_idx)) top_of_bulk_idx else valid_k
  extended_bulk_k <- max(K_extended, top_bulk_safe - 1)
  sig_var_bulk <- sum(eig_vals[1:extended_bulk_k])
  total_sig_vol_bulk <- (R_eff * (sig_var_bulk / m_plus)) + redundant_signal_volume
  # Rebundle gravity strictly into the Observed Generative Rank (K_roots):
  FSIG_extended_bulk <- p_core * total_sig_vol_bulk
  
  # --- EXTENDED MODEL B: KAISER RULE BOUNDARY ---
  extended_kaiser_k <- max(K_extended, valid_k)
  sig_var_kaiser <- sum(eig_vals[1:extended_kaiser_k])
  total_sig_vol_kaiser <- (R_eff * (sig_var_kaiser / m_plus)) + redundant_signal_volume
  # Rebundle gravity strictly into the Observed Generative Rank (K_roots):
  FSIG_extended_kaiser <- p_core * total_sig_vol_kaiser
  
  # ============================================================================
  # EXTRACT BIPOLAR MODULES (TOPOLOGICAL POLES)
  # ============================================================================
  bipolar_modules_out <- NULL
  if (extract_bipolar_modules && !is.null(vectors_out) && K_roots > 0) {
    bipolar_modules_out <- list()
    
    for (k in 1:K_roots) {
      # Bind variable names to the specific eigenvector
      loadings <- setNames(vectors_out[, k], valid_vars)
      
      # Extract the extreme positive anchor
      pos_pole <- sort(loadings[loadings > 0], decreasing = TRUE)
      pos_pole <- head(pos_pole, bipolar_top_n)
      
      # Extract the extreme negative anchor (structural estrangement)
      neg_pole <- sort(loadings[loadings < 0], decreasing = FALSE)
      neg_pole <- head(neg_pole, bipolar_top_n)
      
      bipolar_modules_out[[paste0("Factor_", k)]] <- list(
        Positive_Anchor = pos_pole,
        Negative_Anchor = neg_pole
      )
    }
  }
  
  return(list(
    eigenvalues = eig_vals,
    similarity_matrix = NMI_mat_c,
    retained_features = valid_vars,
    bin_distributions = bin_sample_sizes,
    R_eff = R_eff,
    K_log_gap = K_log_gap,
    K_triple_tap = K_triple_tap,
    triple_tap_multiplier = triple_tap_multiplier,
    K_roots = K_roots,
    K_extended = K_extended,
    top_of_bulk = top_bulk_safe,
    total_signal_volume = total_signal_volume,
    unique_signal_volume = unique_signal_volume,
    redundant_signal_volume = redundant_signal_volume,
    idiosyncratic_volume = idiosyncratic_volume,
    AIG = AIG,
    FSIG_final = FSIG_final,
    structural_topology_profile = STP_final,
    FSIG_extended_bulk = FSIG_extended_bulk,
    FSIG_extended_kaiser = FSIG_extended_kaiser,
    eigenvectors = vectors_out,
    bipolar_modules = bipolar_modules_out,
    processed_data = processed_data_out
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
generate_true_mixed_proxies <- function(s1_continuous, m_proxies,
                                        max_interaction_order = 3, 
                                        max_polynomial_order = 3, 
                                        int_scaling = 1, 
                                        continuous_ratio = 1.0, 
                                        MIN_SPARSITY_RATIO = 0.01, 
                                        MAX_SPARSITY_RATIO = 0.10, 
                                        MIN_ROOT_SPARSITY_RATIO = 0.10, 
                                        MAX_ROOT_SPARSITY_RATIO = 0.50, 
                                        structural_snr = 5.0, 
                                        pure_static_ratio = 0.20, 
                                        pure_root_ratio = 0.10) {
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
      # FIX: Added check.names = FALSE to prevent R from changing ':' to '.'
      design_mat <- cbind(design_mat, as.matrix(as.data.frame(power_list, check.names = FALSE)))
    }
  }

  # --- STANDARDIZE THE ENTIRE DESIGN MATRIX ---
  # Every term (main effects and all interactions/powers) now has Variance = 1.0
  design_mat <- scale(design_mat)
  # -------------------------------------------------

  term_names <- colnames(design_mat)
  interaction_orders <- stringr::str_count(term_names, ":") + 1
  n_terms <- ncol(design_mat)

  # Parse every term in the design matrix to map its dependent roots
  term_dependencies <- lapply(term_names, function(name) {
    matches <- stringr::str_extract_all(name, "X\\d+")[[1]]
    unique(as.integer(gsub("X", "", matches)))
  })

  # --- KNOBS ---
  term_sds <- rep(1.0, n_terms)
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

  # --- DYNAMIC SPARSITY, PURE ROOTS & PURE STATIC IMPLEMENTATION ---
  m_static <- floor(m_proxies * pure_static_ratio)
  m_pure_root <- floor(m_proxies * pure_root_ratio)
  m_mixed <- m_proxies - m_static - m_pure_root

  min_roots_per_proxy <- max(1, round(k * MIN_ROOT_SPARSITY_RATIO))
  max_roots_per_proxy <- max(min_roots_per_proxy, round(k * MAX_ROOT_SPARSITY_RATIO))

  cat(sprintf("      -> Injecting %d Pure Static Proxies (%.1f%% of total dataset)...\n", m_static, pure_static_ratio * 100))
  cat(sprintf("      -> Injecting %d Pure Root Proxies (%.1f%% of total dataset, 1 main effect only)...\n", m_pure_root, pure_root_ratio * 100))
  cat(sprintf("      -> Generating %d Entangled Mixed Proxies...\n", m_mixed))
  cat(sprintf("         * Activating %d to %d Base Roots per proxy (%.0f%% to %.0f%% Root Sparsity)\n", min_roots_per_proxy, max_roots_per_proxy, MIN_ROOT_SPARSITY_RATIO * 100, MAX_ROOT_SPARSITY_RATIO * 100))
  cat(sprintf("         * Randomly sampling %.0f%% to %.0f%% of ALL aligned terms (Main & Higher-Order).\n", MIN_SPARSITY_RATIO * 100, MAX_SPARSITY_RATIO * 100))

  mask <- matrix(0, nrow = n_terms, ncol = m_proxies)
  main_effect_indices <- which(is_main)

  # 1. Assign Pure Root Proxies (Exactly 1 term, must be a main effect)
  if (m_pure_root > 0) {
    for (j in 1:m_pure_root) {
      root_idx <- sample(main_effect_indices, 1)
      mask[root_idx, j] <- 1
    }
  }

  # 2. Assign Entangled Mixed Proxies (Unrestricted Generative Sampling)
  if (m_mixed > 0) {
    start_col <- m_pure_root + 1
    end_col <- m_pure_root + m_mixed
    for (j in start_col:end_col) {

      # Step A: Determine how many base roots this proxy gets
      if (min_roots_per_proxy == max_roots_per_proxy) {
        num_roots <- min_roots_per_proxy
      } else {
        num_roots <- sample(min_roots_per_proxy:max_roots_per_proxy, 1)
      }

      # Step B: Select the specific roots
      selected_roots <- sample(1:k, num_roots)

      # Step C & D: Identify ALL aligned terms (both main effects and higher-order interactions)
      aligned_term_indices <- which(sapply(term_dependencies, function(deps) {
        all(deps %in% selected_roots)
      }))

      # Step E: Sample a random amount BETWEEN Min% and Max% of ALL aligned terms
      if (length(aligned_term_indices) > 0) {
        # Calculate the absolute floors and ceilings
        absolute_min_to_sample <- max(1, round(length(aligned_term_indices) * MIN_SPARSITY_RATIO))
        absolute_max_to_sample <- max(absolute_min_to_sample, round(length(aligned_term_indices) * MAX_SPARSITY_RATIO))

        # Randomly pick how many terms to actually sample (protecting against min == max)
        if (absolute_min_to_sample == absolute_max_to_sample) {
          n_to_sample <- absolute_min_to_sample
        } else {
          n_to_sample <- sample(absolute_min_to_sample:absolute_max_to_sample, 1)
        }

        if (length(aligned_term_indices) == 1) {
          selected_terms <- aligned_term_indices
        } else {
          selected_terms <- sample(aligned_term_indices, n_to_sample)
        }
        mask[selected_terms, j] <- 1
      }
    }
  }

  # Randomly shuffle the columns so the pure static and pure root variables are randomly dispersed
  mask <- mask[, sample(ncol(mask))]

  coeffs <- coeffs * mask
  # ---------------------------------------

  # --- CALCULATE K_eff STRICTLY ON REALIZED TERMS ---
  active_indices <- which(rowSums(abs(coeffs)) > 0)
  active_terms <- length(active_indices)

  if (active_terms > 0) {
    # Apply scaling to active columns to measure theoretical continuous volume
    active_design <- design_mat[, active_indices, drop = FALSE]
    active_orders <- interaction_orders[active_indices]
    active_weights <- ifelse(active_orders == 1, term_sds[active_indices], term_sds[active_indices] * int_scaling)

    muzzled_active_design <- t(t(active_design) * active_weights)

    expanded_cov <- cov(muzzled_active_design)
    expanded_eigen <- pmax(eigen(expanded_cov, symmetric = TRUE)$values, 0)
    p_expanded <- expanded_eigen[expanded_eigen > 1e-9] / sum(expanded_eigen[expanded_eigen > 1e-9])
    K_eff <- exp(-sum(p_expanded * log(p_expanded)))
  } else {
    K_eff <- 0
  }
  # --------------------------------------------------------

  # 3. Generate the Deterministic Core
  deterministic_core <- design_mat %*% coeffs

  # 3.5 SET THRESHOLDS & INJECT STRUCTURAL UNCERTAINTY (ROUND 1 NOISE)
  cat(sprintf("      -> Setting Base Rates and Injecting Structural Uncertainty (Round 1 Noise, SNR = %.2f)...\n", structural_snr))
  true_manifestation <- matrix(0, nrow = n, ncol = m_proxies)

  # Pre-calculate standard deviations for exact noise scaling
  # Signal Var = 1.0, Noise Var = 1/SNR, Total Var = 1.0 + (1/SNR)
  total_sd <- sqrt(1 + (1 / structural_snr))

  for (j in 1:m_proxies) {
    # Step A: Standardize the deterministic core (Forces pure signal variance to exactly 1.0)
    core_vec <- deterministic_core[, j]
    if (sd(core_vec) > 1e-9) {
      core_vec <- as.vector(scale(core_vec))
    } else {
      core_vec <- rep(0, n)
    }

    # Step B: Set the Threshold PRIOR to Structural Noise (For Binary Variables Only)
    is_binary_var <- (j > m_cont)
    if (is_binary_var) {
      # 1. Target a prevalence between 10% and 90%
      target_prev <- runif(1, min = 0.10, max = 0.90)

      # 2. Calculate the exact shift needed to hit this probability
      #    considering the noise that is about to be added.
      shift <- sqrt(2) * total_sd * qnorm(target_prev)

      # 3. Embed the threshold permanently into the structural manifold
      core_vec <- core_vec + shift
    }

    # Step C: Inject Structural Uncertainty (Variance = 1 / SNR)
    structural_noise <- rnorm(n, mean = 0, sd = sqrt(1 / structural_snr))

    # Step D: Finalize the True Probabilistic Manifestation
    true_manifestation[, j] <- core_vec + structural_noise
  }

  # 4A. The True Continuous Signal
  cat(sprintf("      -> Generating %d True Continuous Proxies...\n", m_cont))
  true_cont <- true_manifestation[, 1:m_cont, drop = FALSE]

  # 4B. The True Binary Signal
  if (m_bin > 0) {
    cat(sprintf("      -> Generating %d True Binary Proxies...\n", m_bin))
    signal_bin <- true_manifestation[, (m_cont + 1):m_proxies, drop = FALSE]

    apply_copula_mapping <- function(scores) {
      # Fallback for zero-variance
      if (sd(scores) < 1e-9) return(rep(runif(1, min = 0.50, max = 0.95), length(scores)))

      # Because the threshold shift and structural noise are already permanently baked
      # into the scores, we NO LONGER scale() them here. We simply map them against
      # the theoretical baseline distribution N(0, total_sd) to get the final probability.
      probs <- pnorm(scores, mean = 0, sd = total_sd)

      return(probs)
    }

    prob_mat <- apply(signal_bin, 2, apply_copula_mapping)
    true_bin <- matrix(rbinom(length(prob_mat), 1, prob_mat), nrow = n, ncol = m_bin)

    true_proxies <- cbind(true_cont, true_bin)
    is_continuous <- c(rep(TRUE, m_cont), rep(FALSE, m_bin))
  } else {
    true_proxies <- true_cont
    is_continuous <- rep(TRUE, m_cont)
  }

  # 5. Randomly shuffle the columns
  mix_idx <- sample(m_proxies)
  true_proxies <- true_proxies[, mix_idx, drop = FALSE]
  is_continuous <- is_continuous[mix_idx]

  return(list(
    data_matrix = true_proxies,
    is_continuous = is_continuous,
    active_terms = active_terms,
    K_eff = K_eff
  ))
}

# STEP 2B: Apply independent Measurement Error to the True Data (ROUND 2 NOISE)
apply_measurement_error <- function(true_universe, snr_continuous = 2.0, binary_error_rate = 0.15) {
  cat(sprintf("      -> Applying Measurement Error (Round 2 Noise: Continuous SNR = %.2f, Binary Bit-Flip Rate = %.3f)...\n", snr_continuous, binary_error_rate))

  obs_mat <- true_universe$data_matrix
  is_cont <- true_universe$is_continuous
  n <- nrow(obs_mat)
  m <- ncol(obs_mat)

  for (j in 1:m) {
    if (is_cont[j]) {
      # Add Gaussian Measurement Error
      true_var <- var(obs_mat[, j])
      if (true_var < 1e-9) true_var <- 1e-9

      measurement_error <- rnorm(n, mean = 0, sd = 1)
      error_var <- var(measurement_error)

      scaling_factor <- sqrt(true_var / (error_var * snr_continuous))
      obs_mat[, j] <- obs_mat[, j] + (measurement_error * scaling_factor)

    } else {
      # Add Bit-Flip Measurement Error
      flip_mask <- rbinom(n, 1, binary_error_rate)
      obs_mat[, j] <- abs(obs_mat[, j] - flip_mask)
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

# MACRO-LEVEL GENERATED DATA KNOBS
K_TRUE <- 20             # Intrinsic Generative Rank (r)
N_ROWS <- 10000          # Sample Size (N)
M_PROXIES <- 20000       # Dimensionality (m)
CONTINUOUS_RATIO <- 0.80 # SET TO 1.00 FOR PURE CONTINUOUS, 0.00 FOR PURE BINARY, OR ANYWHERE IN BETWEEN

# DESIGN MATRIX KNOBS
MAX_INTERACTION_ORDER <- 5 # Controls cross-interactions
MAX_POLYNOMIAL_ORDER <- 4  # Controls pure powers

# SPARSITY KNOBS
MIN_ROOT_SPARSITY_RATIO <- 0.50 # Minimum % of base roots a proxy will select
MAX_ROOT_SPARSITY_RATIO <- 0.50 # Maximum % of base roots a proxy will select
MIN_SPARSITY_RATIO <- 0.10      # Minimum % of aligned higher-order terms a proxy will sample
MAX_SPARSITY_RATIO <- 0.10      # Maximum % of aligned higher-order terms a proxy will sample
PURE_ROOT_RATIO <- 0.00         # Exactly X% of proxies will purely measure 1 main effect (Anchor Proxies)
PURE_STATIC_RATIO <- 0.05       # Exactly X% of proxies will be pure, completely idiosyncratic noise

# IDIOSYNCRATIC NOISE KNOBS
STRUCTURAL_SNR <- 10      # Structural Uncertainty: Determines how perfectly the proxy reflects the structural manifold (Idiosyncratic variance)
CONTINUOUS_SNR <- 2       # Pure continuous measurement error (Sensor degradation)
BINARY_ERROR_RATE <- 0.15 # Pure binary measurement error (Bit-flip)

cat(sprintf("Generating Mixed Synthetic Universe: %s rows, %d Proxies, %d Latent Drivers...\n", format(N_ROWS, big.mark=","), M_PROXIES, K_TRUE))

# 1. Generate Generative Space
# SET ROOT ENTANGLEMENT LEVEL HERE:
#   Inf   = Strictly Orthogonal (Identity Matrix)
#   15.0  = Low Entanglement
#   0.75  = Medium Entanglement
#   0.001 = High Entanglement
ALPHAD_PARAM <- Inf

if (ALPHAD_PARAM == Inf) {
  latent_sigma <- diag(K_TRUE)
} else {
  if (!requireNamespace("clusterGeneration", quietly = TRUE)) install.packages("clusterGeneration")

  # Generate a valid, heterogeneously correlated positive-definite matrix
  latent_sigma <- clusterGeneration::rcorrmatrix(K_TRUE, alphad = ALPHAD_PARAM)
}

avg_cor <- mean(abs(latent_sigma[upper.tri(latent_sigma)]))
cat(sprintf(" -> Generated Root Correlation Matrix (Avg Pairwise |r| = %.3f)\n", avg_cor))

# --- CALCULATE EFFECTIVE INTRINSIC GENERATIVE RANK (r_eff) ---
root_eigen <- pmax(eigen(latent_sigma, symmetric = TRUE)$values, 0)
p_root <- root_eigen[root_eigen > 1e-9] / sum(root_eigen[root_eigen > 1e-9])
effective_r <- exp(-sum(p_root * log(p_root)))
cat(sprintf(" -> Effective Intrinsic Generative Rank (r_eff) compressed to: %.2f\n", effective_r))
# ---------------------------------------------------

Z_latent_continuous <- MASS::mvrnorm(n = N_ROWS
                                     , mu = rep(0, K_TRUE)
                                     , Sigma = latent_sigma)

# 2. Generate the Data (Uncoupled Truth and Error)
true_universe <- generate_true_mixed_proxies(
  Z_latent_continuous,
  M_PROXIES,
  max_interaction_order = MAX_INTERACTION_ORDER,
  max_polynomial_order = MAX_POLYNOMIAL_ORDER,
  int_scaling = 1.0,
  continuous_ratio = CONTINUOUS_RATIO,
  MIN_SPARSITY_RATIO = MIN_SPARSITY_RATIO,
  MAX_SPARSITY_RATIO = MAX_SPARSITY_RATIO,
  MIN_ROOT_SPARSITY_RATIO = MIN_ROOT_SPARSITY_RATIO,
  MAX_ROOT_SPARSITY_RATIO = MAX_ROOT_SPARSITY_RATIO,
  structural_snr = STRUCTURAL_SNR,
  pure_static_ratio = PURE_STATIC_RATIO,
  pure_root_ratio = PURE_ROOT_RATIO
)

observed_data <- apply_measurement_error(true_universe, snr_continuous = CONTINUOUS_SNR, binary_error_rate = BINARY_ERROR_RATE)

cat("\nDataset is ready. Starting pipeline...\n\n")

# Print Latent Configurational Ranks
cat(sprintf("\n[***] r  (Intrinsic Generative Rank): %d\n\n", K_TRUE))
cat(sprintf("\n[***] r_eff  (Effective Intrinsic Generative Rank): %.2f\n", effective_r))
cat(sprintf("[***] K_rlzd (Realized Intrinsic Root Expansion Rank)   : %d\n\n", true_universe$active_terms))
cat(sprintf("[***] K_eff (Effective Intrinsic Root Expansion Rank)  : %.2f\n", true_universe$K_eff))

# ==============================================================================
# 4. RUN ENTROPIC SCREE
# ==============================================================================
results <- Entropic.Scree(observed_data
                                    , purge_constants = FALSE
                                    , check_collinearity = FALSE
                                    , interactive_mode = FALSE
)

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
# Helper function to safely pad (if m > N) or truncate (if N > m) to match m_total
safe_pad <- function(eig, m) {
  if (length(eig) > m) return(eig[1:m])
  if (length(eig) < m) return(c(eig, rep(0, m - length(eig))))
  return(eig)
}

pca_padded <- safe_pad(pca_eigenvalues, m_total)
spearman_padded <- safe_pad(spearman_eigenvalues, m_total)
kpca_padded <- safe_pad(kpca_eigenvalues, m_total)

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

  # Extract K_rlzd and K_eff
  K_rlzd <- true_universe$active_terms
  K_eff <- true_universe$K_eff
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
    ggplot2::theme_minimal(base_size = 12) + ggplot2::theme(plot.title
= ggplot2::element_text(face="bold"))

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
    ggplot2::theme_minimal(base_size = 12) + ggplot2::theme(plot.title
= ggplot2::element_text(face="bold"))

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

# ==============================================================================
# 7. MULTI-METHOD COMPARISON SUMMARY TABLE
# ==============================================================================
# Helper to find the naive empirical elbow for baselines (Max Log-Gap, skipping K=1)
get_naive_elbow <- function(eig_vals, max_rank) {
  valid_eigs <- eig_vals[eig_vals > 1e-8 & !is.na(eig_vals)]
  limit <- min(length(valid_eigs), max_rank)
  if (limit < 3) return(NA)
  gaps <- abs(diff(log(valid_eigs[1:limit])))
  # Skip the first gap (index 1 to 2) to force K >= 2, mirroring Entropic Scree logic
  if (length(gaps) >= 2) {
    return(which.max(gaps[-1]) + 1)
  }
  return(NA)
}

# Extract naive elbows
pca_elbow <- get_naive_elbow(pca_eigenvalues, n_rows - 1)
spearman_elbow <- get_naive_elbow(spearman_eigenvalues, n_rows - 1)
kpca_elbow <- get_naive_elbow(kpca_eigenvalues, n_rows - 1)

# Dynamic formatting function for the table output
format_elbow_result <- function(elbow_val, k_true, k_rlzd) {
  if (is.na(elbow_val)) return("Failed (No Gap)")
  if (elbow_val == k_true) return(sprintf("%d (Perfect Match)", elbow_val))
  # If it falls within 10% of K_rlzd, flag it as hitting the expansion ceiling
  if (abs(elbow_val - k_rlzd) <= max(2, k_rlzd * 0.1)) return(sprintf("%d (~K_rlzd)", elbow_val))
  return(sprintf("%d (Spurious)", elbow_val))
}

pca_result_str <- format_elbow_result(pca_elbow, K_TRUE, K_rlzd)
spearman_result_str <- format_elbow_result(spearman_elbow, K_TRUE, K_rlzd)
kpca_result_str <- format_elbow_result(kpca_elbow, K_TRUE, K_rlzd)
ent_result_str <- format_elbow_result(results$K_roots, K_TRUE, K_rlzd)

cat("\n========================================================================================================\n")
cat(" MULTI-METHOD DIMENSIONALITY EXTRACTION COMPARISON\n")
cat("========================================================================================================\n")
cat(sprintf(" -> True Intrinsic Generative Rank (r) : %d\n", K_TRUE))
cat(sprintf(" -> Realized Expansion Rank (K_rlzd)   : %d\n", K_rlzd))
cat("--------------------------------------------------------------------------------------\n")

# Format columns
cat(sprintf(" %-20s | %-14s | %-19s | %-20s\n",
            "Method", "Rank Ceiling", "Kaiser Extraction",
"Auto-Detected Rank"))
cat("----------------------|----------------|---------------------|------------------------\n")

# Row 1: Standard PCA
cat(sprintf(" %-20s | %-14d | %-19d | %-20s\n",
            "Standard PCA",
            (n_rows - 1),
            pca_kaiser,
            pca_result_str))

# Row 2: Spearman Rank PCA
cat(sprintf(" %-20s | %-14d | %-19d | %-20s\n",
            "Spearman Rank PCA",
            (n_rows - 1),
            spearman_kaiser,
            spearman_result_str))

# Row 3: Kernel PCA
cat(sprintf(" %-20s | %-14d | %-19d | %-20s\n",
            "Kernel PCA (RBF)",
            (n_rows - 1),
            kpca_kaiser,
            kpca_result_str))

# Row 4: Entropic Scree
cat(sprintf(" %-20s | %-14d | %-19d | %-20s\n",
            "Entropic Scree",
            (m_total - 1),
            ent_kaiser,
            ent_result_str))
cat("========================================================================================================\n\n")
