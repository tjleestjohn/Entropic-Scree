# ==============================================================================
# Update.Entropic.Scree (Post-Hoc Rank Modification)
# ==============================================================================
Update.Entropic.Scree <- function(scree_obj, new_K_roots = NULL, new_K_extended = NULL, bipolar_top_n = 0.20) {
  
  # 1. Input Validation
  if (!is.list(scree_obj) || is.null(scree_obj$eigenvalues)) {
    stop("Error: 'scree_obj' must be a valid output list from Entropic.Scree().")
  }
  
  eig_vals <- scree_obj$eigenvalues
  m_valid <- length(scree_obj$retained_features)
  
  # Resolve new ranks (default to the object's current ranks if not provided)
  K_roots <- if (!is.null(new_K_roots)) as.integer(new_K_roots) else scree_obj$K_roots
  K_extended <- if (!is.null(new_K_extended)) as.integer(new_K_extended) else scree_obj$K_extended
  
  # Safety bounds
  K_roots <- max(1, min(K_roots, m_valid))
  K_extended <- max(K_roots, min(K_extended, m_valid))
  
  if (new_K_roots > new_K_extended && !is.null(new_K_roots) && !is.null(new_K_extended)) {
    warning("K_roots is greater than K_extended. K_extended has been bounded to K_roots.")
  }

  # 2. Extract Static Core Physics from Object
  R_eff <- scree_obj$R_eff
  m_plus <- sum(eig_vals)
  redundant_signal_volume <- scree_obj$redundant_signal_volume
  vectors_out <- scree_obj$eigenvectors
  valid_vars <- scree_obj$retained_features
  
  # ============================================================================
  # 3. RECALCULATE GRAVITY & STRUCTURAL COMPOSITION
  # ============================================================================
  
  # Capture Total Shared Signal Volume strictly up to the NEW Extended Signal Tail
  signal_variance <- sum(eig_vals[1:K_extended])
  signal_weight <- signal_variance / m_plus
  unique_signal_volume <- R_eff * signal_weight
  total_signal_volume <- unique_signal_volume + redundant_signal_volume
  
  # Rebundle AIG and FSIG strictly into the NEW Observed Generative Rank
  AIG <- total_signal_volume / K_roots
  core_eigenvals <- eig_vals[1:K_roots]
  p_core <- core_eigenvals / sum(core_eigenvals)
  FSIG_final <- p_core * total_signal_volume
  STP_final <- FSIG_final / FSIG_final[1]
  
  # Idiosyncratic Variance
  idiosyncratic_variance <- if (K_extended < m_valid) sum(eig_vals[(K_extended + 1):m_valid]) else 0
  idiosyncratic_weight <- idiosyncratic_variance / m_plus
  idiosyncratic_volume <- R_eff * idiosyncratic_weight
  
  # ============================================================================
  # 4. RE-EXTRACT BIPOLAR MODULES
  # ============================================================================
  bipolar_modules_out <- NULL
  if (!is.null(vectors_out) && K_roots > 0) {
    bipolar_modules_out <- list()
    
    resolved_top_n <- if (bipolar_top_n > 0 && bipolar_top_n < 1) {
      max(1, floor(bipolar_top_n * m_valid))
    } else {
      max(1, as.integer(bipolar_top_n))
    }
    
    for (k in 1:K_roots) {
      loadings <- setNames(vectors_out[, k], valid_vars)
      
      pos_pole <- sort(loadings[loadings > 0], decreasing = TRUE)
      pos_pole <- head(pos_pole, resolved_top_n)
      
      neg_pole <- sort(loadings[loadings < 0], decreasing = FALSE)
      neg_pole <- head(neg_pole, resolved_top_n)
      
      bipolar_modules_out[[paste0("Factor_", k)]] <- list(
        Positive_Anchor = pos_pole,
        Negative_Anchor = neg_pole
      )
    }
  }

  # ============================================================================
  # 5. PRINT UPDATED DASHBOARD
  # ============================================================================
  cat("\n===================================================================================\n")
  cat(" (UPDATED) ENTROPIC SCREE METRICS (User Override)\n")
  cat("===================================================================================\n")
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

  # 6. Mutate and Return Object
  scree_obj$K_roots <- K_roots
  scree_obj$K_extended <- K_extended
  scree_obj$total_signal_volume <- total_signal_volume
  scree_obj$unique_signal_volume <- unique_signal_volume
  scree_obj$idiosyncratic_volume <- idiosyncratic_volume
  scree_obj$AIG <- AIG
  scree_obj$FSIG_final <- FSIG_final
  scree_obj$structural_topology_profile <- STP_final
  scree_obj$bipolar_modules <- bipolar_modules_out
  
  return(scree_obj)
}
