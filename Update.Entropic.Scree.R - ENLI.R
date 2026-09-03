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
  
  # Capture original state for the comparison printout
  orig_K_roots <- scree_obj$K_roots
  orig_K_extended <- scree_obj$K_extended
  orig_AIG <- scree_obj$AIG
  
  orig_total_signal <- scree_obj$total_signal_volume
  orig_unique_signal <- scree_obj$unique_signal_volume
  orig_redundant_signal <- scree_obj$redundant_signal_volume
  orig_idio_vol <- scree_obj$idiosyncratic_volume
  
  orig_pct_signal <- (orig_total_signal / m_valid) * 100
  orig_pct_unique <- (orig_unique_signal / m_valid) * 100
  orig_pct_redundant <- (orig_redundant_signal / m_valid) * 100
  orig_pct_idio <- (orig_idio_vol / m_valid) * 100
  
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
  Tr_Mc <- sum(diag(scree_obj$similarity_matrix))
  
  # ============================================================================
  # 3. RECALCULATE GRAVITY & STRUCTURAL COMPOSITION
  # ============================================================================
  
  # Capture Total Shared Signal Volume strictly up to the NEW Extended Signal Tail
  signal_variance <- sum(eig_vals[1:K_extended])
  signal_weight <- signal_variance / m_plus
  unique_signal_volume <- R_eff * signal_weight
  total_signal_volume <- unique_signal_volume + redundant_signal_volume
  new_pct_signal <- (total_signal_volume / m_valid) * 100
  new_pct_unique <- (unique_signal_volume / m_valid) * 100
  new_pct_redundant <- (redundant_signal_volume / m_valid) * 100
  
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
  new_pct_idio <- (idiosyncratic_volume / m_valid) * 100
  
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
  # 5. PRINT UPDATED DASHBOARDS
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

  cat("===================================================================================\n")
  cat(sprintf(" (UPDATED) TRIPARTITE STRUCTURAL COMPOSITION (Rebundling Based on K_extended = %d)\n", K_extended))
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
  cat(sprintf(" -> %-50s : %.3f%%\n", "% Total Shared Signal", new_pct_signal))
  cat(sprintf(" -> %-40s : %.3f%%\n", "   (% Unique Signal)", new_pct_unique))
  cat(sprintf(" -> %-40s : %.3f%%\n", "   (% Redundant Signal)", new_pct_redundant))
  cat("-----------------------------\n")
  cat(sprintf(" -> %-50s : %.3f%%\n", "% Idiosyncratic Informational Variance", new_pct_idio))
  cat("      (Structural Uncertainty + Independent Measurement Error\n")
  cat("      + Unshared Signal Geometry)\n")
  cat("===================================================================================\n\n")

  cat("===================================================================================\n")
  cat(" (DELTA) METRIC SHIFT COMPARISON\n")
  cat("===================================================================================\n")
  cat(sprintf(" -> %-37s : %d -> %d\n", "Observed Generative Rank (K_roots)", orig_K_roots, K_roots))
  cat(sprintf(" -> %-37s : %d -> %d\n", "Extended Signal Tail (K_extended)", orig_K_extended, K_extended))
  cat("-----------------------------------------------------------------\n")
  cat(sprintf(" -> %-37s : %.3f -> %.3f\n", "Average Info Gravity (AIG)", orig_AIG, AIG))
  cat("-----------------------------------------------------------------\n")
  cat(sprintf(" -> %-37s : %.3f -> %.3f\n", "Total Shared Signal Volume", orig_total_signal, total_signal_volume))
  cat(sprintf(" -> %-37s : %.3f -> %.3f\n", "   (Unique Signal Volume)", orig_unique_signal, unique_signal_volume))
  cat(sprintf(" -> %-37s : %.3f -> %.3f\n", "   (Redundant Signal Volume)", orig_redundant_signal, redundant_signal_volume))
  cat("-----------------------------\n")
  cat(sprintf(" -> %-37s : %.3f -> %.3f\n", "Idiosyncratic Variance Volume", orig_idio_vol, idiosyncratic_volume))
  cat("-----------------------------------------------------------------\n")
  cat(sprintf(" -> %-37s : %.3f%% -> %.3f%%\n", "% Total Shared Signal", orig_pct_signal, new_pct_signal))
  cat(sprintf(" -> %-37s : %.3f%% -> %.3f%%\n", "   (% Unique Signal)", orig_pct_unique, new_pct_unique))
  cat(sprintf(" -> %-37s : %.3f%% -> %.3f%%\n", "   (% Redundant Signal)", orig_pct_redundant, new_pct_redundant))
  cat("-----------------------------\n")
  cat(sprintf(" -> %-37s : %.3f%% -> %.3f%%\n", "% Idiosyncratic Variance", orig_pct_idio, new_pct_idio))
  cat("===================================================================================\n\n")

  cat("===================================================================================\n")
  cat(" ENTROPIC SCREE (v1.0.1) - METHODOLOGICAL REFERENCE & LICENSE\n")
  cat("===================================================================================\n")
  cat(" -> Framework developed by Terrence J. Lee-St. John (Enli)\n")
  cat(" -> Released under the Apache License 2.0 (Open Source)\n")
  cat(" -> For full methods and metric definitions, see:\n")
  cat("    The Entropic Scree (2026) - https://doi.org/10.5281/zenodo.22028087\n")
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
