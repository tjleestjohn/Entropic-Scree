# ENTROPIC SCREE ONLY: REPLICATED ROOT-COUNT STUDY (2026-10-01)
# Open this ONE file in RStudio, select all, and Run. R >= 4.2 required.
# 11 root counts x 3 entanglement conditions x 20 seeds = 660 datasets.
# N = 10,000 records; m = 20,000 generated measurements for EVERY dataset.
# Exact optimized simulation generator; direct published CRAN Entropic.Scree call.
# No PCA baselines, alternative elbow detector, or manual selection is run.
# Missing CRAN packages install automatically (Windows binaries when available).
# Source-only platforms need the standard compiler tools required by CRAN packages.
# Generated data are TEMPORARY and deleted after each verified successful result.
# Keep results/plots/seeds/settings; true-data and centered-matrix saving is OFF.
# Estimated RAM varies by root count; concurrency is chosen from available RAM.
# Default full-size study is NOT automatically run as part of development tests.
# Results go beneath getwd(); change output_dir to a large local drive if needed.
# Avoid an actively synchronized cloud folder for this large study.
# To resume after interruption, rerun with the SAME settings and output_dir.
# Keep the computer awake. Completed checkpoints are not recomputed.

setwd("C:/Users/t.leestjohn/OneDrive - IHLAD/Desktop/Sim")

local({
  # Optional: draw the summary from existing SEED-RESULTS, with no simulation.
  PLOT_ONLY <- FALSE
  PLOT_INPUT <- NULL                # NULL opens a file picker when PLOT_ONLY=TRUE.
  CONFIG <- list(
    roots = 6:26,
    cran_version = "1.0.1",
    n_rows = 5000L,
    m_proxies = 10000L,
    repetitions = 20L,
    seed = 19862026L,
    continuous_ratio = 0.50,
    max_interaction_order = 5L,
    max_polynomial_order = 4L,
    structural_snr = 4,
    continuous_snr = 4,
    binary_error_rate = 0.0528, # Reliability parity: preserves R^2 = 0.80 across continuous & binary
    pure_static_ratio = 0.05,
    pure_root_ratio = 0,
    term_fraction = 0.10,
    output_dir = file.path(getwd(), "Entropic-Scree-CRAN-N10000-M20000"),
    max_workers = NULL,             # NULL: floor(usable CPU threads); optional lower cap.
    reserve_ram_fraction = 0.00,
    minimum_reserve_gb = 1.5,
    memory_safety_factor = 1.20,
    poll_seconds = 3,
    progress_seconds = 30,
    # Generated data are only restart checkpoints, removed after verified success.
    # Optional true-data/centered-matrix files are OFF; if enabled they are temporary too.
    save_true_data = FALSE,
    save_similarity_matrix = FALSE,
    smoke_test = FALSE              # Explicit reduced test only; leave FALSE for paper.
  )
  
  
  ENGINE_CODE <- r"----(
# Generate latent-term mixtures, structural uncertainty and measurement error.
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
      # Preserve term names: check.names = FALSE prevents changing ':' to '.'
      design_mat <- cbind(design_mat, as.matrix(as.data.frame(power_list, check.names = FALSE)))
    }
  }

  # --- STANDARDIZE THE ENTIRE DESIGN MATRIX ---
  # Standardize each main effect, interaction and polynomial term to unit variance.
  # Column-wise standardization is identical to scale(), in bounded blocks.
  term_center <- term_scale <- numeric(ncol(design_mat))
  for (begin in seq.int(1L,ncol(design_mat),by=256L)) {
    ids <- begin:min(ncol(design_mat),begin+255L)
    block <- scale(design_mat[,ids,drop=FALSE])
    design_mat[,ids] <- block
    term_center[ids] <- attr(block,"scaled:center")
    term_scale[ids] <- attr(block,"scaled:scale")
  }
  names(term_center) <- names(term_scale) <- colnames(design_mat)
  attr(design_mat,"scaled:center") <- term_center
  attr(design_mat,"scaled:scale") <- term_scale
  rm(block,s1_df); gc(FALSE)

  term_names <- colnames(design_mat)
  interaction_orders <- stringr::str_count(term_names, ":") + 1
  n_terms <- ncol(design_mat)

  # Parse every term in the design matrix to map its dependent roots
  term_dependencies <- lapply(term_names, function(name) {
    matches <- stringr::str_extract_all(name, "X\\d+")[[1]]
    unique(as.integer(gsub("X", "", matches)))
  })

  # Generating-term weights.
  term_sds <- rep(1.0, n_terms)

  # 2. Build the Weight Matrix
  # Preserve the specified random-number stream without allocating a dense
  # term-by-measurement coefficient matrix. A second pass recovers only the
  # coefficients actually selected, using the saved RNG state.
  coefficient_rng <- .Random.seed
  is_main <- interaction_orders == 1
  is_int <- interaction_orders > 1
  n_main <- sum(is_main)
  n_int <- sum(is_int)

  blueprint_roots <- vector("list", m_proxies)
  for(j in 1:m_proxies) {
    # Main effects
    if (n_main > 0) invisible(rnorm(n_main, mean = 0, sd = term_sds[is_main]))
    # Interactions and Polynomials
    if (n_int > 0) invisible(rnorm(n_int, mean = 0, sd = term_sds[is_int] * int_scaling))
  }

  # Select generating terms and root participation.
  m_static <- floor(m_proxies * pure_static_ratio)
  m_pure_root <- floor(m_proxies * pure_root_ratio)
  m_mixed <- m_proxies - m_static - m_pure_root

  min_roots_per_proxy <- max(1, ceiling(k * MIN_ROOT_SPARSITY_RATIO))
  max_roots_per_proxy <- max(min_roots_per_proxy, ceiling(k * MAX_ROOT_SPARSITY_RATIO))

  cat(sprintf("      -> Injecting %d Pure Static Proxies (%.1f%% of total dataset)...\n", m_static, pure_static_ratio * 100))
  cat(sprintf("      -> Injecting %d Pure Root Proxies (%.1f%% of total dataset, 1 main effect only)...\n", m_pure_root, pure_root_ratio * 100))
  cat(sprintf("      -> Generating %d Entangled Mixed Proxies...\n", m_mixed))
  cat(sprintf("         * Activating %d to %d Base Roots per proxy (%.0f%% to %.0f%% Root Sparsity)\n", min_roots_per_proxy, max_roots_per_proxy, MIN_ROOT_SPARSITY_RATIO * 100, MAX_ROOT_SPARSITY_RATIO * 100))
  cat(sprintf("         * Randomly sampling %.0f%% to %.0f%% of ALL aligned terms (Main & Higher-Order).\n", MIN_SPARSITY_RATIO * 100, MAX_SPARSITY_RATIO * 100))

  selected_by_column <- rep(list(integer()),m_proxies)
  # Exact set containment using integer bit masks (not a statistical heuristic).
  # Retain the general set test when more than 30 roots would overflow integers.
  dependency_bits <- if (k<=30L) vapply(term_dependencies,
    function(deps) as.integer(sum(2^(deps-1L))),integer(1)) else NULL
  aligned_cache <- new.env(hash=TRUE,parent=emptyenv())
  main_effect_indices <- which(is_main)

  # 1. Assign Pure Root Proxies (Exactly 1 term, must be a main effect)
  if (m_pure_root > 0) {
    for (j in 1:m_pure_root) {
      root_idx <- sample(main_effect_indices, 1)
      selected_by_column[[j]] <- root_idx
      blueprint_roots[[j]] <- term_dependencies[[root_idx]]
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
      blueprint_roots[[j]] <- selected_roots

      # Step C & D: Identify ALL aligned terms (both main effects and higher-order interactions)
      key <- paste(sort(selected_roots),collapse=",")
      if (exists(key,aligned_cache,inherits=FALSE)) {
        aligned_term_indices <- get(key,aligned_cache,inherits=FALSE)
      } else {
        aligned_term_indices <- if (!is.null(dependency_bits)) {
          bits <- as.integer(sum(2^(selected_roots-1L)))
          which(bitwAnd(dependency_bits,bits)==dependency_bits)
        } else which(vapply(term_dependencies,function(deps)
          all(deps %in% selected_roots),logical(1)))
        # Bound cache storage: repeated blueprints benefit without retaining
        # thousands of large menus for high-root configurations.
        if (length(aligned_cache)<512L) assign(key,aligned_term_indices,aligned_cache)
      }

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
        selected_by_column[[j]] <- sort(selected_terms)
      }
    }
  }

  # Randomly shuffle the columns so the pure static and pure root variables are randomly dispersed
  mask_permutation <- sample(m_proxies)
  selected_by_column <- selected_by_column[mask_permutation]
  blueprint_roots <- blueprint_roots[mask_permutation]
  after_mask_rng <- .Random.seed
  assign(".Random.seed",coefficient_rng,envir=.GlobalEnv)
  sparse_values <- vector("list",m_proxies)
  for (j in seq_len(m_proxies)) {
    values <- numeric(n_terms)
    if (n_main>0) values[is_main] <- rnorm(n_main,0,term_sds[is_main])
    if (n_int>0) values[is_int] <- rnorm(n_int,0,term_sds[is_int]*int_scaling)
    ids <- selected_by_column[[j]]
    nonzero <- values[ids]!=0
    selected_by_column[[j]] <- ids[nonzero]
    sparse_values[[j]] <- values[ids[nonzero]]
  }
  assign(".Random.seed",after_mask_rng,envir=.GlobalEnv)
  coeffs <- Matrix::sparseMatrix(i=unlist(selected_by_column,use.names=FALSE),
    p=c(0L,cumsum(lengths(selected_by_column))),
    x=unlist(sparse_values,use.names=FALSE),dims=c(n_terms,m_proxies))
  rm(values,aligned_cache); gc(FALSE)

  # Effective rank of the active expanded terms.
  active_indices <- sort(unique(unlist(selected_by_column,use.names=FALSE)))
  active_terms <- length(active_indices)

  if (active_terms > 0) {
    # Apply scaling to active columns to measure theoretical continuous volume
    active_orders <- interaction_orders[active_indices]
    active_weights <- ifelse(active_orders==1,term_sds[active_indices],
                             term_sds[active_indices]*int_scaling)
    centered_active <- design_mat[,active_indices,drop=FALSE]
    for (begin in seq.int(1L,length(active_indices),by=256L)) {
      ids <- begin:min(length(active_indices),begin+255L)
      block <- t(t(centered_active[,ids,drop=FALSE])*active_weights[ids])
      centered_active[,ids] <- scale(block,center=TRUE,scale=FALSE)
    }
    rm(block); gc(FALSE)
    # Complete covariance spectrum via the smaller exact dual Gram problem.
    if (ncol(centered_active) > nrow(centered_active)) {
      expanded_cov <- tcrossprod(centered_active) / (nrow(centered_active) - 1)
      expanded_eigen <- c(pmax(eigen(expanded_cov, symmetric = TRUE,
                                     only.values = TRUE)$values, 0),
                          rep(0, ncol(centered_active) - nrow(centered_active)))
    } else {
      expanded_cov <- crossprod(centered_active) / (nrow(centered_active) - 1)
      expanded_eigen <- pmax(eigen(expanded_cov, symmetric = TRUE,
                                  only.values = TRUE)$values, 0)
    }
    p_expanded <- expanded_eigen[expanded_eigen > 1e-9] / sum(expanded_eigen[expanded_eigen > 1e-9])
    K_eff <- exp(-sum(p_expanded * log(p_expanded)))
    rm(centered_active,expanded_cov); gc(FALSE)
  } else {
    expanded_eigen <- numeric()
    K_eff <- 0
  }

  gc(verbose = FALSE)
  # Exact sparse multiplication, in bounded output blocks. Every selected term
  # is included; there is no projection, sampling or rank truncation.
  deterministic_block <- NULL

  # 3.5 SET THRESHOLDS & INJECT STRUCTURAL UNCERTAINTY (ROUND 1 NOISE)
  cat(sprintf("      -> Setting Base Rates and Injecting Structural Uncertainty (Round 1 Noise, SNR = %.2f)...\n", structural_snr))
  true_manifestation <- matrix(0, nrow = n, ncol = m_proxies)
  core_center <- core_scale <- thresholds <- numeric(m_proxies)

  # Pre-calculate standard deviations for exact noise scaling
  # Signal Var = 1.0, Noise Var = 1/SNR, Total Var = 1.0 + (1/SNR)
  total_sd <- sqrt(1 + (1 / structural_snr))

  for (j in 1:m_proxies) {
    # Step A: Standardize the deterministic core (Forces pure signal variance to exactly 1.0)
    if ((j-1L) %% 128L == 0L) {
      block_ids <- j:min(m_proxies,j+127L)
      deterministic_block <- as.matrix(design_mat %*% coeffs[,block_ids,drop=FALSE])
    }
    core_vec <- deterministic_block[,(j-1L) %% 128L+1L]
    core_center[j] <- mean(core_vec)
    core_scale[j] <- sd(core_vec)
    if (core_scale[j] > 1e-9) {
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
      thresholds[j] <- shift
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

      # Map threshold-shifted, structurally perturbed scores against
      # N(0, total_sd) to obtain probabilities without rescaling the scores.
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
    K_eff = K_eff,
    expanded_eigenvalues = expanded_eigen,
    term_names = term_names,
    term_dependencies = term_dependencies,
    term_center = attr(design_mat, "scaled:center"),
    term_scale = attr(design_mat, "scaled:scale"),
    # Save nonzero coefficients losslessly without writing a mostly-zero matrix.
    coefficients = lapply(mix_idx,function(j)
      list(term_index=selected_by_column[[j]],value=sparse_values[[j]])),
    blueprint_roots = blueprint_roots[mix_idx],
    mask_permutation = mask_permutation,
    measurement_permutation = mix_idx,
    core_center = core_center[mix_idx],
    core_scale = core_scale[mix_idx],
    binary_threshold_shifts = thresholds[mix_idx]
  ))
}

# Apply independent measurement error to the generated values.
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



)----"
RUNNER_CODE <- r"----(
package_ready <- function(package) {
  if (!requireNamespace(package,quietly=TRUE)) return(FALSE)
  if (package=="Entropic.Scree")
    return(!inherits(tryCatch(cran_description(),error=function(e) e),"error"))
  desc <- tryCatch(suppressWarnings(utils::packageDescription(package)),error=function(e) NULL)
  is.list(desc) && length(desc$Version)==1L && !is.na(desc$Version)
}

release_package <- function(package) {
  path <- if (isNamespaceLoaded(package)) getNamespaceInfo(getNamespace(package),"path") else NULL
  attached <- paste0("package:",package)
  if (attached %in% search())
    tryCatch(detach(attached,character.only=TRUE,unload=FALSE),error=function(e)
      stop("Restart R, then rerun the script to repair ",package," automatically.",call.=FALSE))
  if (isNamespaceLoaded(package))
    tryCatch(unloadNamespace(package),error=function(e)
      stop("Restart R, then rerun the script: ",package," is in use and cannot be repaired in this session.",call.=FALSE))
  if (!is.null(path)) {
    prefix <- paste0(tolower(normalizePath(path,winslash="/",mustWork=FALSE)),"/")
    for (dll in getLoadedDLLs()) {
      dll_path <- dll[["path"]]
      if (startsWith(tolower(normalizePath(dll_path,winslash="/",mustWork=FALSE)),prefix))
        tryCatch(library.dynam.unload(tools::file_path_sans_ext(basename(dll_path)),path),error=function(e)
          stop("Restart R, then rerun: the native library for ",package," is still in use.",call.=FALSE))
    }
  }
}

install_cran_version <- function(version,lib,repo) {
  package <- "Entropic.Scree"
  release_package(package)
  binary_type <- .Platform$pkgType
  available <- tryCatch(suppressWarnings(utils::available.packages(repos=repo,type=binary_type)),
                        error=function(e) matrix(character(),0,0))
  if (package %in% rownames(available) && identical(unname(available[package,"Version"]),version)) {
    utils::install.packages(package,lib=lib,repos=repo,type=binary_type,Ncpus=1L)
  } else {
    # Exact archived releases may require the platform's R build tools.
    if (!package_ready("remotes")) {
      release_package("remotes")
      utils::install.packages("remotes",lib=lib,repos=repo,Ncpus=1L)
    }
    if (!package_ready("remotes")) stop("Could not install remotes for the exact package release. Check Internet access and restart R.",call.=FALSE)
    remotes::install_version(package,version=version,lib=lib,repos=repo,
      dependencies=NA,upgrade="never",force=TRUE)
  }
  assert_cran_engine(version)
}

install_required_packages <- function(packages,lib,version) {
  repo <- "https://cloud.r-project.org"
  oldtimeout <- getOption("timeout");options(timeout=max(600,oldtimeout))
  on.exit(options(timeout=oldtimeout),add=TRUE)
  .libPaths(c(lib,.libPaths()))
  others <- setdiff(packages,"Entropic.Scree")
  missing <- others[!vapply(others,package_ready,logical(1))]
  if (length(missing)) {
    cat("Installing or repairing required R packages:",paste(missing,collapse=", "),"\n")
    for (package in missing) release_package(package)
    utils::install.packages(missing,lib=lib,repos=repo,Ncpus=1L)
  }
  target_ok <- tryCatch({assert_cran_engine(version);TRUE},error=function(e) FALSE)
  if (!target_ok) {
    cat("Installing or repairing Entropic.Scree",version,"in the private library...\n")
    tryCatch(install_cran_version(version,lib,repo),error=function(e)
      stop("Automatic installation of Entropic.Scree ",version," did not complete. ",
           conditionMessage(e)," Check Internet access; source installation may require Rtools on Windows or compiler tools on other systems. Saved results are unchanged.",call.=FALSE))
  }
  missing <- packages[!vapply(packages,package_ready,logical(1))]
  if (length(missing)) stop("Required packages still unavailable: ",paste(missing,collapse=", "),
    ". Check the installation messages above, Internet access and R version, then restart R and rerun.",call.=FALSE)
  assert_cran_engine(version)
}


# Entropic Scree rank-recovery simulation and resource-aware execution.
VERSION <- "es-cran-replicated-2.0.0"
GB <- 1024^3
CONDITION_ORDER <- c("orthogonal", "low", "high")

atomic_rds <- function(x, path, compress = FALSE) {
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(x, tmp, compress = compress)
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) stop("Cannot finalize checkpoint: ", path)
}

fingerprint <- function(x) {
  f <- tempfile(); on.exit(unlink(f)); saveRDS(x, f, compress = FALSE)
  unname(tools::md5sum(f))
}

set_threads <- function(n) {
  n <- max(1L, as.integer(n))
  Sys.setenv(OMP_NUM_THREADS = n, OPENBLAS_NUM_THREADS = n,
             MKL_NUM_THREADS = n, VECLIB_MAXIMUM_THREADS = n,
             BLIS_NUM_THREADS = n, RCPP_PARALLEL_NUM_THREADS = n)
  data.table::setDTthreads(n)
  RhpcBLASctl::blas_set_num_threads(n)
  RhpcBLASctl::omp_set_num_threads(n)
}

memory_estimate <- function(k, cfg) {
  n <- as.double(cfg$n_rows); m <- as.double(cfg$m_proxies)
  q <- sum(choose(k,seq_len(min(k,cfg$max_interaction_order))))+k*(cfg$max_polynomial_order-1)
  b <- ceiling(k/2)
  local_terms <- sum(choose(b,seq_len(min(b,cfg$max_interaction_order))))+b*(cfg$max_polynomial_order-1)
  nnz <- m*max(1,round(local_terms*cfg$term_fraction))
  # Account for the optimized generator's live dense designs/Gram workspace,
  # sparse coefficients and generated data.
  generation <- 8*(4*n*q+2*min(n,q)^2+4*n*m)+32*nnz
  # CRAN holds several full m x m arrays and the full eigensolver's workspace.
  # Reserve their peak phase, not the sum of sequential generation+Scree phases.
  scree <- 8*(7*m*m+4*n*m)
  list(bytes=(max(generation,scree)+0.5*GB)*cfg$memory_safety_factor,
       term_pool=q,generation=generation,scree=scree)
}

# Scheduling diagnostics do not change generated data or CRAN calculations.
MEMORY_POLICY <- "whole-job-rss-peak-by-root-v1"

memory_profile <- function(samples) {
  samples <- samples[is.finite(samples) & samples > 0]
  list(samples=length(samples), median=if(length(samples)) stats::median(samples) else NA_real_,
    peak=if(length(samples)) max(samples) else NA_real_)
}

memory_host <- function() {
  list(host=unname(Sys.info()[["nodename"]]),R=R.version.string,
       platform=R.version$platform,BLAS=unname(extSoftVersion()["BLAS"]),
       LAPACK=as.character(La_version()))
}

saved_memory_profiles <- function(jobs) {
  host <- memory_host()
  rows <- lapply(jobs,function(job) {
    path <- file.path(job$dir,"memory-samples.rds")
    obj <- if(file.exists(path)) tryCatch(readRDS(path),error=function(e) NULL) else NULL
    if (is.null(obj) || !job_complete(job) || !isTRUE(obj$full_lifecycle) || !isTRUE(obj$success) ||
        length(obj$policy)!=1L || !obj$policy %in% c(MEMORY_POLICY,"whole-job-rss-p75-v1") || !identical(obj$host,host) ||
        !identical(obj$signature,job$signature) || !is.numeric(obj$rss_bytes)) return(NULL)
    profile <- memory_profile(obj$rss_bytes)
    if (profile$samples<4L || !is.finite(profile$peak)) return(NULL)
    data.frame(configuration=job$id,k=job$k,peak_bytes=profile$peak)
  })
  rows <- Filter(Negate(is.null),rows)
  if (!length(rows)) return(data.frame(configuration=character(),k=integer(),median_bytes=numeric()))
  do.call(rbind,rows)
}

memory_reservation <- function(job, profiles) {
  if (!nrow(profiles)) return(list(bytes=job$memory$bytes,source_k=NA_integer_,profiles=0L))
  # Prefer the same root count. Until it is measured, use the closest measured
  # higher count, or the largest measured count when no higher one is available.
  ks <- sort(unique(profiles$k))
  source_k <- if(any(ks>=job$k)) min(ks[ks>=job$k]) else max(ks)
  values <- profiles$peak_bytes[profiles$k==source_k]
  list(bytes=0.80*max(values),source_k=source_k,profiles=length(values))
}

# Only scheduling functions may differ when reusing an existing experiment.
# All worker calculations, summaries, package identity and scientific settings
# must match; checks compare parsed R code without executing saved source.
scientific_runner <- function(code) {
  # Normalize the condition-list display text when checking resumable results.
  code <- gsub("Groups: orthogonal, low, medium, high;", "Groups: orthogonal, low, high;",code,fixed=TRUE)
  code <- gsub('utils::packageVersion("Entropic.Scree")','cran_description()$Version',code,fixed=TRUE)
  code <- gsub('  job$threads <- refresh_worker_threads(job$threads)\n','',code,fixed=TRUE)
  code <- gsub('  if (!is.null(getOption("simulation.thread_control"))) refresh_worker_threads(1L)\n','',code,fixed=TRUE)
  scheduling <- c("run_sweep","calibration_multiplier","memory_profile",
                  "memory_host","saved_memory_ratios","saved_memory_profiles","memory_reservation","scientific_runner",
                  "resume_batch","scientific_generator","MEMORY_POLICY","CONDITION_ORDER",
                  "cran_description","assert_cran_engine","package_ready",
                  "release_package","install_cran_version","install_required_packages",
                  "launch_threads","worker_thread_limit","worker_count_limit","refresh_worker_threads")
  nodes <- as.list(parse(text=code,keep.source=FALSE))
  nodes[!vapply(nodes,function(node) is.call(node) &&
    identical(node[[1]],as.name("<-")) && is.symbol(node[[2]]) &&
    as.character(node[[2]]) %in% scheduling,logical(1))]
}

scientific_generator <- function(code) {
  nodes <- as.list(parse(text=code,keep.source=FALSE))
  # This standalone helper is not referenced by the generating pipeline.
  nodes[!vapply(nodes,function(node) is.call(node) &&
    identical(node[[1]],as.name("<-")) && is.symbol(node[[2]]) &&
    identical(as.character(node[[2]]),"generate_random_corr_matrix"),logical(1))]
}

resume_batch <- function(root, signature, scientific, engine_code, runner_code, package_identity) {
  expected <- file.path(root,paste0("experiment-",substr(signature,1,12)))
  candidates <- list.dirs(root,recursive=FALSE,full.names=TRUE)
  candidates <- candidates[grepl("^experiment-[a-f0-9]{12}$",basename(candidates))]
  compatible <- lapply(candidates,function(path) tryCatch({
    saved <- readRDS(file.path(path,"SETTINGS.rds"))
    provenance <- readRDS(file.path(path,"CRAN-AND-MATH-PROVENANCE.rds"))
    previous_engine <- paste(readLines(file.path(path,"engine.R"),warn=FALSE),collapse="\n")
    previous_runner <- paste(readLines(file.path(path,"runner.R"),warn=FALSE),collapse="\n")
    if (!identical(saved$scientific,scientific) ||
        !identical(saved$version,VERSION) || !identical(provenance$package,package_identity) ||
        !identical(scientific_generator(previous_engine),scientific_generator(engine_code)) ||
        !identical(scientific_runner(previous_runner),scientific_runner(runner_code))) return(NULL)
    if (!is.character(saved$signature) || length(saved$signature)!=1L ||
        !grepl("^[a-f0-9]{32}$",saved$signature) ||
        !identical(basename(path),paste0("experiment-",substr(saved$signature,1,12)))) return(NULL)
    list(path=path,signature=saved$signature)
  },error=function(e) NULL))
  compatible <- Filter(Negate(is.null),compatible)
  exact <- Filter(function(x) identical(x$path,expected),compatible)
  if (length(exact)) return(exact[[1]])
  if (length(compatible)>1L)
    stop("Multiple scientifically equivalent experiment folders found. Keep the intended experiment in the output folder before resuming; nothing was changed.")
  if (length(compatible)) {
    cat("Resuming scientifically identical results with updated scheduling:\n",compatible[[1]]$path,"\n")
    return(compatible[[1]])
  }
  if (dir.exists(expected)) stop("Existing experiment provenance could not be verified; refusing to overwrite it.")
  list(path=expected,signature=signature)
}


cran_description <- function() {
  package <- "Entropic.Scree"
  if (!requireNamespace(package,quietly=TRUE))
    stop("Entropic.Scree cannot be loaded. Install the required CRAN version, restart R and rerun.",call.=FALSE)
  ns <- getNamespace(package)
  path <- getNamespaceInfo(ns,"path")
  field <- function(x,key) {
    if (!is.list(x)) return(NULL)
    value <- x[[key]]
    if (length(value)!=1L || is.na(value) || !nzchar(as.character(value))) return(NULL)
    as.character(value)
  }
  desc <- tryCatch(suppressWarnings(utils::packageDescription(package,lib.loc=dirname(path))),
                   error=function(e) NULL)
  if (is.null(field(desc,"Version")) || is.null(field(desc,"Package"))) {
    metadata <- tryCatch(suppressWarnings(readRDS(file.path(path,"Meta","package.rds"))),error=function(e) NULL)
    desc <- if (is.list(metadata) && !is.null(metadata$DESCRIPTION))
      as.list(metadata$DESCRIPTION) else NULL
  }
  if (is.null(field(desc,"Version")) || is.null(field(desc,"Package"))) {
    dcf <- tryCatch(suppressWarnings(read.dcf(file.path(path,"DESCRIPTION"))),error=function(e) NULL)
    if (is.matrix(dcf) && nrow(dcf)==1L) desc <- as.list(dcf[1,])
  }
  if (!identical(field(desc,"Package"),package) || is.null(field(desc,"Version")))
    stop("Entropic.Scree loaded, but its installed package metadata cannot be read at: ",path,
         ". Restart R; if this persists, reinstall the required CRAN package version. Saved results are unchanged.",call.=FALSE)
  loaded_version <- as.character(getNamespaceInfo(ns,"spec")[["version"]])
  if (!identical(field(desc,"Version"),loaded_version))
    stop("Entropic.Scree's loaded version differs from its installed metadata. Restart R before rerunning; do not update packages during a run.",call.=FALSE)
  desc
}

assert_cran_engine <- function(version) {
  desc <- cran_description()
  if (!identical(as.character(desc$Version),version) ||
      !identical(as.character(desc$Repository),"CRAN"))
    stop("This study requires published CRAN Entropic.Scree ",version,
         ". Loaded version: ",as.character(desc$Version),
         "; repository: ",if(is.null(desc$Repository)) "not recorded" else as.character(desc$Repository),
         ". Restart R after installing the required version.",call.=FALSE)
  invisible(TRUE)
}


math_backend <- function() {
  list(BLAS=unname(extSoftVersion()["BLAS"]),LAPACK=La_library(),
       LAPACK_version=as.character(La_version()),
       BLAS_threads=RhpcBLASctl::blas_get_num_procs(),R=R.version.string)
}

temporary_data_names <- function() {
  c("observed_data.rds","true_data.rds","centered_NMI.rds","generation.rds","DATA-COMPLETE.rds")
}

retained_result_names <- function() {
  c("Entropic_Scree.rds","seed-result.rds","seed-result.csv","spectrum.csv","FSIG.csv",
    "generation-summary.rds","scree.png","scree.pdf","sessionInfo.txt")
}

job_disk_estimate <- function(job) {
  cfg <- job$cfg; n <- as.double(cfg$n_rows); m <- as.double(cfg$m_proxies)
  b <- ceiling(job$k/2)
  q <- sum(choose(job$k,seq_len(min(job$k,cfg$max_interaction_order)))) +
    job$k*(cfg$max_polynomial_order-1)
  local_terms <- sum(choose(b,seq_len(min(b,cfg$max_interaction_order)))) + b*(cfg$max_polynomial_order-1)
  sparse_terms <- m*max(1,round(local_terms*cfg$term_fraction))
  # Small permanent results: spectra, bin counts, names, diagnostics, CSVs and plots.
  # These are conservative planning estimates, not guaranteed bounds.
  permanent <- 1.3*(0.02*GB + m*(1024+16*max(2,ceiling(n^(1/3)))) + 16*job$k^2)
  temporary <- 1.3*(8*n*m*(1+as.integer(cfg$save_true_data)) +
    8*m*m*as.integer(cfg$save_similarity_matrix) + 16*sparse_terms +
    8*n*job$k + 256*q + 0.02*GB)
  size <- function(files) sum(file.info(file.path(job$dir,files))$size,na.rm=TRUE)
  list(permanent=max(0,permanent-size(retained_result_names())),
       temporary=max(0,temporary-size(temporary_data_names())))
}

disk_estimate <- function(jobs, max_concurrent=1L) {
  if (!length(jobs)) return(0)
  plans <- lapply(jobs,job_disk_estimate)
  retained <- sum(vapply(plans,`[[`,numeric(1),"permanent"))
  active <- sort(vapply(plans,`[[`,numeric(1),"temporary"),decreasing=TRUE)
  retained + sum(head(active,max(1L,as.integer(max_concurrent)))) + 2*GB
}

disk_available <- function(path) {
  tryCatch({
    info <- ps::ps_disk_usage(path)
    if (!"available" %in% names(info)) return(Inf)
    value <- min(info$available)
    if (is.finite(value)) value else Inf
  },error=function(e) Inf)
}

launch_disk_required <- function(outstanding, active) {
  # Space already occupied by interrupted checkpoints is excluded: the live
  # free-space reading already accounts for it. Failed jobs keep their data.
  sum(vapply(outstanding,function(j) job_disk_estimate(j)$permanent,numeric(1))) +
    sum(vapply(active,function(j) job_disk_estimate(j)$temporary,numeric(1))) + 2*GB
}

system_memory <- function() {
  m <- ps::ps_system_memory()
  avail <- if ("avail" %in% names(m)) m[["avail"]] else m[["available"]]
  if (length(avail)!=1L || !is.finite(avail) || avail<=0 || !is.finite(m[["total"]]))
    stop("Cannot read available physical RAM; no unsafe concurrency guess was made.")
  list(total=unname(m[["total"]]),available=unname(avail))
}

worker_thread_limit <- function(logical_cpus) {
  max(1L,as.integer(logical_cpus)-2L)
}

worker_count_limit <- function(usable_threads, requested=NULL) {
  cap <- max(1L,as.integer(floor(usable_threads)))
  if (!is.null(requested)) cap <- min(cap,max(1L,as.integer(requested)))
  cap
}

status_write <- function(dir, stage) {
  writeLines(c(stage,format(Sys.time(),tz="UTC",usetz=TRUE)),file.path(dir,"CURRENT-STAGE.txt"))
  cat(format(Sys.time(),"%H:%M:%S"), stage,"\n"); flush.console()
}

job_complete <- function(job) {
  f <- file.path(job$dir,"COMPLETE.rds")
  if (!file.exists(f)) return(FALSE)
  mark <- tryCatch(readRDS(f),error=function(e) NULL)
  if (is.null(mark) || !identical(mark$signature,job$signature) || !length(mark$files) ||
      !identical(mark$storage_policy,"temporary-data-v1") ||
      !setequal(names(mark$files),retained_result_names())) return(FALSE)
  info <- file.info(file.path(job$dir,names(mark$files)))$size
  all(is.finite(info)) && identical(as.numeric(info),as.numeric(mark$files))
}

compact_generation <- function(meta,job) {
  result <- meta[c("configuration","settings","latent_correlation","root_eigenvalues",
    "root_effective_rank","mean_absolute_root_correlation","initial_rng","final_rng",
    "n_rows","n_measurements","blueprint_size")]
  result$generation <- meta$generation[c("active_terms","K_eff")]
  result$dataset_signature <- job$signature
  result$storage_policy <- "temporary-data-v1"
  result
}

generation_summary <- function(job) {
  small <- file.path(job$dir,"generation-summary.rds")
  if (file.exists(small)) {
    meta <- readRDS(small)
    if (!identical(meta$dataset_signature,job$signature)) stop("Generation summary signature mismatch")
    return(meta)
  }
  large <- file.path(job$dir,"generation.rds")
  if (!file.exists(large)) stop("Generation metadata missing; cannot finalize this result")
  meta <- compact_generation(readRDS(large),job)
  atomic_rds(meta,small)
  meta
}

cleanup_temporary_data <- function(job) {
  if (!job_complete(job)) stop("Temporary data retained: saved result files are not complete")
  dir <- normalizePath(job$dir,winslash="/",mustWork=TRUE)
  candidates <- list.files(dir,all.files=TRUE,no..=TRUE,full.names=FALSE)
  allowed <- temporary_data_names()
  # Only these exact files and their interrupted atomic-save siblings qualify.
  # Never traverse directories, glob-delete unrelated files, or follow an input path.
  candidates <- candidates[vapply(candidates,function(n) n %in% allowed ||
    any(vapply(allowed,function(a) startsWith(n,paste0(a,".tmp-")) &&
      grepl("^[0-9]+$",substring(n,nchar(a)+6L)),logical(1))),logical(1))]
  if (!length(candidates)) return(invisible(TRUE))
  marker <- readRDS(file.path(dir,"COMPLETE.rds"))
  core <- c("Entropic_Scree.rds","seed-result.rds","generation-summary.rds")
  if (!identical(names(marker$checksums),core) ||
      !identical(unname(tools::md5sum(file.path(dir,core))),unname(marker$checksums)))
    stop("Temporary data retained: saved result verification failed")
  obj <- readRDS(file.path(dir,core[1])); row <- readRDS(file.path(dir,core[2]))
  meta <- readRDS(file.path(dir,core[3]))
  if (!identical(obj$dataset_signature,job$signature) ||
      !identical(meta$dataset_signature,job$signature) ||
      nrow(row)!=1L || !identical(row$configuration,job$id) ||
      !isTRUE(all.equal(row$extracted_rank,obj$K_roots)))
    stop("Temporary data retained: numerical result or generation identity mismatch")
  paths <- file.path(dir,candidates)
  resolved <- normalizePath(paths,winslash="/",mustWork=TRUE)
  if (any(dirname(resolved)!=dir) || any(dir.exists(paths)))
    stop("Temporary cleanup target is outside this job or is a directory; nothing deleted")
  sizes <- file.info(paths)$size
  for (path in paths) {
    unlink(path,recursive=FALSE)
    if (file.exists(path)) stop("Results are safe, but cannot delete temporary file: ",path,
      ". Close programs using it and rerun to retry cleanup.")
  }
  writeLines(c(format(Sys.time(),tz="UTC",usetz=TRUE),
    "Verified numerical results, generation summary, tables and plots before cleanup.",
    paste(candidates,format(sizes,scientific=FALSE),"bytes removed")),
    file.path(dir,"TEMPORARY-DATA-REMOVED.txt"))
  invisible(TRUE)
}

numeric_result <- function(observed, job) {
  start <- proc.time()[["elapsed"]]
  pdf(file.path(job$dir,"entropic-native-plots.pdf"),width=12,height=6)
  on.exit(dev.off(),add=TRUE)
  assert_cran_engine(job$cfg$cran_version)
  # Direct exported CRAN call. All scientific defaults, including preparation
  # and both automatic detectors, remain package defaults. Only automation and
  # the OpenMP thread budget are specified. Full spectrum; eigenvectors not needed.
  # The package's optional collinearity check can sample rows when N > m.
  # Set its RNG explicitly so resumed data checkpoints give the same analysis.
  RNGkind("Mersenne-Twister","Inversion","Rejection"); set.seed(job$seed)
  obj <- Entropic.Scree::Entropic.Scree(observed,
    interactive_mode=FALSE,num_cores=job$threads)
  obj$package_version <- as.character(cran_description()$Version)
  obj$math_backend <- math_backend()
  obj$call_arguments <- list(interactive_mode=FALSE,num_cores=job$threads)
  # This is the package's default bin calculation; no second MI/eigen solve.
  obj$num_bins_used <- max(2,ceiling(nrow(observed)^(1/3)))
  if (job$cfg$save_similarity_matrix)
    atomic_rds(obj$similarity_matrix,file.path(job$dir,"centered_NMI.rds"))
  obj$similarity_matrix <- NULL
  obj$selected_rank <- obj$K_roots
  obj$selection_method <- paste0("CRAN_Entropic.Scree_",obj$package_version,"_noninteractive")
  obj$seconds <- proc.time()[["elapsed"]]-start
  obj$dataset_signature <- job$signature
  obj
}

result_row <- function(obj, meta, job) {
  get <- function(key) if(is.null(obj[[key]])) NA_real_ else obj[[key]]
  rank <- obj$K_roots
  stopifnot(length(rank)==1L,is.finite(rank),rank>=1,rank==as.integer(rank))
  data.frame(configuration=job$id,true_rank=job$k,blueprint_roots=ceiling(job$k/2),
    entanglement=job$condition,alphad=job$alpha,repetition=job$repetition,seed=job$seed,
    N=job$cfg$n_rows,m_generated=job$cfg$m_proxies,m_retained=length(obj$retained_features),
    root_effective_rank=meta$root_effective_rank,
    root_mean_abs_correlation=meta$mean_absolute_root_correlation,
    active_terms=meta$generation$active_terms,expanded_effective_rank=meta$generation$K_eff,
    status="complete",extracted_rank=rank,rank_error=rank-job$k,
    absolute_error=abs(rank-job$k),exact_match=rank==job$k,within_one=abs(rank-job$k)<=1,
    selection_method=obj$selection_method,package_version=obj$package_version,
    K_extended=get("K_extended"),K_log_gap=get("K_log_gap"),K_triple_tap=get("K_triple_tap"),
    top_of_bulk=get("top_of_bulk"),bins=get("num_bins_used"),R_eff=get("R_eff"),
    AIG=get("AIG"),SCDR=get("SCDR"),mean_trace=get("mean_trace"),
    total_signal_volume=get("total_signal_volume"),unique_signal_volume=get("unique_signal_volume"),
    redundant_signal_volume=get("redundant_signal_volume"),
    idiosyncratic_volume=get("idiosyncratic_volume"),seconds=obj$seconds,error="")
}

save_plots <- function(obj, job) {
  draw <- function() {
    par(mfrow=c(1,2),mar=c(4.3,4.2,3.6,1),oma=c(0,0,3.8,0))
    vals <- obj$eigenvalues
    for (detail in c(FALSE,TRUE)) {
      last <- if(detail) min(length(vals),max(60L,ceiling(1.25*max(job$k,obj$K_roots)))) else length(vals)
      ids <- seq_len(last)
      plot(ids,vals[ids],type="l",log="x",col="#286e85",lwd=1.3,
        xlab="Eigenvalue index (log scale)",ylab="Working eigenvalue",
        main=if(detail) "Primary-boundary detail" else "Full spectrum")
      abline(v=job$k,col="#c76b28",lty=2,lwd=2)
      abline(v=obj$K_roots,col="#724496",lty=3,lwd=2)
      legend("topright",c(paste("Generating rank",job$k),paste("Automatic rank",obj$K_roots)),
        col=c("#c76b28","#724496"),lty=c(2,3),bty="n",cex=.78)
    }
    mtext(sprintf("Entropic Scree | %s | r = %d | repetition %d | seed %d",
      job$condition,job$k,job$repetition,job$seed),outer=TRUE,line=1.3,font=2)
    mtext(sprintf("N = %s; generated m = %s; retained m = %s | log x, linear y",
      format(job$cfg$n_rows,big.mark=","),format(job$cfg$m_proxies,big.mark=","),
      format(length(obj$retained_features),big.mark=",")),outer=TRUE,line=.1,cex=.85)
  }
  png(file.path(job$dir,"scree.png"),width=2000,height=1000,res=160)
  tryCatch(draw(),finally=dev.off())
  pdf(file.path(job$dir,"scree.pdf"),width=12.5,height=6.25)
  tryCatch(draw(),finally=dev.off())
}

run_job <- function(job) {
  dir.create(job$dir,recursive=TRUE,showWarnings=FALSE)
  set_threads(job$threads)
  gc(reset=TRUE)
  tryCatch({
    unlink(file.path(job$dir,"JOB-ERROR.txt"))
    if (job_complete(job)) {
      cleanup_temporary_data(job)
      status_write(job$dir,"Complete; temporary data removed")
      return(TRUE)
    }
    path <- file.path(job$dir,"Entropic_Scree.rds")
    if (file.exists(path)) {
      status_write(job$dir,"Reusing numerical checkpoint; completing outputs")
      obj <- readRDS(path)
      if (!identical(obj$dataset_signature,job$signature)) stop("Numerical checkpoint signature mismatch")
    } else {
      status_write(job$dir,"Loading published CRAN Entropic.Scree")
      assert_cran_engine(job$cfg$cran_version)
      marker <- file.path(job$dir,"DATA-COMPLETE.rds")
      required <- c("observed_data.rds","generation.rds",if(job$cfg$save_true_data) "true_data.rds")
      if (file.exists(marker)) {
        if (!identical(readRDS(marker)$signature,job$signature)) stop("Data checkpoint signature mismatch")
        if (!all(file.exists(file.path(job$dir,required)))) stop("Completed data checkpoint has missing files")
        status_write(job$dir,"Reusing saved dataset")
      } else {
        status_write(job$dir,"Generating dataset; saving temporary restart checkpoint")
        meta <- generate_dataset(job)
        rm(meta); gc(FALSE)
      }
      observed <- readRDS(file.path(job$dir,"observed_data.rds"))
      status_write(job$dir,"Running Entropic Scree (unchanged automatic detectors)")
      obj <- numeric_result(observed,job)
      atomic_rds(obj,path)
      rm(observed); gc(FALSE)
    }
    meta <- generation_summary(job)
    table <- result_row(obj,meta,job)
    atomic_rds(table,file.path(job$dir,"seed-result.rds"))
    data.table::fwrite(table,file.path(job$dir,"seed-result.csv"))
    data.table::fwrite(data.frame(index=seq_along(obj$eigenvalues),
      eigenvalue=obj$eigenvalues),file.path(job$dir,"spectrum.csv"))
    data.table::fwrite(data.frame(index=seq_along(obj$FSIG_final),FSIG=obj$FSIG_final),
      file.path(job$dir,"FSIG.csv"))
    status_write(job$dir,"Saving full/detail scree plots")
    save_plots(obj,job)
    writeLines(capture.output(sessionInfo()),file.path(job$dir,"sessionInfo.txt"))
    allocation <- gc()
    atomic_rds(list(max_R_heap_bytes=sum(allocation[,6])*1024^2,
      backend=math_backend()),file.path(job$dir,"memory-diagnostics.rds"))
    required <- retained_result_names()
    sizes <- setNames(file.info(file.path(job$dir,required))$size,required)
    if (any(!is.finite(sizes)) || any(sizes<=0)) stop("Output files missing; completion not recorded")
    core <- c("Entropic_Scree.rds","seed-result.rds","generation-summary.rds")
    checksums <- setNames(unname(tools::md5sum(file.path(job$dir,core))),core)
    if (anyNA(checksums)) stop("Result verification failed; temporary data retained")
    # The completion marker lists permanent outputs only. If interrupted during
    # cleanup, restart checks those results and retries deletion without a refit.
    atomic_rds(list(signature=job$signature,files=sizes,checksums=checksums,
      storage_policy="temporary-data-v1",finished=Sys.time()),file.path(job$dir,"COMPLETE.rds"))
    cleanup_temporary_data(job)
    status_write(job$dir,"Complete; temporary data removed")
    TRUE
  },error=function(e) {
    writeLines(conditionMessage(e),file.path(job$dir,"JOB-ERROR.txt"))
    status_write(job$dir,paste("FAILED:",conditionMessage(e)))
    FALSE
  })
}

summarize_recovery <- function(all) {
  all <- as.data.frame(all)
  rows <- list()
  for (condition in CONDITION_ORDER) for (k in sort(unique(all$true_rank))) {
    z <- all[all$entanglement==condition & all$true_rank==k,,drop=FALSE]
    if (!nrow(z)) next
    ok <- z$status=="complete" & is.finite(z$extracted_rank)
    ranks <- z$extracted_rank[ok]; errors <- ranks-k; n <- length(ranks)
    stat <- function(f) if(n) f(ranks) else NA_real_
    rows[[length(rows)+1L]] <- data.frame(entanglement=condition,true_rank=k,
      seeds_planned=nrow(z),seeds_valid=n,seeds_failed=sum(z$status=="failed"),
      seeds_pending=sum(z$status=="pending"),
      all_seeds_valid=n==nrow(z),
      mean_extracted_rank=stat(mean),median_extracted_rank=stat(median),
      min_extracted_rank=stat(min),max_extracted_rank=stat(max),
      sd_extracted_rank=if(n>1) sd(ranks) else NA_real_,
      exact_matches=sum(errors==0),within_one_matches=sum(abs(errors)<=1),
      exact_match_pct=if(n) 100*mean(errors==0) else NA_real_,
      within_one_pct=if(n) 100*mean(abs(errors)<=1) else NA_real_,
      mean_absolute_error=if(n) mean(abs(errors)) else NA_real_,
      mean_signed_error=if(n) mean(errors) else NA_real_)
  }
  do.call(rbind,rows)
}

save_recovery_plot <- function(all,root) {
  good <- all[all$status=="complete" & is.finite(all$extracted_rank),]
  if (!nrow(good)) return(invisible(NULL))
  draw <- function() {
    par(mfrow=c(2,2),mar=c(4.3,4.4,3,1),oma=c(1.4,0,3.5,0))
    xr <- range(all$true_rank)+c(-.5,.5)
    yr <- range(c(all$true_rank,good$extracted_rank))+c(-.5,.5)
    for (condition in CONDITION_ORDER) {
      plot(NA,xlim=xr,ylim=yr,xlab="Generating rank",ylab="Automatic primary rank",
        main=paste0(toupper(substr(condition,1,1)),substring(condition,2)),xaxt="n")
      axis(1,at=sort(unique(all$true_rank)))
      abline(0,1,lty=2,col="#a45d29",lwd=1.5)
      z <- good[good$entanglement==condition,]
      for (k in sort(unique(z$true_rank))) {
        a <- z[z$true_rank==k,]; v <- a$extracted_rank
        segments(k,min(v),k,max(v),col="#92a6ac",lwd=2)
        # Deterministic offsets; no RNG use or perturbation of vertical estimates.
        offset <- (a$repetition-(max(all$repetition)+1)/2)*min(.07,.35/max(all$repetition))
        points(k+offset,v,pch=16,cex=.68,col=adjustcolor("#286e85",.72))
        segments(k-.08,median(v),k+.08,median(v),lwd=2.5,col="#102f38")
      }
    }
    mtext("Entropic Scree: automatic primary-rank recovery",outer=TRUE,line=1.8,font=2,cex=1.15)
    mtext(sprintf("%d / %d datasets have numerical results; N = %s, generated m = %s",
      nrow(good),nrow(all),format(all$N[1],big.mark=","),format(all$m_generated[1],big.mark=",")),
      outer=TRUE,line=.4,cex=.9)
    mtext("Dots = seeds; bars = observed range; dark marks = median; dashed line = exact recovery. Ranges are not confidence intervals.",
      outer=TRUE,side=1,line=.2,cex=.74)
  }
  png(file.path(root,"RANK-RECOVERY.png"),width=1900,height=1550,res=170)
  tryCatch(draw(),finally=dev.off())
  pdf(file.path(root,"RANK-RECOVERY.pdf"),width=11.18,height=9.12)
  tryCatch(draw(),finally=dev.off())
}

write_master <- function(jobs,root,failed=character(),make_plot=FALSE) {
  old_options <- options(width=220)
  on.exit(options(old_options),add=TRUE)
  tables <- lapply(jobs,function(j) {
    f <- file.path(j$dir,"seed-result.rds")
    table <- if(file.exists(f)) tryCatch(readRDS(f),error=function(e) NULL) else NULL
    if (!is.null(table)) {
      table$outputs_complete <- job_complete(j)
      return(table)
    }
    err <- file.path(j$dir,"JOB-ERROR.txt")
    is_failed <- j$id %in% failed || file.exists(err)
    data.frame(configuration=j$id,true_rank=j$k,blueprint_roots=ceiling(j$k/2),
      entanglement=j$condition,alphad=j$alpha,repetition=j$repetition,seed=j$seed,
      N=j$cfg$n_rows,m_generated=j$cfg$m_proxies,status=if(is_failed) "failed" else "pending",
      extracted_rank=NA_real_,error=if(file.exists(err)) paste(readLines(err,warn=FALSE),collapse=" ") else
        if(is_failed) "Worker exited before saving numerical results; see worker.log" else "",
      outputs_complete=FALSE)
  })
  all <- as.data.frame(data.table::rbindlist(tables,fill=TRUE))
  all <- all[order(match(all$entanglement,CONDITION_ORDER),all$true_rank,all$repetition),]
  rownames(all) <- NULL
  data.table::fwrite(all,file.path(root,"SEED-RESULTS.csv"))
  atomic_rds(all,file.path(root,"SEED-RESULTS.rds"))
  summary <- summarize_recovery(all)
  data.table::fwrite(summary,file.path(root,"RANK-SUMMARY.csv"))
  # Companion compact manuscript-oriented table; the numeric CSV keeps full precision.
  fmt <- function(x) ifelse(is.na(x),"",formatC(x,format="f",digits=1))
  compact <- data.frame(Entanglement=summary$entanglement,True_rank=summary$true_rank,
    Seeds=paste0(summary$seeds_valid,"/",summary$seeds_planned),
    Median_rank_and_range=ifelse(summary$seeds_valid>0,
      paste0(fmt(summary$median_extracted_rank)," [",summary$min_extracted_rank,"-",summary$max_extracted_rank,"]"),""),
    Exact_match_percent=fmt(summary$exact_match_pct),Within_one_percent=fmt(summary$within_one_pct),
    Mean_absolute_error=fmt(summary$mean_absolute_error),Failed=summary$seeds_failed,
    Pending=summary$seeds_pending)
  data.table::fwrite(compact,file.path(root,"RANK-SUMMARY-COMPACT.csv"))
  notes <- c("ENTROPIC SCREE: AUTOMATIC PRIMARY-RANK RECOVERY",
    "Groups: orthogonal, low, high; ascending generating rank within group.",
    "Exact = selected K_roots equals generating r. Within one includes exact matches.",
    "Percentages use seeds_valid as denominator; failures and pending seeds are listed separately.",
    "Do not describe an incomplete row as a fully replicated result. Check Seeds and Failed/Pending.",
    sprintf("At %d valid seeds per configuration, percentages change in increments of %g percentage points.",
      unique(summary$seeds_planned),100/unique(summary$seeds_planned)),
    "[min-max] is the empirical spread of extracted ranks, not a confidence interval.",
    "Numerical completeness and output-file completeness are separate in SEED-RESULTS.csv.",
    "No manual overrides or baseline-method comparisons are performed.",
    "RANK-SUMMARY.csv also includes mean, SD, signed error and counts; SEED-RESULTS.csv has every seed.",
    "",capture.output(print(compact,row.names=FALSE,right=FALSE)))
  writeLines(notes,file.path(root,"RANK-SUMMARY.txt"))
  if (make_plot) save_recovery_plot(all,root)
  invisible(summary)
}

generate_dataset <- function(job) {
  cfg <- job$cfg; k <- job$k
  RNGkind("Mersenne-Twister","Inversion","Rejection"); set.seed(job$seed)
  initial_rng <- .Random.seed
  sigma <- if (is.infinite(job$alpha)) diag(k) else
    clusterGeneration::rcorrmatrix(k,alphad=job$alpha)
  root_eigen <- pmax(eigen(sigma,symmetric=TRUE)$values,0)
  prob <- root_eigen[root_eigen>1e-9]/sum(root_eigen[root_eigen>1e-9])
  z <- MASS::mvrnorm(n=cfg$n_rows,mu=rep(0,k),Sigma=sigma)
  truth <- generate_true_mixed_proxies(z,cfg$m_proxies,
    max_interaction_order=cfg$max_interaction_order,
    max_polynomial_order=cfg$max_polynomial_order,int_scaling=1,
    continuous_ratio=cfg$continuous_ratio,
    MIN_SPARSITY_RATIO=cfg$term_fraction,MAX_SPARSITY_RATIO=cfg$term_fraction,
    MIN_ROOT_SPARSITY_RATIO=0.5,MAX_ROOT_SPARSITY_RATIO=0.5,
    structural_snr=cfg$structural_snr,pure_static_ratio=cfg$pure_static_ratio,
    pure_root_ratio=cfg$pure_root_ratio)
  observed <- apply_measurement_error(truth,cfg$continuous_snr,cfg$binary_error_rate)
  # Save generation products before any method runs. Each RDS is finalized atomically.
  atomic_rds(observed,file.path(job$dir,"observed_data.rds"))
  if (cfg$save_true_data) atomic_rds(truth$data_matrix,file.path(job$dir,"true_data.rds"))
  truth$data_matrix <- NULL
  meta <- list(configuration=job[c("k","condition","alpha","seed","repetition")],
    settings=cfg,latent_correlation=sigma,latent_values=z,root_eigenvalues=root_eigen,
    root_effective_rank=exp(-sum(prob*log(prob))),
    mean_absolute_root_correlation=mean(abs(sigma[upper.tri(sigma)])),
    generation=truth,initial_rng=initial_rng,final_rng=.Random.seed,
    n_rows=nrow(observed),n_measurements=ncol(observed),blueprint_size=ceiling(k/2))
  atomic_rds(meta,file.path(job$dir,"generation.rds"))
  # Marker is only written after every generation file succeeds.
  atomic_rds(list(signature=job$signature,seed=job$seed),file.path(job$dir,"DATA-COMPLETE.rds"))
  meta
}


run_sweep <- function(cfg,engine_code,runner_code) {
  if (getRversion()<"4.2.0") stop("Use R 4.2 or later (current R recommended).")
  packages <- c("Entropic.Scree","Matrix","data.table","ggplot2","patchwork",
                "MASS","stringr","clusterGeneration","ps","processx","RhpcBLASctl")
  # A private user-writable library keeps installation independent of administrator rights.
  root <- path.expand(cfg$output_dir)
  if (cfg$smoke_test) {
    cfg$roots <- c(4L,5L); cfg$n_rows <- 100L; cfg$m_proxies <- 50L
    cfg$repetitions <- 2L
    cfg$max_interaction_order <- 3L; cfg$max_polynomial_order <- 2L
    cfg$minimum_reserve_gb <- 0.25
    root <- paste0(root,"-SMOKE-TEST")
    cat("EXPLICIT SMOKE TEST: reduced dimensions, two repetitions; not paper results.\n")
  }
  dir.create(root,recursive=TRUE,showWarnings=FALSE)
  root <- normalizePath(root,winslash="/",mustWork=TRUE)
  lib <- file.path(root,"_R-library",paste0(R.version$major,".",strsplit(R.version$minor,"\\.")[[1]][1]))
  dir.create(lib,recursive=TRUE,showWarnings=FALSE); .libPaths(c(lib,.libPaths()))
  install_required_packages(packages,lib,cfg$cran_version)
  stopifnot(is.numeric(cfg$roots), !anyDuplicated(cfg$roots),
    is.finite(cfg$repetitions), cfg$repetitions==as.integer(cfg$repetitions),
    cfg$n_rows==as.integer(cfg$n_rows), cfg$m_proxies==as.integer(cfg$m_proxies),
    cfg$continuous_ratio>0, cfg$continuous_ratio<=1,
    cfg$term_fraction>0, cfg$term_fraction<=1,
    cfg$pure_static_ratio>=0, cfg$pure_root_ratio>=0,
    cfg$pure_static_ratio+cfg$pure_root_ratio<1,
    cfg$structural_snr>0,cfg$continuous_snr>0,
    cfg$binary_error_rate>=0,cfg$binary_error_rate<=1,
    cfg$progress_seconds>=1,cfg$poll_seconds>0)
  stopifnot(length(cfg$roots)>0,all(cfg$roots>=2),all(cfg$roots==as.integer(cfg$roots)),
    cfg$n_rows>=10,cfg$m_proxies>=4,cfg$repetitions>=1,
    cfg$reserve_ram_fraction>=0,cfg$reserve_ram_fraction<1,cfg$memory_safety_factor>=1)
  
  logical_cpus <- parallel::detectCores(logical=TRUE)
  if (is.na(logical_cpus)) logical_cpus <- 1L
  original_affinity <- tryCatch(ps::ps_get_cpu_affinity(),error=function(e) seq_len(logical_cpus))
  logical_cpus <- length(original_affinity)
  cpu_set <- head(original_affinity,worker_thread_limit(logical_cpus))
  try(ps::ps_set_cpu_affinity(affinity=as.integer(cpu_set)), silent=TRUE)
  on.exit(try(ps::ps_set_cpu_affinity(affinity=as.integer(original_affinity)), silent=TRUE), add=TRUE)
  
  threads <- worker_thread_limit(logical_cpus)
  max_workers <- worker_count_limit(threads,cfg$max_workers)
  set_threads(1L)
  mem <- system_memory()
  reserve <- max(cfg$minimum_reserve_gb*GB,cfg$reserve_ram_fraction*mem$total)
  # Experiment identity excludes machine/performance settings, includes algorithm sources.
  scientific <- cfg[setdiff(names(cfg),c("output_dir","max_workers","threads_per_worker",
    "reserve_ram_fraction","minimum_reserve_gb","memory_safety_factor","poll_seconds","progress_seconds"))]
  package_identity <- list(version=as.character(cran_description()$Version),
    function_hash=fingerprint(list(formals(Entropic.Scree::Entropic.Scree),
                                   body(Entropic.Scree::Entropic.Scree))))
  signature <- fingerprint(list(version=VERSION,settings=scientific,generator=engine_code,
                                 runner=runner_code,cran=package_identity))
  selected_batch <- resume_batch(root,signature,scientific,engine_code,runner_code,package_identity)
  signature <- selected_batch$signature
  batch <- selected_batch$path
  dir.create(batch,showWarnings=FALSE)
  lock_dir <- file.path(batch,"CONTROLLER-LOCK")
  lock <- file.path(lock_dir,"owner.rds")
  if (!dir.create(lock_dir,showWarnings=FALSE)) {
    last <- if(file.exists(lock)) tryCatch(readRDS(lock),error=function(e) NULL) else NULL
    if (is.null(last)) stop("Controller lock has no readable owner; verify no controller is running before removing CONTROLLER-LOCK.")
    if (!identical(last$host,unname(Sys.info()[["nodename"]])))
      stop("Lock belongs to another computer; verify it has stopped before reusing this output directory.")
    alive <- tryCatch({
      h <- ps::ps_handle(last$pid)
      identical(as.numeric(ps::ps_create_time(h)),last$created)
    },error=function(e) FALSE)
    if (alive) stop("This experiment already has a running controller. Do not launch it twice.")
    # Only this experiment's own lock directory is eligible for stale-lock removal.
    stopifnot(identical(normalizePath(dirname(lock_dir),winslash="/"),normalizePath(batch,winslash="/")),
              basename(lock_dir)=="CONTROLLER-LOCK")
    unlink(lock_dir,recursive=TRUE)
    if (!dir.create(lock_dir,showWarnings=FALSE)) stop("Another controller acquired the lock.")
  }
  atomic_rds(list(pid=Sys.getpid(),created=as.numeric(ps::ps_create_time(ps::ps_handle())),
    host=unname(Sys.info()[["nodename"]])),lock)
  running <- list()
  on.exit({
    for (x in running) try(x$process$kill_tree(),silent=TRUE)
    unlink(lock)
    # rmdir through unlink is restricted to the verified fixed child above.
    if (dir.exists(lock_dir) && !length(list.files(lock_dir,all.files=TRUE,no..=TRUE))) unlink(lock_dir,recursive=TRUE)
  },add=TRUE)
  unlink(file.path(batch,"WAITING-FOR-RESOURCES.txt"))
  atomic_rds(list(settings=cfg,scientific=scientific,signature=signature,version=VERSION,
    threads=threads,max_workers=max_workers,total_ram_gb=mem$total/GB,
    available_ram_gb=mem$available/GB,reserve_ram_gb=reserve/GB),file.path(batch,"SETTINGS.rds"))
  writeLines(capture.output(str(cfg)),file.path(batch,"SETTINGS.txt"))
  atomic_rds(list(policy=MEMORY_POLICY,runner_hash=fingerprint(runner_code),
    changed_at=Sys.time(),host=memory_host(),quantile=.50,headroom=0.80,
    theoretical_floor=0,system_reserve_bytes=reserve),
    file.path(batch,"SCHEDULER-PROVENANCE.rds"))
  writeLines(CONDITION_ORDER,file.path(batch,"ACTIVE-CONDITIONS.txt"))
  atomic_rds(list(package=package_identity,defaults=formals(Entropic.Scree::Entropic.Scree),
    description=cran_description(),math=math_backend()),
    file.path(batch,"CRAN-AND-MATH-PROVENANCE.rds"))
  writeLines(capture.output(str(math_backend())),file.path(batch,"MATH-BACKEND.txt"))
  cat("Full exact eigendecomposition uses R's linked BLAS/LAPACK.\n")
  cat("A threaded MKL/OpenBLAS/BLIS R installation is used when linked; stock single-thread BLAS cannot be accelerated by a thread flag alone.\n")
  engine_path <- file.path(batch,"engine.R"); runner_path <- file.path(batch,"runner.R")
  writeLines(engine_code,engine_path); writeLines(runner_code,runner_path)
  worker_path <- file.path(batch,"worker.R")
  writeLines(c('args <- commandArgs(TRUE)', 'job <- readRDS(args[1])',
    '.libPaths(c(job$lib_paths,.libPaths()))','try(ps::ps_set_cpu_affinity(affinity=as.integer(job$cpu_set)), silent=TRUE)','source(job$runner_path)',
    'source(job$engine_path)','quit(save="no", status=if (run_job(job)) 0L else 1L)'),worker_path)
  source(engine_path,local=.GlobalEnv)
  conditions <- c(orthogonal=Inf,low=15,high=.001)
  seed_offsets <- c(orthogonal=0L,low=1L,high=3L)
  jobs <- list()
  for (rep in seq_len(cfg$repetitions)) for (k in cfg$roots) for (i in seq_along(conditions)) {
    seed_double <- cfg$seed+(rep-1)*100000+(k-20)*100+seed_offsets[[names(conditions)[i]]]
    if (!is.finite(seed_double) || seed_double<0 || seed_double>.Machine$integer.max)
      stop("Seed is outside R's integer range.")
    seed <- as.integer(seed_double)
    id <- sprintf("roots_%02d_%s_rep_%02d_seed_%d",k,names(conditions)[i],rep,seed)
    job <- list(id=id,k=as.integer(k),condition=names(conditions)[i],alpha=unname(conditions[i]),
      repetition=rep,seed=seed,threads=if(is.null(threads)) 1L else threads,cfg=cfg,dir=file.path(batch,id),
      engine_path=engine_path,runner_path=runner_path,
      lib_paths=.libPaths(),signature=signature)
    job$memory <- memory_estimate(k,cfg)
    jobs[[id]] <- job
  }
  if (anyDuplicated(vapply(jobs,function(j) j$seed,integer(1)))) stop("Generated seeds are not unique; adjust the seed schedule.")
  manifest <- do.call(rbind,lapply(jobs,function(j) data.frame(configuration=j$id,true_rank=j$k,
    entanglement=j$condition,alpha=j$alpha,seed=j$seed,repetition=j$repetition,
    blueprint_roots=ceiling(j$k/2),term_pool=j$memory$term_pool,
    estimated_peak_gb=j$memory$bytes/GB)))
  manifest <- manifest[order(match(manifest$entanglement,CONDITION_ORDER),manifest$true_rank,manifest$repetition),]
  data.table::fwrite(manifest,file.path(batch,"CONFIGURATIONS.csv"))
  # Largest configuration first for a conservative, useful calibration.
  pending <- names(jobs)[order(-vapply(jobs,function(j) j$memory$bytes,numeric(1)),
    -vapply(jobs,function(j) j$k,integer(1)))]
  completed <- character(); failed <- character()
  for (id in pending) {
    if (job_complete(jobs[[id]])) {
      cleanup_temporary_data(jobs[[id]])
      completed <- c(completed,id)
    }
  }
  pending <- setdiff(pending,completed)
  write_master(jobs,batch,failed)
  cat(sprintf("\n%d root-count/entanglement configurations x %d seeds = %d datasets; %d already complete. Output:\n%s\n",
    length(cfg$roots)*length(conditions),cfg$repetitions,length(jobs),length(completed),batch))
  cat(sprintf("RAM: %.1f GB available / %.1f GB total; reserve %.1f GB.\nCPU: %d logical; %s threads/job; worker cap %d (RAM decides concurrency).\n",
    mem$available/GB,mem$total/GB,reserve/GB,logical_cpus,
    as.character(threads),max_workers))
  if (length(pending)) cat(sprintf("Estimated largest job: %.1f GB. A fresh full job calibrates when no reusable empirical profile is available.\n",
    max(vapply(jobs[pending],function(j) j$memory$bytes,numeric(1)))/GB))
  # Disk space estimate: datasets + centered NMI + coefficients + temporary files.
  disk_required <- disk_estimate(jobs[pending])
  cat(sprintf("Approximate additional disk reservation: %.1f GiB for retained results plus one temporary working dataset.\n",
    disk_required/GB))
  cat("Temporary generated data are deleted after each verified result; concurrent jobs require additional working space.\n")
  du <- tryCatch(ps::ps_disk_usage(batch),error=function(e) NULL)
  if (!is.null(du) && "available" %in% names(du) && is.finite(min(du$available)) && min(du$available)<disk_required)
    stop(sprintf("Insufficient disk space: approximately %.1f GiB additional space needed for results and temporary working files. Choose a larger output drive.",disk_required/GB))
  cat(sprintf("Approximate disk reservation: %.1f GB (not a strict upper bound).\n",disk_required/GB))
  memory_profiles <- saved_memory_profiles(jobs[completed])
  calibrated <- nrow(memory_profiles)>0L
  cat(sprintf("Whole-job RAM calibration: peak x 0.80 by root count, no theoretical floor. %d saved profiles reused.\n",
    nrow(memory_profiles)))
  rscript <- file.path(R.home("bin"),if(.Platform$OS.type=="windows") "Rscript.exe" else "Rscript")
  last_progress <- Sys.time()-3600; resource_rows <- list()
  while (length(pending) || length(running)) {
    # Refresh RSS and collect finished workers. No completed configuration is rerun.
    for (id in names(running)) {
      x <- running[[id]]
      rss <- tryCatch(as.numeric(ps::ps_memory_info(ps::ps_handle(x$process$get_pid()))[["rss"]]),error=function(e) 0)
      if (is.finite(rss) && rss>0) {
        x$rss <- rss; x$peak <- max(x$peak,rss)
        x$samples <- c(x$samples,rss)
        x$sample_times <- c(x$sample_times,as.numeric(Sys.time()))
      }
      if (!x$process$is_alive()) {
        x$process$wait(timeout=1000)
        ok <- identical(x$process$get_exit_status(),0L) && job_complete(jobs[[id]])
        if (ok) completed <- c(completed,id) else failed <- c(failed,id)
        # Sample the entire lifecycle, including design construction, CRAN and exports.
        # Peak/heap are retained as diagnostics, not substituted for the percentile.
        diagnostic <- file.path(jobs[[id]]$dir,"memory-diagnostics.rds")
        heap <- if (file.exists(diagnostic)) tryCatch(readRDS(diagnostic)$max_R_heap_bytes,
          error=function(e) 0) else 0
        profile <- memory_profile(x$samples)
        eligible <- ok && x$full_lifecycle && profile$samples>=4L && is.finite(profile$median)
        if (eligible) {
          memory_profiles <- rbind(memory_profiles,data.frame(configuration=id,
            k=jobs[[id]]$k,peak_bytes=profile$peak))
        }
        calibrated <- nrow(memory_profiles)>0L
        atomic_rds(c(profile,list(policy=MEMORY_POLICY,signature=jobs[[id]]$signature,
          host=memory_host(),full_lifecycle=x$full_lifecycle,success=ok,
          rss_bytes=x$samples,sample_times=x$sample_times,poll_seconds=cfg$poll_seconds,
          theoretical_bytes=jobs[[id]]$memory$bytes)),file.path(jobs[[id]]$dir,"memory-samples.rds"))
        data.table::fwrite(data.frame(configuration=id,root_count=jobs[[id]]$k,policy=MEMORY_POLICY,
          samples=profile$samples,median_rss_gb=profile$median/GB,peak_rss_gb=x$peak/GB,
          full_lifecycle=x$full_lifecycle,used_for_calibration=eligible,
          next_reservation_gb=memory_reservation(jobs[[id]],memory_profiles)$bytes/GB),
          file.path(batch,"RESOURCE-MEDIAN-CALIBRATION.csv"),append=file.exists(file.path(batch,"RESOURCE-MEDIAN-CALIBRATION.csv")))
        resource_rows[[length(resource_rows)+1L]] <- data.frame(configuration=id,
          estimated_gb=x$estimate/GB,sampled_peak_rss_gb=x$peak/GB,
          max_R_heap_gb=heap/GB,next_calibration_multiplier=NA_real_,
          threads=x$threads,seconds=as.numeric(difftime(Sys.time(),x$started,units="secs")),exit_code=x$process$get_exit_status())
        data.table::fwrite(do.call(rbind,resource_rows),file.path(batch,"RESOURCE-USAGE.csv"),
                          append=file.exists(file.path(batch,"RESOURCE-USAGE.csv")))
        resource_rows <- list()
        cat(sprintf("[%d/%d] %s: %s\n",length(completed)+length(failed),length(jobs),id,if(ok) "complete" else "FAILED (see worker.log)"))
        if (!ok && !file.exists(file.path(jobs[[id]]$dir,"JOB-ERROR.txt")))
          writeLines(paste("Worker exit code",x$process$get_exit_status(),"; see worker.log"),file.path(jobs[[id]]$dir,"JOB-ERROR.txt"))
        running[[id]] <- NULL
        write_master(jobs,batch,failed)
      } else running[[id]] <- x
    }
    mem <- system_memory()
    # Live usage may exceed the planning percentile. Recheck free RAM on every
    # admission; do not kill workers or pretend the percentile bounds their peak.
    future_reserved <- sum(vapply(running,function(x)
      max(0,max(x$estimate,0.80*x$rss)-x$rss),numeric(1)))
    capacity <- mem$available-reserve-future_reserved
    slots <- if (calibrated) max_workers else 1L
    disk_message <- NULL
    while (length(pending) && length(running)<slots) {
      estimates <- vapply(jobs[pending],function(j) memory_reservation(j,memory_profiles)$bytes,numeric(1))
      fits <- which(estimates<=capacity)
      if (!length(fits)) break
      idx <- fits[1]; id <- pending[idx]; job <- jobs[[id]]; estimate <- estimates[idx]
      job$threads <- worker_thread_limit(logical_cpus)
      job$cpu_set <- cpu_set
      # Reserve future growth of active temporary files plus remaining results.
      # Completed data are removed by workers before another job takes their slot.
      disk_need <- launch_disk_required(jobs[unique(c(pending,names(running)))],
        jobs[c(names(running),id)])
      disk_free <- disk_available(batch)
      if (disk_free<disk_need) {
        disk_message <- sprintf("Waiting for disk space: %.1f GiB free; %.1f GiB estimated for results and active temporary files. Completed datasets are cleaned automatically. Free space or use a larger output drive, then rerun if no workers remain.",
          disk_free/GB,disk_need/GB)
        writeLines(disk_message,file.path(batch,"WAITING-FOR-DISK.txt"))
        break
      }
      unlink(file.path(batch,"WAITING-FOR-DISK.txt"))
      dir.create(job$dir,recursive=TRUE,showWarnings=FALSE)
      unlink(file.path(job$dir,"JOB-ERROR.txt"))
      job_file <- file.path(job$dir,"job.rds"); atomic_rds(job,job_file)
      env <- Sys.getenv(); env[c("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS",
         "VECLIB_MAXIMUM_THREADS","BLIS_NUM_THREADS")] <- as.character(job$threads)
      full_lifecycle <- !any(file.exists(file.path(job$dir,c("DATA-COMPLETE.rds","Entropic_Scree.rds"))))
      proc <- processx::process$new(rscript,c("--vanilla",worker_path,job_file),
        stdout=file.path(job$dir,"worker.log"),stderr="2>&1",env=env,cleanup_tree=TRUE,
        windows_hide_window=TRUE)
      running[[id]] <- list(process=proc,estimate=estimate,rss=0,peak=0,
        samples=numeric(),sample_times=numeric(),threads=job$threads,started=Sys.time(),
        full_lifecycle=full_lifecycle)
      pending <- pending[-idx]; capacity <- capacity-estimate
      basis <- memory_reservation(job,memory_profiles)
      cat(sprintf("Started %s | reserve %.1f GB | %d threads | %d active | %s\n",id,estimate/GB,job$threads,length(running),
        if(basis$profiles) sprintf("peak x 0.80 from %d completed profiles at r=%d",basis$profiles,basis$source_k)
        else "initial theoretical calibration"))
    }
    if (length(pending) && !length(running)) {
      if (!is.null(disk_message)) {
        cat(disk_message,"\n")
        break
      }
      need <- min(vapply(jobs[pending],function(j) memory_reservation(j,memory_profiles)$bytes,numeric(1)))
      writeLines(sprintf("Pending jobs need at least %.1f GB available RAM plus %.1f GB system reserve. Currently available %.1f GB. Close other programs or use a larger machine, then rerun. No data sizes were changed.",
        need/GB,reserve/GB,mem$available/GB),file.path(batch,"WAITING-FOR-RESOURCES.txt"))
      cat(readLines(file.path(batch,"WAITING-FOR-RESOURCES.txt")),"\n")
      break
    }
    if (as.numeric(difftime(Sys.time(),last_progress,units="secs"))>=cfg$progress_seconds) {
      lines <- c(sprintf("%s | complete %d/%d | failed %d | running %d | queued %d | available RAM %.1f GB",
        format(Sys.time(),"%H:%M:%S"),length(completed),length(jobs),length(failed),length(running),length(pending),mem$available/GB))
      if (length(pending)) {
        next_need <- min(vapply(jobs[pending],function(j) memory_reservation(j,memory_profiles)$bytes,numeric(1)))
        cpu_limit <- worker_thread_limit(logical_cpus)
        lines <- c(lines,sprintf("  Admission budget: %.1f GiB; smallest next reservation %.1f GiB; system reserve %.1f GiB; target threads per worker %d (%s)",
          capacity/GB,next_need/GB,reserve/GB,cpu_limit,
          if(calibrated) "whole-job peak x 0.80" else "awaiting full-job calibration"))
      }
      for(id in names(running)) {
        path <- file.path(jobs[[id]]$dir,"CURRENT-STAGE.txt")
        stage <- if(file.exists(path)) tryCatch(readLines(path,warn=FALSE)[1],error=function(e) "working") else "starting"
        logpath <- file.path(jobs[[id]]$dir,"worker.log")
        # Tail only, not the entire log, to surface the MI/eigen/generation stage.
        detail <- tryCatch((function() {
          con <- file(logpath,"rb"); on.exit(close(con),add=TRUE)
          seek(con,max(0,file.info(logpath)$size-4096),origin="start")
          tail(strsplit(rawToChar(readBin(con,"raw",4096L)),"\n",fixed=TRUE)[[1]],2L)
        })(),error=function(e) character())
        lines <- c(lines,sprintf("  %s - %s (%.1f GiB sampled RSS; %.1f min elapsed)",id,stage,
          running[[id]]$rss/GB,as.numeric(difftime(Sys.time(),running[[id]]$started,units="mins"))),
          if(length(detail)) paste("    ",paste(trimws(detail),collapse=" | ")))
      }
      writeLines(lines,file.path(batch,"PROGRESS.txt")); cat(paste(lines,collapse="\n"),"\n")
      last_progress <- Sys.time()
    }
    if (length(running)) Sys.sleep(cfg$poll_seconds)
  }
  table <- write_master(jobs,batch,failed,make_plot=TRUE)
  writeLines(c(sprintf("Complete: %d / %d; failed: %d; pending: %d",length(completed),length(jobs),length(failed),length(pending)),
    "Rerun the same script with the same scientific settings/output directory to resume.",
    "Completed numerical results are reused. Temporary generated data are deleted after verified success; incomplete jobs keep restart checkpoints.",
    "Within-one recovery includes exact recovery. Percentages use valid numerical results; failed and pending counts remain visible.",
    "Summary rows group all orthogonal, then low, high; increasing root count within each group.",
    "Read RANK-SUMMARY-COMPACT.csv or RANK-SUMMARY.txt; all per-seed values are in SEED-RESULTS.csv."),file.path(batch,"RUN-SUMMARY.txt"))
  writeLines(c(format(Sys.time(),tz="UTC",usetz=TRUE),readLines(file.path(batch,"RUN-SUMMARY.txt"))),
             file.path(batch,"PROGRESS.txt"))
  cat("\n",readLines(file.path(batch,"RUN-SUMMARY.txt")),sep="\n")
  cat(paste(readLines(file.path(batch,"RANK-SUMMARY.txt")),collapse="\n"),"\n")
  cat("\nResults saved to:\n",batch,"\n")
  invisible(batch)
}

)----"
PLOT_CODE <- r"----(
# Summary figures from saved numerical results; no simulation is required.
summary_figure_data <- function(all) {
  all <- as.data.frame(all)
  required <- c("entanglement","true_rank","status","extracted_rank")
  if (!all(required %in% names(all))) stop("Choose the SEED-RESULTS.csv or SEED-RESULTS.rds file, not an aggregate table.")
  if (!nrow(all) || any(!is.finite(all$true_rank)) ||
      any(!all$entanglement %in% CONDITION_ORDER)) stop("Invalid or empty seed-results table.")
  numeric_ok <- all$status=="complete" & is.finite(all$extracted_rank)
  if (any(all$extracted_rank[numeric_ok]<1) ||
      any(all$extracted_rank[numeric_ok]!=round(all$extracted_rank[numeric_ok])))
    stop("Extracted ranks must be positive integers.")
  summarize_recovery(all)
}

save_summary_figure <- function(all,root,title="Entropic Scree: rank-recovery summary") {
  all <- as.data.frame(all)
  s <- summary_figure_data(all)
  roots <- sort(unique(all$true_rank))
  colours <- c(orthogonal="#0072B2",low="#D55E00",high="#9467BD")
  labels <- c("Orthogonal","Low","High")
  offsets <- c(-.10,0,.10)
  valid <- all$status=="complete" & is.finite(all$extracted_rank)
  planned <- unique(s$seeds_planned)
  counts <- if(length(planned)==1L) paste(planned,"seeds planned per configuration") else "See the table for planned seeds per configuration"
  scale_text <- if(all(c("N","m_generated") %in% names(all))) {
    sprintf("N = %s; generated m = %s | ",format(all$N[1],big.mark=","),
      format(all$m_generated[1],big.mark=","))
  } else ""
  median_limits <- range(c(roots,s$min_extracted_rank,s$max_extracted_rank),na.rm=TRUE)
  if (diff(median_limits)==0) median_limits <- median_limits+c(-.5,.5)
  max_error <- max(c(0,s$mean_absolute_error),na.rm=TRUE)
  panels <- list(
    list(key="exact_match_pct",title="A. Exact-match percentage",ylab="Exact matches (%)",ylim=c(0,100)),
    list(key="within_one_pct",title="B. Within-one percentage",ylab="Within +/-1 root (%)",ylim=c(0,100)),
    list(key="median_extracted_rank",title="C. Median extracted rank and range",ylab="Extracted primary rank",ylim=median_limits),
    list(key="mean_absolute_error",title="D. Mean absolute error",ylab="Mean absolute error (roots)",ylim=c(0,max(1,max_error*1.08)))
  )
  draw <- function() {
    # A shared title/key and footer keep the four metric panels uncluttered.
    layout(matrix(c(1,1,2,3,4,5,6,6),nrow=4,byrow=TRUE),heights=c(.95,3.8,3.8,.75))
    par(mar=rep(0,4),family="sans")
    plot.new()
    text(.5,.85,title,font=2,cex=1.3)
    text(.5,.56,sprintf("%s%d / %d numerical results | %s",scale_text,sum(valid),nrow(all),counts),cex=.84)
    legend("bottom",legend=labels,col=unname(colours),pch=16,lty=1:3,lwd=1.8,
      ncol=3,bty="n",cex=.92,x.intersp=.9,y.intersp=1,
      text.width=strwidth("Orthogonal")*1.2)
    for (panel in panels) {
      par(mar=c(4.1,4.8,2.9,1.2),mgp=c(2.7,.8,0),tcl=-.25)
      plot(NA,xlim=range(roots)+c(-.35,.35),ylim=panel$ylim,xaxt="n",yaxt="n",
        xlab="Number of generating roots",ylab=panel$ylab,main=panel$title,
        cex.main=1.02,col.axis="#303840",col.lab="#303840",bty="l",xaxs="i")
      ticks <- if(panel$key %in% c("exact_match_pct","within_one_pct")) seq(0,100,20) else axTicks(2)
      abline(h=ticks,col="#e4e8eb",lwd=.7)
      axis(1,at=roots,col="#8b969e",col.axis="#303840")
      axis(2,at=ticks,las=1,col="#8b969e",col.axis="#303840")
      if (panel$key=="median_extracted_rank") abline(0,1,col="#666666",lty=2,lwd=1)
      for (i in seq_along(CONDITION_ORDER)) {
        condition <- CONDITION_ORDER[i]
        z <- s[s$entanglement==condition,,drop=FALSE]
        z <- z[match(roots,z$true_rank),,drop=FALSE]
        x <- roots+offsets[i]; y <- z[[panel$key]]
        ok <- is.finite(y)
        # Keep NA entries: missing configurations break lines instead of being
        # silently interpolated. Horizontal offsets only separate overlapping keys.
        lines(x,y,col=colours[i],lty=i,lwd=1.8)
        if (panel$key=="median_extracted_rank" && any(ok)) {
          segments(x[ok],z$min_extracted_rank[ok],x[ok],z$max_extracted_rank[ok],
            col=adjustcolor(colours[i],.6),lwd=1.1)
          segments(x[ok]-.028,z$min_extracted_rank[ok],x[ok]+.028,z$min_extracted_rank[ok],
            col=colours[i],lwd=1.1)
          segments(x[ok]-.028,z$max_extracted_rank[ok],x[ok]+.028,z$max_extracted_rank[ok],
            col=colours[i],lwd=1.1)
        }
        complete <- ok & !is.na(z$all_seeds_valid) & z$all_seeds_valid
        partial <- ok & !complete
        points(x[complete],y[complete],pch=16,cex=.8,col=colours[i])
        points(x[partial],y[partial],pch=21,cex=.82,col=colours[i],bg="white",lwd=1.2)
      }
      if (!any(is.finite(s[[panel$key]]))) text(mean(range(roots)),mean(panel$ylim),"No numerical results yet",col="#65727a")
    }
    par(mar=rep(0,4))
    plot.new()
    text(.5,.85,"Within +/-1 includes exact matches. Vertical bars in C show minimum-maximum, not confidence intervals.",cex=.79)
    text(.5,.54,"Filled dots: all planned seeds completed. Open dots: partial results; percentages use completed numerical results.",cex=.79)
    text(.5,.23,"Dashed grey line in C: extracted rank = generating rank. Series are slightly offset horizontally to reveal overlap.",cex=.79)
  }
  dir.create(root,recursive=TRUE,showWarnings=FALSE)
  png(file.path(root,"RANK-SUMMARY.png"),width=2400,height=2050,res=200)
  tryCatch(draw(),finally=dev.off())
  pdf(file.path(root,"RANK-SUMMARY.pdf"),width=12,height=10.25,useDingbats=FALSE)
  tryCatch(draw(),finally=dev.off())
  cat("Four-panel summary saved:\n",file.path(root,"RANK-SUMMARY.png"),"\n",
    file.path(root,"RANK-SUMMARY.pdf"),"\n")
  invisible(s)
}

plot_saved_summary <- function(path=NULL) {
  if (is.null(path)) {
    cat("Choose SEED-RESULTS.csv or SEED-RESULTS.rds from the completed/partial experiment. No simulations will run.\n")
    path <- file.choose()
  }
  if (!file.exists(path)) stop("Saved seed-results file does not exist: ",path)
  if (!basename(path) %in% c("SEED-RESULTS.csv","SEED-RESULTS.rds"))
    stop("Select SEED-RESULTS.csv or SEED-RESULTS.rds.")
  all <- if (tolower(tools::file_ext(path))=="rds") readRDS(path) else
    read.csv(path,stringsAsFactors=FALSE,check.names=FALSE)
  save_summary_figure(all,dirname(normalizePath(path,winslash="/",mustWork=TRUE)))
}

)----"
env <- new.env(parent = globalenv())
eval(parse(text = RUNNER_CODE), envir = env)
eval(parse(text = PLOT_CODE), envir = env)
if (PLOT_ONLY) {
  env$plot_saved_summary(PLOT_INPUT)
} else {
  results_dir <- env$run_sweep(CONFIG, ENGINE_CODE, RUNNER_CODE)
  env$plot_saved_summary(file.path(results_dir,"SEED-RESULTS.rds"))
}
})