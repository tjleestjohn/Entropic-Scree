# ENTROPIC SCREE: ORTHOGONAL, LOW AND HIGH COMPARISONS AND BASELINE BINNING ABLATION
# Open this file in RStudio, select all, and Run. R >= 4.2 is required.
# Main section: r=20 orthogonal dataset, repetition 1, seed 19862026.
# Includes six-panel broad/zoom ablation figures and low/high four-method comparisons.
# Low/high seeds: 19862027/19862029; binning ablation uses the main dataset only.
# N=5,000 records; m=10,000 measurements; 50% continuous and 50% binary.
# All methods within each entanglement context use the same observations.
# Missing dependencies install automatically; Entropic.Scree is pinned to 1.0.1.
# An archived source package can require Rtools/compiler tools on a new machine.
# Results are saved beneath getwd(). Rerunning resumes completed stages.
# The large observed-data checkpoint is deleted after successful completion.
# No true-data matrix or centered NMI matrix is retained on disk.

CONFIG <- list(
  r = 20L,
  entanglement = "orthogonal",   # "orthogonal", "low", or "high"
  repetition = 1L,
  n_rows = 5000L,
  m_proxies = 10000L,
  cran_version = "1.0.1",
  seed_base = 19862026L,
  continuous_ratio = 0.50,
  max_interaction_order = 5L,
  max_polynomial_order = 4L,
  structural_snr = 4,
  continuous_snr = 4,
  binary_error_rate = 0.0528,
  pure_static_ratio = 0.05,
  pure_root_ratio = 0,
  term_fraction = 0.10,
  num_cores = NULL,              # NULL uses detected logical processors minus 2.
  output_dir = file.path(getwd(), "Entropic-Scree-Single-Run")
)

stage_message <- function(...) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ..., "\n", sep="")
  flush.console()
}

ensure_packages <- function(lib, version) {
  if (getRversion() < "4.2.0") stop("Use R 4.2 or later.", call.=FALSE)
  dir.create(lib, recursive=TRUE, showWarnings=FALSE)
  .libPaths(c(lib, .libPaths()))
  old_timeout <- getOption("timeout")
  options(timeout=max(600, old_timeout))
  on.exit(options(timeout=old_timeout), add=TRUE)
  ready <- function(p) {
    if (!requireNamespace(p, quietly=TRUE)) return(FALSE)
    d <- suppressWarnings(utils::packageDescription(p))
    is.list(d) && length(d$Version)==1L && !is.na(d$Version)
  }
  repo <- "https://cloud.r-project.org"
  required <- c("MASS", "Matrix", "data.table", "stringr", "clusterGeneration", "ggplot2", "patchwork")
  missing <- required[!vapply(required, ready, logical(1))]
  if (length(missing)) {
    stage_message("Installing dependencies: ", paste(missing, collapse=", "))
    utils::install.packages(missing, lib=lib, repos=repo, Ncpus=1L)
  }
  correct <- ready("Entropic.Scree") &&
    identical(as.character(utils::packageVersion("Entropic.Scree")), version)
  if (!correct) {
    if (isNamespaceLoaded("Entropic.Scree")) {
      tryCatch(unloadNamespace("Entropic.Scree"), error=function(e)
        stop("Restart R and rerun to install Entropic.Scree ", version, ".", call.=FALSE))
    }
    available <- tryCatch(suppressWarnings(utils::available.packages(repos=repo)),
      error=function(e) matrix(character(), 0, 0))
    stage_message("Installing Entropic.Scree ", version, " in the study library.")
    tryCatch({
      if ("Entropic.Scree" %in% rownames(available) &&
          identical(unname(available["Entropic.Scree", "Version"]), version)) {
        utils::install.packages("Entropic.Scree", lib=lib, repos=repo, Ncpus=1L)
      } else {
        if (!ready("remotes")) utils::install.packages("remotes", lib=lib, repos=repo)
        remotes::install_version("Entropic.Scree", version=version, lib=lib,
          repos=repo, dependencies=NA, upgrade="never", force=TRUE)
      }
    }, error=function(e) stop("Package installation failed: ", conditionMessage(e),
      " Source installation may require Rtools on Windows or compiler tools on other systems.", call.=FALSE))
  }
  missing <- c(required, "Entropic.Scree")[!vapply(c(required, "Entropic.Scree"), ready, logical(1))]
  if (length(missing)) stop("Unavailable packages: ", paste(missing, collapse=", "),
    ". Check the installation messages, Internet access and R version, then rerun.", call.=FALSE)
  if (as.character(utils::packageVersion("Entropic.Scree")) != version)
    stop("Entropic.Scree must be version ", version, ". Restart R and rerun.", call.=FALSE)
}

save_rds <- function(object, path, compress=TRUE) {
  temporary <- paste0(path, ".partial")
  saveRDS(object, temporary, compress=compress)
  if (file.exists(path) && !file.remove(path)) stop("Cannot replace ", path)
  if (!file.rename(temporary, path)) stop("Cannot finalize ", path)
  invisible(path)
}

read_checkpoint <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(readRDS(path), error=function(e) {
    stage_message("Unreadable checkpoint will be recomputed: ", basename(path))
    NULL
  })
}

fingerprint <- function(object) {
  f <- tempfile(fileext=".rds")
  on.exit(unlink(f), add=TRUE)
  saveRDS(object, f, version=2, compress=FALSE)
  unname(tools::md5sum(f))
}

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

  m_cont <- floor(m_proxies * continuous_ratio)
  m_bin <- m_proxies - m_cont

  s1_df <- as.data.frame(s1_continuous)
  colnames(s1_df) <- paste0("X", 1:k)

  eff_cross_order <- min(k, max_interaction_order)
  formula_str <- if(eff_cross_order == 1) as.formula("~ .") else as.formula(paste0("~ .^", eff_cross_order))

  cat(sprintf("      -> Expanding continuous roots up to interaction order %d and polynomial order %d...\n", eff_cross_order, max_polynomial_order))
  design_mat <- model.matrix(formula_str, data = s1_df)[, -1, drop = FALSE]

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
      design_mat <- cbind(design_mat, as.matrix(as.data.frame(power_list, check.names = FALSE)))
    }
  }

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

  term_dependencies <- lapply(term_names, function(name) {
    matches <- stringr::str_extract_all(name, "X\\d+")[[1]]
    unique(as.integer(gsub("X", "", matches)))
  })

  term_sds <- rep(1.0, n_terms)

  coefficient_rng <- .Random.seed
  is_main <- interaction_orders == 1
  is_int <- interaction_orders > 1
  n_main <- sum(is_main)
  n_int <- sum(is_int)

  blueprint_roots <- vector("list", m_proxies)
  for(j in 1:m_proxies) {
    if (n_main > 0) invisible(rnorm(n_main, mean = 0, sd = term_sds[is_main]))
    if (n_int > 0) invisible(rnorm(n_int, mean = 0, sd = term_sds[is_int] * int_scaling))
  }

  m_static <- floor(m_proxies * pure_static_ratio)
  m_pure_root <- floor(m_proxies * pure_root_ratio)
  m_mixed <- m_proxies - m_static - m_pure_root

  min_roots_per_proxy <- max(1, ceiling(k * MIN_ROOT_SPARSITY_RATIO))
  max_roots_per_proxy <- max(min_roots_per_proxy, ceiling(k * MAX_ROOT_SPARSITY_RATIO))

  cat(sprintf("      -> Injecting %d Pure Static Proxies (%.1f%% of total dataset)...\n", m_static, pure_static_ratio * 100))
  cat(sprintf("      -> Injecting %d Pure Root Proxies (%.1f%% of total dataset, 1 main effect only)...\n", m_pure_root, pure_root_ratio * 100))
  cat(sprintf("      -> Generating %d Entangled Mixed Proxies...\n", m_mixed))
  cat(sprintf("         * Activating %d to %d Base Roots per proxy (%.0f%% to %.0f%% Root Participation)\n", min_roots_per_proxy, max_roots_per_proxy, MIN_ROOT_SPARSITY_RATIO * 100, MAX_ROOT_SPARSITY_RATIO * 100))
  cat(sprintf("         * Randomly sampling %.0f%% to %.0f%% of ALL aligned terms (Main & Higher-Order).\n", MIN_SPARSITY_RATIO * 100, MAX_SPARSITY_RATIO * 100))

  selected_by_column <- rep(list(integer()),m_proxies)
  dependency_bits <- if (k<=30L) vapply(term_dependencies,
    function(deps) as.integer(sum(2^(deps-1L))),integer(1)) else NULL
  aligned_cache <- new.env(hash=TRUE,parent=emptyenv())
  main_effect_indices <- which(is_main)

  if (m_pure_root > 0) {
    for (j in 1:m_pure_root) {
      root_idx <- sample(main_effect_indices, 1)
      selected_by_column[[j]] <- root_idx
      blueprint_roots[[j]] <- term_dependencies[[root_idx]]
    }
  }

  if (m_mixed > 0) {
    start_col <- m_pure_root + 1
    end_col <- m_pure_root + m_mixed
    for (j in start_col:end_col) {

      if (min_roots_per_proxy == max_roots_per_proxy) {
        num_roots <- min_roots_per_proxy
      } else {
        num_roots <- sample(min_roots_per_proxy:max_roots_per_proxy, 1)
      }

      selected_roots <- sample(1:k, num_roots)
      blueprint_roots[[j]] <- selected_roots

      key <- paste(sort(selected_roots),collapse=",")
      if (exists(key,aligned_cache,inherits=FALSE)) {
        aligned_term_indices <- get(key,aligned_cache,inherits=FALSE)
      } else {
        aligned_term_indices <- if (!is.null(dependency_bits)) {
          bits <- as.integer(sum(2^(selected_roots-1L)))
          which(bitwAnd(dependency_bits,bits)==dependency_bits)
        } else which(vapply(term_dependencies,function(deps)
          all(deps %in% selected_roots),logical(1)))
        if (length(aligned_cache)<512L) assign(key,aligned_term_indices,aligned_cache)
      }

      if (length(aligned_term_indices) > 0) {
        absolute_min_to_sample <- max(1, round(length(aligned_term_indices) * MIN_SPARSITY_RATIO))
        absolute_max_to_sample <- max(absolute_min_to_sample, round(length(aligned_term_indices) * MAX_SPARSITY_RATIO))

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

  active_indices <- sort(unique(unlist(selected_by_column,use.names=FALSE)))
  active_terms <- length(active_indices)

  if (active_terms > 0) {
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
  deterministic_block <- NULL

  cat(sprintf("      -> Setting Base Rates and Injecting Structural Uncertainty (Round 1 Noise, SNR = %.2f)...\n", structural_snr))
  true_manifestation <- matrix(0, nrow = n, ncol = m_proxies)
  core_center <- core_scale <- thresholds <- numeric(m_proxies)

  total_sd <- sqrt(1 + (1 / structural_snr))

  for (j in 1:m_proxies) {
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

    is_binary_var <- (j > m_cont)
    if (is_binary_var) {
      target_prev <- runif(1, min = 0.10, max = 0.90)

      shift <- sqrt(2) * total_sd * qnorm(target_prev)

      thresholds[j] <- shift
      core_vec <- core_vec + shift
    }

    structural_noise <- rnorm(n, mean = 0, sd = sqrt(1 / structural_snr))

    true_manifestation[, j] <- core_vec + structural_noise
  }

  cat(sprintf("      -> Generating %d True Continuous Proxies...\n", m_cont))
  true_cont <- true_manifestation[, 1:m_cont, drop = FALSE]

  if (m_bin > 0) {
    cat(sprintf("      -> Generating %d True Binary Proxies...\n", m_bin))
    signal_bin <- true_manifestation[, (m_cont + 1):m_proxies, drop = FALSE]

    apply_copula_mapping <- function(scores) {
      if (sd(scores) < 1e-9) return(rep(runif(1, min = 0.50, max = 0.95), length(scores)))

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

apply_measurement_error <- function(true_universe, snr_continuous = 2.0, binary_error_rate = 0.15) {
  cat(sprintf("      -> Applying Measurement Error (Round 2 Noise: Continuous SNR = %.2f, Binary Bit-Flip Rate = %.3f)...\n", snr_continuous, binary_error_rate))

  obs_mat <- true_universe$data_matrix
  is_cont <- true_universe$is_continuous
  n <- nrow(obs_mat)
  m <- ncol(obs_mat)

  for (j in 1:m) {
    if (is_cont[j]) {
      true_var <- var(obs_mat[, j])
      if (true_var < 1e-9) true_var <- 1e-9

      measurement_error <- rnorm(n, mean = 0, sd = 1)
      error_var <- var(measurement_error)

      scaling_factor <- sqrt(true_var / (error_var * snr_continuous))
      obs_mat[, j] <- obs_mat[, j] + (measurement_error * scaling_factor)

    } else {
      flip_mask <- rbinom(n, 1, binary_error_rate)
      obs_mat[, j] <- abs(obs_mat[, j] - flip_mask)
    }
  }

  obs_dt <- data.table::as.data.table(obs_mat)
  data.table::setnames(obs_dt, paste0("V", 1:m))
  return(obs_dt)
}

generate_dataset <- function(cfg, seed) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  initial_rng <- .Random.seed
  alpha <- c(orthogonal=Inf, low=15, high=0.001)[[cfg$entanglement]]
  sigma <- if (is.infinite(alpha)) diag(cfg$r) else
    clusterGeneration::rcorrmatrix(cfg$r, alphad=alpha)
  root_eigen <- pmax(eigen(sigma, symmetric=TRUE)$values, 0)
  p <- root_eigen[root_eigen>1e-9] / sum(root_eigen[root_eigen>1e-9])
  z <- MASS::mvrnorm(n=cfg$n_rows, mu=rep(0, cfg$r), Sigma=sigma)
  truth <- generate_true_mixed_proxies(z, cfg$m_proxies,
    max_interaction_order=cfg$max_interaction_order,
    max_polynomial_order=cfg$max_polynomial_order, int_scaling=1,
    continuous_ratio=cfg$continuous_ratio,
    MIN_SPARSITY_RATIO=cfg$term_fraction, MAX_SPARSITY_RATIO=cfg$term_fraction,
    MIN_ROOT_SPARSITY_RATIO=0.5, MAX_ROOT_SPARSITY_RATIO=0.5,
    structural_snr=cfg$structural_snr, pure_static_ratio=cfg$pure_static_ratio,
    pure_root_ratio=cfg$pure_root_ratio)
  observed <- apply_measurement_error(truth, cfg$continuous_snr, cfg$binary_error_rate)
  meta <- list(seed=seed, initial_rng=initial_rng, final_rng=.Random.seed,
    settings=cfg, latent_correlation=sigma, root_eigenvalues=root_eigen,
    root_effective_rank=exp(-sum(p*log(p))),
    mean_absolute_root_correlation=mean(abs(sigma[upper.tri(sigma)])),
    blueprint_roots=ceiling(cfg$r/2), active_terms=truth$active_terms,
    expanded_effective_rank=truth$K_eff, expanded_eigenvalues=truth$expanded_eigenvalues,
    is_continuous=truth$is_continuous)
  list(observed=observed, metadata=meta)
}

fit_pca <- function(observed, ranks=FALSE) {
  x <- as.matrix(observed)
  if (ranks) x <- apply(x, 2, data.table::frank) # Average ranks for ties.
  if (any(!is.finite(x)) || any(apply(x, 2, sd)==0))
    stop("PCA requires finite, nonconstant generated columns; no columns were silently removed.")
  fit <- stats::prcomp(x, center=TRUE, scale.=TRUE, retx=FALSE)
  values <- fit$sdev^2
  list(eigenvalues=values, threshold=1, threshold_count=sum(values>1),
    threshold_rule="Kaiser: correlation eigenvalue > 1",
    rank_ceiling=min(nrow(x)-1L, ncol(x)), n_records=nrow(x), n_measures=ncol(x))
}

fit_kernel_pca <- function(observed) {
  x <- scale(as.matrix(observed))
  if (any(!is.finite(x))) stop("Kernel PCA requires finite, nonconstant columns.")
  n <- nrow(x); m <- ncol(x)
  norms <- rowSums(x^2)
  distance2 <- outer(norms, norms, "+") - 2*tcrossprod(x)
  gamma <- 1/m
  kernel <- exp(-gamma*distance2)
  rm(x, distance2); gc(FALSE)
  # Exactly H K H, using means instead of dense products with 11'/N.
  centered <- sweep(sweep(kernel, 1, rowMeans(kernel), "-"), 2, colMeans(kernel), "-") + mean(kernel)
  rm(kernel); gc(FALSE)
  values <- pmax(eigen(centered, symmetric=TRUE, only.values=TRUE)$values, 0)
  list(eigenvalues=values, gamma=gamma, threshold=mean(values),
    threshold_count=sum(values>mean(values)),
    threshold_rule="Centered-kernel eigenvalue > mean of all N kernel eigenvalues",
    rank_ceiling=n-1L, n_records=n, n_measures=m)
}

choose_bin_targets <- function(observed) {
  n <- nrow(observed)
  candidates <- which(vapply(observed, function(x) data.table::uniqueN(x)>10, logical(1)))
  set.seed(42)
  selected <- if (length(candidates)) sample(names(observed)[candidates], min(100L, length(candidates))) else character()
  counts <- if (length(selected)) vapply(selected,
    function(j) grDevices::nclass.FD(observed[[j]]), numeric(1)) else numeric()
  fd <- if (length(counts)) round(mean(counts)) else ceiling(n^(2/3))
  if (!is.finite(fd) || fd<2) stop("The FD-derived target is not usable: ", fd)
  list(default=max(2, ceiling(n^(1/3))), dense=ceiling(sqrt(n)), fd=fd,
    fd_sample_seed=42L, fd_measurements=selected, fd_counts=counts,
    rule="Rounded mean nclass.FD over up to 100 sampled columns with >10 unique values")
}

fit_entropic <- function(observed, cfg, seed, cores, bins, native_pdf) {
  stopifnot(as.character(utils::packageVersion("Entropic.Scree"))==cfg$cran_version)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(seed)
  # Only the ablation changes num_bins. Preparation and detector defaults stay unchanged.
  args <- list(data=data.table::copy(observed), interactive_mode=FALSE, num_cores=cores)
  if (!is.null(bins)) args$num_bins <- bins
  grDevices::pdf(native_pdf, width=12, height=6)
  on.exit(grDevices::dev.off(), add=TRUE)
  result <- do.call(Entropic.Scree::Entropic.Scree, args)
  result$similarity_matrix <- NULL
  result$package_version <- as.character(utils::packageVersion("Entropic.Scree"))
  result$call_arguments <- args[setdiff(names(args), "data")]
  result$bin_target <- if (is.null(bins)) max(2, ceiling(nrow(observed)^(1/3))) else bins
  result
}

save_figure <- function(path, draw, width=12, height=8) {
  grDevices::png(paste0(path, ".png"), width=width*160, height=height*160, res=160)
  tryCatch(draw(), finally=grDevices::dev.off())
  grDevices::pdf(paste0(path, ".pdf"), width=width, height=height)
  tryCatch(draw(), finally=grDevices::dev.off())
}

plot_original_comparison <- function(all_results,cfg,meta) {
  pca_eigenvalues <- all_results$standard$eigenvalues
  spearman_eigenvalues <- all_results$rank$eigenvalues
  kpca_eigenvalues <- all_results$kernel$eigenvalues
  pca_kaiser <- all_results$standard$threshold_count
  spearman_kaiser <- all_results$rank$threshold_count
  results <- all_results$entropic_default
  m_total <- cfg$m_proxies; n_rows <- cfg$n_rows
  M_PROXIES <- length(results$retained_features); K_TRUE <- cfg$r
  true_universe <- list(active_terms=meta$active_terms,K_eff=meta$expanded_effective_rank)
  pad <- function(e) c(e,rep(0,max(0,m_total-length(e))))[seq_len(m_total)]
  df_compare <- data.frame(Rank=rep(seq_len(m_total),4),
    Eigenvalue=c(pad(pca_eigenvalues),pad(spearman_eigenvalues),pad(kpca_eigenvalues),pad(results$eigenvalues)),
    Method=factor(rep(c("Standard PCA","Spearman PCA","Kernel PCA","Entropic Scree"),each=m_total),
      levels=c("Standard PCA","Spearman PCA","Kernel PCA","Entropic Scree")))
if (requireNamespace("ggplot2", quietly = TRUE)) {

  pca_y_max <- if(length(pca_eigenvalues) >= 2) pca_eigenvalues[2] * 1.1 else max(pca_eigenvalues)
  spear_y_max <- if(length(spearman_eigenvalues) >= 2) spearman_eigenvalues[2] * 1.1 else max(spearman_eigenvalues)
  kpca_y_max <- if(length(kpca_eigenvalues) >= 2) kpca_eigenvalues[2] * 1.1 else max(kpca_eigenvalues)
  ent_y_max <- if(length(results$eigenvalues) >= 2) results$eigenvalues[2] * 1.1 else max(results$eigenvalues)

  ent_y_min <- min(results$eigenvalues[results$eigenvalues > 1e-8])
  pca_y_min <- 0
  kpca_y_min <- 0

  log10_breaks <- function(x) {
    10^seq(floor(log10(min(x))), ceiling(log10(max(x))))
  }

  K_rlzd <- true_universe$active_terms
  K_eff <- true_universe$K_eff
  ref_size <- 3.0 # Slightly scaled down for 4-panel density

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

  pca_data <- df_compare[df_compare$Method == "Standard PCA", ]
  spear_data <- df_compare[df_compare$Method == "Spearman PCA", ]
  kpca_data <- df_compare[df_compare$Method == "Kernel PCA", ]
  ent_data <- df_compare[df_compare$Method == "Entropic Scree" & df_compare$Rank < M_PROXIES, ]

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
}

plot_comparison <- function(results, cfg, meta, annotated=TRUE, cap=TRUE) {
  if(annotated) return(plot_original_comparison(results,cfg,meta))
  par(mfrow=c(2,2),mar=c(4.2,4.2,3.6,1),oma=c(0,0,2.5,0))
  keys <- c("standard","rank","kernel","entropic_default")
  labels <- c("Standard PCA","Spearman rank PCA","RBF kernel PCA","Entropic Scree")
  colors <- c("firebrick","darkgoldenrod","mediumorchid4","dodgerblue")
  for(i in seq_along(keys)) {
    e <- results[[keys[i]]]$eigenvalues
    ymax <- if(cap && length(e)>1L) 1.1*e[2] else 1.05*max(e)
    plot(seq_along(e),e,type="l",log="x",col=colors[i],lwd=1.2,
      ylim=c(0,max(ymax,.Machine$double.eps)),xlab="Eigenvalue index (log scale)",
      ylab=if(i==4L) "Working eigenvalue" else "Eigenvalue",main=labels[i])
    if(cap && e[1]>ymax) mtext("Leading eigenvalue exceeds display range",side=3,line=.3,cex=.67)
  }
  mtext("Spectral comparison without reference boundaries",outer=TRUE,line=.6,font=2)
}

plot_scree <- function(obj, cfg, annotated=TRUE) {
  e <- obj$eigenvalues
  if(!annotated) {
    par(mfrow=c(1,2),mar=c(4.2,4.4,3,1))
    for(last in c(length(e),min(length(e),75L))) {
      plot(seq_len(last),e[seq_len(last)],type="o",log="y",col="dodgerblue",pch=16,cex=.35,
        xlab="Eigenvalue index",ylab="Working eigenvalue (log scale)",
        main=if(last==length(e)) "Full spectrum" else "First 75 axes")
    }
    return(invisible(NULL))
  }
  primary <- obj$K_roots; extended <- obj$K_extended
  if(primary==extended) {
    lines_geom <- ggplot2::geom_vline(xintercept=primary,color="forestgreen",linetype="dashed",linewidth=1.2)
    subtitle <- sprintf("Converged: Primary Rank & Extended Signal Tail = %d",primary)
  } else {
    lines_geom <- list(
      ggplot2::geom_vline(xintercept=primary,color="#D55E00",linetype="dashed",linewidth=1.2),
      ggplot2::geom_vline(xintercept=extended,color="purple",linetype="twodash",linewidth=1.2))
    subtitle <- sprintf("Primary Rank (Orange) = %d | Extended Signal Tail (Purple) = %d",primary,extended)
  }
  last <- max(2L,sum(e>1e-9))
  zoom_start <- max(1L,min(primary,extended)-5L)
  zoom_end <- min(last,max(primary,extended)+15L)
  zoom_start <- min(zoom_start,zoom_end)
  macro_base <- max(50L,max(primary,extended)*10L)
  if(is.finite(obj$top_of_bulk)) macro_base <- max(macro_base,obj$top_of_bulk+25L)
  macro_end <- min(last,macro_base)
  panel <- function(ids,title,macro=FALSE) {
    d <- data.frame(Rank=ids,Eigenvalue=e[ids])
    p <- ggplot2::ggplot(d,ggplot2::aes(x=Rank,y=Eigenvalue))+
      ggplot2::geom_line(color="dodgerblue",linewidth=1)+
      ggplot2::geom_point(color="dodgerblue",size=2)+lines_geom+
      ggplot2::scale_y_continuous(trans="log10")+
      ggplot2::labs(title=title,x="Eigenvalue Index",y="Log(Eigenvalue)")+
      ggplot2::theme_minimal(base_size=14)+
      ggplot2::theme(plot.title=ggplot2::element_text(face="bold",hjust=.5,size=12))
    if(macro) p <- p+ggplot2::coord_cartesian(ylim=c(NA,1.1*e[2])) else
      p <- p+ggplot2::scale_x_continuous(breaks=function(x) unique(floor(pretty(seq(min(x),max(x))))))
    p
  }
  p <- (panel(seq_len(macro_end),"Macro View",TRUE)+panel(seq.int(zoom_start,zoom_end),"Zoomed View"))+
    patchwork::plot_annotation(title="Entropic Scree Results",subtitle=subtitle,
      theme=ggplot2::theme(plot.title=ggplot2::element_text(face="bold",size=16,hjust=.5),
        plot.subtitle=ggplot2::element_text(size=14,hjust=.5)))
  print(p)
}

plot_ablation <- function(results,cfg) {
  keys <- c("entropic_default","entropic_dense","entropic_fd")
  panels <- lapply(seq_along(keys),function(i) {
    obj <- results[[keys[i]]]; e <- obj$eigenvalues
    bins <- obj$bin_target
    title <- switch(i,
      bquote(atop(bold("Conservative Default"~(N^{1/3})),bold(.(sprintf("[%d Bins]",bins))))),
      bquote(atop(bold("Dense"~(N^{1/2})),bold(.(sprintf("[%d Bins]",bins))))),
      bquote(atop(bold("Freedman-Diaconis"~phantom(N^{1/2})),bold(.(sprintf("[%d Bins]",bins))))))
    ids <- seq_len(min(75L,length(e)))
    ggplot2::ggplot(data.frame(Rank=ids,Eigenvalue=e[ids]),ggplot2::aes(x=Rank,y=Eigenvalue))+
      ggplot2::geom_line(color="dodgerblue",linewidth=1)+
      ggplot2::geom_point(color="dodgerblue",size=1.5)+
      ggplot2::geom_vline(xintercept=cfg$r,color="#D55E00",linetype="dashed",linewidth=1.2)+
      ggplot2::scale_y_continuous(trans="log10")+
      ggplot2::coord_cartesian(ylim=c(NA,e[2]*1.025))+
      ggplot2::labs(title=title,x="Eigenvalue Index",y=if(i==1L) "Log(Eigenvalue)" else "")+
      ggplot2::theme_minimal(base_size=14)+
      ggplot2::theme(plot.title=ggplot2::element_text(face="bold",hjust=.5,size=14),
        panel.grid.minor=ggplot2::element_blank(),axis.title.y=ggplot2::element_text(margin=ggplot2::margin(r=10)))
  })
  p <- patchwork::wrap_plots(panels,ncol=3)+patchwork::plot_annotation(
    title="Eigenspectra Across Discretization Regimes",
    subtitle=sprintf("Same generated observations; dashed line: r = %d\nAutomatic primary selections: %s",cfg$r,
      paste(vapply(keys,function(k) results[[k]]$K_roots,numeric(1)),collapse=", ")),
    theme=ggplot2::theme(plot.title=ggplot2::element_text(face="bold",size=18,hjust=.5),
      plot.subtitle=ggplot2::element_text(size=14,hjust=.5,color="gray30",margin=ggplot2::margin(b=15))))
  print(p)
}

make_outputs <- function(results, meta, targets, cfg, seed, out) {
  keys <- c("standard","rank","kernel","entropic_default")
  labels <- c("Standard PCA","Spearman rank PCA","RBF kernel PCA","Entropic Scree")
  comparison <- do.call(rbind,lapply(seq_along(keys),function(i) {
    o <- results[[keys[i]]]; es <- i==4L
    data.frame(method=labels[i],generating_count=cfg$r,seed=seed,N=cfg$n_rows,
      m_input=cfg$m_proxies,m_analyzed=if(es) length(o$retained_features) else cfg$m_proxies,
      rank_ceiling=if(es) length(o$retained_features)-1L else o$rank_ceiling,
      automatic_primary_boundary=if(es) o$K_roots else NA_real_,
      threshold_count=if(es) sum(o$eigenvalues>mean(o$eigenvalues)) else o$threshold_count,
      threshold=if(es) mean(o$eigenvalues) else o$threshold,
      threshold_rule=if(es) "Working eigenvalue > mean (separate diagnostic)" else o$threshold_rule,
      seconds=o$seconds,stringsAsFactors=FALSE)
  }))
  write.csv(comparison,file.path(out,"COMPARISON-SUMMARY.csv"),row.names=FALSE,na="")
  for(key in names(results)) {
    o <- results[[key]]
    write.csv(data.frame(axis=seq_along(o$eigenvalues),eigenvalue=o$eigenvalues),
      file.path(out,paste0("SPECTRUM-",key,".csv")),row.names=FALSE)
    if(startsWith(key,"entropic_")) {
      write.csv(data.frame(axis=seq_along(o$FSIG_final),FSIG=o$FSIG_final),
        file.path(out,paste0("FSIG-",key,".csv")),row.names=FALSE)
      write.csv(data.frame(measure=o$retained_features),
        file.path(out,paste0("RETAINED-",key,".csv")),row.names=FALSE)
    }
  }
  ablation <- do.call(rbind,lapply(c("entropic_default","entropic_dense","entropic_fd"),function(key) {
    o <- results[[key]]; m <- length(o$retained_features)
    data.frame(rule=key,bin_target=o$bin_target,m_retained=m,
      R_eff=o$R_eff,total_shared_signal=o$total_signal_volume,
      total_shared_pct=100*o$total_signal_volume/m,
      unique_shared_signal=o$unique_signal_volume,redundant_shared_signal=o$redundant_signal_volume,
      idiosyncratic_volume=o$idiosyncratic_volume,primary_boundary=o$K_roots,
      extended_boundary=o$K_extended,K_log_gap=o$K_log_gap,K_triple_tap=o$K_triple_tap,
      triple_tap_multiplier=o$triple_tap_multiplier,top_of_bulk=o$top_of_bulk,
      AIG=o$AIG,FSIG_1=o$FSIG_final[1],seconds=o$seconds)
  }))
  write.csv(ablation,file.path(out,"ABLATION-RESULTS.csv"),row.names=FALSE)
  write.csv(data.frame(measure=targets$fd_measurements,FD_bins=targets$fd_counts),
    file.path(out,"FD-TARGET-SAMPLE.csv"),row.names=FALSE)
  dimensions <- data.frame(generating_count=cfg$r,
    generating_root_effective_rank=meta$root_effective_rank,
    mean_absolute_root_correlation=meta$mean_absolute_root_correlation,
    active_expansion_terms=meta$active_terms,expanded_effective_rank=meta$expanded_effective_rank,
    observed_information_effective_volume=results$entropic_default$R_eff,
    automatic_primary_boundary=results$entropic_default$K_roots)
  write.csv(dimensions,file.path(out,"DIMENSION-CONTEXT.csv"),row.names=FALSE)
  write.csv(data.frame(axis=seq_along(meta$root_eigenvalues),eigenvalue=meta$root_eigenvalues),
    file.path(out,"GENERATING-ROOT-SPECTRUM.csv"),row.names=FALSE)
  write.csv(data.frame(axis=seq_along(meta$expanded_eigenvalues),eigenvalue=meta$expanded_eigenvalues),
    file.path(out,"EXPANDED-DESIGN-SPECTRUM.csv"),row.names=FALSE)
  review_path <- file.path(out,"VISUAL-ASSESSMENT.csv")
  if(!file.exists(review_path)) write.csv(data.frame(method=labels,
    visual_primary_boundary=NA_integer_,distinct_elbow=NA_character_,
    reviewer="",display_used="",notes=""), review_path,row.names=FALSE,na="")
  save_figure(file.path(out,"Comparison_Plot"),function() plot_comparison(results,cfg,meta))
  save_figure(file.path(out,"Comparison_Unannotated"),function() plot_comparison(results,cfg,meta,FALSE))
  save_figure(file.path(out,"Comparison_Full_Range"),function() plot_comparison(results,cfg,meta,FALSE,FALSE))
  save_figure(file.path(out,"ES_Plot"),function() plot_scree(results$entropic_default,cfg),12,6)
  save_figure(file.path(out,"ES_Unannotated"),function() plot_scree(results$entropic_default,cfg,FALSE),12,6)
  save_figure(file.path(out,"Ablation_Spectra"),function() plot_ablation(results,cfg),14,5.5)
  writeLines(c("\\begin{table}[htbp]","\\centering\\small",
    "\\caption{Binning ablation on the same generated dataset.}",
    "\\begin{tabular}{lrrrrrr}","\\hline",
    "Target & Bins & $R_{eff}$ & Shared (\\%) & $K_{roots}$ & AIG & $FSIG_1$ \\\\","\\hline",
    vapply(seq_len(nrow(ablation)),function(i) sprintf(
      "%s & %d & %.2f & %.2f & %d & %.2f & %.2f \\\\",
      c("Cube root","Square root","FD-derived")[i],ablation$bin_target[i],ablation$R_eff[i],
      ablation$total_shared_pct[i],ablation$primary_boundary[i],ablation$AIG[i],ablation$FSIG_1[i]),character(1)),
    "\\hline","\\end{tabular}","\\end{table}"), file.path(out,"ABLATION-TABLE.tex"))
  invisible(list(comparison=comparison,ablation=ablation,dimensions=dimensions))
}

check_reference <- function(cfg, seed, meta, es, out) {
  expected <- list(r=20L,entanglement="orthogonal",repetition=1L,n_rows=5000L,
    m_proxies=10000L,cran_version="1.0.1",seed_base=19862026L,continuous_ratio=.5,
    max_interaction_order=5L,max_polynomial_order=4L,structural_snr=4,continuous_snr=4,
    binary_error_rate=.0528,pure_static_ratio=.05,pure_root_ratio=0,term_fraction=.1)
  applicable <- all(vapply(names(expected),function(n) isTRUE(all.equal(cfg[[n]],expected[[n]])),logical(1)))
  if(!applicable) return(invisible(NULL))
  reference <- c(root_effective_rank=20,active_terms=21759,expanded_effective_rank=1998.65678682651,
    K_roots=20,K_extended=173,R_eff=9974.1144020771,AIG=11.2725828048136,
    total_signal_volume=225.451656096271,m_retained=10000)
  observed <- c(root_effective_rank=meta$root_effective_rank,active_terms=meta$active_terms,
    expanded_effective_rank=meta$expanded_effective_rank,K_roots=es$K_roots,K_extended=es$K_extended,
    R_eff=es$R_eff,AIG=es$AIG,total_signal_volume=es$total_signal_volume,
    m_retained=length(es$retained_features))
  tolerance <- 1e-7*pmax(1,abs(reference))
  matched <- is.finite(observed) & abs(observed-reference)<=tolerance
  df <- data.frame(metric=names(reference),saved_iterative_run=unname(reference),
    reproduced=unname(observed),tolerance=unname(tolerance),matches=unname(matched))
  write.csv(df,file.path(out,"REFERENCE-RUN-CHECK.csv"),row.names=FALSE)
  stage_message(if(all(matched)) "Default run agrees with the saved iterative reference diagnostics."
    else "Reference differences found. Inspect REFERENCE-RUN-CHECK.csv before reporting reproduction.")
  invisible(df)
}

run_single_study <- function(cfg) {
  stopifnot(cfg$r>=2L,cfg$n_rows>=4L,cfg$m_proxies>=4L,cfg$repetition>=1L,
    cfg$entanglement %in% c("orthogonal","low","high"),
    cfg$continuous_ratio>0,cfg$continuous_ratio<1,
    cfg$pure_static_ratio+cfg$pure_root_ratio<1)
  dir.create(cfg$output_dir,recursive=TRUE,showWarnings=FALSE)
  base <- normalizePath(cfg$output_dir,winslash="/",mustWork=TRUE)
  lib <- file.path(base,"_R-library",paste0(R.version$major,".",strsplit(R.version$minor,"\\.")[[1]][1]))
  ensure_packages(lib,cfg$cran_version)
  detected <- parallel::detectCores(logical=TRUE)
  if(!is.finite(detected)) detected <- 3L
  cores <- if(is.null(cfg$num_cores)) max(1L,detected-2L) else max(1L,as.integer(cfg$num_cores))
  seed <- as.integer(cfg$seed_base+(cfg$repetition-1L)*100000L+(cfg$r-20L)*100L+
    c(orthogonal=0L,low=1L,high=3L)[[cfg$entanglement]])
  scientific <- cfg[setdiff(names(cfg),c("output_dir","num_cores"))]
  dependencies <- c("Entropic.Scree","MASS","Matrix","data.table","stringr","clusterGeneration","ggplot2","patchwork")
  versions <- setNames(vapply(dependencies,function(x) as.character(utils::packageVersion(x)),character(1)),dependencies)
  signature <- fingerprint(list(settings=scientific,versions=versions,R=as.character(getRversion()),
    generator=body(generate_true_mixed_proxies),error=body(apply_measurement_error),
    generation=body(generate_dataset),pca=body(fit_pca),kernel=body(fit_kernel_pca),
    bins=body(choose_bin_targets),entropic=body(fit_entropic)))
  out <- file.path(base,sprintf("r%d-%s-seed%d-%s",cfg$r,cfg$entanglement,seed,substr(signature,1,12)))
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  stage_message("Output: ",out)
  stage_message("One dataset; r=",cfg$r,"; N=",cfg$n_rows,"; m=",cfg$m_proxies,
    "; seed=",seed,"; Entropic.Scree ",cfg$cran_version,"; cores=",cores)
  writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
  save_rds(list(settings=cfg,seed=seed,signature=signature,package_versions=versions,
    detector_formals=formals(Entropic.Scree::Entropic.Scree),cores=cores),file.path(out,"SETTINGS.rds"))
  writeLines(capture.output(dput(list(settings=cfg,seed=seed,package_versions=versions))),file.path(out,"SETTINGS.txt"))
  data_path <- file.path(out,"observed_data.TEMPORARY.rds")
  meta_path <- file.path(out,"GENERATION-SUMMARY.rds")
  stage_keys <- c("standard","rank","kernel","entropic_default","entropic_dense","entropic_fd")
  completed <- read_checkpoint(file.path(out,"COMPLETE.rds"))
  if(!is.null(completed) && all(file.exists(file.path(out,names(completed$checksums)))) &&
      identical(unname(tools::md5sum(file.path(out,names(completed$checksums)))),unname(completed$checksums))) {
    if(file.exists(data_path)) file.remove(data_path)
    stage_message("All stages already complete. No simulation or analysis repeated.")
    return(invisible(out))
  }
  observed <- NULL; meta <- read_checkpoint(meta_path)
  obtain_data <- function() {
    if(!is.null(observed)) return(observed)
    if(!is.null(meta)) observed <<- read_checkpoint(data_path)
    if(is.null(observed)) {
      stage_message("Generating observations and nonlinear-design diagnostics. This stage can take time.")
      generated <- generate_dataset(cfg,seed)
      observed <<- generated$observed; meta <<- generated$metadata
      save_rds(observed,data_path,compress=FALSE)
      save_rds(meta,meta_path)
      rm(generated); gc(FALSE)
      stage_message("Dataset generated and saved as a temporary restart checkpoint.")
    } else stage_message("Loaded the generated-data checkpoint; generation is not repeated.")
    observed
  }
  if(is.null(meta)) invisible(obtain_data())
  targets <- read_checkpoint(file.path(out,"BIN-TARGETS.rds"))
  if(is.null(targets)) {
    targets <- choose_bin_targets(obtain_data())
    save_rds(targets,file.path(out,"BIN-TARGETS.rds"))
  }
  stage_message("Binning targets: ",targets$default,", ",targets$dense,", ",targets$fd," (FD-derived).")
  results <- setNames(vector("list",length(stage_keys)),stage_keys)
  for(i in seq_along(stage_keys)) {
    key <- stage_keys[i]; path <- file.path(out,paste0(key,".rds"))
    result <- read_checkpoint(path)
    if(is.null(result)) {
      stage_message("Stage ",i,"/",length(stage_keys),": ",key,". R remains busy until this calculation finishes.")
      started <- proc.time()[["elapsed"]]
      x <- obtain_data()
      result <- switch(key,
        standard=fit_pca(x),rank=fit_pca(x,TRUE),kernel=fit_kernel_pca(x),
        entropic_default=fit_entropic(x,cfg,seed,cores,NULL,file.path(out,"Native_Scree_Default.pdf")),
        entropic_dense=fit_entropic(x,cfg,seed,cores,targets$dense,file.path(out,"Native_Scree_Dense.pdf")),
        entropic_fd=fit_entropic(x,cfg,seed,cores,targets$fd,file.path(out,"Native_Scree_FD.pdf")))
      result$seconds <- proc.time()[["elapsed"]]-started
      save_rds(result,path)
      stage_message("Saved ",key," (",round(result$seconds,1)," seconds).")
      rm(x); gc(FALSE)
    } else stage_message("Stage ",i,"/",length(stage_keys),": reusing completed ",key,".")
    results[[key]] <- result
  }
  stage_message("Saving comparison tables, effective-rank context and figures.")
  tables <- make_outputs(results,meta,targets,cfg,seed,out)
  check_reference(cfg,seed,meta,results$entropic_default,out)
  writeLines(c(
    "All analyses used the same generated observations. Completed method checkpoints are reused on restart.",
    "Entropic Scree: unmodified CRAN 1.0.1; preparation and both detectors use published defaults.",
    "The square-root and FD-derived ablations change only num_bins; they reuse the same observations.",
    "PCA: centered, unit-variance columns and full prcomp spectrum. Rank PCA: average ranks for ties, then the same PCA.",
    "RBF kernel: gamma=1/m on standardized observations; full double-centered kernel eigenvalues, clipped at zero.",
    "Threshold counts are separate from visual elbows. Baselines have no invented automatic elbow or assigned visual rank.",
    "VISUAL-ASSESSMENT.csv is a blank optional review record and is never overwritten on restart.",
    "Comparison_Unannotated/Full_Range and ES_Unannotated omit generating and selected boundary markers.",
    "Comparison plots: logarithmic eigenvalue index, linear eigenvalues. Full_Range preserves the leading eigenvalue.",
    "ES plots: linear index and log eigenvalues, with the original macro/detail windows. Ablation: first 75 indices, linear index and log eigenvalues.",
    "Ablation display upper limit=1.025 times the second eigenvalue; the leading eigenvalue can be outside the panel.",
    "root_effective_rank: population generating-correlation entropy rank, retaining eigenvalues >1e-9.",
    "expanded_effective_rank: active standardized nonlinear-design covariance entropy rank (K_eff).",
    "R_eff: observed working-information-spectrum effective volume. These are different quantities.",
    "The observed-data checkpoint is deleted after completion; no true-data/centered-NMI matrices are saved.",
    "Reproduction uses the same generator/RNG sequence; numerical differences across R/dependency/BLAS versions remain possible.",
    "The initial default is the completed study's orthogonal r=20, repetition=1, seed=19862026, N=5000, m=10000.",
    "REFERENCE-RUN-CHECK.csv compares its saved diagnostics when those default scientific settings are used."
  ),file.path(out,"READ-ME.txt"))
  tracked <- c(paste0(stage_keys,".rds"),"GENERATION-SUMMARY.rds","BIN-TARGETS.rds",
    "COMPARISON-SUMMARY.csv","ABLATION-RESULTS.csv","DIMENSION-CONTEXT.csv","ABLATION-TABLE.tex",
    "FD-TARGET-SAMPLE.csv","GENERATING-ROOT-SPECTRUM.csv","EXPANDED-DESIGN-SPECTRUM.csv",
    paste0("SPECTRUM-",stage_keys,".csv"),
    paste0("FSIG-",stage_keys[4:6],".csv"),paste0("RETAINED-",stage_keys[4:6],".csv"),
    as.vector(outer(c("Comparison_Plot","Comparison_Unannotated","Comparison_Full_Range","ES_Plot","ES_Unannotated","Ablation_Spectra"),
      c(".png",".pdf"),paste0)))
  if(!all(file.exists(file.path(out,tracked)))) stop("An output is missing; temporary data retained for restart.")
  checksums <- setNames(unname(tools::md5sum(file.path(out,tracked))),tracked)
  save_rds(list(signature=signature,checksums=checksums,completed=Sys.time()),file.path(out,"COMPLETE.rds"))
  if(file.exists(data_path) && !file.remove(data_path)) warning("Could not remove temporary data: ",data_path)
  stage_message("Complete. Temporary observations removed. Results: ",out)
  print(tables$comparison, row.names=FALSE)
  print(tables$ablation, row.names=FALSE)
  invisible(out)
}

single_study_output <- run_single_study(CONFIG)


# ============================================================================
# LOW- AND HIGH-ENTANGLEMENT FOUR-METHOD COMPARISONS
# Run this section alone after the main study in the same R session.
# It uses CONFIG and the functions above; it does not rerun the orthogonal study.
# Each context has its own dataset, four method checkpoints and output folder.
# All four methods within a context use exactly the same observed dataset.
# Binning ablation remains in the original study; it is not repeated here.
# ============================================================================

run_entanglement_comparison <- function(cfg, context) {
  stopifnot(length(context)==1L, context %in% c("low","high"))
  cfg$entanglement <- context
  dir.create(cfg$output_dir,recursive=TRUE,showWarnings=FALSE)
  base <- normalizePath(cfg$output_dir,winslash="/",mustWork=TRUE)
  lib <- file.path(base,"_R-library",paste0(R.version$major,".",strsplit(R.version$minor,"\\.")[[1]][1]))
  ensure_packages(lib,cfg$cran_version)
  detected <- parallel::detectCores(logical=TRUE)
  if(!is.finite(detected)) detected <- 3L
  cores <- if(is.null(cfg$num_cores)) max(1L,detected-2L) else max(1L,as.integer(cfg$num_cores))
  seed <- as.integer(cfg$seed_base+(cfg$repetition-1L)*100000L+(cfg$r-20L)*100L+
    c(orthogonal=0L,low=1L,high=3L)[[cfg$entanglement]])
  scientific <- cfg[setdiff(names(cfg),c("output_dir","num_cores"))]
  dependencies <- c("Entropic.Scree","MASS","Matrix","data.table","stringr","clusterGeneration","ggplot2","patchwork")
  versions <- setNames(vapply(dependencies,function(x) as.character(utils::packageVersion(x)),character(1)),dependencies)
  function_text <- function(f) paste(deparse(body(f),width.cutoff=500L),collapse="\n")
  signature <- fingerprint(list(settings=scientific,versions=versions,R=as.character(getRversion()),
    generator=function_text(generate_true_mixed_proxies),error=function_text(apply_measurement_error),
    generation=function_text(generate_dataset),pca=function_text(fit_pca),kernel=function_text(fit_kernel_pca),
    bins=function_text(choose_bin_targets),entropic=function_text(fit_entropic)))
  out <- file.path(base,sprintf("r%d-%s-seed%d-%s",cfg$r,context,seed,substr(signature,1,12)))
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  stage_message("Root entanglement: ",context,"; seed=",seed,"; N=",cfg$n_rows,"; m=",cfg$m_proxies)
  stage_message("Output: ",out)
  data_path <- file.path(out,"observed_data.TEMPORARY.rds")
  meta_path <- file.path(out,"GENERATION-SUMMARY.rds")
  completion_path <- file.path(out,"COMPARISON-COMPLETE.rds")
  stage_keys <- c("standard","rank","kernel","entropic_default")
  plot_names <- c("Comparison_Plot","Comparison_Unannotated","Comparison_Full_Range")
  tracked <- c(paste0(stage_keys,".rds"),"GENERATION-SUMMARY.rds",
    "COMPARISON-SUMMARY.csv","DIMENSION-CONTEXT.csv",
    paste0("SPECTRUM-",stage_keys,".csv"),
    as.vector(outer(plot_names,c(".png",".pdf"),paste0)))
  completed <- read_checkpoint(completion_path)
  copy_main_figures <- function() {
    label <- if(context=="low") "Low" else "High"
    paths <- file.path(base,paste0(label,"_Entanglement_Plot",c(".png",".pdf")))
    for(i in seq_along(paths)) {
      if(!file.copy(file.path(out,paste0("Comparison_Plot",c(".png",".pdf")[i])),
          paths[i],overwrite=TRUE)) stop("Cannot save comparison copy: ",paths[i])
    }
    invisible(paths)
  }
  if(!is.null(completed) && identical(completed$signature,signature) &&
      all(tracked %in% names(completed$checksums)) &&
      all(file.exists(file.path(out,tracked))) &&
      identical(unname(tools::md5sum(file.path(out,tracked))),
                unname(completed$checksums[tracked]))) {
    copy_main_figures()
    if(file.exists(data_path)) file.remove(data_path)
    stage_message(context," comparison already complete. No data generation or model fitting repeated.")
    return(invisible(out))
  }
  writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
  save_rds(list(settings=cfg,seed=seed,signature=signature,package_versions=versions,
    detector_formals=formals(Entropic.Scree::Entropic.Scree),cores=cores),file.path(out,"SETTINGS.rds"))
  writeLines(capture.output(dput(list(settings=cfg,seed=seed,package_versions=versions))),
    file.path(out,"SETTINGS.txt"))
  observed <- NULL
  meta <- read_checkpoint(meta_path)
  obtain_data <- function() {
    if(!is.null(observed)) return(observed)
    if(!is.null(meta)) observed <<- read_checkpoint(data_path)
    if(is.null(observed)) {
      stage_message("Generating the ",context,"-entanglement dataset and design diagnostics.")
      generated <- generate_dataset(cfg,seed)
      observed <<- generated$observed
      meta <<- generated$metadata
      save_rds(observed,data_path,compress=FALSE)
      save_rds(meta,meta_path)
      rm(generated); gc(FALSE)
    } else stage_message("Reusing the generated-data checkpoint.")
    observed
  }
  if(is.null(meta)) invisible(obtain_data())
  results <- setNames(vector("list",length(stage_keys)),stage_keys)
  for(i in seq_along(stage_keys)) {
    key <- stage_keys[i]
    path <- file.path(out,paste0(key,".rds"))
    result <- read_checkpoint(path)
    if(is.null(result)) {
      stage_message(context," stage ",i,"/4: ",key,". Calculation in progress.")
      started <- proc.time()[["elapsed"]]
      x <- obtain_data()
      result <- switch(key,standard=fit_pca(x),rank=fit_pca(x,TRUE),kernel=fit_kernel_pca(x),
        entropic_default=fit_entropic(x,cfg,seed,cores,NULL,file.path(out,"Native_Scree_Default.pdf")))
      result$seconds <- proc.time()[["elapsed"]]-started
      save_rds(result,path)
      stage_message("Saved ",key," (",round(result$seconds,1)," seconds).")
      rm(x); gc(FALSE)
    } else stage_message(context," stage ",i,"/4: reusing completed ",key,".")
    results[[key]] <- result
  }
  labels <- c("Standard PCA","Spearman rank PCA","RBF kernel PCA","Entropic Scree")
  comparison <- do.call(rbind,lapply(seq_along(stage_keys),function(i) {
    o <- results[[stage_keys[i]]]; es <- i==4L
    data.frame(method=labels[i],entanglement=context,generating_count=cfg$r,seed=seed,N=cfg$n_rows,
      m_input=cfg$m_proxies,m_analyzed=if(es) length(o$retained_features) else cfg$m_proxies,
      rank_ceiling=if(es) length(o$retained_features)-1L else o$rank_ceiling,
      automatic_primary_boundary=if(es) o$K_roots else NA_real_,
      threshold_count=if(es) sum(o$eigenvalues>mean(o$eigenvalues)) else o$threshold_count,
      threshold=if(es) mean(o$eigenvalues) else o$threshold,
      threshold_rule=if(es) "Working eigenvalue > mean (separate diagnostic)" else o$threshold_rule,
      seconds=o$seconds,stringsAsFactors=FALSE)
  }))
  write.csv(comparison,file.path(out,"COMPARISON-SUMMARY.csv"),row.names=FALSE,na="")
  for(key in stage_keys) write.csv(data.frame(axis=seq_along(results[[key]]$eigenvalues),
    eigenvalue=results[[key]]$eigenvalues),file.path(out,paste0("SPECTRUM-",key,".csv")),row.names=FALSE)
  dimensions <- data.frame(entanglement=context,generating_count=cfg$r,
    generating_root_effective_rank=meta$root_effective_rank,
    mean_absolute_root_correlation=meta$mean_absolute_root_correlation,
    active_expansion_terms=meta$active_terms,expanded_effective_rank=meta$expanded_effective_rank,
    observed_information_effective_volume=results$entropic_default$R_eff,
    automatic_primary_boundary=results$entropic_default$K_roots)
  write.csv(dimensions,file.path(out,"DIMENSION-CONTEXT.csv"),row.names=FALSE)
  review_path <- file.path(out,"VISUAL-ASSESSMENT.csv")
  if(!file.exists(review_path)) write.csv(data.frame(method=labels,
    visual_primary_boundary=NA_integer_,distinct_elbow=NA_character_,
    reviewer="",display_used="",notes=""),review_path,row.names=FALSE,na="")
  stage_message("Saving ",context,"-entanglement four-panel comparison plots.")
  save_figure(file.path(out,"Comparison_Plot"),function() plot_comparison(results,cfg,meta))
  save_figure(file.path(out,"Comparison_Unannotated"),function() plot_comparison(results,cfg,meta,FALSE))
  save_figure(file.path(out,"Comparison_Full_Range"),function() plot_comparison(results,cfg,meta,FALSE,FALSE))
  copy_main_figures()
  writeLines(c(
    paste("Root entanglement:",context,"; seed:",seed),
    "Four methods share one generated dataset; no binning ablation is run in this section.",
    "Entropic Scree uses the unmodified CRAN package 1.0.1 and its default detectors.",
    "Comparison_Plot uses the original four-panel styling, logarithmic index and linear eigenvalues.",
    "Its y-axis cap is 1.1 times the second eigenvalue; the leading value can be outside the panel.",
    "The orange mean-spectrum/Kaiser lines are threshold counts, not fitted visual elbows.",
    "The automatic Entropic Scree primary boundary is recorded separately in COMPARISON-SUMMARY.csv.",
    "Unannotated and full-range plots are also saved. VISUAL-ASSESSMENT.csv is never overwritten.",
    "Completed method checkpoints are reused on restart. Temporary observations are removed on completion."
  ),file.path(out,"COMPARISON-READ-ME.txt"))
  if(!all(file.exists(file.path(out,tracked)))) stop("A comparison output is missing; temporary data retained.")
  checksums <- setNames(unname(tools::md5sum(file.path(out,tracked))),tracked)
  save_rds(list(signature=signature,checksums=checksums,completed=Sys.time()),completion_path)
  if(file.exists(data_path) && !file.remove(data_path)) warning("Could not remove temporary data: ",data_path)
  stage_message(context," comparison complete: ",out)
  print(comparison,row.names=FALSE)
  invisible(out)
}

run_low_high_comparisons <- function(cfg) {
  needed <- c("ensure_packages","generate_dataset","generate_true_mixed_proxies",
    "apply_measurement_error","fit_pca","fit_kernel_pca","fit_entropic",
    "choose_bin_targets","fingerprint","save_rds","read_checkpoint","stage_message",
    "plot_comparison","save_figure")
  missing <- needed[!vapply(needed,function(n) exists(n,mode="function",inherits=TRUE),logical(1))]
  if(length(missing)) stop("Load the main study's setup and functions first. Missing: ",paste(missing,collapse=", "))
  contexts <- c("low","high")
  index <- data.frame(entanglement=contexts,status="pending",output_directory="",error="")
  dir.create(cfg$output_dir,recursive=TRUE,showWarnings=FALSE)
  index_path <- file.path(cfg$output_dir,"LOW-HIGH-COMPARISON-INDEX.csv")
  for(i in seq_along(contexts)) {
    index$status[i] <- "running"
    write.csv(index,index_path,row.names=FALSE)
    result <- tryCatch(run_entanglement_comparison(cfg,contexts[i]),error=function(e) e)
    if(inherits(result,"error")) {
      index$status[i] <- "failed"
      index$error[i] <- conditionMessage(result)
      stage_message(contexts[i]," comparison stopped: ",index$error[i])
    } else {
      index$status[i] <- "completed"
      index$output_directory[i] <- result
    }
    write.csv(index,index_path,row.names=FALSE)
    gc(FALSE)
  }
  if(any(index$status=="failed")) stop("One or more comparisons failed. See ",index_path,
    ". Rerun this section to resume; completed methods and contexts will be reused.",call.=FALSE)
  stage_message("Low and high comparisons complete. Plot copies: ",normalizePath(cfg$output_dir,winslash="/"))
  invisible(index)
}

plot_saved_ablation <- function(results_dir, output_stem, zoom_last=75L) {
  results_dir <- normalizePath(results_dir, winslash="/", mustWork=TRUE)
  keys <- c("entropic_default", "entropic_dense", "entropic_fd")
  titles <- c("Conservative default", "Square-root target", "Freedman-Diaconis target")
  tab <- read.csv(file.path(results_dir, "ABLATION-RESULTS.csv"))
  stopifnot(!anyDuplicated(tab$rule), all(keys %in% tab$rule))
  tab <- tab[match(keys, tab$rule), , drop=FALSE]
  settings <- readRDS(file.path(results_dir, "SETTINGS.rds"))$settings
  spectra <- lapply(keys, function(key) {
    d <- read.csv(file.path(results_dir, paste0("SPECTRUM-", key, ".csv")))
    stopifnot(identical(as.integer(d$axis), seq_len(nrow(d))),
              all(is.finite(d$eigenvalue)), all(d$eigenvalue > 0),
              all(diff(d$eigenvalue) <= 1e-10))
    d$eigenvalue
  })
  stopifnot(all(lengths(spectra) == tab$m_retained))
  # The saved spectrum is floored at 1e-9. Its mean is NOT necessarily the
  # original centered-matrix trace / m, because negative eigenvalues were lost.
  # Label the available quantity accurately rather than reconstructing a trace.
  working_means <- vapply(spectra, mean, numeric(1))
  floor_counts <- vapply(spectra, function(e) sum(e <= 1e-9), integer(1))
  colors <- c("#555555", "#D55E00", "#7B3294", "#228833")
  types <- c(3, 2, 4, 5)
  full_ids <- lapply(spectra, function(e) which(seq_along(e)>1L & e>1e-9))
  stopifnot(all(lengths(full_ids)>=2L))
  all_y <- range(unlist(Map(function(e,ids) e[ids],spectra,full_ids)))
  draw <- function() {
    old <- par(no.readonly=TRUE)
    on.exit(par(old))
    layout(rbind(1:3,4:6,c(7,7,7)),heights=c(1,1,.23))
    par(mar=c(4.8,4.6,4.4,1.0), oma=c(0,0,3.4,0),
        mgp=c(2.7,.8,0), las=1, cex=.95)
    for(view in c("full", "zoom")) for(j in seq_along(keys)) {
      e <- spectra[[j]]
      ids <- if(view=="full") full_ids[[j]] else seq_len(min(zoom_last,length(e)))
      if(view=="full") {
        yl <- all_y * c(.99,1.01)
      } else {
        # Preserve the original early-spectrum view: crop the dominant first
        # eigenvalue vertically, but include every index in the selected window.
        yl <- c(min(e[ids])*.997, e[2]*1.025)
      }
      plot(ids,e[ids],type="n",log=if(view=="full") "xy" else "y",
           ylim=yl, xlab=if(view=="full") "Eigenvalue index (log scale)" else "Eigenvalue index",
           ylab="Working eigenvalue (log scale)", bty="l", cex.axis=.82)
      abline(h=working_means[j],col=colors[4],lty=types[4],lwd=1.2)
      boundaries <- c(settings$r,tab$primary_boundary[j],tab$extended_boundary[j])
      for(k in seq_along(boundaries)) if(boundaries[k] >= min(ids) && boundaries[k] <= max(ids)) {
        abline(v=boundaries[k],col=colors[k],lty=types[k],lwd=if(k==1) 2.8 else 1.4)
      }
      lines(ids,e[ids],col="#187CB8",lwd=1.5)
      if(view=="zoom") points(ids,e[ids],col="#187CB8",pch=16,cex=.4)
      title(main=paste0(LETTERS[(if(view=="full") 0 else 3)+j], ". ",titles[j],
                        "\n",tab$bin_target[j]," bins"),line=2.0,cex.main=.98)
      mtext(if(view=="full") sprintf("Broad view | primary %d; tail %d",boundaries[2],boundaries[3])
            else sprintf("First %d axes | leading value cropped",length(ids)),
            side=3,line=.45,cex=.76)
      if(view=="full") {
        mtext(sprintf("Axis 1 and %d numerical-floor value(s) omitted",floor_counts[j]),side=1,line=4.0,cex=.66)
      }
    }
    mtext("Discretization ablation: broad spectra and early structural breaks",outer=TRUE,
          side=3,line=1.7,font=2,cex=1.25)
    mtext(sprintf("Same observations: N = %s; generated measures = %s; generating roots = %d",
                  format(settings$n_rows,big.mark=","),format(settings$m_proxies,big.mark=","),settings$r),
          outer=TRUE,side=3,line=.25,cex=.87)
    par(mar=c(0,0,0,0))
    plot.new()
    legend("top",legend=c("Generating count", "Automatic primary boundary", "Extended signal tail",
                           "Mean of saved working spectrum"),col=colors,lty=types,lwd=c(2.8,1.4,1.4,1.2),
           horiz=TRUE,bty="n",cex=.83)
    text(.5,.34,"Broad views omit axis 1 and numerical-floor values; original indices and shared top-row vertical scales are retained.",cex=.76)
    text(.5,.06,"Saved-spectrum means are reference lines, not recovered pre-clipping mean traces. Boundaries outside zoom windows are omitted.",cex=.74)
  }
  png_path <- file.path(results_dir,paste0(output_stem,".png"))
  pdf_path <- file.path(results_dir,paste0(output_stem,".pdf"))
  grDevices::png(png_path,width=4500,height=2850,res=300)
  tryCatch(draw(),finally=grDevices::dev.off())
  grDevices::pdf(pdf_path,width=15,height=9.5,useDingbats=FALSE)
  tryCatch(draw(),finally=grDevices::dev.off())
  audit <- data.frame(rule=keys,bin_target=tab$bin_target,
    primary_boundary=tab$primary_boundary,extended_boundary=tab$extended_boundary,
    mean_saved_working_spectrum=working_means,values_at_numerical_floor=floor_counts,
    broad_first_index=vapply(full_ids,min,integer(1)),
    broad_last_index=vapply(full_ids,max,integer(1)))
  write.csv(audit,file.path(results_dir,paste0(output_stem,"-REFERENCES.csv")),row.names=FALSE)
  print(audit,row.names=FALSE)
  message("Saved PNG: ",png_path,"\nSaved PDF: ",pdf_path)
  invisible(audit)
}

# Replot the saved main-study spectra before running the additional contexts.
plot_saved_ablation(single_study_output, "Ablation_Broad_And_Zoom", 75L)
run_low_high_comparisons(CONFIG)
