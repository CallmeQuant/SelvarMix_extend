# Helper functions from MMCEM2
PartitionData <- function(data) {
  d <- ncol(data)
  idx <- seq(1:nrow(data))
  is_comp <- stats::complete.cases(data)
  is_incomp <- !is_comp
  
  # Complete cases
  data_comp <- data[is_comp, , drop = FALSE]
  idx_comp <- idx[is_comp]
  
  # Incomplete cases
  data_incomp <- data[is_incomp, , drop = FALSE]
  idx_incomp <- idx[is_incomp]
  
  # Empty cases
  is_empty <- apply(data_incomp, 1, function(x){
    sum(is.na(x)) == d
  })
  data_empty <- data_incomp[is_empty, , drop = FALSE]
  idx_empty <- idx_incomp[is_empty]
  
  # Remove empty cases
  data_incomp <- data_incomp[!is_empty, , drop = FALSE]
  idx_incomp <- idx_incomp[!is_empty]
  
  # Output
  out <- list()
  out$orig_row_names <- rownames(data)
  out$orig_col_names <- colnames(data)
  
  out$n_row <- nrow(data)
  out$n_col <- ncol(data)
  
  out$n0 <- nrow(data_comp)
  out$n1 <- nrow(data_incomp)
  out$n2 <- nrow(data_empty)
  
  out$data_comp <- data_comp
  out$data_incomp <- data_incomp
  out$data_empty <- data_empty
  
  out$idx_comp <- idx_comp
  out$idx_incomp <- idx_incomp
  out$idx_empty <- idx_empty
  out$init_order <- c(idx_comp, idx_incomp, idx_empty)
  return(out)
}

CalcCondDist <- function(y, idx_a, idx_b, mu, sigma) {
  # Split outcome.
  y_b <- y[idx_b]
  
  # Split mean.
  mu_a <- mu[idx_a]
  mu_b <- mu[idx_b]
  
  # Split covariance.
  sigma_aa <- sigma[idx_a, idx_a, drop = FALSE]
  sigma_ab <- sigma[idx_a, idx_b, drop = FALSE]
  sigma_bb <- sigma[idx_b, idx_b, drop = FALSE]
  
  # Calculate conditional mean and covariance.
  mu_cond <- mu_a + sigma_ab %*% solve(sigma_bb, y_b - mu_b)
  sigma_cond <- sigma_aa - sigma_ab %*% solve(sigma_bb, t(sigma_ab))
  
  out <- list(
    mu = mu_cond,
    sigma = sigma_cond
  )
  return(out)
}

EvalDensIncompObs <- function(y, means, covs, pi) {
  k <- length(pi)
  is_obs <- !is.na(y)
  
  obs_ele <- y[is_obs]
  obs_dens <- numeric(k)
  
  for (j in 1:k) {
    obs_mean <- as.numeric(means[[j]][is_obs])
    obs_cov <- as.matrix(covs[[j]][is_obs, is_obs, drop = FALSE])
    
    # Check if covariance matrix positive definite
    if (any(is.na(obs_cov)) || any(diag(obs_cov) <= 0) || 
        (nrow(obs_cov) > 1 && det(obs_cov) <= 0)) {
      warning("Observed covariance matrix may be ill-conditioned")
      obs_dens[j] <- 0
    } else {
      obs_dens[j] <- mvnfast::dmvn(X = obs_ele, mu = obs_mean, sigma = obs_cov) * pi[j]
    }
  }
  
  return(obs_dens)
}

Responsibility <- function(split_data, means, covs, pi) {
  n0 <- split_data$n0
  n1 <- split_data$n1
  
  k <- length(pi)

  out <- list()
  out$k <- k

  # Density evaluation for complete obss.
  if (n0 > 0) {
    X_matrix <- as.matrix(split_data$data_comp)
    
    dens_eval0 <- matrix(0, nrow = n0, ncol = k)
    for (j in 1:k) {
      mu_vector <- as.numeric(means[[j]])
      sigma_matrix <- as.matrix(covs[[j]])
      
      dens_eval0[, j] <- mvnfast::dmvn(
        X = X_matrix,
        mu = mu_vector,
        sigma = sigma_matrix) * pi[j]
    }

    # Normalize by row to get responsibilities
    row_sums <- rowSums(dens_eval0)
    gamma0 <- dens_eval0 / row_sums
    
    colnames(dens_eval0) <- colnames(gamma0) <- paste0("k", 1:k)
    rownames(dens_eval0) <- rownames(gamma0) <- 1:n0
    out$dens_eval0 <- dens_eval0
    out$gamma0 <- gamma0
  }

  # Evaluation for Incomplete obss.
  if (n1 > 0) {
    data_incomp <- as.matrix(split_data$data_incomp)
    dens_eval1 <- matrix(0, nrow = n1, ncol = k)
    
    # Process each incomplete obs individually
    for (i in 1:n1) {
      x <- data_incomp[i, ]
      dens <- EvalDensIncompObs(x, means, covs, pi)
      dens_eval1[i, ] <- dens
    }

    # Normalize by row to get responsibilities
    row_sums <- rowSums(dens_eval1)
    gamma1 <- dens_eval1 / row_sums
    
    # Format results
    colnames(dens_eval1) <- colnames(gamma1) <- paste0("k", 1:k)
    rownames(dens_eval1) <- rownames(gamma1) <- 1:n1
    out$dens_eval1 <- dens_eval1
    out$gamma1 <- gamma1
  }
  
  return(out)
}

# Convert mclust parameters to list format for MMCEM2 functions
convert_mclust_to_lists <- function(parameters) {
  G <- length(parameters$pro)
  
  means <- lapply(1:G, function(g) as.numeric(parameters$mean[, g]))
  
  covs <- lapply(1:G, function(g) {
    get_component_covariance(parameters$variance, g)
  })
  
  pi <- parameters$pro
  
  return(list(means = means, covs = covs, pi = pi))
}

EM_impute <- function(data, 
                     G = 3, 
                     modelName = "VVV",
                     max_iter = 1000,
                     tol = 1e-6,
                     init_method = "hc",
                     method = "usual", 
                     S = 10, 
                     verbose = FALSE) {
  
  # Validate method parameter
  method <- match.arg(method, c("usual", "sampling"))
  
  if (is.data.frame(data)) {
    data <- as.matrix(data)
  }
  
  n <- nrow(data)
  p <- ncol(data)
  
  missing_pattern <- is.na(data)
  has_missing <- any(missing_pattern)
  
  if (!has_missing) {
    warning("No missing values found in data")
    return(list(
      imputed_data = data,
      converged = TRUE,
      iterations = 0,
      loglik = NA,
      parameters = NULL,
      method = method
    ))
  }
  
  if (verbose) cat("Initializing parameters using method:", method, "\n")
  
  # Initialize based on method
  if (method == "usual") {
    # Original initialization
    complete_cases <- complete.cases(data)
    if (sum(complete_cases) < G) {
      stop("Not enough complete cases for initialization. Consider reducing G or using a different approach.")
    }
    
    complete_data <- data[complete_cases, , drop = FALSE]
    init_result <- InitParameter(complete_data, 
                                nbClust = G, 
                                init = init_method)
    
    # Initial imputation with means
    data_imputed <- data
    for (j in 1:p) {
      missing_j <- is.na(data[, j])
      if (any(missing_j)) {
        data_imputed[missing_j, j] <- mean(data[!missing_j, j], na.rm = TRUE)
      }
    }
    
    em_control <- emControl(tol = tol, itmax = max_iter)
    z_current <- resp_to_full_data(init_result, data_imputed, complete_cases, G)
    
  } else {
    # MMCEM2 initialization
    split_data <- PartitionData(data)
    comp_data <- split_data$data_comp
    
    if(nrow(comp_data) < G) {
      stop("Number of complete cases is less than the number of clusters")
    }
    
    # Initialize with k-means on complete cases
    # init_kmeans <- kmeans(comp_data, centers = G)
    # means <- lapply(1:G, function(g) as.numeric(colMeans(comp_data[init_kmeans$cluster == g,,drop=F])))
    # covs <- lapply(1:G, function(g) cov(comp_data[init_kmeans$cluster == g,,drop=F]))
    # pi <- table(init_kmeans$cluster)/nrow(comp_data)
    
    init_results <- InitParameter(data = comp_data, 
                             nbClust = G, 
                             init = init_method)
    
    means <- lapply(1:G, function(g) {
      as.numeric(init_results$Mu[, g])
    })
    
    covs <- lapply(1:G, function(g) {
      init_results$SigmaCube[, , g]
    })
    
    pi <- init_results$prop
    
    
    # Initial imputation
    data_imputed <- data
    
    if(split_data$n1 > 0) {
      incomp_data <- split_data$data_incomp
      for(i in 1:nrow(incomp_data)) {
        y <- incomp_data[i,]
        idx_miss <- which(is.na(y))
        idx_obs <- which(!is.na(y))
        
        if(length(idx_obs) > 0) {
          # Use first component for initial imputation
          cond_dist <- CalcCondDist(y, idx_miss, idx_obs, means[[1]], covs[[1]])
          imputed_values <- mvnfast::rmvn(1, mu = cond_dist$mu, sigma = cond_dist$sigma)
          data_imputed[split_data$idx_incomp[i], idx_miss] <- imputed_values
        } else {
          # If no observed values, use marginal mean
          data_imputed[split_data$idx_incomp[i], idx_miss] <- means[[1]][idx_miss]
        }
      }
    }
  }
  
  if (verbose) cat("Starting EM algorithm...\n")
  
  loglik_old <- -Inf
  converged <- FALSE
  iter <- 0
  
  # Main EM loop
  for (iter in 1:max_iter) {
    if (verbose && iter %% 10 == 0) {
      cat(sprintf("Iteration %d\n", iter))
    }
    
    if (method == "usual") {
      me_result <- me(data = data_imputed, 
                      modelName = modelName,
                      z = z_current,
                      control = em_control)
      
      if (is.null(me_result) || is.null(me_result$parameters)) {
        if (verbose) cat("EM algorithm failed - using previous iteration\n")
        break
      }
      
      # Update responsibilities
      estep_result <- estep(modelName = modelName,
                           data = data_imputed,
                           parameters = me_result$parameters)
      
      z_current <- estep_result$z
      
      # Impute missing values using current parameters 
      data_imputed_new <- impute_missing_values(data, me_result$parameters, z_current, missing_pattern)
      
      # Check convergence
      loglik_new <- me_result$loglik
      
    } else {
      split_data <- PartitionData(data) 
      
      # E-step
      # Calculate responsibilities based on observed data
      # resp$gamma0 for complete cases, resp$gamma1 for incomplete cases
      resp <- Responsibility(split_data, means, covs, pi) 

      # Augment incomplete data with MC samples and prepare corresponding weights
      if (split_data$n1 > 0) { # If there are incomplete observations
        # Get original column names once, if needed for reconstruction matching
        # original_colnames <- colnames(data) 

        # Pre-allocate lists to store augmented data and weights for each original incomplete observation
        all_aug_data_list <- vector("list", split_data$n1)
        all_aug_weights_list <- vector("list", split_data$n1)

        for (obs_idx in 1:split_data$n1) { # For each original incomplete observation
          original_data_row_idx <- split_data$idx_incomp[obs_idx] # Get original row index
          y_original_incomplete_row <- data[original_data_row_idx, , drop = FALSE] # Get the original row with NAs
      
          is_na_y <- is.na(y_original_incomplete_row)
          idx_miss <- which(is_na_y)
          idx_obs <- which(!is_na_y)
          
          # Ensure y_obs_values is a vector, even if only one observed variable
          y_obs_values <- y_original_incomplete_row[1, idx_obs, drop = TRUE] 
      
          # Responsibilities for this specific original incomplete observation
          gamma_i_for_this_obs <- resp$gamma1[obs_idx, , drop = TRUE] 
      
          # Store S imputations & G weight vectors for *this* original incomplete observation
          # Each list element will correspond to imputations for one component g
          observation_specific_aug_data_list <- vector("list", G)
          observation_specific_aug_weights_list <- vector("list", G)
      
          for (g_comp in 1:G) { # For each component g = 1,...,G
            current_imps <- matrix(NA, nrow = S, ncol = length(idx_miss))
      
            if (length(idx_miss) > 0) { # Only impute if there are missing values
                if (length(idx_obs) > 0) { # Condition on observed parts if they exist
                  cond_dist <- CalcCondDist(y_original_incomplete_row, idx_miss, idx_obs, means[[g_comp]], covs[[g_comp]])
                  
                  # Check if conditional covariance is positive definite for sampling
                  is_sigma_cond_pd <- TRUE
                  if(any(is.na(cond_dist$sigma)) || any(diag(cond_dist$sigma) <= 1e-8) ) is_sigma_cond_pd <- FALSE
                  if(is_sigma_cond_pd && nrow(cond_dist$sigma) > 1 && det(cond_dist$sigma) <= 1e-8) is_sigma_cond_pd <- FALSE
                  
                  if(is_sigma_cond_pd) {
                    current_imps <- mvnfast::rmvn(S, mu = cond_dist$mu, sigma = cond_dist$sigma)
                  } else {
                    warning(paste("Conditional covariance for original row", original_data_row_idx, "component", g_comp, 
                                  "is ill-conditioned. Using marginal imputation for S samples."))
                    # Fallback: impute S samples from the marginal distribution of component g_comp for missing parts
                    # Ensure covariance for marginal imputation is PD if subsetting
                    marginal_cov_miss <- as.matrix(covs[[g_comp]][idx_miss, idx_miss, drop=FALSE])
                    if(any(is.na(marginal_cov_miss)) || any(diag(marginal_cov_miss) <= 1e-8) ) {
                         for(s_idx in 1:S) current_imps[s_idx,] <- means[[g_comp]][idx_miss] 
                    } else {
                         current_imps <- mvnfast::rmvn(S, mu = means[[g_comp]][idx_miss], sigma = marginal_cov_miss)
                    }
                  }
                } else { # No observed parts, impute S samples from marginal of component g_comp
                  marginal_cov_miss <- as.matrix(covs[[g_comp]][idx_miss, idx_miss, drop=FALSE])
                  if(any(is.na(marginal_cov_miss)) || any(diag(marginal_cov_miss) <= 1e-8) ) {
                       for(s_idx in 1:S) current_imps[s_idx,] <- means[[g_comp]][idx_miss]
                  } else {
                       current_imps <- mvnfast::rmvn(S, mu = means[[g_comp]][idx_miss], sigma = marginal_cov_miss)
                  }
                }
            }
      
            # Reconstruct S full data vectors for these S samples
            reconstructed_s_samples <- matrix(NA, nrow = S, ncol = p)
            colnames(reconstructed_s_samples) <- colnames(data)
      
            if (length(idx_obs) > 0) {
              reconstructed_s_samples[, idx_obs] <- matrix(rep(y_obs_values, S), nrow = S, byrow = TRUE)
            }
            if (length(idx_miss) > 0) {
              reconstructed_s_samples[, idx_miss] <- current_imps
            }
            observation_specific_aug_data_list[[g_comp]] <- reconstructed_s_samples
      
            # Weights for these S samples: sparse vector (0, ..., z_ig/S, ..., 0)
            # gamma_i_for_this_obs[g_comp] is z_ig for this observation and this component g_comp
            current_s_samples_weights <- matrix(0, nrow = S, ncol = G)
            if (S > 0) { # Avoid division by zero if S is passed as 0
                current_s_samples_weights[, g_comp] <- gamma_i_for_this_obs[g_comp] / S 
            }
            observation_specific_aug_weights_list[[g_comp]] <- current_s_samples_weights
          } # End loop over g_comp

          all_aug_data_list[[obs_idx]] <- do.call(rbind, observation_specific_aug_data_list)
          all_aug_weights_list[[obs_idx]] <- do.call(rbind, observation_specific_aug_weights_list)
        } # End loop over obs_idx

        aug_data_combined <- do.call(rbind, all_aug_data_list)
        aug_weights_combined <- do.call(rbind, all_aug_weights_list)
      } else { # No incomplete data
        aug_data_combined <- matrix(nrow = 0, ncol = p)
        # Ensure aug_weights_combined has G columns even if 0 rows
        aug_weights_combined <- matrix(nrow = 0, ncol = G)
      }

      # Combine data and weights for complete and augmented incomplete cases
      full_data <- rbind(split_data$data_comp, aug_data_combined)
      
      weights_for_mstep <- matrix(0, nrow = nrow(full_data), ncol = G)
      if (split_data$n0 > 0) {
        weights_for_mstep[1:split_data$n0, ] <- resp$gamma0
      }
      if (nrow(aug_data_combined) > 0) { # Check if there was any augmented data
        weights_for_mstep[(split_data$n0 + 1):nrow(full_data), ] <- aug_weights_combined
      }

      
      # M-step with mclust
      mstep_result <- mclust::mstep(data = full_data, modelName = modelName, z = weights_for_mstep)
      
      # Update parameters
      means <- lapply(1:G, function(g) as.numeric(mstep_result$parameters$mean[,g]))
      covs <- lapply(1:G, function(g) mstep_result$parameters$variance$sigma[,,g])
      pi <- mstep_result$parameters$pro
      
      # Calculate new log-likelihood
      new_resp <- Responsibility(split_data, means, covs, pi)
      loglik_new <- 0
      if(split_data$n0 > 0) {
        loglik_new <- loglik_new + sum(log(rowSums(new_resp$dens_eval0)))
      }
      if(split_data$n1 > 0) {
        loglik_new <- loglik_new + sum(log(rowSums(new_resp$dens_eval1)))
      }
      
      # Update imputed data for next iteration
      data_imputed_new <- data_imputed
      if(split_data$n1 > 0) {
        incomp_data <- split_data$data_incomp
        for(i in 1:nrow(incomp_data)) {
          y <- incomp_data[i,]
          idx_miss <- which(is.na(y))
          idx_obs <- which(!is.na(y))
          
          gamma_i <- new_resp$gamma1[i,]
          best_cluster <- which.max(gamma_i)
          
          cond_dist <- CalcCondDist(y, idx_miss, idx_obs, means[[best_cluster]], covs[[best_cluster]])
          imputed_values <- mvnfast::rmvn(1, mu = cond_dist$mu, sigma = cond_dist$sigma)
          data_imputed_new[split_data$idx_incomp[i], idx_miss] <- imputed_values
        }
      }
      
      # Store mstep_result for final output
      me_result <- list(
        parameters = mstep_result$parameters,
        loglik = loglik_new
      )
    }
    
    # Check convergence
    if (abs(loglik_new - loglik_old)/(1 + abs(loglik_new)) < tol) {
      converged <- TRUE
      if (verbose) cat(sprintf("Converged at iteration %d\n", iter))
      break
    }
    
    # Update for next iteration
    data_imputed <- data_imputed_new
    loglik_old <- loglik_new
  }
  
  if (!converged && verbose) {
    cat("Algorithm did not converge within maximum iterations\n")
  }
  
  # Final imputation
  if (method == "usual") {
    final_imputed <- impute_missing_values(data, me_result$parameters, z_current, missing_pattern)
    final_responsibilities <- z_current
  } else {
    # Final MMCEM2 imputation
    split_data <- PartitionData(data)
    final_resp <- Responsibility(split_data, means, covs, pi)
    
    final_imputed <- data
    if(split_data$n1 > 0) {
      incomp_data <- split_data$data_incomp
      for(i in 1:nrow(incomp_data)) {
        y <- incomp_data[i,]
        idx_miss <- which(is.na(y))
        idx_obs <- which(!is.na(y))
        
        gamma_i <- final_resp$gamma1[i,]
        best_cluster <- which.max(gamma_i)
        
        cond_dist <- CalcCondDist(y, idx_miss, idx_obs, means[[best_cluster]], covs[[best_cluster]])
        imputed_values <- mvnfast::rmvn(1, mu = cond_dist$mu, sigma = cond_dist$sigma)
        final_imputed[split_data$idx_incomp[i], idx_miss] <- imputed_values
      }
    }
    
    # Convert final responsibilities to matrix format
    final_responsibilities <- matrix(0, nrow = n, ncol = G)
    if(split_data$n0 > 0) {
      final_responsibilities[split_data$idx_comp, ] <- final_resp$gamma0
    }
    if(split_data$n1 > 0) {
      final_responsibilities[split_data$idx_incomp, ] <- final_resp$gamma1
    }
  }
  
  return(list(
    imputed_data = final_imputed,
    converged = converged,
    iterations = iter,
    loglik = loglik_new,
    parameters = me_result$parameters,
    responsibilities = final_responsibilities,
    method = method
  ))
}

impute_missing_values <- function(data, parameters, z, missing_pattern) {
  n <- nrow(data)
  p <- ncol(data)
  G <- length(parameters$pro)
  
  data_imputed <- data
  
  for (i in 1:n) {
    missing_i <- missing_pattern[i, ]
    
    if (any(missing_i)) {
      # Get observed values for this obs
      obs_i <- !missing_i
      
      if (any(obs_i)) {
        # Conditional imputation based on observed values
        imputed_values <- conditional_imputation(
          obs_values = data[i, obs_i],
          obs_indices = which(obs_i),
          miss_indices = which(missing_i),
          parameters = parameters,
          z_i = z[i, ]
        )
        
        data_imputed[i, missing_i] <- imputed_values
      } else {
        marginal_mean <- compute_marginal_mean(parameters)
        data_imputed[i, missing_i] <- marginal_mean[missing_i]
      }
    }
  }
  
  return(data_imputed)
}

conditional_imputation <- function(obs_values, obs_indices, miss_indices, parameters, z_i) {
  G <- length(parameters$pro)
  p_miss <- length(miss_indices)
  imputed_values <- numeric(p_miss)
  epsilon_pd <- sqrt(.Machine$double.eps)
  
  for (g in 1:G) {
    if (z_i[g] > 1e-10) {
      mu_g <- parameters$mean[, g]
      sigma_g <- get_component_covariance(parameters$variance, g)
      
      mu_obs <- mu_g[obs_indices]
      mu_miss <- mu_g[miss_indices]
      sigma_obs_obs <- sigma_g[obs_indices, obs_indices, drop = FALSE]
      sigma_miss_obs <- sigma_g[miss_indices, obs_indices, drop = FALSE]
      
      eigv <- eigen(sigma_obs_obs, symmetric = TRUE, only.values = TRUE)$values
      min_eig <- min(eigv)
      
      if (min_eig > epsilon_pd) {
        sigma_obs_inv <- solve(sigma_obs_obs)
      } else {
        # Add minimal ridge to ensure PD
        ridge <- max(0, -min_eig) + epsilon_pd
        sigma_obs_inv <- solve(sigma_obs_obs + diag(ridge, nrow(sigma_obs_obs)))
      }
      
      conditional_mean <- mu_miss + sigma_miss_obs %*% sigma_obs_inv %*% (obs_values - mu_obs)
      imputed_values <- imputed_values + z_i[g] * as.vector(conditional_mean)
    }
  }
  
  return(imputed_values)
}

# Helper function to extract covariance matrix from mclust variance structure
get_component_covariance <- function(variance_struct, g) {
  modelName <- variance_struct$modelName
  
  # For EEE model - equal covariance matrices
  if (modelName == "EEE") {
    return(variance_struct$Sigma)
  }
  
  # For spherical models (EII, VII)
  if (modelName == "EII") {
    p <- variance_struct$d
    return(variance_struct$sigmasq * diag(p))
  }
  
  if (modelName == "VII") {
    p <- variance_struct$d
    return(variance_struct$sigmasq[g] * diag(p))
  }
  
  # For models with full covariance matrices
  if (modelName %in% c("VVV", "EVE", "VEV", "EEV", "VVE", "EVV", "VEE")) {
    if (!is.null(variance_struct$sigma)) {
      return(variance_struct$sigma[,,g])
    }
  }
  
  # For diagonal models
  if (modelName %in% c("EEI", "VEI", "EVI", "VVI")) {
    if (!is.null(variance_struct$scale)) {
      return(diag(variance_struct$scale[,g]))
    }
  }
  
  # Try to reconstruct from eigenvalue decomposition if available
  if (!is.null(variance_struct$d) && !is.null(variance_struct$shape) && !is.null(variance_struct$orientation)) {
    # Todo: implementing the eigenvalue decomposition
    # Currently fall back to identity matrix
    warning("Could not extract covariance matrix, using identity matrix")
    return(diag(variance_struct$d))
  }
  stop(paste("Cannot extract covariance matrix for model", modelName, "and variance structure"))
}

resp_to_full_data <- function(init_result, data_imputed, complete_cases, G) {
  n_full <- nrow(data_imputed)
  z_full <- matrix(0, nrow = n_full, ncol = G)

  z_full[complete_cases, ] <- init_result$Z

  incomplete_cases <- !complete_cases
  
  if (any(incomplete_cases)) {
    temp_params <- list(
      pro = init_result$prop,
      mean = init_result$Mu,
      variance = list(
        modelName = "VVV",  
        d = ncol(data_imputed),
        G = G,
        sigma = init_result$SigmaCube
      )
    )
    
    incomplete_indices <- which(incomplete_cases)
    N_incomplete <- length(incomplete_indices)
    
    dens <- matrix(NA, N_incomplete, G)
    
    for (k in 1:G) {
      mu_k <- temp_params$mean[, k]
      sigma_k <- temp_params$variance$sigma[, , k]
      epsilon_pd <- sqrt(.Machine$double.eps)
    
      eigv <- eigen(sigma_k, symmetric = TRUE, only.values = TRUE)$values
      if (min(eigv) <= epsilon_pd) {
        ridge <- max(0, -min(eigv)) + epsilon_pd
        sigma_k <- sigma_k + diag(ridge, ncol(sigma_k))
      }
    
      for (idx in 1:N_incomplete) {
        i <- incomplete_indices[idx]
        dens[idx, k] <- mvtnorm::dmvnorm(data_imputed[i, ], 
                                        mean = mu_k, 
                                        sigma = sigma_k, 
                                        log = TRUE)
      }
    }
  
    denspro <- sweep(dens, 2, log(temp_params$pro), "+")
    
    z_max <- apply(denspro, 1, max)
    
    log_sum_exp <- z_max + log(rowSums(exp(denspro - z_max)))
    z_incomplete <- exp(denspro - log_sum_exp)
    
    for (idx in 1:N_incomplete) {
      if (any(is.na(z_incomplete[idx, ])) || sum(z_incomplete[idx, ]) == 0) {
        z_incomplete[idx, ] <- rep(1/G, G)
      }
    }
    z_full[incomplete_indices, ] <- z_incomplete
  }
  
  return(z_full)
}

compute_marginal_mean <- function(parameters) {
  G <- length(parameters$pro)
  p <- nrow(parameters$mean)
  
  marginal_mean <- numeric(p)
  for (g in 1:G) {
    marginal_mean <- marginal_mean + parameters$pro[g] * parameters$mean[, g]
  }
  
  return(marginal_mean)
}

map_mixall_to_mclust <- function(mixall_model) {
  mixall_to_mclust_map <- list(
    "gaussian_pk_sjk"       = "VVI",  
    "gaussian_pk_sj"        = "VII", 
    "gaussian_pk_sk"        = "VEI",  
    "gaussian_pk_s"         = "EII", 
    "gaussian_p_sjk"        = "EVI", 
    "gaussian_p_sj"         = "VII", 
    "gaussian_p_sk"         = "EEI",
    "gaussian_p_s"          = "EII")
  
  model_lower <- tolower(mixall_model)
  
  if (model_lower %in% names(mixall_to_mclust_map)) {
    return(mixall_to_mclust_map[[model_lower]])
  }
  
  if (grepl("gaussian.*pk.*sjk", model_lower)) return("VVI")
  if (grepl("gaussian.*pk.*sj", model_lower)) return("VII") 
  if (grepl("gaussian.*pk.*sk", model_lower)) return("VEI")
  if (grepl("gaussian.*pk.*s$", model_lower)) return("EII")
  if (grepl("gaussian.*p.*sjk", model_lower)) return("EVI")
  if (grepl("gaussian.*p.*sj", model_lower)) return("VII")
  if (grepl("gaussian.*p.*sk", model_lower)) return("EEI")
  if (grepl("gaussian.*p.*s$", model_lower)) return("EII")
  
  warning(paste("Unknown MixAll model:", mixall_model, 
                "- using VVI as default"))
  return("VVI")
}


map_rmixmod_to_mclust <- function(rmixmod_model) {
  # Rmixmod to mclust mapping
  rmixmod_to_mclust_map <- list(
    "Gaussian_p_L_I"        = "EII", 
    "Gaussian_p_Lk_I"       = "VII",  
    "Gaussian_p_L_B"        = "EEE",  
    "Gaussian_p_Lk_B"       = "VEE",  
    "Gaussian_p_L_Bk"       = "EEV",  
    "Gaussian_p_Lk_Bk"      = "VEV",  
    "Gaussian_p_L_C"        = "EVI",  
    "Gaussian_p_Lk_C"       = "VVI",  
    "Gaussian_p_L_D_Ak_D"   = "EEI",  
    "Gaussian_p_Lk_D_Ak_D"  = "VEI",  
    "Gaussian_p_L_Dk_A_Dk"  = "EVE",  
    "Gaussian_p_Lk_Dk_A_Dk" = "VVE",
    "Gaussian_p_L_Dk_Ak_Dk" = "EVV",  
    "Gaussian_p_Lk_Dk_Ak_Dk"= "VVV", 
    "Gaussian_pk_L_I"       = "EII",  
    "Gaussian_pk_Lk_I"      = "VII",
    "Gaussian_pk_L_B"       = "EEE", 
    "Gaussian_pk_Lk_B"      = "VEE",
    "Gaussian_pk_L_Bk"      = "EEV",
    "Gaussian_pk_Lk_Bk"     = "VEV",
    "Gaussian_pk_L_C"       = "EVI",
    "Gaussian_pk_Lk_C"      = "VVI",
    "Gaussian_pk_L_D_Ak_D"  = "EEI",
    "Gaussian_pk_Lk_D_Ak_D" = "VEI",
    "Gaussian_pk_L_Dk_A_Dk" = "EVE",
    "Gaussian_pk_Lk_Dk_A_Dk"= "VVE", 
    "Gaussian_pk_L_Dk_Ak_Dk"= "EVV",
    "Gaussian_pk_Lk_Dk_Ak_Dk"="VVV"
  )

  model_clean <- gsub("^Gaussian_", "", rmixmod_model)
  model_clean <- paste0("Gaussian_", model_clean)
  
  if (model_clean %in% names(rmixmod_to_mclust_map)) {
    return(rmixmod_to_mclust_map[[model_clean]])
  }
  
  if (rmixmod_model %in% names(rmixmod_to_mclust_map)) {
    return(rmixmod_to_mclust_map[[rmixmod_model]])
  }
  
  # Pattern matching
  model_lower <- tolower(rmixmod_model)
  if (grepl("spherical|p_l_i", model_lower)) return("EII")
  if (grepl("diagonal|p_lk_c", model_lower)) return("VVI")
  if (grepl("general|lk_dk_ak_dk", model_lower)) return("VVV")
  
  warning(paste("Unknown Rmixmod model:", rmixmod_model, 
                "- using VVV as default"))
  return("VVV")
} 

extract_rmixmod_model <- function(model_string) {
  if (grepl("mixmodGaussianModel", model_string)) {
    
    listmodels_pattern <- 'listModels\\s*=\\s*c\\s*\\(\\s*"([^"]+)"'
    listmodels_match <- regmatches(model_string, 
                                   regexpr(listmodels_pattern, model_string, perl = TRUE))
    
    if (length(listmodels_match) > 0) {
      model_name <- gsub(listmodels_pattern, '\\1', listmodels_match, perl = TRUE)
      return(model_name)
    }
    
    family_match <- regmatches(model_string, 
                               regexpr('family\\s*=\\s*"([^"]+)"', model_string, perl = TRUE))
    
    if (length(family_match) > 0) {
      family <- gsub('family\\s*=\\s*"([^"]+)"', '\\1', family_match, perl = TRUE)
      
      family_to_rmixmod <- list(
        "spherical" = "Gaussian_pk_L_I",
        "diagonal" = "Gaussian_pk_Lk_C", 
        "general" = "Gaussian_pk_Lk_Dk_Ak_Dk",
        "gaussian_pk_l_i" = "Gaussian_pk_L_I",
        "gaussian_pk_lk_c" = "Gaussian_pk_Lk_C",
        "gaussian_pk_lk_dk_ak_dk" = "Gaussian_pk_Lk_Dk_Ak_Dk"
      )
      
      family_lower <- tolower(family)
      
      if (family_lower %in% names(family_to_rmixmod)) {
        return(family_to_rmixmod[[family_lower]])
      } else {
        # If it's already in Rmixmod format, return as is (with proper casing)
        if (grepl("^gaussian_p[k]?_l[k]?_", family_lower)) {
          # Convert to proper Rmixmod format
          parts <- strsplit(family_lower, "_")[[1]]
          proper_case <- paste(
            "Gaussian",
            ifelse(parts[2] == "pk", "pk", "p"),
            paste(toupper(substring(parts[3:length(parts)], 1, 1)), 
                  substring(parts[3:length(parts)], 2), 
                  sep = "", collapse = "_"),
            sep = "_"
          )
          return(proper_case)
        }
        
        warning(paste("Unknown family:", family, "- using default Gaussian_pk_Lk_C"))
        return("Gaussian_pk_Lk_C")
      }
    }
  }
  if (grepl("^Gaussian_p[k]?_", model_string)) {
    return(model_string)
  }
  
  warning(paste("Could not parse model string:", model_string, "- using default"))
  return("Gaussian_pk_Lk_C")
}

map_to_mclust <- function(model_name, package_source = NULL) {
  if (is.null(package_source)) {
    model_lower <- tolower(model_name)
    
    if (grepl("gaussian_p[k]?_s[jk]?[jk]?", model_lower)) {
      package_source <- "mixall"
    }
    else if (grepl("gaussian_p[k]?_l[k]?_[bicd]", model_lower)) {
      package_source <- "rmixmod"
    }
    else if (grepl("^[evi]{3}$", model_lower)) {
      package_source <- "mclust"
      return(toupper(model_name)) 
    }
    else {
      package_source <- "mixall"
    }
  }
  
  switch(tolower(package_source),
         "mixall" = map_mixall_to_mclust(model_name),
         "rmixmod" = map_rmixmod_to_mclust(model_name), 
         "mclust" = toupper(model_name),
         {
           warning(paste("Unknown package source:", package_source))
           map_mixall_to_mclust(model_name)
         }
  )
}

get_mclust_models <- function() {
  c("EII", "VII", "EEI", "VEI", "EVI", "VVI", 
    "EEE", "EVE", "VEE", "EEV", "VEV", "EVV", "VVV", "X")
}

is_valid_mclust_model <- function(model_name) {
  toupper(model_name) %in% get_mclust_models()
}

map_and_validate_model <- function(model_name, 
                                   package_source = NULL, 
                                   validate = TRUE) {
  mclust_model <- map_to_mclust(model_name, package_source)
  
  if (validate && !is_valid_mclust_model(mclust_model)) {
    warning(paste("Mapped model", mclust_model, "is not a valid mclust model. Using VVV as default"))
    mclust_model <- "VVV"
  }
  
  return(mclust_model)
}

test_model_mappings <- function() {
  cat("Testing MixAll mappings:\n")
  mixall_models <- c("gaussian_pk_sjk", "gaussian_pk_sj", "gaussian_p_s",)
  for (model in mixall_models) {
    mapped <- map_mixall_to_mclust(model)
    cat(sprintf("  %s -> %s\n", model, mapped))
  }
  
  cat("\nTesting Rmixmod mappings:\n") 
  rmixmod_models <- c("Gaussian_p_L_I", "Gaussian_pk_Lk_C", "Gaussian_p_Lk_Dk_Ak_Dk", "spherical")
  for (model in rmixmod_models) {
    mapped <- map_rmixmod_to_mclust(model)  
    cat(sprintf("  %s -> %s\n", model, mapped))
  }
  
  cat("\nTesting universal mapper:\n")
  test_models <- c("gaussian_pk_sjk", "Gaussian_p_L_I", "VVI")
  for (model in test_models) {
    mapped <- map_to_mclust(model)
    cat(sprintf("  %s -> %s\n", model, mapped))
  }
}