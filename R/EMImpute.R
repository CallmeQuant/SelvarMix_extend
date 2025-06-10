EM_impute <- function(data, 
                     G = 3, 
                     modelName = "VVV",
                     max_iter = 1000,
                     tol = 1e-6,
                     init_method = "hc",
                     method = "usual", 
                     S = 50, 
                     verbose = FALSE,
                     use_glasso = FALSE,     
                     lambda_omega_0 = 50,     
                     n_samples = 100,           
                     burn_in_ratio = 0.1) {      
  
  # Validate method parameter
  method <- match.arg(method, c("usual", "sampling", "sampling_usual"))
  
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
  
  if (verbose) cat("Initializing parameters using method:", init_method, "\n")
  
  # Initialize based on method
  if (method %in% c("usual", "sampling_usual")) {
    complete_cases <- complete.cases(data)
    if (sum(complete_cases) < G) {
      warning("Not enough complete cases - using mean-imputed data for initialization")
      init_data <- data
      for (j in 1:p) {
        missing_j <- is.na(data[, j])
        if (any(missing_j)) {
          init_data[missing_j, j] <- mean(data[!missing_j, j], na.rm = TRUE)
        }
      }
      complete_cases <- rep(TRUE, n)
    } else {
      init_data <- data[complete_cases, , drop = FALSE]
    }
    
    init_result <- InitParameterRobust(init_data, 
                                 nbClust = G, 
                                 init = init_method,
                                 use_glasso = use_glasso,
                                 lambda_omega_0 = lambda_omega_0)
    
    # Initial imputation with means
    data_imputed <- data
    for (j in 1:p) {
      missing_j <- is.na(data[, j])
      if (any(missing_j)) {
        data_imputed[missing_j, j] <- mean(data[!missing_j, j], na.rm = TRUE)
      }
    }
    
    em_control <- mclust::emControl(tol = tol, itmax = max_iter)
    z_current <- resp_to_full_data(init_result, data_imputed, complete_cases, G)
    
  } else {  # MMCEM2 method
    split_data <- PartitionData(data)
    comp_data <- split_data$data_comp
    
    if(nrow(comp_data) < G) {
      stop("Number of complete cases is less than the number of clusters")
    }
    
    init_results <- InitParameterRobust(comp_data, 
                                  nbClust = G, 
                                  init = init_method,
                                  use_glasso = use_glasso,
                                  lambda_omega_0 = lambda_omega_0)
    
    means <- lapply(1:G, function(g) {
      as.numeric(init_results$Mu[, g])
    })
    
    covs <- lapply(1:G, function(g) {
      regularize_cov(init_results$SigmaCube[, , g])
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
          # Use most probable component for initial imputation
          g <- which.max(pi)
          cond_dist <- CalcCondDist(y, idx_miss, idx_obs, means[[g]], covs[[g]])
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
      me_result <- mclust::me(data = data_imputed, 
                      modelName = modelName,
                      z = z_current,
                      control = em_control)
      
      if (is.null(me_result) || is.null(me_result$parameters)) {
        if (verbose) cat("EM algorithm failed - using previous iteration\n")
        break
      }
      
      # Update responsibilities
      estep_result <- mclust::estep(modelName = modelName,
                           data = data_imputed,
                           parameters = me_result$parameters)
      
      z_current <- estep_result$z
      
      # Impute missing values
      data_imputed_new <- impute_missing_values(
        data, 
        me_result$parameters, 
        z_current, 
        missing_pattern
      )
      
      # Compute observed-data log-likelihood
      loglik_new <- compute_observed_loglik(data, me_result$parameters)
      
    } else if (method == "sampling_usual") {
      me_result <- mclust::me(data = data_imputed, 
                      modelName = modelName,
                      z = z_current,
                      control = em_control)
      
      if (is.null(me_result) || is.null(me_result$parameters)) {
        if (verbose) cat("EM algorithm failed - using previous iteration\n")
        break
      }
      
      # Update responsibilities
      estep_result <- mclust::estep(modelName = modelName,
                           data = data_imputed,
                           parameters = me_result$parameters)
      
      z_current <- estep_result$z
      
      # Sampling-based imputation
      data_imputed_new <- impute_with_sampling(
        data, 
        me_result$parameters, 
        z_current, 
        missing_pattern,
        n_samples = n_samples,
        burn_in_ratio = burn_in_ratio
      )
      
      # Compute observed-data log-likelihood
      loglik_new <- compute_observed_loglik(data, me_result$parameters)
      
    } else {  # MMCEM2 method
      split_data <- PartitionData(data) 
      
      # E-step: Calculate responsibilities
      resp <- Responsibility(split_data, means, covs, pi) 

      # Augment incomplete data with MC samples
      if (split_data$n1 > 0) {
        all_aug_data_list <- vector("list", split_data$n1)
        all_aug_weights_list <- vector("list", split_data$n1)

        for (obs_idx in 1:split_data$n1) {
          original_data_row_idx <- split_data$idx_incomp[obs_idx]
          y_original_incomplete_row <- data[original_data_row_idx, , drop = FALSE]
          is_na_y <- is.na(y_original_incomplete_row)
          idx_miss <- which(is_na_y)
          idx_obs <- which(!is_na_y)
          y_obs_values <- y_original_incomplete_row[1, idx_obs, drop = TRUE] 
          gamma_i <- resp$gamma1[obs_idx, ]
          
          observation_aug_data_list <- vector("list", G)
          observation_aug_weights_list <- vector("list", G)
      
          for (g_comp in 1:G) {
            current_imps <- matrix(NA, nrow = S, ncol = length(idx_miss))
            
            if (length(idx_miss) > 0) {
              if (length(idx_obs) > 0) {
                cond_dist <- CalcCondDist(y_original_incomplete_row, idx_miss, idx_obs, means[[g_comp]], covs[[g_comp]])
                
                # Check covariance PD status
                eigv <- eigen(cond_dist$sigma, symmetric = TRUE, only.values = TRUE)$values
                is_sigma_cond_pd <- min(eigv) > 1e-8
                
                if(is_sigma_cond_pd) {
                  current_imps <- mvnfast::rmvn(S, mu = cond_dist$mu, sigma = cond_dist$sigma)
                } else {
                  # Use regularized covariance
                  cond_dist$sigma <- regularize_cov(cond_dist$sigma)
                  current_imps <- mvnfast::rmvn(S, mu = cond_dist$mu, sigma = cond_dist$sigma)
                }
              } else {
                # Marginal imputation
                marginal_cov_miss <- covs[[g_comp]][idx_miss, idx_miss, drop = FALSE]
                marginal_cov_miss <- regularize_cov(marginal_cov_miss)
                current_imps <- mvnfast::rmvn(S, mu = means[[g_comp]][idx_miss], sigma = marginal_cov_miss)
              }
            }
            
            # Reconstruct full data vectors
            reconstructed_samples <- matrix(NA, nrow = S, ncol = p)
            if (length(idx_obs) > 0) {
              reconstructed_samples[, idx_obs] <- matrix(rep(y_obs_values, S), nrow = S, byrow = TRUE)
            }
            if (length(idx_miss) > 0) {
              reconstructed_samples[, idx_miss] <- current_imps
            }
            observation_aug_data_list[[g_comp]] <- reconstructed_samples
            
            # Weights for these samples
            current_weights <- matrix(0, nrow = S, ncol = G)
            current_weights[, g_comp] <- gamma_i[g_comp] / S
            observation_aug_weights_list[[g_comp]] <- current_weights
          }
          
          all_aug_data_list[[obs_idx]] <- do.call(rbind, observation_aug_data_list)
          all_aug_weights_list[[obs_idx]] <- do.call(rbind, observation_aug_weights_list)
        }
        
        aug_data_combined <- do.call(rbind, all_aug_data_list)
        aug_weights_combined <- do.call(rbind, all_aug_weights_list)
      } else {
        aug_data_combined <- matrix(nrow = 0, ncol = p)
        aug_weights_combined <- matrix(nrow = 0, ncol = G)
      }

      # Combine data and weights
      full_data <- rbind(split_data$data_comp, aug_data_combined)
      weights_for_mstep <- matrix(0, nrow = nrow(full_data), ncol = G)
      if (split_data$n0 > 0) {
        weights_for_mstep[1:split_data$n0, ] <- resp$gamma0
      }
      if (nrow(aug_data_combined) > 0) {
        weights_for_mstep[(split_data$n0 + 1):nrow(full_data), ] <- aug_weights_combined
      }
      
      # M-step with mclust
      mstep_result <- mclust::mstep(data = full_data, modelName = modelName, z = weights_for_mstep)
      
      # Update parameters
      means <- lapply(1:G, function(g) as.numeric(mstep_result$parameters$mean[,g]))
      covs <- lapply(1:G, function(g) {
        cov_mat <- mstep_result$parameters$variance$sigma[,,g]
        regularize_cov(cov_mat)
      })
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
      
      # Update imputed data
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
      
      # Store mstep_result
      me_result <- list(
        parameters = mstep_result$parameters,
        loglik = loglik_new
      )
    }
    
    # Check convergence
    if (method %in% c("usual", "sampling_usual")) {
      delta_ll <- (loglik_new - loglik_old) / max(1, abs(loglik_old))
    } else {
      delta_ll <- (loglik_new - loglik_old) / max(1, abs(loglik_new))
    }
    
    if (abs(delta_ll) < tol) {
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
  } else if (method == "sampling_usual") {
    final_imputed <- impute_with_sampling(
      data, 
      me_result$parameters, 
      z_current, 
      missing_pattern,
      n_samples = n_samples,
      burn_in_ratio = burn_in_ratio
    )
    final_responsibilities <- z_current
  } else {  # MMCEM2 method
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
    gamma0[is.na(gamma0) | is.nan(gamma0)] <- 1/k 
    
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
    gamma1[is.na(gamma1) | is.nan(gamma1)] <- 1/k
    
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

InitParameterRobust <- function(data, nbClust, init = c("kmeans", "hc"), n.start = 25,
                          use_glasso = FALSE, lambda_omega_0 = 50) {
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  init <- match.arg(init)
  
  # Initialization with dimensionality reduction
  robust_init <- function(data, G) {
    # Try PCA-based k-means first
    tryCatch({
      if (p > 1) {
        pc <- prcomp(data, scale. = TRUE)
        data_red <- pc$x[, 1:min(2, ncol(pc$x)), drop = FALSE]
      } else {
        data_red <- data
      }
      km <- kmeans(data_red, centers = G, nstart = n.start)
      return(mclust::unmap(km$cluster))
    }, error = function(e) {
      # Fallback to hierarchical clustering
      tryCatch({
        hc <- hclust(dist(data))
        return(mclust::unmap(cutree(hc, k = G)))
      }, error = function(e) {
        # Final fallback: random assignment
        memb <- sample(1:G, n, replace = TRUE)
        return(mclust::unmap(memb))
      })
    })
  }
  
  if (init == "kmeans") {
    z_mat <- robust_init(data, nbClust)
  } else {
    tryCatch({
      hcobj <- mclust::hc(data, modelName = "VVV", use = "SVD")
      z_mat <- mclust::unmap(mclust::hclass(hcobj, nbClust))
    }, error = function(e) {
      z_mat <- robust_init(data, nbClust)
    })
  }
  if (use_glasso){
    n_k <- pmax(colSums(z_mat), 1)
    prop <- n_k / n
    Mu <- t(z_mat) %*% data / n_k
  
    # Compute regularized covariance
    SigmaCube <- array(0, dim = c(p, p, nbClust))
    for (g in 1:nbClust) {
      cluster_data <- data[z_mat[, g] > 0.5, , drop = FALSE]
      if (nrow(cluster_data) > 0) {
        SigmaCube[, , g] <- regularize_cov(cov(cluster_data))
      } else {
        SigmaCube[, , g] <- diag(p)
      }
    }
    return(list(
     prop = prop,
     Mu = Mu,
     SigmaCube = SigmaCube,
     Z = z_mat))
  }
  else{
    n_k <- pmax(colSums(z_mat), 1)
    prop <- n_k / n  
    temp <- mclust::covw(data, z_mat, normalize = FALSE) 
    sigma0 <- temp$S                              
    omega0 <- array(NA_real_, dim = dim(sigma0))    
    for (k in 1:nbClust) {
      eigv <- eigen(sigma0[,,k], symmetric = TRUE, only.values = TRUE)$values
      # If invertible, invert it
      # Otherwise, use glassoFast to estimate the precision matrix
      if (min(eigv) > sqrt(.Machine$double.eps)) {             
        omega0[,,k] <- solve(sigma0[,,k])
      } else {                                       
        gl <- glassoFast::glassoFast(
                S   = sigma0[,,k],
                rho = lambda_omega_0 * sqrt(log(p)/n_k[k]),
                thr = 1e-4, maxIt = 1000)
        omega0[,,k] <- gl$wi
        sigma0[,,k] <- gl$w
      }
    }
    Mu <- t(z_mat) %*% data / n_k
    return(list(prop      = prop,
         Mu        = Mu,
         SigmaCube = sigma0,    
         Z         = z_mat))}
}

# Regularization function
regularize_cov <- function(sigma, epsilon = 1e-8) {
  eig <- eigen(sigma, symmetric = TRUE, only.values = TRUE)
  min_eig <- min(eig$values)
  if (min_eig < epsilon) {
    ridge <- epsilon - min_eig
    return(sigma + diag(ridge, nrow(sigma)))
  }
  return(sigma)
}

# Compute observed-data log-likelihood
compute_observed_loglik <- function(data, params) {
  loglik <- 0
  for (i in 1:n) {
    obs <- !is.na(data[i, ])
    if (any(obs)) {
      log_probs <- sapply(1:G, function(g) {
        mu_g <- params$mean[obs, g]
        sigma_g <- get_component_covariance(params$variance, g)[obs, obs, drop = FALSE]
        sigma_g <- regularize_cov(sigma_g)
        log_dens <- mvtnorm::dmvnorm(
          data[i, obs], 
          mean = mu_g, 
          sigma = sigma_g, 
          log = TRUE
        )
        log(params$pro[g]) + log_dens
      })
      max_log <- max(log_probs)
      loglik <- loglik + max_log + log(sum(exp(log_probs - max_log)))
    }
  }
  return(loglik)
}
  
# Sampling-based imputation function
impute_with_sampling <- function(data, parameters, z, missing_pattern, 
                                  n_samples = 100, burn_in_ratio = 0.1) {
  n <- nrow(data)
  p <- ncol(data)
  G <- length(parameters$pro)
  burn_in <- max(1, floor(n_samples * burn_in_ratio))
  effective_samples <- n_samples - burn_in
  
  # Validate sample counts
  if (effective_samples < 1) {
    stop("n_samples too small after burn-in. Increase n_samples or decrease burn_in_ratio.")
  }
  
  # Identify unique missing patterns
  pattern_str <- apply(missing_pattern, 1, function(x) paste(which(x), collapse = ","))
  unique_patterns <- unique(pattern_str)
  
  data_imputed <- data
  
  for (pat in unique_patterns) {
    pat_indices <- which(pattern_str == pat)
    if (length(pat_indices) == 0) next
    
    obs_indices <- !missing_pattern[pat_indices[1], ]
    mis_indices <- missing_pattern[pat_indices[1], ]
    num_missing <- sum(mis_indices)
    num_obs <- length(pat_indices)
    
    if (any(mis_indices)) {
      if (any(obs_indices)) {
        imp_vals <- matrix(0, nrow = num_obs, ncol = num_missing)
        
        for (i in 1:num_obs) {
          idx <- pat_indices[i]
          obs_data <- data[idx, obs_indices]
          z_i <- z[idx, ]
          
          # Clean responsibilities vector
          z_i[is.na(z_i)] <- 0
          z_i[z_i < 0] <- 0
          if (sum(z_i) < .Machine$double.eps) {
            z_i <- rep(1/G, G)  # Fallback to uniform
          } else {
            z_i <- z_i / sum(z_i)  # Normalize
          }
          
          # Identify valid components
          valid_g <- which(z_i > 1e-6)
          if (length(valid_g) == 0) valid_g <- 1:G
          
          # Draw samples
          all_samples <- matrix(0, nrow = n_samples, ncol = num_missing)
          for (s in 1:n_samples) {
            # Sample component safely
            g <- if (length(valid_g) == 1) {
              valid_g
            } else {
              sample(valid_g, 1, prob = z_i[valid_g])
            }
            
            mu_g <- parameters$mean[, g]
            sigma_g <- get_component_covariance(parameters$variance, g)
            
            mu_obs <- mu_g[obs_indices]
            mu_miss <- mu_g[mis_indices]
            sigma_obs_obs <- sigma_g[obs_indices, obs_indices, drop = FALSE]
            sigma_miss_obs <- sigma_g[mis_indices, obs_indices, drop = FALSE]
            sigma_miss_miss <- sigma_g[mis_indices, mis_indices, drop = FALSE]
            
            # Regularize and compute conditional distribution
            sigma_obs_obs <- regularize_cov(sigma_obs_obs)
            sigma_obs_inv <- solve(sigma_obs_obs)
            
            cond_mean <- mu_miss + sigma_miss_obs %*% sigma_obs_inv %*% (obs_data - mu_obs)
            cond_cov <- sigma_miss_miss - sigma_miss_obs %*% sigma_obs_inv %*% t(sigma_miss_obs)
            cond_cov <- regularize_cov(cond_cov)
            
            # Draw sample safely
            all_samples[s, ] <- tryCatch(
              mvtnorm::rmvnorm(1, cond_mean, cond_cov),
              error = function(e) cond_mean  # Fallback to mean
            )
          }
          
          # Apply burn-in and compute robust estimate
          post_burn_samples <- all_samples[(burn_in + 1):n_samples, , drop = FALSE]
          imp_vals[i, ] <- apply(post_burn_samples, 2, median)
        }
        data_imputed[pat_indices, mis_indices] <- imp_vals
      } else {
        # Completely missing rows
        marg_mean <- rowSums(parameters$mean %*% diag(parameters$pro))
        data_imputed[pat_indices, mis_indices] <- matrix(
          marg_mean[mis_indices], 
          nrow = num_obs, 
          ncol = num_missing, 
          byrow = TRUE
        )
      }
    }
  }
  return(data_imputed)
}

# Conditional mean imputation function
impute_missing_values <- function(data, parameters, z, missing_pattern) {
  n <- nrow(data)
  p <- ncol(data)
  G <- length(parameters$pro)
  data_imputed <- data
  
  # Identify unique missing patterns
  pattern_str <- apply(missing_pattern, 1, function(x) paste(which(x), collapse = ","))
  unique_patterns <- unique(pattern_str)
  
  for (pat in unique_patterns) {
    pat_indices <- which(pattern_str == pat)
    obs_indices <- !missing_pattern[pat_indices[1], ]
    mis_indices <- missing_pattern[pat_indices[1], ]
    num_missing <- sum(mis_indices)
    num_obs <- length(pat_indices)
    
    if (any(mis_indices)) {
      if (any(obs_indices)) {
        # Batch impute for all observations with same pattern
        imp_vals <- matrix(0, nrow = num_obs, ncol = num_missing)
        obs_data <- data[pat_indices, obs_indices, drop = FALSE]
        
        for (g in 1:G) {
          # Only process components with significant responsibility
          comp_responsibilities <- z[pat_indices, g]
          if (max(comp_responsibilities) > 1e-3) {
            mu_g <- parameters$mean[, g]
            sigma_g <- get_component_covariance(parameters$variance, g)
            
            mu_obs <- mu_g[obs_indices]
            mu_miss <- mu_g[mis_indices]
            sigma_obs_obs <- sigma_g[obs_indices, obs_indices, drop = FALSE]
            sigma_miss_obs <- sigma_g[mis_indices, obs_indices, drop = FALSE]
            
            # Regularize
            sigma_obs_obs <- regularize_cov(sigma_obs_obs)
            sigma_obs_inv <- solve(sigma_obs_obs)
            
            # Compute conditional means
            obs_diff <- t(t(obs_data) - mu_obs)
            cond_mean <- mu_miss + 
              (sigma_miss_obs %*% sigma_obs_inv) %*% 
              t(obs_diff)
            
            # Weight by responsibilities
            weights <- matrix(comp_responsibilities, nrow = num_obs, ncol = num_missing)
            weighted_cond <- t(cond_mean) * weights
            
            imp_vals <- imp_vals + weighted_cond
          }
        }
        data_imputed[pat_indices, mis_indices] <- imp_vals
      } else {
        # Completely missing rows
        marg_mean <- rowSums(parameters$mean %*% diag(parameters$pro))
        data_imputed[pat_indices, mis_indices] <- matrix(
          marg_mean[mis_indices], 
          nrow = num_obs, 
          ncol = num_missing, 
          byrow = TRUE
        )
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
  
  if (any(!is.finite(z_i)) || abs(sum(z_i[is.finite(z_i)]) - 1) > 1e-6 
      || sum(z_i[is.finite(z_i)]) == 0) {
    warning("Responsibilities in conditional_imputation are NA, non-finite.")
    valid_entries <- is.finite(z_i) & z_i >= 0
    if (sum(z_i[valid_entries]) > 0) {
        z_i[!valid_entries] <- 0
        z_i <- z_i / sum(z_i)
    } else {
        z_i <- rep(1/G, G)
    }
  }
  
  if (abs(sum(z_i) - 1) > 1e-6) { 
      if(sum(z_i) > .Machine$double.eps) { 
          z_i <- z_i / sum(z_i)
      } else {
          z_i <- rep(1/G, G) 
      }
  }

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

get_component_covariance <- function(variance_struct, g) {
  # Ensure variance_struct is not NULL and is a list
  if (is.null(variance_struct) || !is.list(variance_struct)) {
    stop("variance_struct is NULL or not a list.")
  }

  # Ensure modelName, d (dimension), and G (num components) are present and valid
  if (is.null(variance_struct$modelName) || !is.character(variance_struct$modelName) || length(variance_struct$modelName) != 1) {
    stop("variance_struct must contain a single string 'modelName'.")
  }
  if (is.null(variance_struct$d) || !is.numeric(variance_struct$d) || length(variance_struct$d) != 1 || variance_struct$d < 0 || floor(variance_struct$d) != variance_struct$d) {
    stop("variance_struct must contain 'd' (dimension) as a single non-negative integer.")
  }
  # G=0 is not typical for a fitted model's parameters. If it occurs, it implies no components.
  # If G=0, then g cannot be valid.
  if (is.null(variance_struct$G) || !is.numeric(variance_struct$G) || length(variance_struct$G) != 1 || variance_struct$G < 1 || floor(variance_struct$G) != variance_struct$G) {
    stop("variance_struct must contain 'G' (number of components) as a single positive integer.")
  }
  
  modelName <- variance_struct$modelName
  d <- variance_struct$d
  G <- variance_struct$G

  # Validate g (g > G is now safe as G is confirmed to be a positive integer)
  if (!is.numeric(g) || length(g) != 1 || g < 1 || g > G ) { 
    stop(paste("Invalid component index 'g'=", g, ". Must be an integer between 1 and G (", G, ").", sep=""))
  }

  # --- Univariate Models ---
  if (d == 1) {
    if (modelName == "E") { # Equal variance
      if (!is.null(variance_struct$sigmasq) && is.numeric(variance_struct$sigmasq) && length(variance_struct$sigmasq) == 1) {
        return(matrix(variance_struct$sigmasq))
      } else {
        stop("For univariate 'E' model, variance_struct$sigmasq should be a single numeric value.")
      }
    }
    if (modelName == "V") { # Variable variance
      if (!is.null(variance_struct$sigmasq) && is.numeric(variance_struct$sigmasq) && length(variance_struct$sigmasq) == G) {
        return(matrix(variance_struct$sigmasq[g]))
      } else {
        stop(paste("For univariate 'V' model, variance_struct$sigmasq should be a numeric vector of length G (",G,").",sep=""))
      }
    }
    # If d=1 but modelName is a multivariate one, it will fall through. 
    # This is generally okay as mclust would fit "E" or "V" for d=1.
    # If a multivariate model name is forced on d=1 data, the following sections should handle it if d=1.
  }

  # --- Multivariate Models ---

  # Spherical models: Sigma_k = lambda_k * I
  if (modelName == "EII") { # Spherical, equal volume (lambda * I)
    if (!is.null(variance_struct$sigmasq) && is.numeric(variance_struct$sigmasq) && length(variance_struct$sigmasq) == 1) {
      return(diag(x = variance_struct$sigmasq, nrow = d, ncol = d))
    } else {
      stop("For 'EII' model, variance_struct$sigmasq should be a single numeric value.")
    }
  }
  if (modelName == "VII") { # Spherical, unequal volume (lambda_k * I)
    if (!is.null(variance_struct$sigmasq) && is.numeric(variance_struct$sigmasq) && length(variance_struct$sigmasq) == G) {
      return(diag(x = variance_struct$sigmasq[g], nrow = d, ncol = d))
    } else {
      stop(paste("For 'VII' model, variance_struct$sigmasq should be a numeric vector of length G (",G,").",sep=""))
    }
  }
  
  if (modelName %in% c("EEI", "VEI", "EVI", "VVI")) {
    # Path 1: Use variance_struct$sigma directly (preferred for diagonal models)
    # According to mclust documentation, variance_struct$sigma is a (d x d x G) array for these.
    if (!is.null(variance_struct$sigma) && is.array(variance_struct$sigma) && 
        length(dim(variance_struct$sigma)) == 3 &&
        dim(variance_struct$sigma)[1] == d && dim(variance_struct$sigma)[2] == d && dim(variance_struct$sigma)[3] == G) {
      return(variance_struct$sigma[,,g, drop = FALSE])
    } else if (modelName == "EEI" && !is.null(variance_struct$Sigma) && is.matrix(variance_struct$Sigma) && 
               all(dim(variance_struct$Sigma) == c(d,d))) {
      # Path 1b: EEI might store a single common covariance matrix variance_struct$Sigma
       return(variance_struct$Sigma)
    } else if (!is.null(variance_struct$scale) && !is.null(variance_struct$shape)) {
        # Path 2: Fallback to reconstruction from scale and shape if sigma array is not as expected.
        # This path is less common for these models if mclust output is standard.
        
        current_scale_val <- NA # Renamed from current_scale to avoid conflict with base::scale
        if (length(variance_struct$scale) == 1) { # Common scale (E models like EEI, EVI)
            current_scale_val <- variance_struct$scale
        } else if (length(variance_struct$scale) == G) { # Varying scale (V models like VEI, VVI)
            current_scale_val <- variance_struct$scale[g]
        } else {
            stop(paste("For diagonal model", modelName, ", 'scale' has unexpected length. Expected 1 or G (",G,"). Actual:", length(variance_struct$scale)))
        }

        current_shape_diag_elements <- NULL
        if (is.vector(variance_struct$shape) && length(variance_struct$shape) == d) { # Common shape (E models like EEI, VEI)
            current_shape_diag_elements <- variance_struct$shape
        } else if (is.matrix(variance_struct$shape) && nrow(variance_struct$shape) == d && ncol(variance_struct$shape) == G) { # Varying shape (V models like EVI, VVI)
            current_shape_diag_elements <- variance_struct$shape[,g]
        } else {
            stop(paste("For diagonal model", modelName, ", 'shape' has unexpected dimensions. Expected vector of length d (",d,") or d x G matrix (",d,"x",G,"). Actual dims:", paste(dim(variance_struct$shape),collapse="x"), "or length:", length(variance_struct$shape)))
        }
        
        # The elements of 'shape' are normalized such that prod(shape[,k]) = 1 (for d>0).
        # The diagonal elements of Sigma_k are lambda_k * shape_elements_k,
        # where lambda_k is variance_struct$scale[k]^d (volume).
        lambda_k_volume <- current_scale_val^d
        if (d == 0 && current_scale_val == 0) lambda_k_volume <- 0 
        else if (d == 0 && current_scale_val != 0) lambda_k_volume <- 1

        diag_elements <- lambda_k_volume * current_shape_diag_elements
        return(diag(diag_elements, nrow = d, ncol = d))
    } else {
      stop(paste("For diagonal model '", modelName, "', expected variance_struct$sigma (dxdxG array), or variance_struct$Sigma (dxd matrix for EEI), or variance_struct$scale & shape were not found or structured as expected.", sep=""))
    }
  }

  # EEE: Sigma_k = Sigma (common volume, shape, orientation)
  if (modelName == "EEE") {
    if (!is.null(variance_struct$Sigma) && is.matrix(variance_struct$Sigma) && all(dim(variance_struct$Sigma) == c(d,d))) {
      return(variance_struct$Sigma)
    } else if (!is.null(variance_struct$cholSigma) && is.matrix(variance_struct$cholSigma) && all(dim(variance_struct$cholSigma) == c(d,d))) {
      return(crossprod(variance_struct$cholSigma))
    } else {
      stop("For 'EEE' model, variance_struct$Sigma (d x d matrix) or variance_struct$cholSigma was not found or structured as expected.")
    }
  }

  # VVV: Sigma_k = lambda_k * D_k * A_k * D_k^T (varying volume, shape, orientation)
  # mclust stores these directly in variance_struct$sigma[,,g]
  if (modelName == "VVV") {
    if (!is.null(variance_struct$sigma) && is.array(variance_struct$sigma) && 
        length(dim(variance_struct$sigma)) == 3 &&
        dim(variance_struct$sigma)[1] == d && dim(variance_struct$sigma)[2] == d && dim(variance_struct$sigma)[3] == G) {
      return(variance_struct$sigma[,,g, drop = FALSE])
    } else if (!is.null(variance_struct$cholsigma) && is.array(variance_struct$cholsigma) &&
               length(dim(variance_struct$cholsigma)) == 3 &&
               dim(variance_struct$cholsigma)[1] == d && dim(variance_struct$cholsigma)[2] == d && dim(variance_struct$cholsigma)[3] == G) {
      return(crossprod(variance_struct$cholsigma[,,g, drop = FALSE]))
    } else {
      stop(paste("For 'VVV' model, variance_struct$sigma (dxdxG array) or variance_struct$cholsigma (dxdxG array) was not found or structured as expected. Dims of sigma:", paste(dim(variance_struct$sigma), collapse="x")))
    }
  }
  
  # Models requiring reconstruction from scale, shape, and orientation: EEV, EVE, EVV, VEE, VEV, VVE
  # These models have specific combinations of equal/varying scale, shape, orientation.
  
  current_scale_val <- NA
  # Determine scale for component g
  if (modelName %in% c("EEV", "EVE", "EVV")) { # Equal volume
    if (!is.null(variance_struct$scale) && is.numeric(variance_struct$scale) && length(variance_struct$scale) == 1) {
      current_scale_val <- variance_struct$scale
    } else {
      stop(paste("For model '", modelName, "', variance_struct$scale should be a single numeric value (equal volume). Found length:", length(variance_struct$scale)))
    }
  } else if (modelName %in% c("VEE", "VEV", "VVE")) { # Varying volume
    if (!is.null(variance_struct$scale) && is.numeric(variance_struct$scale) && length(variance_struct$scale) == G) {
      current_scale_val <- variance_struct$scale[g]
    } else {
      stop(paste("For model '", modelName, "', variance_struct$scale should be a numeric vector of length G (",G,"). Found length:", length(variance_struct$scale)))
    }
  } else {
    # This block should only be reached if modelName is one of EEV, EVE, EVV, VEE, VEV, VVE
    # but the scale was not assigned, which implies a logic error or unhandled model.
    # However, the previous specific model blocks (EII, VII, EEI..., EEE, VVV) should catch most models.
    # This path is for the remaining ellipsoidal models.
  }

  current_shape_matrix <- NULL
  # Determine shape matrix (diag of normalized eigenvalues) for component g
  if (modelName %in% c("EVE", "VEE")) { # Equal shape
    if (!is.null(variance_struct$shape) && is.numeric(variance_struct$shape) && length(variance_struct$shape) == d) {
      if(d > 0 && abs(prod(variance_struct$shape) - 1) > 1e-6) warning(paste("Product of shape eigenvalues for model", modelName, "is not 1. Product:", prod(variance_struct$shape)))
      current_shape_matrix <- diag(variance_struct$shape, nrow = d, ncol = d)
    } else {
      stop(paste("For model '", modelName, "', variance_struct$shape should be a numeric vector of length d (",d,") (equal shape). Found length:", length(variance_struct$shape)))
    }
  } else if (modelName %in% c("EEV", "VEV", "EVV", "VVE")) { # Varying shape
     if (!is.null(variance_struct$shape) && is.matrix(variance_struct$shape) && 
        nrow(variance_struct$shape) == d && ncol(variance_struct$shape) == G) {
      if(d > 0 && abs(prod(variance_struct$shape[,g]) - 1) > 1e-6) warning(paste("Product of shape eigenvalues for model", modelName, "component", g, "is not 1. Product:", prod(variance_struct$shape[,g])))
      current_shape_matrix <- diag(variance_struct$shape[,g], nrow = d, ncol = d)
    } else {
      stop(paste("For model '", modelName, "', variance_struct$shape should be a d x G matrix (d=",d,", G=",G,"). Actual dims:", paste(dim(variance_struct$shape), collapse="x")))
    }
  }

  current_orientation_matrix <- NULL
  # Determine orientation matrix for component g
  if (modelName %in% c("EEV", "VEE")) { # Equal orientation
    if (!is.null(variance_struct$orientation) && is.matrix(variance_struct$orientation) && 
        all(dim(variance_struct$orientation) == c(d,d))) {
      current_orientation_matrix <- variance_struct$orientation
    } else {
      stop(paste("For model '", modelName, "', variance_struct$orientation should be a d x d matrix (equal orientation). Actual dims:", paste(dim(variance_struct$orientation), collapse="x")))
    }
  } else if (modelName %in% c("EVE", "VEV", "EVV", "VVE")) { # Varying orientation
    if (!is.null(variance_struct$orientation) && is.array(variance_struct$orientation) &&
        length(dim(variance_struct$orientation)) == 3 &&
        dim(variance_struct$orientation)[1] == d && dim(variance_struct$orientation)[2] == d && 
        dim(variance_struct$orientation)[3] == G) {
      current_orientation_matrix <- variance_struct$orientation[,,g, drop = FALSE]
    } else {
      stop(paste("For model '", modelName, "', variance_struct$orientation should be a d x d x G array (d=",d,", G=",G,"). Actual dims:", paste(dim(variance_struct$orientation), collapse="x")))
    }
  }
  
  # Reconstruct covariance if all parts for ellipsoidal models (EEV, EVE, EVV, VEE, VEV, VVE) are found
  if (!is.na(current_scale_val) && !is.null(current_shape_matrix) && !is.null(current_orientation_matrix)) {
    lambda_k_volume <- current_scale_val^d
    if (d == 0 && current_scale_val == 0) lambda_k_volume <- 0 
    else if (d == 0 && current_scale_val != 0) lambda_k_volume <- 1 

    if (d > 0) {
      # Sigma_k = D_k * (lambda_k * A_k) * D_k^T
      # where lambda_k is volume, D_k is orientation, A_k is diagonal matrix of normalized eigenvalues (shape)
      cov_matrix <- current_orientation_matrix %*% (lambda_k_volume * current_shape_matrix) %*% t(current_orientation_matrix)
      return(cov_matrix)
    } else if (d == 0) { # d == 0
      return(matrix(numeric(0), nrow=0, ncol=0)) 
    }
  }
  stop(paste("Unsupported or unrecognized modelName '", modelName, 
             "' or variance structure not correctly specified for reconstruction after checking all known model types.", sep=""))
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