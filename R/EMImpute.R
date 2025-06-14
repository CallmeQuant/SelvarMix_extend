EM_impute <- function(data, 
                     G = 3, 
                     modelName = "VVV",
                     max_iter = 1000,
                     tol = 1e-6,
                     init_method = "hc",
                     use_glasso = FALSE,
                     lambda_omega_0 = 50,
                     sampling=FALSE,
                     n_samples = 100,
                     verbose = FALSE) {
  # Validate inputs
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("data must be a data frame or matrix")
  }
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  
  # Identify missing data patterns
  missing_pattern <- is.na(data)
  has_missing <- any(missing_pattern)
  if (!has_missing) {
    warning("No missing values found in data")
    return(list(
      imputed_data = data,
      converged = TRUE,
      iterations = 0,
      loglik = NA,
      parameters = NULL
    ))
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
  
  # Initialize parameters
  complete_cases <- complete.cases(data)
  init_data <- data
  if (sum(complete_cases) < G) {
    if (verbose) message("Insufficient complete cases - using mean-imputed data for initialization")
    for (j in 1:p) {
      missing_j <- is.na(data[, j])
      if (any(missing_j)) {
        init_data[missing_j, j] <- mean(data[, j], na.rm = TRUE)
      }
    }
    complete_cases <- rep(TRUE, n)
  } else {
    init_data <- data[complete_cases, , drop = FALSE]
  }
  
  init_result <- InitParameterRobust(
    data = init_data,
    nbClust = G,
    init = init_method,
    use_glasso = use_glasso,
    lambda_omega_0 = lambda_omega_0
  )
  
  # Create initial parameters structure
  parameters <- list(
    pro = init_result$prop,
    mean = init_result$Mu,
    variance = list(
      modelName = modelName,
      d = p,
      G = G,
      sigma = init_result$SigmaCube
    )
  )
  
  # Initial imputation
  data_imputed <- data
  for (j in 1:p) {
    missing_j <- is.na(data[, j])
    if (any(missing_j)) {
      data_imputed[missing_j, j] <- mean(data[, j], na.rm = TRUE)
    }
  }
  
  # Get initial responsibilities
  z_current <- resp_to_full_data(init_result, data_imputed, complete_cases, G)
  
  # EM control
  em_control <- mclust::emControl(tol = tol, itmax = max_iter)
  
  # EM algorithm
  converged <- FALSE
  loglik_old <- -Inf
  iter <- 0
  
  for (iter in 1:max_iter) {
    if (verbose && iter %% 10 == 0) {
      cat(sprintf("Iteration %d\n", iter))
    }
    
    # prior <- mclust::priorControl(
    #                   functionName = "defaultPrior",
    #                   scale = prior_scale,
    #                   dof = p + 2)
    
    # M-step using mclust::me
    me_result <- mclust::me(
      data = data_imputed, 
      modelName = modelName,
      z = z_current,
      control = em_control
    )
    
    if (is.null(me_result) || is.null(me_result$parameters)) {
      if (verbose) cat("EM algorithm failed - using previous iteration\n")
      break
    }
    
    # E-step using mclust::estep
    estep_result <- mclust::estep(
      modelName = modelName,
      data = data_imputed,
      parameters = me_result$parameters
    )
    z_current <- estep_result$z
    
    # Compute observed-data log-likelihood
    loglik_new <- compute_observed_loglik(data, me_result$parameters)
    
    # Check convergence
    if (iter > 1) {
      delta_ll <- (loglik_new - loglik_old) / max(1, abs(loglik_old))
      if (abs(delta_ll) < tol) {
        converged <- TRUE
        if (verbose) cat(sprintf("Converged at iteration %d\n", iter))
        break
      }
      if (delta_ll < -1e-4 && verbose) {
        warning(sprintf("Log-likelihood decreased at iteration %d (delta = %.2e)", iter, delta_ll))
      }
    }
    loglik_old <- loglik_new
    
    # Impute missing values using pattern-based approach
    data_imputed_new <- impute_missing_values(
      data = data,
      parameters = me_result$parameters,
      z = z_current,
      missing_pattern = missing_pattern
    )
    
    # Update for next iteration
    data_imputed <- data_imputed_new
  }
  
  # Final imputation
  if (sampling) {
    final_imputed <- impute_with_sampling(
      data = data,
      parameters = me_result$parameters,
      z = z_current,
      missing_pattern = missing_pattern,
      n_samples = n_samples 
    )
  } else {
    final_imputed <- impute_missing_values(
      data = data,
      parameters = me_result$parameters,
      z = z_current,
      missing_pattern = missing_pattern
    )
  }

  return(list(
    imputed_data = final_imputed,
    converged = converged,
    iterations = iter,
    loglik = loglik_new,
    parameters = me_result$parameters,
    responsibilities = z_current
  ))
}

regularize_cov <- function(sigma, epsilon = sqrt(.Machine$double.eps)) {
  eig <- eigen(sigma, symmetric = TRUE, only.values = TRUE)
  min_eig <- min(eig$values)
  if (min_eig < epsilon) {
    ridge <- epsilon - min_eig
    return(sigma + diag(ridge, nrow(sigma)))
  }
  return(sigma)
}

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
            
            # Regularize covariance
            sigma_obs_obs <- regularize_cov(sigma_obs_obs)
            sigma_obs_inv <- solve(sigma_obs_obs)
            
            # Compute conditional means
            obs_diff <- t(t(obs_data) - mu_obs) 
            cond_mean <- mu_miss + 
              (sigma_miss_obs %*% sigma_obs_inv) %*% 
              t(obs_diff)  # num_missing x num_obs 
            
            # Weight by responsibilities (proper dimension alignment)
            weights <- matrix(comp_responsibilities, nrow = num_obs, ncol = num_missing)
            weighted_cond <- t(cond_mean) * weights  # num_obs x num_missing
            
            imp_vals <- imp_vals + weighted_cond
          }
        }
        data_imputed[pat_indices, mis_indices] <- imp_vals
      } else {
        # Completely missing rows
        marg_mean <- compute_marginal_mean(parameters)
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

impute_with_sampling <- function(data, parameters, z, missing_pattern, 
                                 n_samples = 100, burn_in_ratio = 0.1) {
  n <- nrow(data)
  p <- ncol(data)
  G <- length(parameters$pro)
  data_imputed <- data
  burn_in <- max(1, floor(n_samples * burn_in_ratio))
  effective_samples <- n_samples - burn_in
  
  if (effective_samples < 1) {
    stop("n_samples too small after burn-in. Increase n_samples or decrease burn_in_ratio.")
  }
  
  pattern_str <- apply(missing_pattern, 1, function(x) paste(which(x), collapse = ","))
  unique_patterns <- unique(pattern_str)
  
  for (pat in unique_patterns) {
    pat_indices <- which(pattern_str == pat)
    if (length(pat_indices) == 0) next
    
    obs_indices <- !missing_pattern[pat_indices[1], ]
    mis_indices <- missing_pattern[pat_indices[1], ]
    num_missing <- sum(mis_indices)
    num_obs <- length(pat_indices)
    
    if (any(mis_indices) && any(obs_indices)) {
      imp_vals <- matrix(0, nrow = num_obs, ncol = num_missing)
      
      for (i in 1:num_obs) {
        idx <- pat_indices[i]
        obs_data <- data[idx, obs_indices]
        z_i <- z[idx, ]
        
        z_i[is.na(z_i)] <- 0
        z_i[z_i < 0] <- 0
        if (sum(z_i) < .Machine$double.eps) {
          z_i <- rep(1/G, G)  # Fallback to uniform
        } else {
          z_i <- z_i / sum(z_i) 
        }
      
        valid_g <- which(z_i > 1e-6)
        if (length(valid_g) == 0) valid_g <- 1:G
        
        all_samples <- matrix(0, n_samples, num_missing)
        for (s in 1:n_samples) {
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
    }
  }
  return(data_imputed)
}

compute_marginal_mean <- function(parameters) {
  G <- length(parameters$pro)
  p <- nrow(parameters$mean)
  marginal_mean <- rowSums(parameters$mean %*% diag(parameters$pro))
  return(marginal_mean)
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
    # Mu <- t(sweep(t(z_mat) %*% data, 1, n_k, "/"))
    Mu <- t(z_mat) %*% data / n_k # (p x K)
    return(list(prop      = prop,
         Mu        = Mu,
         SigmaCube = sigma0,    
         Z         = z_mat))}
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
      sigma_k <- regularize_cov(sigma_k)
      
      for (idx in 1:N_incomplete) {
        i <- incomplete_indices[idx]
        obs <- !is.na(data_imputed[i, ])
        if (any(obs)) {
          dens[idx, k] <- mvtnorm::dmvnorm(
            data_imputed[i, obs], 
            mean = mu_k[obs], 
            sigma = sigma_k[obs, obs, drop = FALSE], 
            log = TRUE
          )
        } else {
          dens[idx, k] <- 0
        }
      }
    }
    
    denspro <- sweep(dens, 2, log(temp_params$pro), "+")
    z_max <- apply(denspro, 1, max)
    log_sum_exp <- z_max + log(rowSums(exp(denspro - z_max)))
    z_incomplete <- exp(denspro - log_sum_exp)
    
    # Handle numerical issues
    z_incomplete[is.na(z_incomplete) | !is.finite(z_incomplete)] <- 1/G
    z_incomplete <- z_incomplete / rowSums(z_incomplete)
    
    z_full[incomplete_indices, ] <- z_incomplete
  }
  
  return(z_full)
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

map_to_mclust <- function(model_name, package_source = NULL) {
  if (toupper(model_name) %in% c("EII", "VII", "EEI", "VEI", "EVI", "VVI", 
                                 "EEE", "EVE", "VEE", "EEV", "VEV", "EVV", "VVV")) {
    return(toupper(model_name))
  }
  
  # Convert to lowercase for case-insensitive
  model_lower <- tolower(model_name)
  
  # MixAll model mapping
  mixall_mapping <- list(
    "gaussian_pk_sjk" = "VVI",
    "gaussian_pk_sj" = "VII",
    "gaussian_pk_sk" = "VEI",
    "gaussian_pk_s" = "EII",
    "gaussian_p_sjk" = "EVI",
    "gaussian_p_sj" = "VII",
    "gaussian_p_sk" = "EEI",
    "gaussian_p_s" = "EII"
  )
  
  # Rmixmod model mapping
  rmixmod_mapping <- list(
    "gaussian_p_l_i" = "EII",
    "gaussian_p_lk_i" = "VII",
    "gaussian_p_l_b" = "EEE",
    "gaussian_p_lk_b" = "VEE",
    "gaussian_p_l_bk" = "EEV",
    "gaussian_p_lk_bk" = "VEV",
    "gaussian_p_l_c" = "EVI",
    "gaussian_p_lk_c" = "VVI",
    "gaussian_p_l_d_ak_d" = "EEI",
    "gaussian_p_lk_d_ak_d" = "VEI",
    "gaussian_p_l_dk_a_dk" = "EVE",
    "gaussian_p_lk_dk_a_dk" = "VVE",
    "gaussian_p_l_dk_ak_dk" = "EVV",
    "gaussian_p_lk_dk_ak_dk" = "VVV",
    "gaussian_pk_l_i" = "EII",
    "gaussian_pk_lk_i" = "VII",
    "gaussian_pk_l_b" = "EEE",
    "gaussian_pk_lk_b" = "VEE",
    "gaussian_pk_l_bk" = "EEV",
    "gaussian_pk_lk_bk" = "VEV",
    "gaussian_pk_l_c" = "EVI",
    "gaussian_pk_lk_c" = "VVI",
    "gaussian_pk_l_d_ak_d" = "EEI",
    "gaussian_pk_lk_d_ak_d" = "VEI",
    "gaussian_pk_l_dk_a_dk" = "EVE",
    "gaussian_pk_lk_dk_a_dk" = "VVE",
    "gaussian_pk_l_dk_ak_dk" = "EVV",
    "gaussian_pk_lk_dk_ak_dk" = "VVV"
  )
  
  # Try MixAll mapping
  if (model_lower %in% names(mixall_mapping)) {
    return(mixall_mapping[[model_lower]])
  }
  
  # Try Rmixmod mapping
  if (model_lower %in% names(rmixmod_mapping)) {
    return(rmixmod_mapping[[model_lower]])
  }
  
  patterns <- list(
    c("pk_sjk", "VVI"), c("pk_sj", "VII"), c("pk_sk", "VEI"), 
    c("pk_s$", "EII"), c("p_sjk", "EVI"), c("p_sj", "VII"), 
    c("p_sk", "EEI"), c("p_s$", "EII"), c("lk_c$", "VVI"), 
    c("l_i$", "EII"), c("lk_i$", "VII"), c("l_b$", "EEE"),
    c("lk_b$", "VEE"), c("l_c$", "EVI"), c("lk_dk_ak_dk$", "VVV"),
    c("l_dk_ak_dk$", "VVV")  # Added pattern for "l_dk_ak_dk"
  )
  
  for (p in patterns) {
    if (grepl(p[1], model_lower, ignore.case = TRUE)) {
      return(p[2])
    }
  }
  
  # Default to VVV
  warning(paste("Unknown model name:", model_name, "- using VVV as default"))
  return("VVV")
}

get_component_covariance <- function(variance_struct, g) {
  if (!is.list(variance_struct)) stop("variance_struct must be a list")
  
  modelName <- variance_struct$modelName
  d <- variance_struct$d
  G <- variance_struct$G
  
  if (is.null(modelName) || !is.character(modelName) || length(modelName) != 1) 
    stop("Missing or invalid modelName")
  if (is.null(d) || d < 1 || d != as.integer(d)) 
    stop("Invalid dimension d")
  if (is.null(G) || G < 1 || G != as.integer(G)) 
    stop("Invalid component count G")
  if (g < 1 || g > G) 
    stop("Component index g out of range")
  
  if (!is.null(variance_struct$sigma) && 
      is.array(variance_struct$sigma) &&
      length(dim(variance_struct$sigma)) == 3) {
    dims <- dim(variance_struct$sigma)
    if (dims[1] == d && dims[2] == d && dims[3] >= g) {
      return(variance_struct$sigma[,,g])
    }
  }
  
  # univariate models
  if (d == 1) {
    if (modelName %in% c("E", "EII", "EEI")) {
      if (length(variance_struct$sigmasq) == 1) 
        return(matrix(variance_struct$sigmasq))
    }
    else if (modelName %in% c("V", "VII", "VVI")) {
      if (length(variance_struct$sigmasq) >= g) 
        return(matrix(variance_struct$sigmasq[g]))
    }
    stop("Univariate model parameters incomplete")
  }
  
  # Multivariate model
  switch(modelName,
         "EII" = {
           if (length(variance_struct$sigmasq) == 1)
             return(diag(variance_struct$sigmasq, d))
         },
         "VII" = {
           if (length(variance_struct$sigmasq) >= g)
             return(diag(variance_struct$sigmasq[g], d))
         },
         "EEI" =, "VEI" =, "EVI" =, "VVI" = {
           if (!is.null(variance_struct$scale) && !is.null(variance_struct$shape)) {
             scale_val <- if (modelName %in% c("EEI", "EVI")) 
               variance_struct$scale else variance_struct$scale[g]
             shape_vec <- if (modelName %in% c("EEI", "VEI")) 
               variance_struct$shape else variance_struct$shape[,g]
             return(diag(scale_val * shape_vec, d))
           }
         },
         "EEE" = {
           if (!is.null(variance_struct$Sigma))
             return(variance_struct$Sigma)
         },
         "VVV" = {
         },
         # Ellipsoidal models
         "EEV" =, "EVE" =, "EVV" =, "VEE" =, "VEV" =, "VVE" = {
           if (!is.null(variance_struct$scale) && 
               !is.null(variance_struct$shape) && 
               !is.null(variance_struct$orientation)) {
             # Corrected scale extraction
             scale_val <- if (modelName %in% c("EEV", "EVE", "EVV")) 
               variance_struct$scale else variance_struct$scale[g]
             shape_vec <- if (modelName %in% c("EVE", "VEE")) 
               variance_struct$shape else variance_struct$shape[,g]
             orientation <- if (modelName %in% c("EEV", "VEE")) 
               variance_struct$orientation else variance_struct$orientation[,,g]
             
             # Correct covariance reconstruction
             shape_mat <- diag(shape_vec, d)
             volume <- scale_val
             return(volume * orientation %*% shape_mat %*% t(orientation))
           }
         }
  )
  
  # Final fallback: identity matrix with warning
  warning("Using identity covariance fallback for component ", g)
  diag(1, d)
}

extract_model_name <- function(model_input) {
  # Direct mclust model names
  mclust_models <- c("EII", "VII", "EEI", "VEI", "EVI", "VVI", 
                     "EEE", "EVE", "VEE", "EEV", "VEV", "EVV", "VVV")
  if (toupper(model_input) %in% mclust_models) {
    return(toupper(model_input))
  }
  
  # Handle Rmixmod strings
  if (grepl("mixmodGaussianModel", model_input, ignore.case = TRUE)) {
    # Extract model from listModels parameter
    if (grepl('listModels\\s*=', model_input, ignore.case = TRUE)) {
      match <- regmatches(model_input, 
                          regexec('listModels\\s*=\\s*c\\(?\\s*"([^"]+)"', 
                                  model_input, ignore.case = TRUE))[[1]]
      if (length(match) >= 2) {
        return(tolower(match[2]))
      }
    }
    
    # Extract model from family parameter
    if (grepl('family\\s*=', model_input, ignore.case = TRUE)) {
      match <- regmatches(model_input, 
                          regexec('family\\s*=\\s*"([^"]+)"', 
                                  model_input, ignore.case = TRUE))[[1]]
      if (length(match) >= 2) {
        family <- tolower(match[2])
        family_map <- c(
          "spherical" = "gaussian_pk_l_i",
          "diagonal"  = "gaussian_pk_lk_c",
          "general"   = "gaussian_pk_lk_dk_ak_dk"
        )
        if (family %in% names(family_map)) {
          return(family_map[family])
        }
        return(family)
      }
    }
    
    # Default if parsing fails
    warning("Couldn't parse Rmixmod model: ", model_input)
    return("gaussian_pk_lk_c")
  }
  
  # Handle MixAll strings
  if (grepl("gaussian_p", model_input, ignore.case = TRUE)) {
    return(tolower(model_input))
  }
  
  # Default for unknown models
  warning("Unrecognized model format: ", model_input)
  "gaussian_pk_lk_c"
}

is_rmixmod_model <- function(model_string) {
  grepl("mixmodGaussianModel", model_string, ignore.case = TRUE)
}