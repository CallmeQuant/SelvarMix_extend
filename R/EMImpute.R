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
  
  # Create regularization wrapper for covariance matrices
  regularize_cov <- function(sigma, epsilon = sqrt(.Machine$double.eps)) {
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
