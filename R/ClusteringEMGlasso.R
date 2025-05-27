source("R/spectral_distance.R")
ClusteringEMGlassoWeighted <- function(data,
                                       nbcluster,
                                       lambda,
                                       rho,
                                       group_shrinkage_method = c("common", 
                                       "weighted_by_W0", 
                                       "weighted_by_dist_to_I", 
                                       "weighted_by_dist_to_diag_W0",
                                       "laplacian_spectral"),
                                       distance_method = "Euclidean",
                                       lambda_omega_0 = 50, 
                                       epsilon_weighted_by_W0 = sqrt(.Machine$double.eps),
                                       penalize_diag = FALSE,
                                       laplacian_target_type = c("identity", "diag_Omega_hat"), # For laplacian_spectral
                                       adj_threshold = 1e-4,          # For laplacian_spectral
                                       laplacian_norm_type = c("symmetric", "unsymmetric"), # For laplacian_spectral
                                       initialize = c("kmeans", "hc"),
                                       nbcores = 1,
                                       n.start = 250)
{
  compute_Pk <- function(omega0_cube,
                         method           = "common",
                         distance_method = "Euclidean",
                         laplacian_target_type_pk = "diag_Omega_hat",
                         adj_threshold_pk       = 1e-4,
                         laplacian_norm_type_pk = "symmetric",
                         eps_w0                 = sqrt(.Machine$double.eps),
                         penalize_diag_pk       = FALSE){
    p <- dim(omega0_cube)[1]; K <- dim(omega0_cube)[3]

    out <- switch(method,
      common = array(1, dim = c(p, p, K)),

      weighted_by_W0 = {
        res <- 1 / (eps_w0 + abs(omega0_cube))
        res
      },  

      weighted_by_dist_to_I = {
        d <- apply(omega0_cube, 3L,
                  function(omega) shapes::distcov(omega, diag(p), method = distance_method))
        d[d < eps_w0] <- eps_w0
        array(rep(1/d, each = p*p), dim = c(p, p, K))
      },

      weighted_by_dist_to_diag_W0 = {
        d <- apply(omega0_cube, 3L,
                  function(omega){
                    diag_omega <- diag(diag(omega))
                    shapes::distcov(omega, diag_omega, method = distance_method)
                  })
        d[d < eps_w0] <- eps_w0
        array(rep(1/d, each = p*p), dim = c(p, p, K))
      },

      laplacian_spectral = {
        Pk_out_cube <- apply(omega0_cube, 3L,
                  function(omega){
                    spectral_distance(omega,
                                      epsilon = eps_w0,
                                      laplacian_target_type = laplacian_target_type_pk,
                                      adj_threshold = adj_threshold_pk,
                                      laplacian_norm_type = laplacian_norm_type_pk)
                  })

        dim(Pk_out_cube) <- c(p, p, K)
        Pk_out_cube
      },
      stop("unknown group_shrinkage_method")
    )

    if (!penalize_diag)
      for (k in seq_len(K)) diag(out[,,k]) <- 0

    out[!is.finite(out)] <- 1
    out
                       }
  #  Input Validation and Setup 
  data <- as.matrix(data)
  n <- as.integer(dim(data)[1])
  p <- as.integer(dim(data)[2])
  nbcluster <- as.integer(nbcluster)
  group_shrinkage_method <- match.arg(group_shrinkage_method)
  laplacian_target_type  <- match.arg(laplacian_target_type)
  adj_threshold = as.double(adj_threshold)
  laplacian_norm_type = match.arg(laplacian_norm_type)
  initialize <- match.arg(initialize)
  n.start <- as.integer(n.start)
  lambda_omega_0 <- as.double(lambda_omega_0)

  if (nbcores < 1) nbcores <- 1
  num_jobs <- length(lambda) * length(rho)
  if (num_jobs == 0) {
      warning("Lambda or Rho grid is empty. No models will be fit.")
      return(array(0, dim = c(0, p, length(nbcluster))))
  }
  if (num_jobs < nbcores) {
      nbcores <- num_jobs
  }

  is_windows <- Sys.info()["sysname"] == "Windows"
  cl <- NULL 

  # Setup cluster if needed
  if (nbcores > 1 && is_windows) {
      cl <- makeCluster(nbcores)
      clusterEvalQ(cl, {
          library(glassoFast)
          library(mclust)
          library(Rcpp)
          library(RcppArmadillo)
      })
      # Export Rcpp function if not automatically found via ::
      # clusterExport(cl, "rcppClusteringEMGlassoWeighted")
  }

  #  Parameter Initialization 
  message("Initializing parameters...")
  if (length(nbcluster) == 1) {
    junk_list <- list(InitParameter(data, nbcluster, init = initialize, n.start = n.start,
                                lambda_omega_0 = lambda_omega_0))
  } else {
    wrapper.init.parameter <- function(k) {
      library(glassoFast) 
      return(InitParameter(data, k, init = initialize, n.start = n.start, lambda_omega_0 = lambda_omega_0))
    }
    if (nbcores > 1 && is_windows) {
      clusterExport(cl = cl, varlist = c("InitParameter", "data", "initialize", "n.start", "lambda_omega_0"), envir = environment())
      clusterEvalQ(cl, {
                        library(glassoFast)
                        library(mclust)} )
      junk_list <- clusterApply(cl, x = as.integer(nbcluster), fun = wrapper.init.parameter)
    } else {
      junk_list <- mclapply(X = as.integer(nbcluster),
                            FUN = wrapper.init.parameter,
                            mc.cores = nbcores,
                            mc.preschedule = TRUE,
                            mc.cleanup = TRUE)
    }
  }
  message("Initialization complete.")


  #  Grid Setup 
  pen.grid <- as.matrix(expand.grid(lambda, rho))
  VarRole <- array(0, dim = c(nrow(pen.grid), p, length(nbcluster)),
                   dimnames = list(paste("L", pen.grid[,1], "R", pen.grid[,2]),
                                   colnames(data),
                                   paste("K", nbcluster)))

  #  Main Loop over K 
  for (k_idx in 1:length(nbcluster)) {
    K_current <- nbcluster[k_idx]
    message(sprintf("Processing K = %d...", K_current))

    # Get initial parameters list for this K
    P_init <- junk_list[[k_idx]]
    # Validate the structure returned by InitParameter
    if (!is.list(P_init) || length(P_init) != 5 ||
        !is.numeric(P_init[[1]]) || !is.matrix(P_init[[2]]) || !is.array(P_init[[3]]) ||
        !is.array(P_init[[4]]) || !is.matrix(P_init[[5]])) {
        warning(paste("InitParameter returned invalid list structure for K =", K_current, ". Skipping this K."))
        next # Skip to next K
    }
     # Check dimensions
    if(ncol(P_init[[2]]) != K_current || dim(P_init[[3]])[3] != K_current ||
       dim(P_init[[4]])[3] != K_current || ncol(P_init[[5]]) != K_current) {
        warning(paste("InitParameter returned list with inconsistent dimensions for K =", K_current, ". Skipping this K."))
        next
    }


    #  Calculate Pk_cube (once per K) 
    message(sprintf("Calculating Pk cube for K = %d using method: %s...", K_current, group_shrinkage_method))
    initial_omega_0 <- P_init[[4]]

    Pk_cube <- compute_Pk(omega0_cube         = initial_omega_0,
                          method           = group_shrinkage_method,
                          distance_method = distance_method,
                          laplacian_target_type_pk= laplacian_target_type,
                          adj_threshold_pk      = adj_threshold,
                          laplacian_norm_type_pk= laplacian_norm_type,
                          eps_w0                = epsilon_weighted_by_W0,
                          penalize_diag_pk      = penalize_diag)

    # Pk_cube <- compute_Pk(initial_omega_0,
    #                         method        = group_shrinkage_method,
    #                         distance_method = distance_method,
    #                         eps_w0          = epsilon_weighted_by_W0,
    #                         penalize_diag   = penalize_diag)

    # # Robustly calculate Pk_cube (similar logic to previous response)
    # Pk_cube <- tryCatch({
    #   # Check if initial_omega_0 is valid before using it
    #   if(is.null(initial_omega_0) || !is.array(initial_omega_0) || any(!is.finite(initial_omega_0))) {
    #       stop("Initial Omega_0 from InitParameter is invalid.") # Stop if init Omega is bad
    #   }
    #   switch (
    #     group_shrinkage_method,
    #     "common" = array(1, dim = c(p, p, K_current)),
    #     "weighted_by_W0" = {
    #         Pk_temp <- 1 / (epsilon_weighted_by_W0 + abs(initial_omega_0))
    #         Pk_temp[!is.finite(Pk_temp)] <- 1 # Handle potential Inf/NaN
    #         Pk_temp
    #     },
    #     "weighted_by_dist_to_I" = {
    #         dist_vals <- apply(initial_omega_0, 3, function(omega) {
    #             omega[!is.finite(omega)] <- 0
    #             omega <- omega + diag(epsilon_weighted_by_W0, p)
    #             if (!matrixcalc::is.positive.definite(omega, tol=epsilon_weighted_by_W0)) return(1) # Default dist
    #             res <- try(shapes::distcov(S1 = omega, S2 = diag(p), method = distance_method), silent=TRUE)
    #             if (inherits(res, "try-error")) return(1) # Default dist
    #             return(res)
    #         })
    #         dist_vals[dist_vals < epsilon_weighted_by_W0] <- epsilon_weighted_by_W0
    #         1 / array(rep(dist_vals, each = p ^ 2), dim = c(p, p, K_current))
    #     },
    #     "weighted_by_dist_to_diag_W0" = {
    #          dist_vals <- apply(initial_omega_0, 3, function(omega) {
    #             omega[!is.finite(omega)] <- 0
    #             omega <- omega + diag(epsilon_weighted_by_W0, p)
    #              if (!matrixcalc::is.positive.definite(omega, tol=epsilon_weighted_by_W0)) return(1)
    #             diag_omega <- diag(diag(omega))
    #             diag_omega[diag_omega <= epsilon_weighted_by_W0] <- epsilon_weighted_by_W0
    #             res <- try(shapes::distcov(S1 = omega, S2 = diag_omega, method = distance_method), silent=TRUE)
    #              if (inherits(res, "try-error")) return(1)
    #              return(res)
    #         })
    #         dist_vals[dist_vals < epsilon_weighted_by_W0] <- epsilon_weighted_by_W0
    #         1 / array(rep(dist_vals, each = p ^ 2), dim = c(p, p, K_current))
    #     }
    #   ) # end switch
    # }, error = function(e) {
    #     warning("Error calculating Pk_cube for K=", K_current, ": ", e$message, ". Defaulting to 'common'.")
    #     array(1, dim = c(p, p, K_current)) # Default to common if error
    # })

    # # Ensure Pk diagonal is zero if not penalizing diagonal
    # if (!penalize_diag) {
    #   for (k_inner in 1:K_current) { diag(Pk_cube[, , k_inner]) <- 0 }
    # }
    # Pk_cube[!is.finite(Pk_cube)] <- 1.0 # Ensure finite weights

    message("  Pk cube calculation complete.")


    #  Prepare InputList for Rcpp 
    # The list needs elements in the order expected by the C++ constructor
    # 0: Data, 1: prop, 2: Mu, 3: Covariance, 4: Precision, 5: ProbCond
    InputList_rcpp <- list(
        data,          # 0: Original data matrix
        P_init[[1]],   # 1: Initial proportions (pi_k)
        P_init[[2]],   # 2: Initial means (mu)
        P_init[[3]],   # 3: Initial covariance cube (Sigma_k)
        P_init[[4]],   # 4: Initial precision cube (Omega_k)
        P_init[[5]]    # 5: Initial posterior probabilities (z_ik)
    )

    #  Define Wrapper for Parallel Execution 
    wrapper.clusteringEMGlassoWeighted <- function(prm) {
      # prm is a row from pen.grid: prm[1] = lambda, prm[2] = rho
      result <- tryCatch(
          rcppClusteringEMGlassoWeighted(InputList_local, prm[1], prm[2], Pk_cube_local),
          error = function(e) {
              warning(sprintf("Rcpp function failed for K=%d, lambda=%.4f, rho=%.4f: %s",
                              K_local, prm[1], prm[2], e$message))
              rep(NA_integer_, p_local) 
          }
      )
      # Post-call validation
      if (!is.integer(result) || length(result) != p_local) {
           warning(sprintf("Rcpp function returned unexpected result for K=%d, lambda=%.4f, rho=%.4f. Returning NAs.",
                           K_local, prm[1], prm[2]))
           result <- rep(NA_integer_, p_local)
      }
      return(result)
    }

    #  Parallel Execution 
    message(sprintf("  Running penalized EM grid for K = %d (%d parameter pairs)...", K_current, nrow(pen.grid)))
    parallel.varrole_k <- list() # Results for the current K

    if (nbcores > 1 && is_windows) {
      #  Windows Parallel (parApply) 
      clusterExport(cl = cl,
                    varlist = c("InputList_rcpp", "Pk_cube", "p", "K_current", "rcppClusteringEMGlassoWeighted"),
                    envir = environment())
      # Define local variables for the wrapper function scope
      InputList_local <- InputList_rcpp
      Pk_cube_local <- Pk_cube
      p_local <- p
      K_local <- K_current
      clusterExport(cl=cl, varlist=c("InputList_local", "Pk_cube_local", "p_local", "K_local"), envir=environment())

      parallel.varrole_k_matrix <- parApply(cl,
                                            X = pen.grid,
                                            MARGIN = 1,
                                            FUN = wrapper.clusteringEMGlassoWeighted)
      # Transpose and handle errors
      if (is.matrix(parallel.varrole_k_matrix) && nrow(parallel.varrole_k_matrix) == p) {
          parallel.varrole_k <- t(parallel.varrole_k_matrix)
      } else {
          warning("parApply did not return a matrix[p, n_grid]. Results might be inconsistent.")
          # Coerce/fill with NA
          results_mat <- matrix(NA_integer_, nrow=nrow(pen.grid), ncol=p)
          if(is.list(parallel.varrole_k_matrix) && length(parallel.varrole_k_matrix) == nrow(pen.grid)) {
              for(i in 1:length(parallel.varrole_k_matrix)) {
                  if(is.integer(parallel.varrole_k_matrix[[i]]) && length(parallel.varrole_k_matrix[[i]]) == p) {
                      results_mat[i,] <- parallel.varrole_k_matrix[[i]]
                  }
              }
          } # Add other recovery logic if needed
          parallel.varrole_k <- results_mat
      }

    } else {
      #  Non-Windows Parallel (mclapply) or Single Core 
      InputList_local <- InputList_rcpp
      Pk_cube_local <- Pk_cube
      p_local <- p
      K_local <- K_current
      parallel.varrole_k_list <- mclapply(X = as.list(data.frame(t(pen.grid))),
                                          FUN = wrapper.clusteringEMGlassoWeighted,
                                          mc.cores = nbcores,
                                          mc.preschedule = TRUE,
                                          mc.cleanup = TRUE)
      # Combine list results into a matrix [n_grid x p]
      parallel.varrole_k <- do.call(rbind, parallel.varrole_k_list)
       # Handle errors
      if (!is.matrix(parallel.varrole_k) || ncol(parallel.varrole_k) != p) {
          warning("mclapply results inconsistent. Filling with NAs.")
          results_mat <- matrix(NA_integer_, nrow=nrow(pen.grid), ncol=p)
          if(is.list(parallel.varrole_k_list)) {
               for(i in 1:length(parallel.varrole_k_list)) {
                  if(is.integer(parallel.varrole_k_list[[i]]) && length(parallel.varrole_k_list[[i]]) == p) {
                      results_mat[i,] <- parallel.varrole_k_list[[i]]
                  } else {
                      results_mat[i,] <- rep(NA_integer_, p)
                  }
               }
          }
          parallel.varrole_k <- results_mat
      }
    } # End parallel execution choice
    message(sprintf("  Grid processing complete for K = %d.", K_current))

    # Store results, handling NAs
    var.role.k <- parallel.varrole_k
    var.role.k[is.na(var.role.k)] <- 0 # Replace NAs (from errors) with 0 score
    if(nrow(var.role.k) == nrow(VarRole[,,k_idx]) && ncol(var.role.k) == ncol(VarRole[,,k_idx])) {
        VarRole[, , k_idx] <- var.role.k
    } else {
        warning(paste("Dimension mismatch storing results K =", K_current))
    }

  } # End loop over K

  #  Cleanup 
  on.exit(if (!is.null(cl)) parallel::stopCluster(cl), add = TRUE)

  message("Variable ranking complete.")
  return(VarRole)
}



