SortvarClust <- function(x,
                         nbcluster,
                         type="lasso", 
                         lambda=seq(20, 100, by = 10),
                         rho=seq(0.1, 1, length=5),
                         group_shrinkage_method = "common", 
                         distance_method = "Euclidean",
                         lambda_omega_0 = 50,
                         epsilon_weighted_by_W0 = sqrt(.Machine$double.eps),
                         laplacian_target_type = c("identity", "diag_Omega_hat"),
                         adj_threshold = 1e-4,
                         laplacian_norm_type = c("symmetric", "unsymmetric"),
                         penalize_diag = FALSE,
                         initialize = "hc",
                         nbcores=min(2, parallel::detectCores(logical = FALSE)), 
                         n.start = 250,
                         scale = FALSE,
                         ...)
{
  # --- Input Checks ---
  if(missing(x)){ stop("x is missing !") }
  if(!is.matrix(x) && !is.data.frame(x)) stop(sQuote("x"), " must be a matrix or data frame")
  x <- data.matrix(x) 
  if(any(!is.finite(x))) stop("Input data 'x' contains non-finite values (NA, NaN, Inf). Consider imputation or filtering.")

  if(missing(nbcluster)){ stop("nbcluster is missing!") }
  if(any(!is.wholenumber(nbcluster)) || any(nbcluster < 1)) stop(sQuote("nbcluster"), " must contain only positive integers!")

  if(!is.vector(lambda) || length(lambda) < 1) stop(sQuote("lambda"), " must be a vector with length >= 1")
  if (any(lambda < 0)) stop("lambda must be non-negative!")

  if(!is.vector(rho) || length(rho) < 1) stop(sQuote("rho"), " must be a vector with length >= 1")
  if(any(rho < 0)) stop("rho must be non-negative!")

  if(!is.wholenumber(nbcores) || nbcores < 1) stop(sQuote("nbcores"), " must be an integer > 0")

  # Validate group_shrinkage_method
  group_shrinkage_method <- tryCatch(
      match.arg(group_shrinkage_method, c("common", "weighted_by_W0", "weighted_by_dist_to_I", "weighted_by_dist_to_diag_W0", "laplacian_spectral")),
      error = function(e) stop("Invalid 'group_shrinkage_method'. Choose from 'common', 'weighted_by_W0', 'weighted_by_dist_to_I', 'weighted_by_dist_to_diag_W0', 'laplacian_spectral'.")
  )

  # Validate laplacian_method
  laplacian_target_type <- match.arg(laplacian_target_type)
  laplacian_norm_type  <- match.arg(laplacian_norm_type)
  # Validate intialize 
  initialize <- tryCatch(
      match.arg(initialize, c("kmeans", "hc")),
      error = function(e) stop("Invalid 'initialize'. Choose from 'kmeans', 'hc'.")
  )


  # Add checks for other new parameters if needed

  # --- Scaling ---
  if (scale){
    x <- scale(x, center = TRUE, scale = TRUE) # Scale data once
  }
  # Check if scaling produced NaNs (e.g., constant columns)
  if (any(!is.finite(x))) {
      stop("Scaling resulted in non-finite values. Check for constant columns in the input data.")
  }
  p <- as.integer(ncol(x))

  # --- Variable Ranking ---
  OrderVariable <- matrix(NA, nrow = length(nbcluster), ncol = p)
  rownames(OrderVariable) <- paste("K", nbcluster, sep="=")
  colnames(OrderVariable) <- colnames(x) # Assign variable names if available

  if(tolower(type) == "lasso")
  {
    VarRole <- ClusteringEMGlassoWeighted(
        data = x,
        nbcluster = nbcluster,
        lambda = lambda,
        rho = rho,
        group_shrinkage_method = group_shrinkage_method, 
        distance_method = distance_method,           
        lambda_omega_0 = lambda_omega_0,           
        epsilon_weighted_by_W0 = epsilon_weighted_by_W0, 
        laplacian_target_type = laplacian_target_type,
        adj_threshold = adj_threshold,
        laplacian_norm_type = laplacian_norm_type,
        penalize_diag = penalize_diag,    
        initialize = initialize,
        nbcores = nbcores,
        n.start = n.start
        )

    # Check if VarRole is valid before proceeding
    if (!is.array(VarRole) || length(dim(VarRole)) != 3) {
        stop("ClusteringEMGlassoWeighted did not return a valid 3D array.")
    }

    # Calculate scores and rank
    MatrixScores <- matrix(0, nrow = length(nbcluster), ncol = p)
    for (k_idx in 1:length(nbcluster)) {
      # Sum scores across the lambda/rho grid for this K
      MatrixScores[k_idx, ] <- colSums(VarRole[, , k_idx], na.rm = TRUE) # Sum scores, treating NA as 0
      # Get the order of variables based on scores (highest score first)
      OrderVariable[k_idx, ] <- order(MatrixScores[k_idx, ], decreasing = TRUE)
    }

  } else if(tolower(type) == "likelihood") {
      warning("Likelihood-based ordering not implemented.")
      # for(k in 1:length(nbcluster))
      #   OrderVariable[k,] <- orderlikC(x, nbcluster[k], nbcores) # Placeholder
  } else {
      stop("Invalid 'type' specified. Currently only 'lasso' is supported.")
  }

  return(OrderVariable) # Return matrix of ranked variable indices
}

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol