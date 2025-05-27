InitParameter <- function(data,
                          nbClust,
                          init = c("kmeans", "hc"),
                          n.start = 250,
                          lambda_omega_0 = 50,
                          epsilon_pd = sqrt(.Machine$double.eps))
{
  data    <- as.matrix(data)
  n <- as.integer(dim(data)[1])
  p <- as.integer(dim(data)[2])
  n.start <- as.integer(n.start)
  nbClust <- as.integer(nbClust)
  # small.pen <- as.double(small.pen)
  init    <- match.arg(init)

  # initial partition 
  if (init == "kmeans") {
    km <- kmeans(data, centers = nbClust,
                 nstart = n.start, iter.max = 1000)
    memb <- km$cluster
    z_mat <- mclust::unmap(memb, G = nbClust)  
  } else {            
    hcobj <- mclust::hc(data, modelName = "VVV", use = "SVD")           
    z_mat <- mclust::unmap(mclust::hclass(hcobj, nbClust))
  }
  n_k   <- colSums(z_mat)
  prop  <- pmax(n_k, 1e-6) / n  

  # weighted sample covariances 
  temp <- mclust::covw(data, z_mat, normalize = FALSE) 

  sigma0 <- temp$S                              
  omega0 <- array(NA_real_, dim = dim(sigma0))      # to be filled

  for (k in 1:nbClust) {
    eigv <- eigen(sigma0[,,k], symmetric = TRUE, only.values = TRUE)$values
    # If invertible, invert it
    # Otherwise, use glassoFast to estimate the precision matrix
    if (min(eigv) > epsilon_pd) {             
      omega0[,,k] <- solve(sigma0[,,k])
    } else {                                       
      gl <- glassoFast::glassoFast(
              S   = sigma0[,,k],
              rho = 2 * lambda_omega_0 / n_k[k],
              thr = 1e-4, maxIt = 1000)
      omega0[,,k] <- gl$wi
      sigma0[,,k] <- gl$w
    }
  }

  # initial means 
  Mu <- t( sweep(t(z_mat) %*% data, 1, pmax(n_k,1), "/") )   # p × K

  # result list 
  list(prop      = prop,
       Mu        = Mu,
       SigmaCube = sigma0,    
       OmegaCube = omega0,  
       Z         = z_mat)       # n × K responsibilities
}




# InitParameter <-
#   function(data,
#            nbClust,
#            n.start = 250,
#            small.pen = 0.5)
#   {
#     # data should ideally be scaled *before* calling this function if desired.
#     # Avoid scaling inside if the main function already does it.
#     # data <- as.matrix(scale(data, TRUE, FALSE)) # Scaling moved outside
#     data <- as.matrix(data)
#     n <- as.integer(dim(data)[1])
#     p <- as.integer(dim(data)[2])
#     n.start <- as.integer(n.start)
#     nbClust <- as.integer(nbClust)
#     small.pen <- as.double(small.pen)

#     if (nbClust <= 0) stop("nbClust must be positive.")
#     if (n < nbClust) stop("Number of observations less than number of clusters.")

#     Mu <- matrix(0, p, nbClust)

#     # K-means for initial partition
#     clust <- tryCatch(
#         kmeans(data, nbClust, iter.max = 1000, nstart = n.start), # Removed algorithm choice for simplicity
#         error = function(e) {
#             warning("kmeans failed during initialization: ", e$message, ". Attempting random assignment.")
#             # Fallback: random assignment
#             list(cluster = sample(1:nbClust, n, replace = TRUE), size = tabulate(sample(1:nbClust, n, replace = TRUE), nbins=nbClust))
#         }
#     )
#     memb <- clust$cluster
#     prop <- clust$size / n
#     prop[prop == 0] <- 1e-6 # Avoid zero proportions
#     prop <- prop / sum(prop)

#     # Calculate initial means and empirical covariances
#     S <- array(0, dim = c(p, p, nbClust))
#     min_obs_per_cluster = p + 1 # Minimum needed for stable cov estimate

#     for (k in 1:nbClust)
#     {
#       obs_in_k <- which(memb == k)
#       n_k <- length(obs_in_k)

#       if (n_k > 0) {
#           Mu[, k] <- colMeans(data[obs_in_k, , drop = FALSE])
#       } else {
#           # Handle empty cluster - assign random mean (e.g., from overall mean + noise)
#           Mu[, k] <- colMeans(data) + rnorm(p, 0, sd(as.vector(data))/10) # Random perturbation
#           warning("Initialization: Cluster ", k, " is empty. Assigning random mean.")
#       }

#       if (n_k >= min_obs_per_cluster) {
#         S[, , k] <- cov(data[obs_in_k, , drop = FALSE])
#         # Add small ridge for stability before glassoFast
#         S[, , k] <- S[, , k] + diag(sqrt(.Machine$double.eps), p)
#       } else {
#         # Handle clusters with too few points for cov calculation
#         warning("Initialization: Cluster ", k, " has ", n_k, " points (< p+1=", min_obs_per_cluster,"). Using pooled covariance estimate.")
#         # Use pooled covariance or identity matrix as fallback
#         S[, , k] <- cov(data) + diag(sqrt(.Machine$double.eps), p)
#         if(any(!is.finite(S[,,k]))) S[,,k] <- diag(p) # Final fallback
#       }
#        if(any(!is.finite(S[,,k]))) S[,,k] <- diag(p) # Ensure S is finite
#     }

#     # Estimate initial W (Sigma) and Wi (Omega) using glassoFast
#     W <- Wi <- array(0, dim = c(p, p, nbClust))
#     for (k in 1:nbClust)
#     {
#       # Ensure S is positive definite enough for glassoFast
#       eig_S <- try(eigen(S[, , k], symmetric = TRUE, only.values = TRUE)$values, silent = TRUE)
#       if (inherits(eig_S, "try-error") || min(eig_S) <= sqrt(.Machine$double.eps)) {
#           S[, , k] <- S[, , k] + diag(sqrt(.Machine$double.eps) * max(1, abs(diag(S[,,k]))), p) # Add adaptive ridge
#       }

#       gg <- tryCatch(
#           glassoFast::glassoFast(S = S[, , k], rho = small.pen, thr = 1e-4, maxIt = 1000),
#           error = function(e) {
#               warning("glassoFast failed for initial estimate (k=", k,"): ", e$message, ". Using identity matrices.")
#               list(w = diag(p), wi = diag(p)) # Fallback
#           }
#       )
#       # Check output validity
#       if(is.null(gg) || !is.list(gg) || is.null(gg$w) || is.null(gg$wi) || any(!is.finite(gg$w)) || any(!is.finite(gg$wi))) {
#            warning("glassoFast returned invalid result for initial estimate (k=", k,"). Using identity matrices.")
#            W[, , k] <- diag(p)
#            Wi[, , k] <- diag(p)
#       } else {
#           W[, , k] <- gg$w
#           Wi[, , k] <- gg$wi
#       }
#     }

#     # Prepare the list in the order expected by Rcpp
#     # 0: Data, 1: prop, 2: Mu, 3: Covariance, 4: Precision, 5: ProbCond (initial z)
#     P <- list()
#     # P$X <- data # Rcpp function takes data separately now
#     P[[1]] <- prop # Index 1 (R) -> Element 1 (List) -> Arg 1 (Rcpp)
#     P[[2]] <- Mu
#     P[[3]] <- W   # CovarianceMatrix
#     P[[4]] <- Wi  # PrecisionMatrix
#     P[[5]] <- mclust::unmap(memb, G = nbClust) # Convert cluster vector to n x K matrix


#     return(P)
#   }
