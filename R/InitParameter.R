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
  
  init    <- match.arg(init)
  
  # Initial partition 
  if (init == "kmeans") {
    km <- kmeans(data, centers = nbClust,
                 nstart = n.start, iter.max = 1000)
    memb <- km$cluster
    z_mat <- mclust::unmap(memb, G = nbClust)  
  } else {
    # Try SVD first
    tryCatch({
      hcobj <- mclust::hc(data, modelName = "VVV", use = "SVD")
      z_mat <<- mclust::unmap(mclust::hclass(hcobj, nbClust))
    }, error = function(e) {
      warning("mclust::hc with SVD failed: ", e$message, 
              ". Trying alternative approaches.")
      
      # Try PCR
      tryCatch({
        hcobj <- mclust::hc(data, modelName = "VVV", use = "PCR")
        z_mat <<- mclust::unmap(mclust::hclass(hcobj, nbClust))
      }, error = function(e) {
        warning("mclust::hc with PCR failed: ", e$message, 
                ". Trying alternative approaches.")
        
        # Try VARS
        tryCatch({
          hcobj <- mclust::hc(data, modelName = "VVV", use = "VARS")
          z_mat <<- mclust::unmap(mclust::hclass(hcobj, nbClust))
        }, error = function(e2) {
          warning("mclust::hc with VARS failed: ", e2$message, 
                  ". Trying standard hierarchical clustering.")
          
          # Try standard hierarchical clustering
          tryCatch({
            if (n > 1000) {
              sample_idx <- sample(1:n, min(1000, n))
              hc_result <- hclust(dist(data[sample_idx, ]), method = "ward.D2")
              cluster_centers <- cutree(hc_result, k = nbClust)
              
              centers <- matrix(0, nbClust, p)
              for (k in 1:nbClust) {
                if (sum(cluster_centers == k) > 0) {
                  centers[k, ] <- colMeans(data[sample_idx[cluster_centers == k], 
                                              , drop = FALSE])
                } else {
                  centers[k, ] <- data[sample(sample_idx, 1), ]
                }
              }
              distances <- as.matrix(dist(rbind(centers, data)))
              dist_to_centers <- distances[(nbClust + 1):(nbClust + n), 1:nbClust]
              memb <- apply(dist_to_centers, 1, which.min)
            } else {
              hc_result <- hclust(dist(data), method = "ward.D2")
              memb <- cutree(hc_result, k = nbClust)
            }
            z_mat <<- mclust::unmap(memb, G = nbClust)
          }, error = function(e3) {
            warning("All hierarchical clustering methods failed: ", e3$message,
                    ". Using k-means as final fallback.")
            
            # Final fallback: k-means
            tryCatch({
              km <- kmeans(data, centers = nbClust, nstart = n.start, iter.max = 1000)
              memb <- km$cluster
              z_mat <<- mclust::unmap(memb, G = nbClust)
            }, error = function(e4) {
              warning("K-means fallback failed. Using random assignment.")
              memb <- sample(1:nbClust, n, replace = TRUE)
              z_mat <<- mclust::unmap(memb, G = nbClust)
            })
          })
        })
      })
    })
  }
  n_k   <- colSums(z_mat)
  prop  <- pmax(n_k, 1e-6) / n  
  temp <- mclust::covw(data, z_mat, normalize = FALSE) 
  sigma0 <- temp$S                              
  omega0 <- array(NA_real_, dim = dim(sigma0))    
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
  Mu <- t( sweep(t(z_mat) %*% data, 1, pmax(n_k,1), "/") )   # p × K
  list(prop      = prop,
       Mu        = Mu,
       SigmaCube = sigma0,    
       OmegaCube = omega0,  
       Z         = z_mat)       # n × K responsibilities
}
