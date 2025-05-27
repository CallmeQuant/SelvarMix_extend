orderlikC <- function(x, g, nbcores) {
  
  loglikwrapper <- function(j, x, g) {
    library(MixAll)
    quick_precise_strategy <- clusterStrategy(
                                              nbTry = 1,
                                              nbInit = 5,
                                              initMethod = "class",
                                              initAlgo = "SEM",
                                              nbInitIteration = 25,
                                              initEpsilon = 1e-4,
                                              nbShortRun = 5,
                                              shortRunAlgo = "EM",
                                              nbShortIteration = 120,
                                              shortEpsilon = 1e-6,
                                              longRunAlgo = "EM",
                                              nbLongIteration = 150,
                                              longEpsilon = 1e-7
                                            )
    result <- clusterDiagGaussian(data = x[, j, drop = FALSE],
                                  strategy = quick_precise_strategy,
                                  nbCluster = g,
                                  models = "gaussian_pk_sjk")
    return(result@lnLikelihood)
  }
  
  return(order(simplify2array(mclapply(X = 1:ncol(x), 
                                       FUN = loglikwrapper, 
                                       x = x,
                                       g = g,
                                       mc.preschedule = TRUE, 
                                       mc.silent = TRUE, 
                                       mc.cleanup = TRUE, 
                                       mc.cores = min(nbcores, ncol(x)))), 
               decreasing = TRUE))
}