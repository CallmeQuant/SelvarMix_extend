ModelSelectionClust <- function(VariableSelectRes,
                                data,
                                rmodel,
                                imodel,
                                nbcores) {
  
  mylist.size <- length(VariableSelectRes)
  if (mylist.size == 1) {
    junk <- try(rcppCrit(data, VariableSelectRes[[1]], rmodel, imodel), silent = TRUE)
    # junk <- list(junk)
    # if (inherits(junk, "try-error")) stop("Model evaluation failed")
  } else {
    wrapper.rcppCrit <- function(idx) {
      mylist <- VariableSelectRes[[idx]]
      res <- rcppCrit(data, mylist, rmodel, imodel)
      return(res)
    }
    
    if (mylist.size < nbcores) 
      nbcores <- mylist.size
    
    if (.Platform$OS.type == "windows") {
      cl <- makeCluster(nbcores)
      common.objects <- c("data", "VariableSelectRes", "rmodel", "imodel")
      clusterExport(cl = cl, varlist = common.objects, envir = environment())
      junk <- clusterApply(cl, 
                           x = seq_len(mylist.size), 
                           fun = wrapper.rcppCrit)
      stopCluster(cl)
    } else {
      junk <- mclapply(X = seq_len(mylist.size), 
                       FUN = wrapper.rcppCrit,
                       mc.cores = nbcores,
                       mc.silent = TRUE,
                       mc.preschedule = TRUE,
                       mc.cleanup = TRUE)
    }
  } 
  
  # Processing the results
  if ((mylist.size == 1) && (class(junk) != "try-error")) {
    bestModel <- junk
  } else { 
    lmax <- -Inf
    for (idx in seq_along(junk)) {
      if (inherits(junk[[idx]], "try-error")) {
          cat("Error in model", idx, ":", conditionMessage(attr(junk[[idx]], "condition")), "\n")
          next
        }
        
      if (is.null(junk[[idx]]$criterionValue) || is.na(junk[[idx]]$criterionValue)) {
        cat("Model", idx, "has null/NA criterion value\n")
        next
      }
        
      if (junk[[idx]]$criterionValue > lmax) {
        bestModel <- junk[[idx]]
        lmax <- bestModel$criterionValue 
      }
      }
      }
      
  if (is.null(bestModel)) {
    stop("No valid models found. All model evaluations failed.")
  }
  
  if (length(bestModel$R) == 0) {
    bestModel$R <- NULL
    bestModel$W <- c(bestModel$U, bestModel$W)
    bestModel$U <- NULL
  }
  
  if (length(bestModel$W) == 0)
    bestModel$W <- NULL
  
  return(bestModel)
}
