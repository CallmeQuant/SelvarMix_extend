selvarmixext <- function(bestModel) {
  if (length(bestModel) == 2) {
    output <- vector(mode = "list", length = 2)
    for (el in 1:2) {
      object <- list(S = bestModel[[el]]$S, 
                     R = bestModel[[el]]$R, 
                     U = bestModel[[el]]$U,
                     W = bestModel[[el]]$W,  
                     criterionValue = bestModel[[el]]$criterionValue, 
                     criterion = bestModel[[el]]$criterion, 
                     model = bestModel[[el]]$model, 
                     nbCluster = bestModel[[el]]$nbCluster, 
                     partition = bestModel[[el]]$partition, 
                     proba = bestModel[[el]]$proba,
                     missingValues = if (!is.null(bestModel[[el]]$missingValues)) bestModel[[el]]$missingValues else NULL)
      class(object) <- "selvarmixext"
      output[[el]] <- object
    }
  } else {
    output <- list(S = bestModel[[1]]$S, 
                   R = bestModel[[1]]$R, 
                   U = bestModel[[1]]$U,
                   W = bestModel[[1]]$W,  
                   criterionValue = bestModel[[1]]$criterionValue, 
                   criterion = bestModel[[1]]$criterion, 
                   model = bestModel[[1]]$model, 
                   nbCluster = bestModel[[1]]$nbCluster, 
                   partition = bestModel[[1]]$partition, 
                   proba = bestModel[[1]]$proba,
                   missingValues = if (!is.null(bestModel[[1]]$missingValues)) bestModel[[1]]$missingValues else NULL)
    class(output) <- "selvarmixext"
  }
  return(output)
}
