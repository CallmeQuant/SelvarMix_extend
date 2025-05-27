summary.selvarmixext <- function(object, ...) {
  x <- object
  if (class(x) != "selvarmixext") {
    stop(paste(sQuote("x"), sep = ""), " must be of class ", paste(dQuote("selvarmixext"), sep = ""), sep = "")
  }
  
  if (length(x) == 2) {
    for (i in 1:2) {
      cat("Criterion:", x[[i]]$criterion, "\n")
      cat("Criterion value:", x[[i]]$criterionValue, "\n")
      cat("Number of clusters:", x[[i]]$nbcluster, "\n")
      cat("Gaussian mixture model:", x[[i]]$model, "\n")
      if (!is.null(x[[i]]$imputedData)) {
        cat("Handling of missing values: Likelihood-based approach\n")
      }
      cat("Regression covariance model:", x[[i]]$rmodel, "\n")
      cat("Independent covariance model:", x[[i]]$imodel, "\n")
      cat("The SRUW model:\n")
      cat(" S:", x[[i]]$S, "\n")
      cat(" R:", x[[i]]$R, "\n")
      cat(" U:", x[[i]]$U, "\n")
      cat(" W:", x[[i]]$W, "\n")
    }
  } else {
    if (is.null(x$error)) {
      cat("Criterion:", x$criterion, "\n")
      cat("Criterion value:", x$criterionValue, "\n")
      cat("Number of clusters:", x$nbcluster, "\n")
      cat("Gaussian mixture model:", x$model, "\n")
      if (!is.null(x$imputedData)) {
        cat("Handling of missing values: Likelihood-based approach\n")
      }
      cat("Regression covariance model:", x$rmodel, "\n")
      cat("Independent covariance model:", x$imodel, "\n")
      cat("The SRUW model:\n")
      cat(" S:", x$S, "\n")
      cat(" R:", x$R, "\n")
      cat(" U:", x$U, "\n")
      cat(" W:", x$W, "\n")
    } else {
      cat("Criterion:", x$criterion, "\n")
      cat("Criterion value:", x$criterionValue, "\n")
      cat("Number of clusters:", x$nbcluster, "\n")
      cat("Gaussian mixture model:", x$model, "\n")
      cat("Prediction error:", round(x$error, 2), "\n")
      if (!is.null(x$imputedData)) {
        cat("Handling of missing values: Likelihood-based approach\n")
      }
      cat("Regression covariance model:", x$rmodel, "\n")
      cat("Independent covariance model:", x$imodel, "\n")
      cat("The SRUW model:\n")
      cat(" S:", x$S, "\n")
      cat(" R:", x$R, "\n")
      cat(" U:", x$U, "\n")
      cat(" W:", x$W, "\n")
    }
  }
}
