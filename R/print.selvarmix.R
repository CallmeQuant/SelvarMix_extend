print.selvarmixext <- function(object) {
  x <- object
  if (class(x) != "selvarmixext") {
    stop(paste(sQuote("x"), sep = ""), " must be of class ", paste(dQuote("selvarmixext"), sep = ""), sep = "")
  }
  
  if (length(x) == 2) {
    for (i in 1:2) {
      print(x[[i]]$parameters)
      if (!is.null(x[[i]]$missingValues)) {
        cat("Handling of missing values: Likelihood-based approach\n")
      }
      cat("Regression parameters:\n")
      print(x[[i]]$regparameters)
    }
  } else {
    print(x$parameters)
    if (!is.null(x$missingValues)) {
      cat("Handling of missing values: Likelihood-based approach\n")
    }
    cat("Regression parameters:\n")
    print(x$regparameters)
  }
}
