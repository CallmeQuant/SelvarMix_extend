CheckInputsC <- function(x, nbcluster, lambda, rho, type, hsize, criterion, models, rmodel, imodel, nbcores) {
  
  # Check if x is missing
  if (missing(x)) {
    stop("x is missing!")
  }
  
  # Check if x is a matrix or a data frame
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop(paste(sQuote("x"), "must be a matrix or data frame!"))
  }
  
  # Check if nbcluster is missing
  if (missing(nbcluster)) {
    stop("nbcluster is missing!")
  }
  
  # Check if nbcluster contains only positive integers
  if (sum(!is.wholenumber(nbcluster))) {
    stop("nbcluster must contain only integer values!")
  }
  
  if (sum(nbcluster < 1)) {
    stop(paste(sQuote("nbcluster"), "must be an integer greater than 0!"))
  }
  
  # Check if lambda is missing
  if (missing(lambda)) {
    stop("lambda is missing!")
  }
  
  # # Check if lambda is a vector of length >= 2
  # if (!is.vector(lambda) || length(lambda) <= 1) {
  #   stop(paste(sQuote("lambda"), "must be a vector with length >= 2!"))
  # }
  
  # if (sum(lambda <= 0)) {
  #   stop("lambda must be greater than 0!")
  # }
  
  # Check if rho is missing
  if (missing(rho)) {
    stop("rho is missing!")
  }
  
  # # Check if rho is a vector
  # if (!is.vector(rho)) {
  #   stop(paste(sQuote("rho"), "must be a vector!"))
  # }
  
  # if (sum(rho <= 0)) {
  #   stop("rho must be greater than 0!")
  # }
  
  # Check if type is valid and matches available types
  if (missing(type)) {
    stop("type is missing!")
  }
  
  if (!all(type %in% c("lasso", "likelihood")) || length(type) != 1) {
    stop(paste(type, "is not a valid type! Must be either 'lasso' or 'likelihood'."))
  }
  
  # Check if hsize is a valid positive integer <= number of columns of x
  if (!is.wholenumber(hsize) || sum(hsize < 1) || hsize > ncol(x)) {
    stop(paste(sQuote("hsize"), "must be a positive integer <= ncol(x)!"))
  }
  
  # Check if criterion is valid
  if (missing(criterion)) {
    criterion <<- "BIC"
  }
  
  if (!all(criterion %in% c("BIC", "ICL"))) {
    stop(paste(criterion, "is not a valid criterion! Must be either 'BIC' or 'ICL'."))
  }
  
  # Check if models match MixAll/Rmixmod/Mclust requirements 
  is_mixall_model <- models %in% c(
    "gaussian_pk_sjk", "gaussian_pk_sj", "gaussian_pk_sk", "gaussian_pk_s",
    "gaussian_p_sjk", "gaussian_p_sj", "gaussian_p_sk", "gaussian_p_s"
  )

  is_rmixmod_model <- grepl("^mixmodGaussianModel\\(", models)

  valid_mclust_models <- c(
    # Spherical models
    "EII", "VII", 
    
    # Diagonal models
    "EEI", "VEI", "EVI", "VVI", 
    
    # Ellipsoidal models
    "EEE", "VEE", "EVE", "VVE", 
    "EEV", "VEV", "EVV", "VVV"
  )

  is_mclust_model <- models %in% valid_mclust_models

  if (!is_mixall_model && !is_rmixmod_model && !is_mclust_model) {
    stop(paste("models must be either:",
              "- a valid MixAll model (e.g., 'gaussian_pk_sjk')",
              "- a valid Rmixmod model specification (e.g., 'mixmodGaussianModel(family=\"general\")')",
              "- a valid Mclust model (one of the following):",
              paste(valid_mclust_models, collapse=", "),
              sep = "\n"))
  }

  # For Rmixmod models, validate the family parameter if present
  if (is_rmixmod_model) {
    valid_families <- c("general", "diagonal", "spherical", "all")
    # Extract family parameter if present
    family <- "general" # Default 
    if (grepl('family="([^"]+)"', models, perl=TRUE)) {
      family <- gsub('.*family="([^"]+)".*', "\\1", models)
    }
    
    if (!family %in% valid_families) {
      stop(paste("For Rmixmod models, family must be one of:", paste(valid_families, collapse=", ")))
    }
  }

  # For Mclust models, validate against known model names
  if (is_mclust_model) {
    if (!models %in% valid_mclust_models) {
      stop(paste("For Mclust models, use one of the following:",
                paste(valid_mclust_models, collapse=", ")))
    }
  }
    
  # Check if rmodel is valid
  if (sum(!rmodel %in% c("LI", "LB", "LC"))) {
    stop(cat(rmodel[which(!(rmodel %in% c("LI", "LB", "LC")))], "is not a valid rmodel name!\n"))
  }
  
  # Check if imodel is valid
  if (sum(!imodel %in% c("LI", "LB"))) {
    stop(cat(imodel[which(!(imodel %in% c("LI", "LB")))], "is not a valid imodel name!\n"))
  }
  
  # Check if nbcores is a valid positive integer
  if (!is.wholenumber(nbcores) || nbcores < 1) {
    stop(paste(sQuote("nbcores"), "must be an integer > 0"))
  }
  
  # Check for missing values in x
  if (any(is.na(x))) {
    cat("Warning: x contains missing values. Missing values will be handled using a likelihood-based approach\n")
  }
}
