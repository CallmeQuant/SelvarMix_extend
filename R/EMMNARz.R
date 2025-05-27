EMClustMNARz <- function(x,
                         K,
                         mecha     = "MNARz",
                         criterion = "BIC",
                         diag      = TRUE,
                         rmax      = 100,
                         init      = NULL,
                         tol       = 1e-4,
                         is_mnar   = NULL) {
  if (is.null(is_mnar)) {
    warning("is_mnar must be provided for mixed mechanism")
    is_mnar <- rep(TRUE, ncol(x))
  }

  if (mecha == "MNARz"){
    result <- tryCatch(
        EMGaussian(x, K, mecha, diag, rmax, init, tol),
        error = function(e) stop("EMGaussian failed: ", e$message)
      )
  }
  else if (mecha == "mixed"){
    result <- tryCatch(
        EMGaussianMixed(x, K, mecha, is_mnar, diag, rmax, init, tol),
        error = function(e) stop("EMGaussian failed: ", e$message)
      )
  }
  else {
    stop("Unknown mechanism '", mecha, "'")
  }

  loglik_final <- tail(result$loglik_vec, 1)
  if (!is.finite(loglik_final)) {
    stop("EMGaussian returned NA or infinite log‑likelihood; check covariance/mechanism settings")
  }

  n   <- nrow(x)
  d   <- ncol(x)
  tik <- result$tik

  ## partition
  partition <- apply(tik, 1, which.max)

  ## compute BIC and ICL
  if (mecha == "MNARz"){
    num_params <- (K - 1)       +  # mixing proportions
                  K * d         +  # cluster means
                  if (diag) K * d else K * d * (d + 1) / 2 + # covariances
                  K * d         # cluster-specific means
  } else if (mecha == "mixed") {
    num_params <-
                  (K - 1)         +           
                  K * d           +           
                  if(diag) K*d else K*d*(d+1)/2  
                  K * d           +     
                  K * sum(is_mnar) 
  }    

  bic  <- -2 * loglik_final + num_params * log(n)
  entropy <- -sum(tik * log(tik + 1e-10))
  icl     <- bic + 2 * entropy

  crit_list <- list(BIC = bic, ICL = icl)
  if (!criterion %in% names(crit_list)) {
    warning("Unknown criterion “", criterion, "”, defaulting to BIC")
    criterion <- "BIC"
  }
  if (mecha == "MNARz"){
    obj_ret <- list(
    loglik_obs     = loglik_final,
    partition      = partition,
    imputedData    = result$imputedData,
    criterionValue = crit_list,
    parameters     = list(
      pik   = result$pik,
      mu    = result$mu,
      sigma = result$sigma,
      alpha = result$alpha
    ),
    proba = tik) 
  }
  else {
    obj_ret <- list(
    loglik_obs     = loglik_final,
    partition      = partition,
    imputedData    = result$imputedData,
    criterionValue = crit_list,
    parameters     = list(
      pik   = result$pik,
      mu    = result$mu,
      sigma = result$sigma,
      alpha = result$alpha,
      beta  = result$beta        
    ),
    proba = tik)
  }
  return(obj_ret)
}
