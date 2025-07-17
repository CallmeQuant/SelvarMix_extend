EMClustMNARz <- function(x,
                         K,
                         mecha     = "MNARz",
                         criterion = "BIC",
                         diag      = TRUE,
                         rmax      = 100,
                         init      = NULL,
                         tol       = 1e-4,
                         is_mnar   = NULL) {
  
  if (mecha == "mixed" && is.null(is_mnar)) {
    stop("`is_mnar` must be provided as a logical vector for `mecha = 'mixed'`.")
  }

  if (mecha == "MNARz"){
    result <- tryCatch(
        EMGaussian(x, K, mecha, diag, rmax, init, tol),
        error = function(e) stop("EMGaussian for MNARz failed: ", e$message)
      )
  } else if (mecha == "mixed"){
    result <- tryCatch(
        EMGaussianMixed(x, K, mecha, is_mnar, diag, rmax, init, tol),
        error = function(e) stop("EMGaussianMixed for mixed mechanism failed: ", e$message)
      )
  } else {
    stop("Unknown mechanism '", mecha, "'. Supported mechanisms are 'MNARz' and 'mixed'.")
  }

  loglik_final <- tail(result$loglik_vec, 1)
  if (!is.finite(loglik_final)) {
    stop("EM algorithm returned a non-finite log-likelihood. Check model parameters or data scaling.")
  }

  n   <- nrow(x)
  d   <- ncol(x)
  tik <- result$tik

  # Generate partition from posterior probabilities
  partition <- apply(tik, 1, which.max)

  # Compute number of parameters for BIC/ICL
  if (mecha == "MNARz"){
    # pi_k, mu_k, sigma_k, alpha_k (one alpha per cluster)
    num_params <- (K - 1)       +  # mixing proportions
                  K * d         +  # cluster means
                  (if (diag) K * d else K * d * (d + 1) / 2) + # covariances
                  K             # MNARz: one alpha_k per cluster
  } else if (mecha == "mixed") {
    # pi_k, mu_k, sigma_k, alpha_kj, beta_kj (for MNAR vars)
    num_params <- (K - 1)         +
                  K * d           +
                  (if(diag) K*d else K*d*(d+1)/2) +
                  K * d           +     # alpha_kj for all variables
                  K * sum(is_mnar)      # beta_kj only for MNAR variables
  }

  bic  <- -2 * loglik_final + num_params * log(n)
  # Add a small epsilon to log to avoid log(0)
  entropy <- -sum(tik * log(tik + 1e-20), na.rm = TRUE)
  icl     <- bic + 2 * entropy

  crit_list <- list(BIC = bic, ICL = icl)
  if (!criterion %in% names(crit_list)) {
    warning("Unknown criterion '", criterion, "', defaulting to BIC.")
    criterion <- "BIC"
  }
  
  # Structure the final output object
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
      beta  = if (!is.null(result$beta)) result$beta else NULL
    ),
    proba = tik
  )
  
  return(obj_ret)
}
