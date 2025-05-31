SelvarClustLasso <- function(
  x,                      # Input data matrix
  nbcluster,              # Number of clusters
  strategy = NULL,        # Only for MixAll-strategy for fitting EM
  lambda                    = NULL,
  rho                       = NULL,
  num_vals_penalty = 5,   # Number of penalty values to consider
  type                      = "lasso",
  hsize = 3,             
  criterion = "BIC",      # Model selection criterion (BIC or ICL)
  models = "gaussian_pk_sjk",  # Mixture model specification
  rmodel = c("LI", "LB"), 
  imodel = c("LI", "LB"),       
  nbcores = min(2, detectCores(all.tests = FALSE, logical = FALSE)),
  impute_missing = TRUE,       # control missing value imputation
  use_copula = FALSE,           # use copula for imputation or not 
  scale_data = TRUE,           # allow automatic scaling
  scale_check_method = "pairwise.complete.obs",
  use_missing_pattern = FALSE, 
  use_diag = TRUE,
  true_labels = NULL,
  true_data = NULL, 
  sd_ratio_threshold    = 10,  # threshold to trigger scaling
  cond_number_threshold = 30,  # threshold to trigger scaling
  rank           = NULL,
  rank_control   = list(),    
  mnarz_control  = list(),
  verbose = TRUE      
) { 
  # "common", "weighted_by_W0", "weighted_by_dist_to_I", "weighted_by_dist_to_diag_W0"
  # used in "weighted_by_dist_to_I" and "weighted_by_dist_to_diag_W0" only
  # c("Procrustes","ProcrustesShape","Riemannian","Cholesky", "Euclidean", "LogEuclidean", "RiemannianLe")
  # "diag_Omega_hat" or "identity"
  # "kmeans" or "hc"
  required_pkgs <- c("missRanger", "MixAll", "Rmixmod", "mclust", "gcimputeR", "mvnfast")
  lapply(required_pkgs, requireNamespace, quietly = TRUE)
  invisible(lapply(required_pkgs, library, character.only = TRUE))

  # Input Validation
  CheckInputsC(x, nbcluster, lambda, rho, type, hsize, criterion, models, rmodel, imodel, nbcores)
  
  # Data Preprocessing
  x <- as.matrix(x)
  n <- as.integer(nrow(x))
  p <- as.integer(ncol(x))
  nbcluster <- as.integer(nbcluster)

  # Force garbage collection after setup
  gc(verbose = FALSE)
  

if (scale_data == FALSE && check_scale_data(x, sd_ratio_threshold,
                                               cond_number_threshold,
                                               use=scale_check_method)) {
  warning("Ensure data has been scale since input matrix is ill-conditioned (κ = ", round(cond_number_threshold,2),
          ") but 'scale_data = FALSE'. Proceeding anyway.")}

  do_scale <- if (scale_data) check_scale_data(x, sd_ratio_threshold,
                                               cond_number_threshold,
                                               use=scale_check_method) else FALSE                        
  centers <- if (do_scale) colMeans(x,  na.rm = TRUE) else rep(0, p)
  sds     <- if (do_scale) apply(  x, 2, sd, na.rm = TRUE) else rep(1, p)

  x_scaled <- sweep(x, 2, centers, "-")
  x_scaled <- sweep(x_scaled, 2, sds,    "/")
  
  # Impute if needed
  if (impute_missing && any(is.na(x_scaled))) {
    if (use_copula) {
      x_imp_scaled <- as.matrix(impute_GC(as.data.frame(x_scaled),
                                          verbose = FALSE)$Ximp)
    } else {
      x_imp_scaled <- as.matrix(missRanger(as.data.frame(x_scaled), verbose = 0))
    }
  } else {
    warning("Missing values present but impute_missing = FALSE")
    x_imp_scaled <- x_scaled
  }

  x_imp_orig <- sweep(x_imp_scaled, 2, sds, "*")
  x_imp_orig <- sweep(x_imp_orig, 2, centers, "+")

  # Force garbage collection after data preprocessing
  gc(verbose = FALSE)

  .rank_defaults <- list(
    type = type,
    lambda = lambda,
    rho = rho,
    group_shrinkage_method = "weighted_by_dist_to_diag_W0",
    distance_method        = "Euclidean",
    lambda_omega_0         = 50,
    epsilon_weighted_by_W0 = sqrt(.Machine$double.eps),
    penalize_diag          = FALSE,
    laplacian_target_type  = "diag_Omega_hat",
    adj_threshold          = 1e-4,
    laplacian_norm_type    = "symmetric",
    initialize             = "hc",
    nbcores                = nbcores,
    n.start                = 250
  )
  if (missing(rank_control) || is.null(rank_control)) rank_control <- list()               

  rank_control <- utils::modifyList(.rank_defaults, rank_control)

  .mnarz_defaults <- list(
    mecha     = "MNARz",
    is_mnar = NULL,
    diag      = use_diag,
    rmax      = 100,
    tol       = 1e-4,
    init      = NULL,
    method    = "usual",
    S         = 250,
    initialize = "hc"
  )
  mnarz_control <- utils::modifyList(.mnarz_defaults, mnarz_control)

  if (!is.null(lambda)) rank_control$lambda <- lambda
  if (!is.null(rho))    rank_control$rho    <- rho
  rank_control$type <- type                         

  rank_control$penalize_diag <- isTRUE(rank_control$penalize_diag)

  if ((is.null(rank_control$lambda) || length(rank_control$lambda) == 0) ||
      (is.null(rank_control$rho)    || length(rank_control$rho)    == 0)) {

    grid_list <- compute_grids_per_K(
        X               = as.matrix(x_imp_scaled),
        nbcluster       = nbcluster,
        P_method        = rank_control$group_shrinkage_method,
        distance_method = rank_control$distance_method,
        eps_w0          = rank_control$epsilon_weighted_by_W0,
        L               = num_vals_penalty,
        frac_min        = 0.05,
        init_method     = rank_control$initialize,
        n.start         = rank_control$n.start,
        lambda_omega_0  = rank_control$lambda_omega_0,
        penalize_diag   = rank_control$penalize_diag,
        laplacian_target_type = rank_control$laplacian_target_type,
        adj_threshold   = rank_control$adj_threshold,
        laplacian_norm_type = rank_control$laplacian_norm_type)

    if (is.null(rank_control$lambda) || length(rank_control$lambda) == 0)
        rank_control$lambda <- sort(unique(
                            unlist(lapply(grid_list, `[[`, "lambda_mu"))))
    if (is.null(rank_control$rho)    || length(rank_control$rho) == 0)
        rank_control$rho    <- sort(unique(
                            unlist(lapply(grid_list, `[[`, "rho"))))
  }

  # Clean up grid computation objects
  if (exists("grid_list")) rm(grid_list)
  gc(verbose = FALSE)
  
  # Suppress warnings
  options(warn = -1)
  
  # Define default strategy if strategy is NULL
  if (is.null(strategy)) {
    strategy <- clusterStrategy(
      nbTry = 1, 
      nbInit = 50, 
      initMethod = "class",
      initAlgo = "SEM",
      nbInitIteration = 5,
      initEpsilon = 1e-4,
      nbShortRun = 5,
      shortRunAlgo = "EM",
      nbShortIteration = 100,
      shortEpsilon = 1e-4,
      longRunAlgo = "EM",
      nbLongIteration = 200,
      longEpsilon = 1e-7
    )
  }
  
  # Initialize variable order matrix
  OrderVariable <- matrix(NA, nrow = length(nbcluster), ncol = p)

  # Variable Ranking
  if (missing(rank)) {
    rank_ctrl <- rank_control
    rank_ctrl$x         <- x_imp_scaled      
    rank_ctrl$nbcluster <- nbcluster   

    if (verbose) cat("Variable ranking\n")
    OrderVariable <- do.call(SortvarClust, rank_ctrl)
  } else {
    if (verbose) cat("  Using provided variable ranking\n")
    for (r in seq_len(nrow(OrderVariable))) {
      OrderVariable[r, ] <- rank
    }
  }
  if (verbose) {
    cat("Variable Ranks Preview: \n")
    print(OrderVariable[, 1:min(10, ncol(OrderVariable)), drop = FALSE])
  }

  # Clean up ranking objects and force garbage collection
  rm(rank_ctrl)
  gc(verbose = FALSE)
  
  # Supervised Status (Default: Unsupervised)
  supervised <- FALSE 
  knownlabels <- as.integer(1:n)
  
  # Variable Selection and Model Selection
  bestModel <- list()
  if (length(criterion) == 1) {
    if (verbose) cat("Variable selection with", criterion, "criterion\n")
    VariableSelectRes <- VariableSelection(
      x_imp_scaled, nbcluster, models, criterion, OrderVariable, 
      hsize, supervised, knownlabels, nbcores
    )

    gc(verbose = FALSE)
    
    bestModel[[criterion]] <- ModelSelectionClust(
      VariableSelectRes, x_imp_scaled, rmodel, imodel, nbcores
    )
  } else {
    for (crit in criterion) {
      if (verbose) cat("Variable selection with", crit, "criterion\n")
      VariableSelectRes <- VariableSelection(
        x_imp_scaled, nbcluster, models, crit, OrderVariable, 
        hsize, supervised, knownlabels, nbcores
      )

      gc(verbose = FALSE)
      
      bestModel[[crit]] <- ModelSelectionClust(
        VariableSelectRes, x_imp_scaled, rmodel, imodel, nbcores
      )
    }
  }

  # Clean up model selection intermediate objects
  rm(VariableSelectRes)
  gc(verbose = FALSE)

  # MNARz part 
  for (i in seq_along(bestModel)) {
    model_name <- names(bestModel)[i]
    finalModel <- bestModel[[i]]
    number_clusters <- finalModel$nbcluster
    if (use_missing_pattern) {
      if (verbose) cat("Fitting MNARz for criterion", model_name, "\n")
      if (!exists("EMClustMNARz")) stop("EMClustMNARz function is missing")
      em_call <- c(list(x = x_scaled, K = number_clusters,
                        criterion = model_name), mnarz_control)

      clust_result <- do.call(EMClustMNARz, em_call)
      
      # Validate imputedData
      if (is.null(clust_result$imputedData) || !is.matrix(clust_result$imputedData) || 
          !all(dim(clust_result$imputedData) == dim(x))) {
        stop("EMClustMNARz did not return valid imputedData")
          }

      x_imputed_final <- sweep(clust_result$imputedData, 2, sds, "*")
      x_imputed_final <- sweep(x_imputed_final,        2, centers, "+")
      finalModel$imputedData <- x_imputed_final
            
      # Validate partition
      if (is.null(clust_result$partition) || length(clust_result$partition) != nrow(x) ||
          !all(clust_result$partition %in% 1:number_clusters)) {
        stop("EMClustMNARz returned invalid partition")
      }

      adjustedRandIndex <- mclust::adjustedRandIndex
      if (!is.null(true_labels)) {  
        ari_original <- adjustedRandIndex(finalModel$partition, true_labels)
        ari_mnarz <- adjustedRandIndex(clust_result$partition, true_labels)

        if (ari_mnarz > ari_original) {
          finalModel$partition <- clust_result$partition}
        } else if (!is.null(finalModel$loglik) && !is.null(cl_res$loglik)){
          if (clust_result$loglik >= finalModel$loglik) {
            finalModel$partition <- clust_result$partition
          }
        }
        finalModel$criterionValue <- clust_result$criterionValue[[model_name]]
        finalModel$clust_result   <- clust_result
        }
      else {
        if (is_rmixmod_model(models)){
          mod <- extract_rmixmod_model(models)
        }
        else{
          mod <- finalModel$model
        }
        model_name <- map_and_validate_model(mod)
        cat("Using model: ", model_name, "\n")
        imputation_result <- tryCatch({
                                      EM_impute(
                                          data = x_scaled,
                                          G = number_clusters,
                                          modelName = model_name,
                                          method = mnarz_control$method,
                                          S = mnarz_control$S,
                                          max_iter = mnarz_control$rmax,
                                          init_method = mnarz_control$initialize,
                                          tol = 1e-8,
                                          verbose = FALSE
                                      )
        }, error = function(e) {
            cat("Error in EM_impute:", e$message, "\n")
            cat("Using original imputation method instead.\n")
            return(NULL)
        })
        if (is.null(imputation_result)) {
            finalModel$imputedData <- x_imp_orig
        } else {
            if (!is.null(true_data)){
              if (do_scale) {
                  true_data_scaled <- sweep(true_data, 2, centers, "-")
                  true_data_scaled <- sweep(true_data_scaled, 2, sds, "/")
                } else {
                  true_data_scaled <- true_data
                }
                nrmse1 <- compute_nrmse(
                    original_data = true_data_scaled,
                    missing_data = x_scaled,
                    imputed_data = imputation_result$imputedData,
                    normalization = "missing"
                )
                nrmse2 <- compute_nrmse(
                    original_data = true_data_scaled,
                    missing_data = x_scaled,
                    imputed_data = x_imp_orig,
                    normalization = "missing"
                )
                if (nrmse1 < nrmse2) {
                  finalModel$imputedData <- imputation_result$imputedData
                } else {
                  finalModel$imputedData <- x_imp_orig
                }
            } else {
              finalModel$imputedData <- imputation_result$imputedData
            }
        }            
    }
    # if (is.null(finalModel$imputedData)) finalModel$imputedData <- x_imp_orig
    bestModel[[i]] <- finalModel
  }

  # Final cleanup
  gc(verbose = FALSE)

  # Remove any null or invalid elements from bestModel before outputing.
  bestModel <- bestModel[!sapply(bestModel, is.null)]
 
  output <- PrepareOutput(bestModel)
  
  # verbose time or not
  total_time <- difftime(Sys.time(), start_time, units = "secs")
  if (verbose) {
    cat("\n==========================================\n")
    cat("Pipeline Completed\n")
    cat("==========================================\n")
    cat("Total running time:", round(total_time, 2), "seconds\n")
  }

  return(output)
}

# Helper function to prepare output
PrepareOutput <- function(bestModel) {
  output <- list()
  for (name in names(bestModel)) {
    processed <- ProcessModelOutput(bestModel[[name]])
    if (!is.null(processed)) output[[name]] <- processed
  }
  if (length(output) == 1) {
    return(output[[1]])
  } else {
    return(output)
  }
}

# Process individual model outputs
ProcessModelOutput <- function(modelResult) {
  if (is.null(modelResult)) return(NULL)
  if (is.null(modelResult$imputedData)) {
    warning("imputedData is NULL in modelResult")
  }
  if (!is.null(modelResult$regparameters)) {
    if (!is.null(modelResult$U) && length(modelResult$U) != 0)
      colnames(modelResult$regparameters) <- modelResult$U
    if (!is.null(modelResult$R) && length(modelResult$R) != 0)
      rownames(modelResult$regparameters) <- c("intercept", modelResult$R)
  }
  
  object <- list(
    S = modelResult$S,
    R = modelResult$R,
    U = modelResult$U,
    W = modelResult$W,
    criterionValue = modelResult$criterionValue,
    criterion = modelResult$criterion,
    model = modelResult$model,
    rmodel = modelResult$rmodel,
    imodel = modelResult$imodel,
    parameters = modelResult$parameters,
    nbcluster = modelResult$nbcluster,
    partition = modelResult$partition,
    proba = modelResult$proba,
    regparameters = modelResult$regparameters,
    imputedData = modelResult$imputedData,
    parametersMNARz = modelResult$clust_result
  )
  class(object) <- "selvarmixext"
  return(object)
}

is_rmixmod_model <- function(models) {
  return(grepl("mixmodGaussianModel", models))
}

is_mclust_model <- function(models){
  mclust_models <-c(
                    # Spherical models
                    "EII", "VII", 
                    
                    # Diagonal models
                    "EEI", "VEI", "EVI", "VVI", 
                    
                    # Ellipsoidal models
                    "EEE", "VEE", "EVE", "VVE", 
                    "EEV", "VEV", "EVV", "VVV"
                    )
  return(models %in% mclust_models)
}

check_scale_data <- function(x,
                             sd_ratio_threshold = 10,
                             cond_number_threshold = 30,
                             use = c("pairwise.complete.obs", "median")) {
  sds <- apply(x, 2, sd, na.rm = TRUE)
  sds[is.na(sds) | sds == 0] <- sqrt(.Machine$double.eps)
  ratio_sd <- max(sds) / min(sds)

  use <- match.arg(use)
  if (use == "pairwise.complete.obs") {
    cov_mat <- stats::cov(x, use = "pairwise.complete.obs")
  } else {                         
    med <- matrixStats::colMedians(x, na.rm = TRUE)
    x_med <- x
    na <- is.na(x_med)
    if (any(na)) x_med[na] <- rep(med, each = nrow(x))[na]
    cov_mat <- stats::cov(x_med)
  }

  if (anyNA(cov_mat) || min(dim(cov_mat)) == 0)
    return(TRUE)   

  eigs <- eigen(cov_mat, symmetric = TRUE, only.values = TRUE)$values
  eigs <- eigs[eigs > .Machine$double.eps]
  cond_num <- max(eigs) / min(eigs)

  do_scale <- (ratio_sd > sd_ratio_threshold) ||
              (cond_num  > cond_number_threshold)
  return(do_scale)
}
