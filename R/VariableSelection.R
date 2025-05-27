VariableSelection <- function(data,
                              nbcluster,
                              models,
                              criterion,
                              OrderVariable,
                              hsize,
                              supervised,
                              z,
                              nbcores) {
  # Load required libraries
  requireNamespace("MixAll", quietly = TRUE)
  requireNamespace("Rmixmod", quietly = TRUE)
  requireNamespace("mclust", quietly = TRUE)
  requireNamespace("parallel", quietly = TRUE)

  data <- as.matrix(data)
  nbcluster <- as.integer(nbcluster)
  criterion <- as.character(criterion)
  p <- ncol(data)
  
  # Define framework and model_name
  mclust_models <-c(
                    # Spherical models
                    "EII", "VII", 
                    
                    # Diagonal models
                    "EEI", "VEI", "EVI", "VVI", 
                    
                    # Ellipsoidal models
                    "EEE", "VEE", "EVE", "VVE", 
                    "EEV", "VEV", "EVV", "VVV"
                    )
  if (models %in% mclust_models) {
    framework <- "Mclust"
    model_name <- models
  } else if (grepl("^mixmodGaussianModel\\(", models)) {
    framework <- "Rmixmod"
    model_name <- models
  } else {
    framework <- "MixAll"
    model_name <- models
  }
  
  # Wrapper function for parallel processing
  wrapper.selectVar <- function(current_nbcluster, current_ordervar) {
    result <- tryCatch({
      rcppSelectS(data,
                  current_ordervar,
                  current_nbcluster,
                  framework, model_name,
                  hsize, criterion,
                  as.integer(z), supervised) 
    }, error = function(e) {
       list(error = paste("Error in rcppSelectS call:", e$message),
            S = integer(0),
            W = integer(0),
            U = integer(0)) 
    })

    if (!is.null(result$error) || is.null(result$S)) {
       cat("Run failed (rcppSelectS):", result$error, "\n")
       return(list(S=integer(0), W=integer(0), U=1:p, error=result$error, nbcluster=current_nbcluster))
    }
    
    # Select W variables using rcppSelectW
    OrderAux <- setdiff(current_ordervar, result$S)
    result$W <- tryCatch(
      rcppSelectW(data, OrderAux, result$S, hsize),
      error = function(e) {
        cat("Error in rcppSelectW:", e$message, "\n")
        integer(0)
      }
    )
    
    # Determine U variables
    result$U <- setdiff(seq_len(ncol(data)), union(result$S, result$W))

    result$nbcluster_run <- current_nbcluster
    return(result)
  }

  max_cores <- parallel::detectCores(logical = FALSE) 
  # nbcores <- min(max(1L, as.integer(nbcores)), length(nbcluster), max_cores)
  nbcores <- min(max(1L, as.integer(nbcores)), max_cores)
  # Prepare input for processing
  input_list <- lapply(seq_along(nbcluster), function(i) {
    list(
      nbcluster = nbcluster[i],
      ordervar = if (is.matrix(OrderVariable) && nrow(OrderVariable) >= length(nbcluster)) OrderVariable[i, ] else if (is.matrix(OrderVariable)) OrderVariable[1, ] else OrderVariable
    )
  })
  # --- Parallel Execution ---
  results <- list()

  # if (nbcores > 1 && length(input_list) > 1) {
  if (nbcores > 1) {
    if (.Platform$OS.type == "windows") {
      cl <- tryCatch({
        parallel::makeCluster(nbcores)
      }, error = function(e) {
        warning("Failed to create cluster, running sequentially. Error: ", e$message, call. = FALSE)
        NULL
      })

      if (!is.null(cl)) {
        on.exit(parallel::stopCluster(cl), add = TRUE)
        tryCatch({
          parallel::clusterExport(cl, varlist = c("data", "framework", "model_name", "hsize", "criterion", "supervised", "z",
          "rcppSelectS", "rcppSelectW", "wrapper.selectVar"), envir = environment())

          # Load packages on workers
          parallel::clusterEvalQ(cl, {
             library(mclust)
             library(Rmixmod) 
             library(MixAll) 
           })

          results <- parallel::parLapply(cl, input_list, function(x) {
            tryCatch({
              wrapper.selectVar(x$nbcluster, x$ordervar)
            }, error = function(e) {
               list(error = paste("Error in parallel wrapper.selectVar:", e$message), nbcluster=x$nbcluster)
            })
          })

        }, error = function(e) {
           warning("Error during parallel execution setup/run on Windows, running sequentially. Error: ", e$message, call. = FALSE)
           parallel::stopCluster(cl)
           cl <- NULL
        })
      }
      if(is.null(cl)) nbcores <- 1

    } else { # Non-Windows
      results <- tryCatch({
        parallel::mclapply(input_list,
                         FUN = function(x) {
                            tryCatch({
                              wrapper.selectVar(x$nbcluster, x$ordervar)
                            }, error = function(e) {
                               list(error = paste("Error in parallel wrapper.selectVar:", e$message), nbcluster=x$nbcluster)
                            })
                         },
                         mc.cores = nbcores,
                         mc.silent = TRUE 
                         )
      }, error = function(e) {
         warning("Error during mclapply execution, running sequentially. Error: ", e$message, call. = FALSE)
         list()
      })
       if(length(results) == 0) nbcores <- 1 
    }
  } else {
     nbcores <- 1
  }


  # --- Sequential Execution---
  if (nbcores <= 1) {
    cat("Running sequentially...\n")
    results <- lapply(input_list, function(x) {
       tryCatch({
         wrapper.selectVar(x$nbcluster, x$ordervar)
       }, error = function(e) {
          list(error = paste("Error in sequential wrapper.selectVar:", e$message), nbcluster=x$nbcluster)
       })
    })
  }

  # --- Process Results ---
  VariableSelectRes <- list()
  valid_results_count <- 0

  for (ll in seq_along(results)) {
    current_result <- results[[ll]]
    current_nbcluster_val <- if(!is.null(current_result$nbcluster_run)) current_result$nbcluster_run else input_list[[ll]]$nbcluster

    if (!is.null(current_result$error)) {
      cat("Run for", current_nbcluster_val, "clusters failed:", current_result$error, "\n")
      next
    }
    if (is.null(current_result$S) || !is.numeric(current_result$S)) {
       cat("Run for", current_nbcluster_val, "clusters produced invalid 'S' component. Skipping.\n")
       next
    }

    # Add the result to the list
    VariableSelectRes[[length(VariableSelectRes) + 1]] <- list(
      S = current_result$S,
      W = if (!is.null(current_result$W)) current_result$W else integer(0), 
      U = if (!is.null(current_result$U)) current_result$U else integer(0),
      criterionValue = if (!is.null(current_result$criterionValue)) current_result$criterionValue else NA_real_,
      criterion = if (!is.null(current_result$criterion)) current_result$criterion else criterion, 
      model = if (!is.null(current_result$model)) current_result$model else model_name,     
      nbcluster = current_nbcluster_val,
      parameters = current_result$parameters,
      partition = current_result$partition, 
      proba = current_result$proba,      
      missingValues = current_result$missingValues
    )
    valid_results_count <- valid_results_count + 1
  }

  # Check if all runs failed
  if (valid_results_count == 0) {
    stop("All model fitting runs during variable selection failed. Please check data and parameters.", call. = FALSE)
  }

  # Sort results by nbcluster
  if (length(VariableSelectRes) > 1) {
     cluster_order <- order(sapply(VariableSelectRes, `[[`, "nbcluster"))
     VariableSelectRes <- VariableSelectRes[cluster_order]
  }


  return(VariableSelectRes)
}




#   if (framework == "Mclust"){
#     # cat("Using sequential processing for MClust.\n")
#     # results <- lapply(input_list, function(x) {
#     #                   tryCatch(
#     #                     wrapper.selectVar(x$nbcluster, x$ordervar),
#     #                     error = function(e) list(error = paste("Error in wrapper.selectVar:", e$message))
#     #                   )
#     #                 })
#     results <- lapply(input_list, function(x) {
#                     wrapper.selectVar(x$nbcluster, x$ordervar)
#                 })
#   } else { 
#     junk <- tryCatch({
#       # cat("Using parallel processing with", nbcores, "cores.\n")
#       if (.Platform$OS.type == "windows") {
#         cl <- makeCluster(nbcores)
#         common.objects <- c("data", "model_name", "framework", "hsize", "criterion", "supervised", "z", 
#                         "rcppSelectS", "rcppSelectW")
#         clusterExport(cl, varlist = common.objects, envir = environment())
        
#         clusterEvalQ(cl, {
#           library(MixAll)
#           library(Rmixmod)
#           library(mclust)
#         })
        
#         results <- parLapply(cl, input_list, function(x) {
#           tryCatch(
#             wrapper.selectVar(x$nbcluster, x$ordervar),
#             error = function(e) {
#               list(error = paste("Error in wrapper.selectVar:", e$message))
#             }
#           )
#         })
        
#         stopCluster(cl)
#         results
#       } else {
#         mclapply(
#           X = input_list,
#           FUN = function(x) {
#             tryCatch(
#               wrapper.selectVar(x$nbcluster, x$ordervar),
#               error = function(e) {
#                 list(error = paste("Error in wrapper.selectVar:", e$message))
#               }
#             )
#           },
#           mc.cores = nbcores,
#           mc.silent = FALSE,
#           mc.preschedule = TRUE,
#           mc.cleanup = TRUE
#         )
#       }
#     }, error = function(e) {
#       cat("Error in parallel processing:", e$message, "\n")
#       list(list(error = paste("Error in parallel processing:", e$message)))
#     })
#   }
  
#   # Initialize the result list
#   VariableSelectRes <- list()
#   idx <- 1
    
#   # Iterate over results
#   for (ll in seq_along(results)) {
#     if (!is.null(results[[ll]]$error)) {
#       cat("Error in run", ll, ":", results[[ll]]$error, "\n")
#       next
#     }
#     VariableSelectRes[[idx]] <- list(
#       S = results[[ll]]$S,
#       W = results[[ll]]$W,
#       U = results[[ll]]$U,
#       criterionValue = if (!is.null(results[[ll]]$criterionValue)) results[[ll]]$criterionValue else NA_real_,
#       criterion = if (!is.null(results[[ll]]$criterion)) results[[ll]]$criterion else criterion,
#       model = if (!is.null(results[[ll]]$model)) results[[ll]]$model else model_name,
#       nbcluster = if (!is.null(results[[ll]]$nbcluster)) results[[ll]]$nbcluster else input_list[[ll]]$nbcluster,
#       parameters = results[[ll]]$parameters,
#       partition = results[[ll]]$partition,
#       proba = results[[ll]]$proba,
#       missingValues = results[[ll]]$missingValues
#     )
#     idx <- idx + 1
#   }
  
#   # Check if all results failed
#   if (length(VariableSelectRes) == 0) {
#     cat("No successful results. All model fittings failed.\n")
#     stop("All model fittings failed. Please check your data and parameters.")
#   }
#   return(VariableSelectRes)
# }