.model_maps <- new.env(parent = emptyenv())

# Initialize model mappings on package load
.onLoad <- function(libname, pkgname) {
  
  # Rmixmod models
  rmixmod_models <- list()
  
  rmixmod_individual <- c(
    "Gaussian_p_L_I", "Gaussian_p_Lk_I", "Gaussian_p_L_B", "Gaussian_p_Lk_B",
    "Gaussian_p_L_Bk", "Gaussian_p_Lk_Bk", "Gaussian_p_L_C", "Gaussian_p_Lk_C",
    "Gaussian_p_L_D_Ak_D", "Gaussian_p_Lk_D_Ak_D", "Gaussian_p_L_Dk_A_Dk",
    "Gaussian_p_Lk_Dk_A_Dk", "Gaussian_p_L_Ck", "Gaussian_p_Lk_Ck",
    "Gaussian_pk_L_I", "Gaussian_pk_Lk_I", "Gaussian_pk_L_B", "Gaussian_pk_Lk_B",
    "Gaussian_pk_L_Bk", "Gaussian_pk_Lk_Bk", "Gaussian_pk_L_C", "Gaussian_pk_Lk_C",
    "Gaussian_pk_L_D_Ak_D", "Gaussian_pk_Lk_D_Ak_D", "Gaussian_pk_L_Dk_A_Dk",
    "Gaussian_pk_Lk_Dk_A_Dk", "Gaussian_pk_L_Ck", "Gaussian_pk_Lk_Ck"
  )
  
  for (model in rmixmod_individual) {
    tryCatch({
      rmixmod_models[[model]] <- Rmixmod::mixmodGaussianModel(listModels = model)
    }, error = function(e) {
      warning(paste("Failed to create Rmixmod model:", model, "-", e$message))
    })
  }
  
  # Family-based models
  rmixmod_families <- list(
    "all" = "all",
    "spherical" = "spherical", 
    "diagonal" = "diagonal",
    "general" = "general"
  )
  
  for (family_name in names(rmixmod_families)) {
    tryCatch({
      rmixmod_models[[family_name]] <- Rmixmod::mixmodGaussianModel(family = rmixmod_families[[family_name]])
    }, error = function(e) {
      warning(paste("Failed to create Rmixmod family model:", family_name, "-", e$message))
    })
  }
  
  # Mclust models
  mclust_models <- list(
    univariate = c("E", "V"),
    multivariate = c(
      "EII", "VII", "EEI", "VEI", "EVI", "VVI",
      "EEE", "EVE", "VEE", "VVE", "EEV", "VEV", "EVV", "VVV"
    ),
    all = c("E", "V", "EII", "VII", "EEI", "VEI", "EVI", "VVI",
            "EEE", "EVE", "VEE", "VVE", "EEV", "VEV", "EVV", "VVV")
  )
  
  # MixAll models
  mixall_models <- list(
    diagonal = c(
      "gaussian_pk_sjk", "gaussian_pk_sk", "gaussian_pk_skj", "gaussian_pk_s",
      "gaussian_p_sjk", "gaussian_p_sk", "gaussian_p_skj", "gaussian_p_s"
    ),
    all = c(
      "gaussian_pk_sjk", "gaussian_pk_sk", "gaussian_pk_skj", "gaussian_pk_s",
      "gaussian_p_sjk", "gaussian_p_sk", "gaussian_p_skj", "gaussian_p_s"
    )
  )
  
  # Store in global environment
  .model_maps$Rmixmod <- rmixmod_models
  .model_maps$Mclust <- mclust_models
  .model_maps$MixAll <- mixall_models
  .model_maps$RmixmodLookup <- names(rmixmod_models)
}