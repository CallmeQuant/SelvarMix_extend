compute_nrmse <- function(original_data, 
                          imputed_data,
                          missing_data, 
                          normalization = "missing") {
  if (is.data.frame(missing_data) || is.matrix(missing_data))
    {missing_indices <- is.na(missing_data)}
  else {missing_indices <- missing_data}
  
  original_missing_values <- original_data[missing_indices]
  imputed_missing_values <- imputed_data[missing_indices]
  
  rmse <- sqrt(mean((original_missing_values - imputed_missing_values)^2))
  
  # Normalization factor
  if (normalization == "sd") {
    # Normalize by the standard deviation of the entire original data
    norm_factor <- sd(as.numeric(unlist(original_data)), na.rm = TRUE)
  } else if (normalization == "missing") {
    # Normalize by the square root of the mean squared actual values at missing positions
    norm_factor <- sqrt(mean(original_missing_values^2))
  } else {
    stop("Invalid normalization method. Use 'sd' or 'missing'.")
  }
  
  nrmse <- rmse / norm_factor
  
  return(nrmse)
}

compute_weighted_nrmse <- function(original_data, 
                                   missing_data, 
                                   imputed_data, 
                                   true_labels, 
                                   normalization = "missing") {
  true_labels <- as.factor(true_labels)
  clusters <- levels(true_labels)
  num_clusters <- length(clusters)
  
  total_nrmse <- 0
  total_weight <- 0
  
  for (cluster in clusters) {
    cluster_indices <- which(true_labels == cluster)
    
    missing_indices <- which(is.na(missing_data[cluster_indices, ]), arr.ind = TRUE)
    if (nrow(missing_indices) == 0) next  # Skip if no missing data in this cluster
    
    original_missing_values <- original_data[cluster_indices, 
    ][cbind(missing_indices[,1],
            missing_indices[,2])]
    imputed_missing_values <- imputed_data[cluster_indices, 
    ][cbind(missing_indices[,1], 
            missing_indices[,2])]
    
    rmse <- sqrt(mean((original_missing_values - imputed_missing_values)^2))
    
    # Normalization factor
    if (normalization == "sd") {
      norm_factor <- sd(as.numeric(original_data), na.rm = TRUE)
    } else if (normalization == "missing") {
      norm_factor <- sqrt(mean(original_missing_values^2))
    } else {
      stop("Invalid normalization method. Use 'sd' or 'missing'.")
    }
    
    nrmse <- rmse / norm_factor
    
    # Weight by the proportion of missing data in each cluster
    weight <- length(original_missing_values)
    
    total_nrmse <- total_nrmse + (nrmse * weight)
    total_weight <- total_weight + weight
    
  }
  weighted_nrmse <- total_nrmse / total_weight
  
  return(weighted_nrmse)
}

# Clustering-Integrated Imputation Error (CIIE)
compute_ciie <- function(original_data, 
                         missing_data, 
                         imputed_data, 
                         true_labels,
                         predicted_labels, 
                         normalization = "missing", 
                         similarity = "ARI",
                         alpha = 0.4,  # Weight for imputation accuracy
                         beta = 0.6   # Weight for clustering quality
) {
  if ((alpha + beta) != 1) {
    stop("Weights alpha and beta must sum to 1.")
  }
  
  nrmse <- compute_nrmse(original_data, imputed_data, missing_data, normalization)
  
  # invert so that lower error has higher score
  nrmse_score <- 1 - nrmse  
  nrmse_score <- ifelse(nrmse_score < 0, 0, nrmse_score)
  
  # Compute clustering similarity metrics
  similarity_metrics <- clustComp(true_labels, predicted_labels)
  
  # Extract desired similarity metrics
  ari <- similarity_metrics$ARI
  ami <- similarity_metrics$AMI
  nvi <- similarity_metrics$NVI
  nid <- similarity_metrics$NID
  nmi <- similarity_metrics$NMI
  
  # Invert NVI and NID to convert distance metrics to similarity metrics
  nvi_sim <- 1 - nvi
  nid_sim <- 1 - nid
  
  # Ensure all similarity metrics are within [0,1]
  similarity_metrics_processed <- c(
    ARI = max(min(ari, 1), 0),
    AMI = max(min(ami, 1), 0),
    NMI = max(min(nmi, 1), 0),
    NVI = max(min(nvi_sim, 1), 0),
    NID = max(min(nid_sim, 1), 0)
  )
  
  # Compute the average similarity score
  if (similarity == "all"){
    similarity_score <- mean(similarity_metrics_processed)}
  else {
    similarity_score <- similarity_metrics_processed[[similarity]]
  }
  
  # Compute CIIE as a weighted sum of NRMSE score and similarity score
  ciie <- (alpha * nrmse_score) + (beta * similarity_score)
  
  return(ciie)
}