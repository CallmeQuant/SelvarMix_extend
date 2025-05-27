spectral_distance <- function(Omega_hat_k0,
                              epsilon=1e-6,
                              laplacian_target_type = c("identity", "diag_Omega_hat"),
                              adj_threshold = 1e-4,
                              laplacian_norm_type = c("symmetric", "unsymmetric")) {

      laplacian_target_type <- match.arg(laplacian_target_type)
      laplacian_norm_type   <- match.arg(laplacian_norm_type)
      # P_k,ij = 1 / (||spec(L_k) - spec(L_target)||_2 + epsilon)
      # Helper function to get adjacency and Laplacian
      get_laplacian <- function(M, threshold, type) {
        adj_M <- abs(M) > threshold
        diag(adj_M) <- FALSE # No self-loops for graph Laplacian
        graph_M <- igraph::graph_from_adjacency_matrix(adj_M, mode = "undirected")
      
        if (type == "unsymmetric") {
          L <- igraph::laplacian_matrix(graph_M, normalization = "unnormalized", 
                                        sparse = FALSE)
        } else if (type == "symmetric") {
          L <- igraph::laplacian_matrix(graph_M, normalization = "symmetric", 
                                        sparse = FALSE)
          # igraph's normalized Laplacian is I - D^-1/2 A D^-1/2.
          # Ensure it's D^-1/2 L D^-1/2 if needed, 
          # or use properties that L_sym has eigenvalues in [0,2]
          # For D^-1/2 (D-A) D^-1/2 = I - D^-1/2 A D^-1/2, this is correct.
        } else {
          stop("Invalid laplacian_norm_type specified.")
        }
        return(L)
      }
      p <- nrow(Omega_hat_k0)
      
      # Laplacian for Omega_hat_k0
      L_k <- get_laplacian(Omega_hat_k0, threshold = adj_threshold, type = laplacian_norm_type)
      spec_L_k <- sort(eigen(L_k, symmetric = TRUE, only.values = TRUE)$values)
      
      # Laplacian for target matrix
      target_matrix_for_laplacian <- NULL
      if (laplacian_target_type == "identity") {
        target_matrix_for_laplacian <- diag(0, p) # Adjacency matrix for empty graph
      } else if (laplacian_target_type == "diag_Omega_hat") {
        # A graph based on diag(Omega_hat_k0) would have no off-diagonal edges
        target_matrix_for_laplacian <- diag(0, p) # Adjacency matrix for empty graph
      } else {
        stop("Invalid laplacian_target_type specified.")
      }
      L_target <- get_laplacian(target_matrix_for_laplacian, threshold = adj_threshold, type = laplacian_norm_type)
      spec_L_target <- sort(eigen(L_target, symmetric = TRUE, only.values = TRUE)$values)
      
      # Calculate Euclidean distance between sorted spectra
      dist_val <- if(length(spec_L_k) == length(spec_L_target) && length(spec_L_k) > 0) {
                    sqrt(sum((spec_L_k - spec_L_target)^2))
                 } else if (length(spec_L_k) == 0 && length(spec_L_target) == 0) {
                    0 
                 } else {
                    warning("Spectra length mismatch or empty spectra in spectral_distance. Defaulting dist_val to a large value (1/epsilon).")
                    1/epsilon 
                 }
      # Ensure distance is not too small
      dist_val <- max(dist_val, epsilon) 
      Pk <- matrix(1 / dist_val, p, p)
      return(Pk)
    }