compute_grids_per_K <- function(
  X, nbcluster, P_method,
  distance_method     = "Euclidean",
  eps_w0              = sqrt(.Machine$double.eps),
  L                   = 5,
  frac_min            = 0.05,
  init_method         = c("kmeans", "hc"),
  n.start             = 250,
  lambda_omega_0      = 50,
  penalize_diag       = FALSE,
  laplacian_target_type = c("identity", "diag_Omega_hat"),
  adj_threshold       = 1e-4,
  laplacian_norm_type = c("symmetric", "unsymmetric")
) {
  stopifnot(L >= 2, frac_min > 0, frac_min < 1)
  init_method <- match.arg(init_method)
  laplacian_target_type <- match.arg(laplacian_target_type)
  laplacian_norm_type  <- match.arg(laplacian_norm_type)

  build_Pk <- function(Omega_cube) {
    dims <- dim(Omega_cube)
    p <- dims[1]; K <- dims[3]

    Pk <- switch(P_method,
      common = array(1, dim = dims),

      weighted_by_W0 = 1 / (eps_w0 + abs(Omega_cube)),

      weighted_by_dist_to_I = {
        d <- apply(Omega_cube, 3, function(Om)
          shapes::distcov(Om, diag(p), method = distance_method)
        )
        d[d < eps_w0] <- eps_w0
        array(rep(1/d, each = p*p), dim = dims)
      },

      weighted_by_dist_to_diag_W0 = {
        d <- apply(Omega_cube, 3, function(Om) {
          D <- diag(diag(Om))
          shapes::distcov(Om, D, method = distance_method)
        })
        d[d < eps_w0] <- eps_w0
        array(rep(1/d, each = p*p), dim = dims)
      },

      laplacian_spectral = {
        Pk_out <- array(0, dim = dims)
        for (k in seq_len(dims[3])) {
          Pk_out[,,k] <- spectral_distance(
            Omega_hat_k0         = Omega_cube[,,k],
            epsilon              = eps_w0,
            laplacian_target_type = laplacian_target_type,
            adj_threshold        = adj_threshold,
            laplacian_norm_type  = laplacian_norm_type
          )
        }
        Pk_out
      },

      stop("unknown P_method: ", P_method)
    )

    if (!penalize_diag) {
      for (k in seq_len(dims[3])) diag(Pk[,,k]) <- 0
    }
    Pk[!is.finite(Pk)] <- 1
    Pk
  }

  out <- vector("list", length(nbcluster))
  names(out) <- paste0("K", nbcluster)

  for (i in seq_along(nbcluster)) {
    K <- nbcluster[i]
    init <- InitParameter(
      X,               # data
      K,               # number of clusters
      init_method,     # init
      n.start,         # n.start
      lambda_omega_0   # lambda_omega_0
    )

    Z0     <- init$Z
    nk0    <- colSums(Z0)
    Sigma0 <- init$SigmaCube
    Omega0 <- init$OmegaCube

    Pk_cube <- build_Pk(Omega0)

    lam_max <- max(abs(t(Z0) %*% X))
    if (!is.finite(lam_max) || lam_max == 0) lam_max <- 1

    rho_max <- 0
    for (k in seq_len(K)) {
      Sk  <- Sigma0[,,k]
      tmp <- nk0[k] * abs(Sk) / Pk_cube[,,k]
      diag(tmp) <- 0
      rho_max <- max(rho_max, tmp[is.finite(tmp)])
    }
    if (!is.finite(rho_max) || rho_max == 0) rho_max <- 1

    geo <- function(maxv) {
      minv <- maxv * frac_min
      maxv * (minv / maxv) ^ (seq(0, L - 1) / (L - 1))
    }

    out[[i]] <- list(
      lambda_mu = geo(lam_max),
      rho       = geo(rho_max)
    )
  }

  out
}
