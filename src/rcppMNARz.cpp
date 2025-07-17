#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Helper function to mimic R's which()
int which_index(const std::vector<int>& vec, int value) {
  auto it = std::find(vec.begin(), vec.end(), value);
  if(it != vec.end()) {
    return std::distance(vec.begin(), it);
  }
  return -1; // Not found
}

// Calculates the observed log-likelihood and posterior probabilities
// [[Rcpp::export]]
List LoglikelihoodObsGaussian(NumericMatrix YNA, List mu, List sigma, NumericMatrix alpha, NumericVector prop_pi) {
  int n = YNA.nrow();
  int d = YNA.ncol();
  int K = prop_pi.size();
  
  NumericMatrix C(n, d);
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < d; j++) {
      C(i, j) = R_IsNA(YNA(i, j)) ? 1.0 : 0.0;
    }
  }
  
  // Create unique patterns of missingness
  std::vector<NumericVector> uniquePatterns;
  std::vector<IntegerVector> patternIndices;
  
  for(int i = 0; i < n; i++) {
    NumericVector currentPattern = C(i, _);
    bool found = false;
    
    for(size_t p = 0; p < uniquePatterns.size(); p++) {
      bool match = true;
      for(int j = 0; j < d; j++) {
        if(currentPattern[j] != uniquePatterns[p][j]) {
          match = false;
          break;
        }
      }
      
      if(match) {
        patternIndices[p].push_back(i);
        found = true;
        break;
      }
    }
    
    if(!found) {
      uniquePatterns.push_back(currentPattern);
      IntegerVector newIndices;
      newIndices.push_back(i);
      patternIndices.push_back(newIndices);
    }
  }
  
  NumericMatrix logprobcond(n, K);
  
  // Process each pattern
  for(size_t p = 0; p < uniquePatterns.size(); p++) {
    NumericVector pattern = uniquePatterns[p];
    IntegerVector indices = patternIndices[p];
    
    // Get observed variables for this pattern
    std::vector<int> var_obs;
    for(int j = 0; j < d; j++) {
      if(pattern[j] == 0) var_obs.push_back(j);
    }
    
    int d_obs = var_obs.size();
    
    // For each observation with this pattern
    for(int idx = 0; idx < indices.size(); idx++) {
      int i = indices[idx];
      
      // Extract observed values
      NumericVector y_obs(d_obs);
      for(int j = 0; j < d_obs; j++) {
        y_obs[j] = YNA(i, var_obs[j]);
      }
      
      // For each cluster
      for(int k = 0; k < K; k++) {
        // Extract cluster parameters for observed variables
        NumericVector mu_obs(d_obs);
        NumericMatrix sigma_obs(d_obs, d_obs);
        
        NumericVector mu_k = as<NumericVector>(mu[k]);
        NumericMatrix sigma_k = as<NumericMatrix>(sigma[k]);
        
        for(int j = 0; j < d_obs; j++) {
          mu_obs[j] = mu_k[var_obs[j]];
          for(int l = 0; l < d_obs; l++) {
            sigma_obs(j, l) = sigma_k(var_obs[j], var_obs[l]);
          }
        }
        
        // Calculate log density of observed data (term1_log)
        double log_density = 0.0;
        if(d_obs > 0) {
          // Convert to arma for multivariate normal density calculation
          // arma::vec y_obs_arma = as<arma::vec>(y_obs);
          arma::vec y_obs_arma(d_obs);
          for(int j = 0; j < d_obs; j++) {
            if(R_IsNA(y_obs[j])) {
              y_obs_arma[j] = 0.0;
            } else {
              y_obs_arma[j] = y_obs[j];
            }
          }

          arma::vec mu_obs_arma = as<arma::vec>(mu_obs);
          arma::mat sigma_obs_arma = as<arma::mat>(sigma_obs);
          
          // Check if matrix is positive definite
          if(!sigma_obs_arma.is_sympd()) {
            sigma_obs_arma = symmatu(sigma_obs_arma);  // Force symmetry
            sigma_obs_arma.diag() += 1e-6;            // Add small value to diagonal
          }
          
          // Calculate multivariate normal density
          double sign = 0.0;
          double logdet = 0.0;
          arma::log_det(logdet, sign, sigma_obs_arma);
          
          arma::vec diff = y_obs_arma - mu_obs_arma;
          arma::mat invSigma = inv_sympd(sigma_obs_arma);
          double mahalanobis = as_scalar(diff.t() * invSigma * diff);
          
          log_density = -0.5 * (d_obs * log(2 * M_PI) + logdet + mahalanobis);
        }
        
        // Calculate term1_log
        double term1_log = log(prop_pi[k]) + log_density;
        
        // Calculate term2_log (missingness mechanism)
        double term2_log = 0.0;
        for(int j = 0; j < d; j++) {
          double pnorm_alpha = R::pnorm(alpha(k, j), 0.0, 1.0, 1, 0);
          pnorm_alpha = std::max(std::min(pnorm_alpha, 1.0 - 1e-12), 1e-12);
          
          if(pattern[j] == 0) {
            // Observed value
            if(pnorm_alpha >= 0.9999) {
              // return List::create(
              //   Named("error") = "Degenerescence (mechanism): loglikelihood is not computable"
              // );
              Rcpp::stop("Degenerate in missing‐data mechanism: log‐likelihood not computable");
            }
        
            term2_log += std::log1p(-pnorm_alpha);
          } else {
            // Missing value
            term2_log += std::log(pnorm_alpha);
          }
        }
        
        // Combine terms
        logprobcond(i, k) = term1_log + term2_log;
      }
    }
  }
  
  // Normalize to avoid numerical issues
  NumericVector normval(n);
  for(int i = 0; i < n; i++) {
    normval[i] = max(logprobcond(i, _));
    for(int k = 0; k < K; k++) {
      logprobcond(i, k) = exp(logprobcond(i, k) - normval[i]);
    }
  }
  
  // Calculate final log-likelihood and posterior probabilities
  double loglik_obs = sum(normval) + sum(log(rowSums(logprobcond)));
  if (!std::isfinite(loglik_obs)) {
    Rcpp::stop("loglik_obs became NaN or infinite; numerical degeneracy detected");
  }
  NumericMatrix tik = clone(logprobcond);
  for(int i = 0; i < n; i++) {
    double row_sum = sum(tik(i, _));
    for(int k = 0; k < K; k++) {
      tik(i, k) /= row_sum;
    }
  }
  
  return List::create(
    Named("loglik_obs") = loglik_obs,
    Named("tik") = tik
  );
}

// Initialize parameters for EM algorithm
// [[Rcpp::export]]
List InitEMGaussian(NumericMatrix YNA, int K, std::string mecha, bool diag, Nullable<List> init, Nullable<int> samplesize) {
  int n = YNA.nrow();
  int d = YNA.ncol();
  
  // Default initialization
  if(init.isNull()) {
    // Declare and initialize Z_init
    NumericMatrix Z_init(n, K);
    
    // Random class assignments using Rcpp's R interface
    Function sample = Environment::base_env()["sample"];
    // Create proper R arguments
    IntegerVector values = seq_len(K);
    IntegerVector assignments = as<IntegerVector>(sample(values, n, true)) - 1;
    
    // Assign one-hot encoding with safety check
    for(int i = 0; i < n; i++) {
      if(assignments[i] < 0 || assignments[i] >= K) {
        assignments[i] = 0; // Fallback to first cluster if invalid
      }
      Z_init(i, assignments[i]) = 1.0;
    }
    
    // Calculate class proportions
    NumericVector prop_pi(K);
    for(int k = 0; k < K; k++) {
      prop_pi[k] = sum(Z_init(_, k)) / n;
    }
    
    // Calculate means and covariances for each class
    List mu_init(K);
    List sigma_init(K);
    
    for(int k = 0; k < K; k++) {
      // Find observations in class k
      IntegerVector class_obs;
      for(int i = 0; i < n; i++) {
        if(Z_init(i, k) == 1.0) class_obs.push_back(i);
      }
      
      // Calculate means
      NumericVector mu_k(d);
      for(int j = 0; j < d; j++) {
        double sum = 0.0;
        int count = 0;
        for(int idx = 0; idx < class_obs.size(); idx++) {
          int i = class_obs[idx];
          if(!R_IsNA(YNA(i, j))) {
            sum += YNA(i, j);
            count++;
          }
        }
        mu_k[j] = count > 0 ? sum / count : 0.0;
      }
      mu_init[k] = mu_k;
      
      // Calculate covariances
      NumericMatrix sigma_k(d, d);
      if(diag) {
        // Diagonal covariance
        for(int j = 0; j < d; j++) {
          double sum_sq_diff = 0.0;
          int count = 0;
          for(int idx = 0; idx < class_obs.size(); idx++) {
            int i = class_obs[idx];
            if(!R_IsNA(YNA(i, j))) {
              double diff = YNA(i, j) - mu_k[j];
              sum_sq_diff += diff * diff;
              count++;
            }
          }
          sigma_k(j, j) = count > 1 ? sum_sq_diff / count : 1.0;
        }
      } else {
        // Full covariance
        for(int j1 = 0; j1 < d; j1++) {
          for(int j2 = 0; j2 <= j1; j2++) {
            double sum_cross = 0.0;
            int count = 0;
            for(int idx = 0; idx < class_obs.size(); idx++) {
              int i = class_obs[idx];
              if(!R_IsNA(YNA(i, j1)) && !R_IsNA(YNA(i, j2))) {
                double diff1 = YNA(i, j1) - mu_k[j1];
                double diff2 = YNA(i, j2) - mu_k[j2];
                sum_cross += diff1 * diff2;
                count++;
              }
            }
            double cov_val = count > 1 ? sum_cross / count : (j1 == j2 ? 1.0 : 0.0);
            sigma_k(j1, j2) = cov_val;
            sigma_k(j2, j1) = cov_val; // Symmetry
          }
        }
      }
      sigma_init[k] = sigma_k;
    }
    
    // Initialize alpha for missingness mechanism
    NumericMatrix alpha_init(K, d);
    if(mecha == "MCAR") {
      double miss_rate = 0.0;
      int total = n * d;
      for(int i = 0; i < n; i++) {
        for(int j = 0; j < d; j++) {
          if(R_IsNA(YNA(i, j))) miss_rate += 1.0;
        }
      }
      miss_rate /= total;
      
      double alpha_val = R::qnorm(miss_rate, 0.0, 1.0, 1, 0);
      for(int k = 0; k < K; k++) {
        for(int j = 0; j < d; j++) {
          alpha_init(k, j) = alpha_val;
        }
      }
    } else {
      for(int k = 0; k < K; k++) {
        for(int j = 0; j < d; j++) {
          alpha_init(k, j) = 0.0;
        }
      }
    }
    
    return List::create(
      Named("pi_init") = prop_pi,
      Named("mu_init") = mu_init,
      Named("sigma_init") = sigma_init,
      Named("alpha_init") = alpha_init
    );
  } else {
    List init_list = as<List>(init);
    NumericMatrix alpha_user;
    if (init_list.containsElementNamed("alpha")) {
      alpha_user = as<NumericMatrix>(init_list["alpha"]);
    } else {
      alpha_user = NumericMatrix(K, d);
      
      double miss_rate = 0.0;
      int total_cells = n * d;
      for(int i = 0; i < total_cells; ++i) if(R_IsNA(YNA[i])) miss_rate++;
      miss_rate /= total_cells;

      double alpha_val = R::qnorm(miss_rate, 0.0, 1.0, 1, 0);
      alpha_user.fill(alpha_val);
    }
    return List::create(
      _["pi_init"]    = as<NumericVector>(init_list["pik"]),
      _["mu_init"]    = as<List>(init_list["mu"]),
      _["sigma_init"] = as<List>(init_list["sigma"]),
      _["alpha_init"] = alpha_user
    );
  } 
}

// [[Rcpp::export]]
NumericMatrix MechanismEMGLM(NumericMatrix YNA, NumericMatrix tik, std::string mecha) {
  int n = YNA.nrow();
  int d = YNA.ncol();
  int K = tik.ncol();
  
  // Create missingness indicator matrix C
  NumericMatrix C(n, d);
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < d; j++) {
      C(i, j) = R_IsNA(YNA(i, j)) ? 1.0 : 0.0;
    }
  }
  
  NumericMatrix alpha_new(K, d);
  
  // Get the binomial family function with probit link from R
  Environment stats = Environment::namespace_env("stats");
  Function glm = stats["glm"];
  Function binomial = stats["binomial"];
  RObject family = binomial(Rcpp::Named("link", "probit"));
  
  if(mecha == "MNARz") {
    // Flatten C into a vector for MNARz
    NumericVector C_flat(n * d);
    for(int i = 0; i < n; i++) {
      for(int j = 0; j < d; j++) {
        C_flat[i * d + j] = C(i, j);
      }
    }
    
    for(int k = 0; k < K; k++) {
      // Prepare weights by repeating tik[,k] for each variable
      NumericVector weights(n * d);
      for(int i = 0; i < n; i++) {
        for(int j = 0; j < d; j++) {
          weights[i * d + j] = tik(i, k);
        }
      }
      
      // Call R's glm function for probit regression
      List glm_fit = glm(
        _["formula"] = Formula("y ~ 1"),
        _["family"] = family,
        _["data"] = DataFrame::create(_["y"] = C_flat),
        _["weights"] = weights
      );
      
      NumericVector coef = glm_fit["coefficients"];
      double alpha_k = coef[0];  // Intercept from GLM fit
      
      // Assign the same alpha to all variables in cluster k
      for(int j = 0; j < d; j++) {
        alpha_new(k, j) = alpha_k;
      }
    }
  } else if(mecha == "MNARzj") {
    // Fit a separate GLM for each variable j
    for(int j = 0; j < d; j++) {
      NumericVector y = C(_, j);  // Missingness indicator for variable j
      
      for(int k = 0; k < K; k++) {
        NumericVector weights = tik(_, k);  // Weights for cluster k
        
        // Call R's glm function for probit regression
        List glm_fit = glm(
          _["formula"] = Formula("y ~ 1"),
          _["family"] = family,
          _["data"] = DataFrame::create(_["y"] = y),
          _["weights"] = weights
        );
        
        NumericVector coef = glm_fit["coefficients"];
        alpha_new(k, j) = coef[0];  // Intercept for this variable and cluster
      }
    }
  }
  
  return alpha_new;
}

// [[Rcpp::export]]
List EMGaussian(NumericMatrix YNA, int K, std::string mecha, bool diag, int rmax, 
                Nullable<List> init = R_NilValue,
                double tol = 0.0001, 
                Nullable<int> samplesize = R_NilValue) {
  int n = YNA.nrow();
  int d = YNA.ncol();
  
  // Missing data indicator matrix
  NumericMatrix C(n, d);
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < d; j++) {
      C(i, j) = R_IsNA(YNA(i, j)) ? 1.0 : 0.0;
    }
  }
  
  // Initialize parameters
  List init_params;
  if (init.isNull()) {
    init_params = InitEMGaussian(YNA, K, mecha, diag, R_NilValue, samplesize);
  } else {
    init_params = InitEMGaussian(YNA, K, mecha, diag, init, samplesize);
  }
  
  NumericVector pi_new = init_params["pi_init"];
  List mu_new = init_params["mu_init"];
  List sigma_new = init_params["sigma_init"];
  NumericMatrix alpha_new = init_params["alpha_init"];
  
  // Ensure covariance matrices are positive definite
  for(int k = 0; k < K; k++) {
    // arma::mat sigma_arma = as<arma::mat>(sigma_new[k]);
    // if(!sigma_arma.is_sympd()) {
    //   sigma_arma = symmatu(sigma_arma);
    //   sigma_arma.diag() += 1e-6;
    //   sigma_new[k] = wrap(sigma_arma);
    // }
    // Implement SPD approximation
    arma::mat S = as<arma::mat>(sigma_new[k]);
    S = 0.5 * (S + S.t());
    
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, S);

    double delta = arma::datum::eps * arma::trace(S) / S.n_rows;
    for(unsigned i = 0; i < eigval.n_elem; ++i)
      if(eigval[i] < delta) eigval[i] = delta;
    arma::mat S_pd = eigvec * arma::diagmat(eigval) * eigvec.t();
    sigma_new[k] = wrap(S_pd);
  }
  
  // Initial log-likelihood
  List loglik_res = LoglikelihoodObsGaussian(YNA, mu_new, sigma_new, alpha_new, pi_new);
  if(loglik_res.containsElementNamed("error")) return loglik_res;
  
  double loglik_prev = as<double>(loglik_res["loglik_obs"]);
  double prec = -std::numeric_limits<double>::infinity();
  NumericVector loglik_vec = NumericVector::create(loglik_prev);
  NumericMatrix tik = as<NumericMatrix>(loglik_res["tik"]);
  
  int r = 0;

  // For later on imputation
  List mu_tilde_list(n);
  List sigma_tilde_list(n);

  // Initialize lists to prevent null entries
  for(int i = 0; i < n; i++) {
    mu_tilde_list[i] = List(K);
    sigma_tilde_list[i] = List(K);
    
    for(int k = 0; k < K; k++) {
      List mu_tilde_i = mu_tilde_list[i];
      List sigma_tilde_i = sigma_tilde_list[i];
      mu_tilde_i[k] = NumericVector(d);
      sigma_tilde_i[k] = NumericMatrix(d, d);
      mu_tilde_list[i] = mu_tilde_i;
      sigma_tilde_list[i] = sigma_tilde_i;
    }
  }
  
  // EM loop
  while(r < rmax && (loglik_prev - prec > tol)) {
    r++;
    
    // E-step: Compute conditional expectations and variances
    for(int i = 0; i < n; i++) {
      // Process E-step
      std::vector<int> indobs, indmiss;
      for(int j = 0; j < d; j++) {
        if(C(i, j) == 0) indobs.push_back(j);
        else indmiss.push_back(j);
      }
      
      for(int k = 0; k < K; k++) {
        NumericVector mu_k = as<NumericVector>(mu_new[k]);
        NumericMatrix sigma_k = as<NumericMatrix>(sigma_new[k]);
        NumericVector mu_tilde_ik(d);
        NumericMatrix sigma_tilde_ik(d, d);
        
        if(indmiss.empty()) {
          for(int j = 0; j < d; j++) mu_tilde_ik[j] = YNA(i, j);
          
          // Set sigma_tilde_ik to zeros for the fully observed case
          for(int j1 = 0; j1 < d; j1++) {
            for(int j2 = 0; j2 < d; j2++) {
              sigma_tilde_ik(j1, j2) = 0.0;
            }
          }
        } else if(!indobs.empty()) {
          int d_obs = indobs.size();
          int d_miss = indmiss.size();
          
          arma::mat sigma_obs_obs(d_obs, d_obs);
          arma::mat sigma_miss_obs(d_miss, d_obs);
          arma::mat sigma_obs_miss(d_obs, d_miss);
          arma::mat sigma_miss_miss(d_miss, d_miss);
          
          for(int j1 = 0; j1 < d_obs; j1++) {
            for(int j2 = 0; j2 < d_obs; j2++) {
              sigma_obs_obs(j1, j2) = sigma_k(indobs[j1], indobs[j2]);
            }
          }
          
          for(int j1 = 0; j1 < d_miss; j1++) {
            for(int j2 = 0; j2 < d_obs; j2++) {
              sigma_miss_obs(j1, j2) = sigma_k(indmiss[j1], indobs[j2]);
              sigma_obs_miss(j2, j1) = sigma_k(indobs[j2], indmiss[j1]);
            }
          }
          
          for(int j1 = 0; j1 < d_miss; j1++) {
            for(int j2 = 0; j2 < d_miss; j2++) {
              sigma_miss_miss(j1, j2) = sigma_k(indmiss[j1], indmiss[j2]);
            }
          }
          
          arma::vec y_obs(d_obs), mu_obs(d_obs), mu_miss(d_miss);
          for(int j = 0; j < d_obs; j++) {
            y_obs[j] = YNA(i, indobs[j]);
            mu_obs[j] = mu_k[indobs[j]];
          }
          for(int j = 0; j < d_miss; j++) mu_miss[j] = mu_k[indmiss[j]];
          
          arma::mat inv_sigma_obs_obs = inv_sympd(sigma_obs_obs);
          arma::vec diff = y_obs - mu_obs;
          arma::vec mu_cond = mu_miss + sigma_miss_obs * inv_sigma_obs_obs * diff;
          arma::mat sigma_cond = sigma_miss_miss - sigma_miss_obs * inv_sigma_obs_obs * sigma_obs_miss;
          
          // Fill mu_tilde_ik
          for(int j = 0; j < d; j++) {
            if(C(i, j) == 0) {
              mu_tilde_ik[j] = YNA(i, j);
            } else {
              int idx = which_index(indmiss, j);
              if(idx >= 0) {
                mu_tilde_ik[j] = mu_cond[idx];
              } else {
                mu_tilde_ik[j] = mu_k[j]; // Fallback
              }
            }
          }
          
          // Initialize sigma_tilde_ik with zeros
          for(int j1 = 0; j1 < d; j1++) {
            for(int j2 = 0; j2 < d; j2++) {
              sigma_tilde_ik(j1, j2) = 0.0;
            }
          }
          
          // Fill sigma_tilde_ik for missing values
          for(int j1 = 0; j1 < d_miss; j1++) {
            for(int j2 = 0; j2 < d_miss; j2++) {
              sigma_tilde_ik(indmiss[j1], indmiss[j2]) = sigma_cond(j1, j2);
            }
          }
        } else {
          // All values are missing
          for(int j = 0; j < d; j++) mu_tilde_ik[j] = mu_k[j];
          for(int j1 = 0; j1 < d; j1++) {
            for(int j2 = 0; j2 < d; j2++) {
              sigma_tilde_ik(j1, j2) = sigma_k(j1, j2);
            }
          }
        }
        
        // Store the results
        // List mu_tilde_i = as<List>(mu_tilde_list[i]);
        // List sigma_tilde_i = as<List>(sigma_tilde_list[i]);
        // mu_tilde_i[k] = mu_tilde_ik;
        // sigma_tilde_i[k] = sigma_tilde_ik;
        // mu_tilde_list[i] = mu_tilde_i;
        // sigma_tilde_list[i] = sigma_tilde_i;
        List mu_tilde_i = clone(as<List>(mu_tilde_list[i]));
        List sigma_tilde_i = clone(as<List>(sigma_tilde_list[i]));
        mu_tilde_i[k] = clone(mu_tilde_ik);
        sigma_tilde_i[k] = clone(sigma_tilde_ik);
        mu_tilde_list[i] = mu_tilde_i;
        sigma_tilde_list[i] = sigma_tilde_i;
      }
    }
    
    // Update tik
    loglik_res = LoglikelihoodObsGaussian(YNA, mu_new, sigma_new, alpha_new, pi_new);
    if(loglik_res.containsElementNamed("error")) return loglik_res;
    tik = as<NumericMatrix>(loglik_res["tik"]);
    
    // M-step
    // Update pi
    for(int k = 0; k < K; k++) pi_new[k] = sum(tik(_, k)) / n;
    
    // Update mu and sigma
    for(int k = 0; k < K; k++) {
      NumericVector mu_new_k(d, 0.0);
      NumericMatrix sigma_new_k(d, d);
      std::fill(sigma_new_k.begin(), sigma_new_k.end(), 0.0);
      double tik_sum = sum(tik(_, k));
      
      for(int i = 0; i < n; i++) {
        List mu_tilde_i = as<List>(mu_tilde_list[i]);
        List sigma_tilde_i = as<List>(sigma_tilde_list[i]);
        NumericVector mu_tilde_ik = as<NumericVector>(mu_tilde_i[k]);
        NumericMatrix sigma_tilde_ik = as<NumericMatrix>(sigma_tilde_i[k]);
        
        for(int j = 0; j < d; j++) {
          mu_new_k[j] += tik(i, k) * mu_tilde_ik[j];
        }
        
        for(int j1 = 0; j1 < d; j1++) {
          for(int j2 = 0; j2 < d; j2++) {
            double diff1 = mu_tilde_ik[j1] - as<NumericVector>(mu_new[k])[j1];
            double diff2 = mu_tilde_ik[j2] - as<NumericVector>(mu_new[k])[j2];
            sigma_new_k(j1, j2) += tik(i, k) * (diff1 * diff2 + sigma_tilde_ik(j1, j2));
          }
        }
      }
      
      for(int j = 0; j < d; j++) mu_new_k[j] /= tik_sum;
      for(int j1 = 0; j1 < d; j1++) {
        for(int j2 = 0; j2 < d; j2++) {
          sigma_new_k(j1, j2) /= tik_sum;
          if(diag && j1 != j2) sigma_new_k(j1, j2) = 0.0;
        }
      }
      
      arma::mat sigma_arma = as<arma::mat>(sigma_new_k);
      if(!sigma_arma.is_sympd()) {
        sigma_arma = symmatu(sigma_arma);
        sigma_arma.diag() += 1e-6;
      }
      mu_new[k] = mu_new_k;
      sigma_new[k] = wrap(sigma_arma);
    }
    
    // Update alpha
    if(mecha != "MCAR") {
      alpha_new = MechanismEMGLM(YNA, tik, mecha);
    }
    
    // Environment base = Environment::namespace_env("base");
    // Function cat = base["cat"];
    // cat("Sucessfully update alpha.\n");
    
    // Update log-likelihood
    prec = loglik_prev;
    loglik_res = LoglikelihoodObsGaussian(YNA, mu_new, sigma_new, alpha_new, pi_new);
    if(loglik_res.containsElementNamed("error")) return loglik_res;
    loglik_prev = as<double>(loglik_res["loglik_obs"]);
    loglik_vec.push_back(loglik_prev);
  }
  
  // Compute imputed data
  NumericMatrix imputedData = clone(YNA);
  for(int i = 0; i < n; i++) {
    std::vector<int> indmiss;
    for(int j = 0; j < d; j++) if(C(i, j) == 1) indmiss.push_back(j);
    if(!indmiss.empty()) {
      for(size_t idx = 0; idx < indmiss.size(); idx++) {
        int j = indmiss[idx];
        double imputed_value = 0.0;
        for(int k = 0; k < K; k++) {
          List mu_tilde_i = as<List>(mu_tilde_list[i]);
          NumericVector mu_tilde_ik = as<NumericVector>(mu_tilde_i[k]);
          imputed_value += tik(i, k) * mu_tilde_ik[j];
        }
        imputedData(i, j) = imputed_value;
      }
    }
  }
  
  // Return results
  return List::create(
    Named("pik") = pi_new,
    Named("mu") = mu_new,
    Named("sigma") = sigma_new,
    Named("alpha") = alpha_new,
    Named("loglik_vec") = loglik_vec,
    Named("tik") = tik,
    Named("imputedData") = imputedData,
    Named("error") = "No error"
  );
}
