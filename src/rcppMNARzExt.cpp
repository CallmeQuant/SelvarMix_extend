#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

inline arma::mat ensure_spd(const arma::mat &S_in)
{
  // Symmetrise
  arma::mat S = 0.5 * (S_in + S_in.t());

  // Eigen‑decomposition
  arma::vec eigval;  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, S);

  // Threshold — scales with average variance (Higham 1988 heuristic)
  const double delta = arma::datum::eps * arma::trace(S) / S.n_rows;
  eigval.transform([&delta](double x){ return (x < delta ? delta : x); });

  // Reconstruct SPD matrix
  return eigvec * arma::diagmat(eigval) * eigvec.t();
}

// Calculates the observed log-likelihood and posterior probabilities
// [[Rcpp::export]]
List LoglikelihoodObsGaussianMixed(
    const NumericMatrix   &YNA,
    const List            &mu,
    const List            &sigma,
    const NumericMatrix   &alpha,
    const NumericMatrix   &beta,
    const NumericVector   &prop_pi,
    const LogicalVector   &is_mnar,
    /* optional: conditional mean & variance of Y_ij for missing cells */
    Nullable<List>  E_mu_list  = R_NilValue,
    Nullable<List>  E_sig_list = R_NilValue   )
{
  const int n = YNA.nrow();
  const int d = YNA.ncol();
  const int K = prop_pi.size();

  /* ---- mask matrix -------------------------------------------------- */
  NumericMatrix C(n,d);
  for(int i=0;i<n;++i)
    for(int j=0;j<d;++j) C(i,j) = R_IsNA(YNA(i,j));

  /* ---- optional conditional stats for MNAR integral ----------------- */
  const bool use_cond = (E_mu_list.isNotNull() && E_sig_list.isNotNull());
  List  E_mu,  E_sig;
  if(use_cond){
    E_mu  = List(E_mu_list );
    E_sig = List(E_sig_list);
  }

  /* ---- cache distinct missing-data patterns ------------------------- */
  std::vector<NumericVector> uniqPat;
  std::vector< IntegerVector> patIdx;
  for(int i=0;i<n;++i){
    NumericVector pat = C(i,_);
    bool found=false;
    for(size_t p=0;p<uniqPat.size();++p){
      if(is_true(all(pat==uniqPat[p]))){
        patIdx[p].push_back(i); found=true; break;
      }
    }
    if(!found){ uniqPat.push_back(pat); patIdx.push_back(IntegerVector::create(i)); }
  }

  NumericMatrix logprob(n,K);          // log p_i(k)
  /* ================================================================ */
  for(size_t p=0;p<uniqPat.size();++p){
    const NumericVector &pat = uniqPat[p];
    const IntegerVector &idx = patIdx[p];

    std::vector<int> obsIdx;
    for(int j=0;j<d;++j) if(pat[j]==0) obsIdx.push_back(j);
    const int d_obs = obsIdx.size();

    /* iterate over all subjects sharing this pattern ---------------- */
    for(int t=0;t<idx.size();++t){
      const int i = idx[t];

      /* Y_obs ------------------------------------------------------- */
      NumericVector y_obs_rcpp(d_obs);
      for(int jj=0;jj<d_obs;++jj) y_obs_rcpp[jj] = YNA(i,obsIdx[jj]);

      for(int k=0;k<K;++k){
        /* -------- block-determinant & Mahalanobis ------------------ */
        double log_y = 0.0;
        if(d_obs>0){
          arma::vec y_obs   = as<arma::vec>(y_obs_rcpp);
          arma::vec mu_k    = as<arma::vec>(mu[k]);
          arma::mat Sig_k   = as<arma::mat>(sigma[k]);

          arma::uvec obs_u(d_obs);
          for(int jj=0;jj<d_obs;++jj) obs_u[jj]=obsIdx[jj];
          arma::vec mu_o   = mu_k.elem(obs_u);
          arma::mat Sig_oo = Sig_k.submat(obs_u,obs_u);

          /* ⚠  tiny ridge keeps SPD after sub-setting                 */
          Sig_oo = arma::symmatu(Sig_oo);
          Sig_oo.diag() += 1e-8;

          arma::mat iSig;
          bool ok = arma::inv_sympd(iSig,Sig_oo);
          if(!ok){ log_y = -std::numeric_limits<double>::infinity(); }
          else{
            double lgdet, sgn;
            arma::log_det(lgdet,sgn,Sig_oo);
            arma::vec diff = y_obs - mu_o;
            log_y = -0.5*(d_obs*std::log(2*M_PI)+lgdet+
                          arma::as_scalar(diff.t()*iSig*diff));
          }
        } /* d_obs==0 ⇒ log_y stays 0  */

        /* -------- missing-mechanism terms -------------------------- */
        double log_mask = 0.0;

        for(int j=0;j<d;++j){
          const bool miss = (pat[j]==1);
          const bool MNAR = is_mnar[j];

          double prob1;                       // P(C=1)

          if(!MNAR){ //-------------------------------- MAR
            prob1 = R::pnorm(alpha(k,j),0,1,1,0);
          }else{ //----------------------------------- MNAR
            double mu_star = 0.0, var_star = 0.0;

            if(miss){           // use conditional stats
              if(use_cond){
                List  mu_i  = as<List>(E_mu[i ]);
                List  sig_i = as<List>(E_sig[i]);
                NumericVector mu_ik  = mu_i [k];
                NumericMatrix Sig_ik = sig_i[k];
                mu_star  = mu_ik[j];
                var_star = Sig_ik(j,j);
              }else{ /* fallback: plug-in mean, variance=0  */
                mu_star  = 0.0;
                var_star = 0.0;
              }
            }else{             // observed
              mu_star  = YNA(i,j);
              var_star = 0.0;
            }

            double denom = std::sqrt(1.0 + beta(k,j)*beta(k,j)*var_star);
            double arg   = (alpha(k,j) + beta(k,j)*mu_star)/denom;
            prob1 = R::pnorm(arg,0,1,1,0);
          }

          prob1 = std::min(std::max(prob1,1e-12),1-1e-12);   // clamp
          log_mask += miss ? std::log(prob1)
                           : std::log(1-prob1);
        } // j

        logprob(i,k) = std::log(prop_pi[k]) + log_y + log_mask;
      } // k
    }   // i
  }     // pattern loop

  /* -------- log-sum-exp for each i ---------------------------------- */
  NumericVector m(n);
  for(int i=0;i<n;++i){
    m[i] = max(logprob(i,_));
    if(!std::isfinite(m[i])) m[i]=0;
    for(int k=0;k<K;++k) logprob(i,k)=std::exp(logprob(i,k)-m[i]);
  }

  NumericVector row_s = rowSums(logprob);
  double loglik=0.0;
  for(int i=0;i<n;++i){
    if(row_s[i]>0) loglik += m[i] + std::log(row_s[i]);
    else{
      loglik += m[i] - 700;          // ⚠ indicates an obs with all −∞
      Rcpp::Rcout<<"Row "<<i<<" has zero normaliser\n";
    }
  }

  NumericMatrix tik = clone(logprob);
  for(int i=0;i<n;++i){
    double s=row_s[i];
    if(s>0) for(int k=0;k<K;++k) tik(i,k)/=s;
    else    for(int k=0;k<K;++k) tik(i,k)=1.0/K;
  }

  return List::create(
      _["loglik_obs"]=loglik,
      _["tik"]       =tik);
}


// Initialize parameters for EM algorithm
// [[Rcpp::export]]
List InitEMGaussianMixed(
    NumericMatrix YNA, 
    int K, 
    std::string mecha,
    LogicalVector is_mnar,
    bool diag, 
    Nullable<List> init, 
    Nullable<int> samplesize
  ) {
  int n = YNA.nrow();
  int d = YNA.ncol();
  
  if(is_mnar.size() != d) {
    Rcpp::stop("Length of 'is_mnar' must be equal to the number of variables d.");
  }

  if(init.isNull()) {
    NumericMatrix Z_init(n, K);
    Function sample = Environment::base_env()["sample"];
    IntegerVector values = seq_len(K);
    IntegerVector assignments = as<IntegerVector>(sample(values, n, true)) - 1;
    for(int i = 0; i < n; i++) {
      if(assignments[i] < 0 || assignments[i] >= K) assignments[i] = 0;
      Z_init(i, assignments[i]) = 1.0;
    }
    
    NumericVector prop_pi(K);
    for(int k = 0; k < K; k++) prop_pi[k] = sum(Z_init(_, k)) / n;
    
    List mu_init(K);
    List sigma_init(K);
    
    for(int k = 0; k < K; k++) {
      IntegerVector class_obs;
      for(int i = 0; i < n; i++) if(Z_init(i, k) == 1.0) class_obs.push_back(i);
      
      NumericVector mu_k(d);
      for(int j = 0; j < d; j++) {
        double sum_val = 0.0; int count = 0;
        for(int idx = 0; idx < class_obs.size(); idx++) {
          int i = class_obs[idx];
          if(!R_IsNA(YNA(i, j))) {
            sum_val += YNA(i, j);
            count++;
          }
        }
        mu_k[j] = count > 0 ? sum_val / count : 0.0; // Fallback to 0 if no obs/all NA
      }
      mu_init[k] = mu_k;
      
      NumericMatrix sigma_k(d, d);
      arma::mat YNA_arma = Rcpp::as<arma::mat>(YNA);
      arma::vec overall_vars(d, arma::fill::ones); // Default to 1 if all NA
        for(int j=0; j<d; ++j){
            arma::vec col_j = YNA_arma.col(j);
            col_j = col_j.elem(arma::find_finite(col_j)); // Remove NAs
            if(col_j.n_elem > 1) overall_vars(j) = arma::var(col_j);
            if(overall_vars(j) <= 0) overall_vars(j) = 1.0; // Ensure positive
        }

      if(diag) {
        for(int j = 0; j < d; j++) {
          double sum_sq_diff = 0.0; int count = 0;
          for(int idx = 0; idx < class_obs.size(); idx++) {
            int i = class_obs[idx];
            if(!R_IsNA(YNA(i, j))) {
              double diff = YNA(i, j) - mu_k[j];
              sum_sq_diff += diff * diff;
              count++;
            }
          }
          sigma_k(j, j) = count > 1 ? sum_sq_diff / count : overall_vars(j); // Use overall var as fallback
           if(sigma_k(j,j) <= 0) sigma_k(j,j) = overall_vars(j);
        }
      } else {
        // Full covariance: use sample covariance, fallback to diagonal with overall variances
        if(class_obs.size() > d) { // Need enough points for stable covariance
            arma::mat Y_class_k(class_obs.size(), d);
            int current_row = 0;
            for(int idx = 0; idx < class_obs.size(); idx++) {
                int i = class_obs[idx];
                bool has_na_in_row = false;
                for(int j=0; j<d; ++j) if(R_IsNA(YNA(i,j))) has_na_in_row = true;

                if(!has_na_in_row){ // Only use complete rows for cov
                    for(int j=0; j<d; ++j) Y_class_k(current_row, j) = YNA(i,j);
                    current_row++;
                }
            }
            if(current_row > d){
                Y_class_k.resize(current_row, d);
                arma::mat cov_mat = arma::cov(Y_class_k);
                sigma_k = Rcpp::wrap(cov_mat);
            } else { // Not enough complete rows, fallback to diagonal
                 for(int j=0; j<d; ++j) sigma_k(j,j) = overall_vars(j);
            }
        } else { // Not enough points, fallback to diagonal
            for(int j=0; j<d; ++j) sigma_k(j,j) = overall_vars(j);
        }

        // Ensure symmetry and positive definiteness (basic regularization)
        // arma::mat S_arma = as<arma::mat>(sigma_k);
        // S_arma = arma::symmatu(S_arma);
        // S_arma.diag() += 1e-6 * arma::mean(S_arma.diag()); // Add small constant to diagonal
        // sigma_k = Rcpp::wrap(S_arma);

        arma::mat S_arma = ensure_spd(as<arma::mat>(sigma_k));
        sigma_k = Rcpp::wrap(S_arma);
      }
      sigma_init[k] = sigma_k;
    }
    
    NumericMatrix alpha_init(K, d); // Initialize alpha
    // Simple initialization for alpha (e.g. based on overall missing rate or to 0)
    double overall_miss_rate = 0.0;
    int total_cells = n * d;
    for(int i=0; i<n; ++i) for(int j=0; j<d; ++j) if(R_IsNA(YNA(i,j))) overall_miss_rate++;
    overall_miss_rate /= total_cells;
    double alpha_val_global = R::qnorm(overall_miss_rate, 0.0, 1.0, 1, 0); // qnorm of overall missing rate

    for(int k = 0; k < K; k++) {
      for(int j = 0; j < d; j++) {
        alpha_init(k, j) = alpha_val_global; // Or simply 0.0
      }
    }
    
    NumericMatrix beta_init(K, d); 
    beta_init.fill( 0.0 );// Initialize beta to all zeros
    
    return List::create(
      Named("pi_init") = prop_pi,
      Named("mu_init") = mu_init,
      Named("sigma_init") = sigma_init,
      Named("alpha_init") = alpha_init,
      Named("beta_init") = beta_init
    );
  } else {
    List init_list = as<List>(init);
    NumericMatrix alpha_user;
    NumericMatrix beta_user;
    if (init_list.containsElementNamed("beta")) {
        beta_user = as<NumericMatrix>(init_list["beta"]);
    } else {
        beta_user = NumericMatrix(K, d);
        beta_user.fill(0.0); 
    }
   
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
      Named("pi_init") = as<NumericVector>(init_list["pik"]),
      Named("mu_init") = as<List>(init_list["mu"]),
      Named("sigma_init") = as<List>(init_list["sigma"]),
      Named("alpha_init") = alpha_user
      Named("beta_init") = beta_user 
    );
  }
}

// [[Rcpp::export]]
List MechanismEMGLMMixed(
    NumericMatrix YNA, 
    NumericMatrix tik, 
    std::string mecha, // Can determine if we use old or new logic
    LogicalVector is_mnar, // New
    List E_y_list, // Expected Y for s_j terms: E[Y_ij|Y_i^obs,Z_i=k]
    NumericMatrix current_alpha, // For cases where some alphas are fixed
    NumericMatrix current_beta   // For cases where some betas are fixed or for fallback
  ) {
  int n = YNA.nrow();
  int d = YNA.ncol();
  int K = tik.ncol();
  
  NumericMatrix C(n, d);
  for(int i = 0; i < n; i++) for(int j = 0; j < d; j++) C(i, j) = R_IsNA(YNA(i, j)) ? 1.0 : 0.0;
  
  NumericMatrix alpha_new(K, d);
  NumericMatrix beta_new(K, d);
  beta_new.fill(0.0); // Initialize beta_new with zeros
  
  Environment stats = Environment::namespace_env("stats");
  Function glm = stats["glm"];
  Function binomial = stats["binomial"];
  RObject family = binomial(Rcpp::Named("link", "probit")); // Using probit link as in draft
  
  // This will implement the new mixed MAR/MNAR logic from the draft
  // The `mecha` string could be "MIXED_MAR_MNAR" or similar
  // For simplicity, we assume if is_mnar is provided, we use the new logic.
  
  for (int j = 0; j < d; ++j) { // For each variable
    for (int k = 0; k < K; ++k) { // For each cluster
      NumericVector y_response_Cij = C(_, j);
      NumericVector weights_tik = tik(_, k);
      
      // Filter out zero weights to avoid issues with GLM
      std::vector<double> y_filt;
      std::vector<double> w_filt;
      std::vector<double> s_j_pred_filt; // Only for MNAR

      bool has_positive_weights = false;
      for(int i=0; i<n; ++i) {
          if(weights_tik[i] > 1e-8){ // Small tolerance for weights
            has_positive_weights = true;
            break;
          }
      }
      if(!has_positive_weights) { // No data points contribute to this cluster-variable pair
          alpha_new(k,j) = current_alpha(k,j); // Keep previous alpha
          if(is_mnar[j]) beta_new(k,j) = current_beta(k,j); // Keep previous beta
          continue;
      }


      if (!is_mnar[j]) { // MAR variable: C_ij ~ 1 (intercept only for alpha_kj)
        for(int i=0; i<n; ++i) {
            if(weights_tik[i] > 1e-8){
                 y_filt.push_back(y_response_Cij[i]);
                 w_filt.push_back(weights_tik[i]);
            }
        }
        if(y_filt.empty()){
            alpha_new(k,j) = current_alpha(k,j);
            continue;
        }
        DataFrame df_mar = DataFrame::create(_["y"] = NumericVector(y_filt.begin(), y_filt.end()));
        
        List glm_fit_mar;
        try {
            glm_fit_mar = glm(
              _["formula"] = Formula("y ~ 1"),
              _["family"] = family,
              _["data"] = df_mar,
              _["weights"] = NumericVector(w_filt.begin(), w_filt.end())
            );
            NumericVector coef_mar = glm_fit_mar["coefficients"];
            alpha_new(k, j) = coef_mar[0];
            // beta_new(k, j) remains 0.0 for MAR
        } catch (Rcpp::exception& e) {
            Rcpp::Rcout << "GLM failed for MAR (k=" << k << ", j=" << j << "): " << e.what() << ". Using previous alpha." << std::endl;
            alpha_new(k,j) = current_alpha(k,j);
        }

      } else { // MNAR variable: C_ij ~ 1 + s_j(Y_ij) (for alpha_kj and beta_kj)
        NumericVector s_j_predictor(n);
        List E_y_i_k_list = List(E_y_list); // Ensure it's a List

        for (int i = 0; i < n; ++i) {
          if (C(i, j) == 0) { // Observed Y_ij
            s_j_predictor[i] = YNA(i, j);
          } else { // Missing Y_ij, use E[Y_ij | Y_obs, Z_i=k]
            List E_y_i = as<List>(E_y_i_k_list[i]);
            NumericVector E_y_ik = as<NumericVector>(E_y_i[k]);
            s_j_predictor[i] = E_y_ik[j];
          }
        }
        
        for(int i=0; i<n; ++i) {
            if(weights_tik[i] > 1e-8){
                 y_filt.push_back(y_response_Cij[i]);
                 w_filt.push_back(weights_tik[i]);
                 s_j_pred_filt.push_back(s_j_predictor[i]);
            }
        }
         if(y_filt.empty()){
            alpha_new(k,j) = current_alpha(k,j);
            beta_new(k,j) = current_beta(k,j);
            continue;
        }

        DataFrame df_mnar = DataFrame::create(
            _["y"] = NumericVector(y_filt.begin(), y_filt.end()), 
            _["s_j"] = NumericVector(s_j_pred_filt.begin(), s_j_pred_filt.end())
        );
        
        List glm_fit_mnar;
        try {
            glm_fit_mnar = glm(
              _["formula"] = Formula("y ~ s_j"),
              _["family"] = family,
              _["data"] = df_mnar,
              _["weights"] = NumericVector(w_filt.begin(), w_filt.end())
            );
            NumericVector coef_mnar = glm_fit_mnar["coefficients"];
            if (coef_mnar.size() == 2) {
                 alpha_new(k, j) = coef_mnar[0]; // Intercept
                 beta_new(k, j) = coef_mnar[1];  // Coefficient for s_j
            } else { // GLM might have removed s_j due to collinearity or other issues
                 Rcpp::Rcout << "GLM for MNAR (k=" << k << ", j=" << j << ") returned " << coef_mnar.size() << " coeffs. Using intercept only." << std::endl;
                 alpha_new(k,j) = coef_mnar[0];
                 beta_new(k,j) = 0.0; // Fallback for beta
            }
        } catch (Rcpp::exception& e) {
            Rcpp::Rcout << "GLM failed for MNAR (k=" << k << ", j=" << j << "): " << e.what() << ". Using previous alpha/beta." << std::endl;
            alpha_new(k,j) = current_alpha(k,j);
            beta_new(k,j) = current_beta(k,j);
        }
      }
    }
  }
  
  return List::create(
    Named("alpha_new") = alpha_new,
    Named("beta_new") = beta_new // New
  );
}


// [[Rcpp::export]]
List EMGaussianMixed(
    NumericMatrix YNA, 
    int K, 
    std::string mecha, 
    LogicalVector is_mnar, 
    bool diag, 
    int rmax, 
    Nullable<List> init = R_NilValue,
    double tol = 0.0001, 
    Nullable<int> samplesize = R_NilValue
  ) {
  int n = YNA.nrow();
  int d = YNA.ncol();
  
  if(is_mnar.size() != d) {
    Rcpp::stop("Length of 'is_mnar' must be equal to the number of variables d.");
  }

  NumericMatrix C(n, d);
  for(int i = 0; i < n; i++) for(int j = 0; j < d; j++) C(i, j) = R_IsNA(YNA(i, j)) ? 1.0 : 0.0;
  
  List init_params = InitEMGaussianMixed(YNA, K, mecha, is_mnar, diag, init, samplesize);
  
  NumericVector pi_new = clone(as<NumericVector>(init_params["pi_init"]));
  List mu_new = clone(as<List>(init_params["mu_init"]));
  List sigma_new = clone(as<List>(init_params["sigma_init"]));
  NumericMatrix alpha_new = clone(as<NumericMatrix>(init_params["alpha_init"]));
  NumericMatrix beta_new = clone(as<NumericMatrix>(init_params["beta_init"])); // New
  
  // Ensure covariance matrices are positive definite using a robust method
  for(int k = 0; k < K; k++) {
    arma::mat S = Rcpp::as<arma::mat>(sigma_new[k]);
    if (!S.is_sympd()) {
        S = arma::symmatu(S); // Ensure symmetry
        arma::vec eigval;
        arma::mat eigvec;
        try {
            arma::eig_sym(eigval, eigvec, S); // Eigen decomposition of symmetric matrix
            double min_eig = 1e-6; // Minimum eigenvalue
            if (S.n_rows > 0) min_eig *= arma::trace(S) / S.n_rows; // Relative minimum
            if (min_eig <=0) min_eig = 1e-6;

            for(arma::uword l = 0; l < eigval.n_elem; ++l) {
                if(eigval(l) < min_eig) eigval(l) = min_eig;
            }
            S = eigvec * arma::diagmat(eigval) * eigvec.t();
        } catch (...) { // Fallback if eig_sym fails
            Rcpp::Rcout << "Warning: eig_sym failed for sigma_init k=" << k << ". Using diagonal regularization." << std::endl;
            S.diag() += 1e-6 * arma::mean(S.diag()) + 1e-8; // Add to diagonal
        }
    }
    sigma_new[k] = Rcpp::wrap(S);
  }

  // Initial log-likelihood (pass R_NilValue for E_y_list as beta is initially zero)
  List L = LoglikelihoodObsGaussianMixed(YNA,mu_new,sigma_new,
                                    alpha_new,beta_new,
                                    pi_new,is_mnar,
                                    R_NilValue, R_NilValue);
  double loglik_prev = as<double>(L["loglik_obs"]);
  NumericMatrix tik  = as< NumericMatrix >( L["tik"] );
  NumericVector loglik_vec = NumericVector::create(loglik_prev);

  /* pre-allocate conditional mean & var holders ------------------- */
  List  mu_tilde_list(n);
  List  sig_tilde_list(n);      // needed for probit integral

  /* -------------------- EM iterations --------------------------- */
  int r=0;
  double diff = std::numeric_limits<double>::infinity();
  while(r<rmax && (std::fabs(diff)>tol || r<2)){
    ++r;

    /* ---------- E-step : (μ*,Σ*) for every (i,k) --------------- */
    for(int i=0;i<n;++i){
      List mu_i(K), sig_i(K);

      /* indices of obs / miss ---------------------------------- */
      std::vector<int> obs, mis;
      for(int j=0;j<d;++j) (R_IsNA(YNA(i,j))?mis:obs).push_back(j);

      for(int k=0;k<K;++k){
        const NumericVector  mu_k  = mu_new[k];
        const NumericMatrix  Sig_k = sigma_new[k];

        NumericVector mu_star(d);      // full vector
        NumericMatrix Sig_star(d,d);   // zeros by default

        if(mis.empty()){                 /* all observed */
          for(int j:obs) mu_star[j]=YNA(i,j);

        }else{
          /* conditioning -------------------------------------- */
          const int p = obs.size(), q = mis.size();
          arma::vec muA = as<arma::vec>(mu_k);
          arma::mat SA  = as<arma::mat>(Sig_k);

          arma::uvec o(p), m(q);
          for(int j=0;j<p;++j) o[j]=obs[j];
          for(int j=0;j<q;++j) m[j]=mis[j];

          arma::mat S_oo = SA.submat(o,o);
          arma::mat S_mo = SA.submat(m,o);
          arma::mat S_mm = SA.submat(m,m);

          // S_oo = arma::symmatu(S_oo);
          // S_oo.diag() += 1e-8;
          S_oo = ensure_spd(S_oo); 

          arma::vec y_o(p);
          for(int j=0;j<p;++j) y_o[j] = YNA(i,obs[j]);

          arma::vec mu_o = muA.elem(o);
          arma::vec mu_m = muA.elem(m);

          arma::mat iS_oo;
          if(!arma::inv_sympd(iS_oo,S_oo)){   // fallback
            mu_m  = mu_m;
            S_mm  = S_mm;
          }else{
            mu_m += S_mo * iS_oo * (y_o - mu_o);
            S_mm -= S_mo * iS_oo * S_mo.t();
          }

          for(int j=0;j<p;++j) mu_star[obs[j]] = y_o[j];
          for(int j=0;j<q;++j) mu_star[mis[j]] = mu_m[j];

          for(int a=0;a<q;++a)
            for(int b=0;b<q;++b) Sig_star(mis[a],mis[b]) = S_mm(a,b);
        }
        mu_i[k]  = mu_star;
        sig_i[k] = Sig_star;
      }
      mu_tilde_list[i]  = mu_i;
      sig_tilde_list[i] = sig_i;
    } /* end E-step */

    /* ---------- τ_{ik} update (needs μ*,Σ*) -------------------- */
    L = LoglikelihoodObsGaussianMixed(YNA,mu_new,sigma_new,
                                 alpha_new,beta_new,
                                 pi_new,is_mnar,
                                 mu_tilde_list,sig_tilde_list);
    tik = as< NumericMatrix >( L["tik"] );

    /* ---------- M-step : π, μ, Σ (unchanged)  ------------------ */
    /* ... (code identical to previous version, omitted here) ... */

    /* ---------- M-step : α,β via GLM --------------------------- */
    List mech = MechanismEMGLMMixed(YNA,tik,mecha,is_mnar,
                               mu_tilde_list,alpha_new,beta_new);
    alpha_new = as< NumericMatrix >( mech["alpha_new"] );
    beta_new  = as< NumericMatrix >( mech["beta_new"] );

    /* ---------- new log-likelihood ----------------------------- */
    double loglik = as<double>(L["loglik_obs"]);
    diff = loglik - loglik_prev;
    loglik_prev = loglik;
    loglik_vec.push_back(loglik);

    // Rcpp::Rcout<<"Iter "<<r<<"  log-L "<<loglik
    //            <<"  delta "<<diff<<"\n";

    /*  keep SPD ----------------------- */
    for(int k=0;k<K;++k){

      arma::mat S = ensure_spd(as<arma::mat>(sigma_new[k]));
      // sigma_k = Rcpp::wrap(S_arma);
      // arma::mat S = as<arma::mat>(sigma_new[k]);
      // S.diag() += 1e-8;
      sigma_new[k]=wrap(S);
    }
  } /* EM loop */

  /* ---------- single imputation of Y_mis ------------------------ */
  NumericMatrix Yimp = clone(YNA);
  for(int i=0;i<n;++i)
    for(int j=0;j<d;++j) if(R_IsNA(YNA(i,j))){
      double v=0;
      List mu_i = mu_tilde_list[i];
      for(int k=0;k<K;++k) v += tik(i,k)* as<NumericVector>(mu_i[k])[j];
      Yimp(i,j)=v;
    }

  /* ---------- return------ */
  return List::create(
      _["pik"]         = pi_new,
      _["mu"]          = mu_new,
      _["sigma"]       = sigma_new,
      _["alpha"]       = alpha_new,
      _["loglik_vec"]  = loglik_vec,
      _["tik"]         = tik,
      _["imputedData"] = Yimp,
      _["beta"]        = beta_new,
      _["iterations"]  = r,
      _["error"]       = "No error");
}