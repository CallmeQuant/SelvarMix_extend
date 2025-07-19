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
    Nullable<List>  E_mu_list  = R_NilValue,
    Nullable<List>  E_sig_list = R_NilValue   )
{
  const int n = YNA.nrow();
  const int d = YNA.ncol();
  const int K = prop_pi.size();

  NumericMatrix C(n,d);
  for(int i=0;i<n;++i)
    for(int j=0;j<d;++j) C(i,j) = R_IsNA(YNA(i,j));

  const bool use_cond = (E_mu_list.isNotNull() && E_sig_list.isNotNull());
  List  E_mu,  E_sig;
  if(use_cond){
    E_mu  = List(E_mu_list );
    E_sig = List(E_sig_list);
  }

  // cache distinct missing-data patterns 
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
  
  for(size_t p=0;p<uniqPat.size();++p){
    const NumericVector &pat = uniqPat[p];
    const IntegerVector &idx = patIdx[p];

    std::vector<int> obsIdx;
    for(int j=0;j<d;++j) if(pat[j]==0) obsIdx.push_back(j);
    const int d_obs = obsIdx.size();

    for(int t=0;t<idx.size();++t){
      const int i = idx[t];

      NumericVector y_obs_rcpp(d_obs);
      for(int jj=0;jj<d_obs;++jj) y_obs_rcpp[jj] = YNA(i,obsIdx[jj]);

      for(int k=0;k<K;++k){
        // Block-diagonal and determinant calculations
        double log_y = 0.0;
        if(d_obs>0){
          arma::vec y_obs   = as<arma::vec>(y_obs_rcpp);
          arma::vec mu_k    = as<arma::vec>(mu[k]);
          arma::mat Sig_k   = as<arma::mat>(sigma[k]);

          arma::uvec obs_u(d_obs);
          for(int jj=0;jj<d_obs;++jj) obs_u[jj]=obsIdx[jj];
          arma::vec mu_o   = mu_k.elem(obs_u);
          arma::mat Sig_oo = Sig_k.submat(obs_u,obs_u);

          // Ensure Sig_oo is symmetric and positive definite
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

        // Missing mechanism probabilities
        // Updated the logic to handle both MAR and MNAR cases
        double log_mask = 0.0;

        // Inside the k-cluster loop:
        double log_mask = 0.0;
        for(int j=0; j<d; ++j) {
          const bool miss = (pat[j]==1);
          const bool MNAR = is_mnar[j];
          
          if(MNAR) {
            // Pure MNARz: depends only on cluster
            double p_mnar = R::pnorm(alpha(k,j), 0.0, 1.0, 1, 0);
            p_mnar = std::max(std::min(p_mnar, 1.0-1e-12), 1e-12);
            log_mask += miss ? std::log(p_mnar) : std::log1p(-p_mnar);
          }
          else {
            // MAR: depends on expected value
            double predictor;
            if(miss) {
              // Missing: use conditional expectation
              if(use_cond) {
                List mu_i = as<List>(E_mu[i]);
                NumericVector mu_ik = mu_i[k];
                predictor = mu_ik[j];
              } else {
                predictor = 0.0; // Fallback
              }
            } else {
              // Observed: use actual value
              predictor = YNA(i,j);
            }
            double linear_pred = alpha(k,j) + beta(k,j) * predictor;
            double p_mar = R::pnorm(linear_pred, 0.0, 1.0, 1, 0);
            p_mar = std::max(std::min(p_mar, 1.0-1e-12), 1e-12);
            log_mask += miss ? std::log(p_mar) : std::log1p(-p_mar);
          }
        }
        logprob(i,k) = std::log(prop_pi[k]) + log_y + log_mask;
      } // k
    }   // i
  }     // pattern loop

  // Normalize log probabilities
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
      loglik += m[i] - 700;          
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
    
    // Initialize alpha and beta based on mechanism
    NumericMatrix alpha_init(K, d);
    NumericMatrix beta_init(K, d);
    beta_init.fill(0.0); 

    // Calculate overall missing rates per variable
    NumericVector var_miss_rate(d);
    for(int j=0; j<d; ++j) {
      int miss_count = 0;
      for(int i=0; i<n; ++i) if(R_IsNA(YNA(i,j))) miss_count++;
      var_miss_rate[j] = std::max(1e-6, std::min(1.0-1e-6, static_cast<double>(miss_count)/n));
    }

    for(int k=0; k<K; ++k) {
      IntegerVector class_obs;
      for(int i=0; i<n; ++i) if(Z_init(i,k) == 1.0) class_obs.push_back(i);
      int nk = class_obs.size();

      // Handle empty clusters
      if (nk == 0) {
        for (int j = 0; j < d; ++j) alpha_init(k, j) = 0.0;
        continue;
      }

      for(int j=0; j<d; ++j) {
        if(is_mnar[j]) {
          // MNARz: Cluster-specific initialization
          int miss_count = 0;
          for(int idx=0; idx<nk; ++idx) {
            int i = class_obs[idx];
            if(R_IsNA(YNA(i,j))) miss_count++;
          }
          double miss_rate = std::max(1e-6, std::min(1.0-1e-6, static_cast<double>(miss_count)/nk));
          alpha_init(k,j) = R::qnorm(miss_rate, 0.0, 1.0, 1, 0);
        } else {
          // MAR: Variable-specific initialization (same for all clusters)
          alpha_init(k,j) = R::qnorm(var_miss_rate[j], 0.0, 1.0, 1, 0);
        }
      }
    }
    
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
      Named("alpha_init") = alpha_user,
      Named("beta_init") = beta_user 
    );
  }
}

// [[Rcpp::export]]
List MechanismEMGLMMixed(
    NumericMatrix YNA, 
    NumericMatrix tik, 
    std::string mecha,
    LogicalVector is_mnar,
    List E_y_list,
    NumericMatrix current_alpha,
    NumericMatrix current_beta
) {
  int n = YNA.nrow();
  int d = YNA.ncol();
  int K = tik.ncol();
  
  NumericMatrix C(n, d);
  for(int i = 0; i < n; i++) 
    for(int j = 0; j < d; j++) 
      C(i, j) = R_IsNA(YNA(i, j)) ? 1.0 : 0.0;
  
  NumericMatrix alpha_new(K, d);
  NumericMatrix beta_new(K, d);
  beta_new.fill(0.0); 
  
  Environment stats = Environment::namespace_env("stats");
  Function glm = stats["glm"];
  Function binomial = stats["binomial"];
  RObject family = binomial(Rcpp::Named("link", "probit"));
  
  // Separate handling for MAR and MNARz
  for (int j = 0; j < d; ++j) {
    for (int k = 0; k < K; ++k) {
      NumericVector y_response = C(_, j);
      NumericVector weights = tik(_, k);
      
      std::vector<double> y_filt;
      std::vector<double> w_filt;
      std::vector<double> predictors; // Only for MAR

      // Filter observations with meaningful weights
      for(int i = 0; i < n; ++i) {
        if(weights[i] > 1e-8) {
          y_filt.push_back(y_response[i]);
          w_filt.push_back(weights[i]);
          
          // MAR: Use expected value from E-step
          if (!is_mnar[j]) {
            List E_y_i = as<List>(E_y_list[i]);
            NumericVector E_y_ik = as<NumericVector>(E_y_i[k]);
            predictors.push_back(E_y_ik[j]);
          }
        }
      }
      
      if(y_filt.empty()) {
        alpha_new(k, j) = current_alpha(k, j);
        beta_new(k, j) = current_beta(k, j);
        continue;
      }

      try {
        if (is_mnar[j]) {
          // MNARz: Intercept-only model (cluster-specific)
          DataFrame df = DataFrame::create(_["y"] = NumericVector(y_filt.begin(), y_filt.end()));
          List glm_fit = glm(
            _["formula"] = Formula("y ~ 1"),
            _["family"] = family,
            _["data"] = df,
            _["weights"] = NumericVector(w_filt.begin(), w_filt.end())
          );
          
          NumericVector coef = glm_fit["coefficients"];
          alpha_new(k, j) = coef[0];
          // beta_new remains 0 (no data dependence)
        } 
        else {
          // MAR: Self-masking model (using expected value)
          DataFrame df = DataFrame::create(
            _["y"] = NumericVector(y_filt.begin(), y_filt.end()),
            _["x"] = NumericVector(predictors.begin(), predictors.end())
          );
          
          List glm_fit = glm(
            _["formula"] = Formula("y ~ x"),
            _["family"] = family,
            _["data"] = df,
            _["weights"] = NumericVector(w_filt.begin(), w_filt.end())
          );
          
          NumericVector coef = glm_fit["coefficients"];
          alpha_new(k, j) = coef[0];
          beta_new(k, j) = coef[1];
        }
      } 
      catch (...) {
        alpha_new(k, j) = current_alpha(k, j);
        beta_new(k, j) = current_beta(k, j);
      }
    }
  }
  
  return List::create(
    _["alpha_new"] = alpha_new,
    _["beta_new"] = beta_new
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