#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <Rmath.h>
#include <Rdefines.h>


using namespace std;
using namespace Rcpp;
using namespace arma;

#include "Mixture.h"
#include "Function.h"

Mixture::Mixture(Rcpp::List     InputList,
    double         lambda_mu,
    double         lambda_omega,
    arma::cube     Pk_in)
{
Xd   = Rcpp::as<arma::mat>(InputList[0]);
prop = Rcpp::as<arma::rowvec>(InputList[1]);
n        = Xd.n_rows;
p        = Xd.n_cols;
nbClust  = prop.n_elem;

Mu = Rcpp::as<arma::mat>(InputList[2]);

Rcpp::NumericVector vecS (InputList[3]);
CovarianceMatrix = arma::cube(vecS.begin(), p, p, nbClust, false);
EmpiricalCovariance = CovarianceMatrix;             

Rcpp::NumericVector vecW (InputList[4]);
PrecisionMatrix = arma::cube(vecW.begin(), p, p, nbClust, false);

ProbCond = Rcpp::as<arma::mat>(InputList[5]);

lambda = lambda_mu;
rho    = lambda_omega;

Pk_cube = Pk_in;   // deep copy

if (   Pk_cube.n_rows   != static_cast<uword>(p)
|| Pk_cube.n_cols   != static_cast<uword>(p)
|| Pk_cube.n_slices != static_cast<uword>(nbClust) )
Rcpp::stop("Pk_in has wrong dimension (must be p×p×K)");

for (uword k = 0; k < Pk_cube.n_slices; ++k)
Pk_cube.slice(k).diag().zeros();
}

Mixture::Mixture(Rcpp::List InputList)
{
    Xd   = Rcpp::as<arma::mat>(InputList[0]);
    prop = Rcpp::as<arma::rowvec>(InputList[1]);
    n        = Xd.n_rows;
    p        = Xd.n_cols;
    nbClust  = prop.n_elem;

    Mu = Rcpp::as<arma::mat>(InputList[2]);
    Rcpp::NumericVector vecS (InputList[3]);
    CovarianceMatrix = arma::cube(vecS.begin(), p, p, nbClust, false);
    EmpiricalCovariance = CovarianceMatrix;      

    Rcpp::NumericVector vecW (InputList[4]);
    PrecisionMatrix = arma::cube(vecW.begin(), p, p, nbClust, false);

    ProbCond = Rcpp::as<arma::mat>(InputList[5]);

    lambda = 0.0;        
    rho    = 0.0;

    Pk_cube.ones(p, p, nbClust);
    for (int k = 0; k < nbClust; ++k)
        Pk_cube.slice(k).diag().zeros();  // no penalty on diagonal
}

// --- Log Determinant ---
double Mixture::safe_log_det(const mat& M, bool& success) const
{
    double log_det_val = -std::numeric_limits<double>::infinity();
    success = false;
    // Check if matrix seems reasonable before attempting log_det
    if (M.is_finite() && M.is_square()) {
        mat R;
        // Use Cholesky decomposition to check for PD and get log_det
        if (chol(R, M)) { // chol returns true on success
            // log(det(M)) = log(det(R'*R)) = log(det(R')*det(R)) = 2 * log(det(R))
            // det(R) is product of diagonal elements since R is upper triangular
            // log(det(R)) = sum(log(diag(R)))
            log_det_val = 2.0 * accu(log(R.diag()));
            if (std::isfinite(log_det_val)) {
                 success = true;
            } else {
                // Handle case where diagonal elements of R might be zero or negative (shouldn't happen if PD)
                log_det_val = -std::numeric_limits<double>::infinity();
            }
        }
        // else: matrix is not positive definite
    }
    return log_det_val;
}

// --- Penalized Log-Likelihood (Weighted Penalty) ---
double Mixture::PenLogLik(void)
{
    mat lD = zeros<mat>(n, nbClust);
    double SLogDet_Sigma = 0.0; // Need log(det(Sigma)) for dmvnorm
    bool log_det_ok;

    for(int k = 0; k < nbClust; k++)
    {
        // Calculate log(det(Sigma_k)) from Omega_k for dmvnorm
        // log(det(Sigma)) = log(det(Omega^-1)) = -log(det(Omega))
        double log_det_Omega = this->safe_log_det(PrecisionMatrix.slice(k), log_det_ok);
        if (!log_det_ok) {
             // Handle non-PD Omega - assign very low likelihood
             lD.col(k).fill(-std::numeric_limits<double>::infinity());
             continue; // Skip likelihood calculation for this cluster
        }
        SLogDet_Sigma = -log_det_Omega; // log(det(Sigma))

        for(int i = 0; i < n; i++) {
            // Assuming ldcppmvt calculates log density using Omega and log(det(Sigma))
            // Note: The original code passed -log(det(Omega)) which is log(det(Sigma))
            lD(i,k) = log(prop(k)) + ldcppmvt(trans(Xd.row(i)), Mu.col(k), PrecisionMatrix.slice(k), SLogDet_Sigma);
        }
    }

    // Calculate observed log-likelihood using log-sum-exp trick
    double loglik_obs = 0.0;
    for (int i = 0; i < n; ++i) {
        rowvec log_dens_i = lD.row(i);
        if (!log_dens_i.is_finite()) {
             // Handle case where all densities are -Inf for point i
             loglik_obs += -std::numeric_limits<double>::infinity(); // Or some large negative value
             continue;
        }
        double max_log_d = log_dens_i.max();
        loglik_obs += max_log_d + log(accu(exp(log_dens_i - max_log_d)));
    }

    // Calculate Penalty Terms
    double mu_penalty = lambda * accu(abs(Mu)); // L1 penalty on means

    double precision_penalty = 0.0; // Weighted L1 penalty on Omega
    for(int k = 0; k < nbClust; k++) {
        // Element-wise product of weights and precision matrix, then sum absolute values
        precision_penalty += accu(abs(Pk_cube.slice(k) % PrecisionMatrix.slice(k)));
    }

    double penloglik = loglik_obs - mu_penalty - (rho * precision_penalty); // Subtract both penalties

    // Return large negative number if calculation failed
    return std::isfinite(penloglik) ? penloglik : -std::numeric_limits<double>::infinity();
};

// --- E-Step: Calculate Posterior Probabilities ---
void Mixture::GetProbCond(void){
    mat T = zeros<mat>(n, nbClust);
    mat lD = zeros<mat>(n, nbClust);
    double SLogDet_Sigma = 0.0;
    bool log_det_ok;

    for(int k = 0; k < nbClust; k++)
    {
        double log_det_Omega = this->safe_log_det(PrecisionMatrix.slice(k), log_det_ok);
        if (!log_det_ok) {
             lD.col(k).fill(-std::numeric_limits<double>::infinity());
             continue;
        }
        SLogDet_Sigma = -log_det_Omega;

        for(int i = 0; i < n; i++)
        {
            lD(i,k) = log(prop(k)) + ldcppmvt(trans(Xd.row(i)), Mu.col(k), PrecisionMatrix.slice(k), SLogDet_Sigma);
        }
    }

    // Normalize probabilities using log-sum-exp for numerical stability
    for(int i = 0; i < n; i++) {
        rowvec log_p = lD.row(i);
         if (!log_p.is_finite()) {
             T.row(i).fill(1.0 / nbClust); // Assign equal probability if all -Inf
             continue;
         }
        double max_log_p = log_p.max();
        rowvec p_norm = exp(log_p - max_log_p);
        double sum_p_norm = accu(p_norm);
        if (sum_p_norm > std::numeric_limits<double>::epsilon()) {
            T.row(i) = p_norm / sum_p_norm;
        } else {
            T.row(i).fill(1.0 / nbClust); // Handle sum being zero
        }
    }
    ProbCond = T;
};

// --- M-Step: Update Proportions ---
void Mixture::GetClassesSizes(void){
    prop = mean(ProbCond, 0);
    // Ensure proportions sum to 1 and are non-negative
    prop.elem( find(prop < 0) ).zeros();
    double sumProp = accu(prop);
    if (sumProp > std::numeric_limits<double>::epsilon()) {
        prop = prop / sumProp;
    } else {
        // Handle case where all proportions are zero (e.g., empty clusters)
        prop.fill(1.0 / nbClust);
    }
}

// --- M-Step: Update Means (Penalized - Zhou et al. 2009 logic) ---
// This implementation matches the user-provided C++ code for the mean update.
void Mixture::UpdateMeans(void)
{
    long double T_num, Tabs_check; // Use long double for precision in sums
    rowvec ProbCond_col_sums = sum(ProbCond, 0);
    // Avoid division by zero if a cluster is empty
    ProbCond_col_sums.elem( find(ProbCond_col_sums < std::numeric_limits<double>::epsilon()) ).fill(std::numeric_limits<double>::epsilon());

    mat MusNew = zeros<mat>(p, nbClust);
    mat W; // Current precision matrix Omega_k
    colvec MuPrec, MuNew; // Previous and new mean vectors

    for(int k = 0; k < nbClust; k++){
        W = PrecisionMatrix.slice(k); // Get current Omega_k
        MuPrec = Mu.col(k);          // Get previous mu_k
        MuNew.zeros(p);             // Initialize new mu_k for this cluster

        // Check for valid diagonal elements in W
        vec W_diag = W.diag();
         if (any(W_diag <= std::numeric_limits<double>::epsilon())) {
             // If Omega_k,jj is zero or negative, cannot perform update. Keep old mean.
             // Rprintf("Warning: Non-positive diagonal in Omega_%d. Skipping mean update.\n", k+1);
             MusNew.col(k) = MuPrec;
             continue; // Skip to next cluster
         }


        // Coordinate ascent loop for each variable j
        for(int j = 0; j < p; j++){
            // Calculate Q_kj term (from Zhou et al. notation, slightly adapted)
            // Q_kj = sum_i z_ik * [ sum_{v!=j} (x_iv - mu_kv) * Omega_kvj + x_ij * Omega_kjj ]
            Tabs_check = 0.;
            for(int i = 0; i < n; i++) {
                double term_i = dot(Xd.row(i) - trans(MuPrec), trans(W.col(j))) + (MuPrec(j) * W(j,j));
                // Original C++ code term: dot(Xd.row(i) - trans(MuPrec), trans(W.col(j))) + (MuPrec(j)*W(j,j))
                // This calculates: (x_i - mu_k)^T * Omega_k_col_j + mu_kj * Omega_kjj
                // Which simplifies to: sum_{v=1..p} (x_iv - mu_kv)*Omega_kvj + mu_kj*Omega_kjj
                // Rearranging: sum_{v!=j} (x_iv - mu_kv)*Omega_kvj + (x_ij - mu_kj)*Omega_kjj + mu_kj*Omega_kjj
                // Rearranging: sum_{v!=j} (x_iv - mu_kv)*Omega_kvj + x_ij*Omega_kjj
                // This matches the definition of Q_kj needed for the check.
                Tabs_check += ProbCond(i,k) * term_i;
            }

            // Check threshold condition
            if(std::abs(Tabs_check) <= lambda) {
                 MuNew(j) = 0.0;
            } else {
                // Calculate numerator for soft-thresholding
                // T_num = Q_kj (as calculated in Tabs_check)
                T_num = Tabs_check;

                // Calculate denominator: nk * Omega_kjj
                double denominator = ProbCond_col_sums(k) * W(j,j);

                // Apply soft-thresholding update
                // MuNew(j) = S_lambda(T_num) / denominator
                if(T_num > 0.0) {
                    MuNew(j) = (T_num - lambda) / denominator;
                } else {
                    MuNew(j) = (T_num + lambda) / denominator;
                }
            }
        } // end loop j
        MusNew.col(k) = MuNew; // Store updated mean for cluster k
    } // end loop k

    Mu = MusNew; // Update the class member Mu
};

// --- M-Step: Calculate Empirical Covariance ---
void Mixture::GetEmpiricalCovariance(void){
    rowvec ProbCond_cols_sums = sum(ProbCond, 0);
    ProbCond_cols_sums.elem( find(ProbCond_cols_sums < std::numeric_limits<double>::epsilon()) ).fill(std::numeric_limits<double>::epsilon());

    mat S = zeros<mat>(p,p);
    mat centered_Xd = Xd; // Temporary matrix for centered data

    for(int k = 0; k < nbClust; k++)
    {
        S.zeros(p,p);
        // Center data efficiently for cluster k
        centered_Xd = Xd; // Reset
        centered_Xd.each_row() -= trans(Mu.col(k));

        // Calculate weighted sum: sum z_ik * (x_i - mu_k)(x_i - mu_k)^T
        for(int i = 0; i < n; i++) {
            S += ProbCond(i,k) * (trans(centered_Xd.row(i)) * centered_Xd.row(i));
        }

        // Calculate weighted sample covariance S_k = Sum / n_k
        EmpiricalCovariance.slice(k) = S / ProbCond_cols_sums(k);
    }
}

// --- M-Step: Update Covariance/Precision Matrices (Weighted Penalty) ---
void Mixture::UpdateCovarianceMatrices(void){
    // Get the glassoFast function from the R environment
    Environment Rcpp_pkg = Environment::namespace_env("Rcpp"); // Get Rcpp environment
    Function Rcpp_List = Rcpp_pkg["List"]; // Get Rcpp::List function

    Environment glassoFast_pkg("package:glassoFast");
    Function RglassoFast = glassoFast_pkg["glassoFast"];

    rowvec ProbCond_cols_sums = sum(ProbCond, 0);
    ProbCond_cols_sums.elem( find(ProbCond_cols_sums < std::numeric_limits<double>::epsilon()) ).fill(std::numeric_limits<double>::epsilon());

    for(int k = 0; k < nbClust; k++)
    {
        long double nk = ProbCond_cols_sums(k);
        if (nk <= 0) continue; // Skip empty clusters

        // Calculate the cluster-specific penalty matrix for glassoFast
        // rho_matrix = (2 * rho / nk) * Pk
        mat rho_k_matrix = (2.0 * rho / nk) * Pk_cube.slice(k);
        // Ensure penalty is non-negative
        rho_k_matrix.elem( find(rho_k_matrix < 0) ).zeros();

        // Call R's glassoFast function
        try {
            List GlassoResult = RglassoFast(
                Named("S")   = EmpiricalCovariance.slice(k),
                Named("rho") = rho_k_matrix,
                Named("wi.init") = PrecisionMatrix.slice(k),
                Named("w.init")  = CovarianceMatrix.slice(k),
                Named("penalize.diagonal") = false,
                Named("thr")  = 1e-4,
                Named("maxIt")= 1000);

            mat W_new = as<mat>(GlassoResult["w"]);
            mat Wi_new = as<mat>(GlassoResult["wi"]);

            // Check for valid output
            if (W_new.is_finite() && Wi_new.is_finite()) {
                CovarianceMatrix.slice(k) = W_new;
                PrecisionMatrix.slice(k) = Wi_new;
            } else {
                // Keep previous estimate if glassoFast failed or returned non-finite
                 // Rprintf("Warning: glassoFast returned non-finite result for cluster %d. Keeping previous estimate.\n", k+1);
            }

        } catch (Rcpp::exception& e) {
             // Rprintf("Warning: glassoFast call failed for cluster %d: %s. Keeping previous estimate.\n", k+1, e.what());
             // Keep previous estimate on error
        } catch (...) {
             // Rprintf("Warning: glassoFast call failed for cluster %d with unknown error. Keeping previous estimate.\n", k+1);
             // Keep previous estimate on error
        }
    } // End loop k
}

// --- Get Variable Role (Checks if means are zero) ---
rowvec Mixture::VarRole(void){
    rowvec MuSum = zeros<rowvec>(p);
    rowvec alive = ones<rowvec>(p);

    // Sum absolute values of means across clusters for each variable
    MuSum = trans(sum(abs(Mu), 1));

    // Mark variable as inactive (0) if sum is close to zero
    for(int j = 0; j < p; ++j) {
        if(MuSum(j) < std::numeric_limits<double>::epsilon()) { // Use epsilon for robust comparison
            alive(j) = 0;
        }
    }
    return alive;
}

// --- Rcpp Export Function (Weighted Version) ---
//[[Rcpp::export]]
IntegerVector rcppClusteringEMGlassoWeighted(List InputList, double l, double r, arma::cube Pk_in, double tol = 1e-3, int max_iter = 250){
    // Create Mixture object using the constructor that accepts Pk_cube
    Mixture MyMixture(InputList, l, r, Pk_in);

    double PenLogLik_1 = 0.0, PenLogLik_0 = -std::numeric_limits<double>::infinity();
    int itr = 0;
    double relative_diff = std::numeric_limits<double>::infinity();

    // Calculate initial penalized log-likelihood
    PenLogLik_1 = MyMixture.PenLogLik();
    if (!std::isfinite(PenLogLik_1)) {
        // Rprintf("Warning: Initial penalized log-likelihood is non-finite (lambda=%.4f, rho=%.4f).\n", l, r);
        // Return vector of zeros if initialization fails
        IntegerVector result(MyMixture.getDimP(), 0);
        return result;
    }

    // EM Algorithm Loop
    while(relative_diff > tol && itr < max_iter)
    {
        PenLogLik_0 = PenLogLik_1;

        // --- E-step ---
        MyMixture.GetProbCond();
        if (!MyMixture.ProbOK()) break;
        // --- M-step ---
        MyMixture.GetClassesSizes(); // Update pi
        if (!MyMixture.ProbOK()) break;

        MyMixture.UpdateMeans(); // Update mu (penalized)
        if (!MyMixture.MuOK()) break;

        MyMixture.GetEmpiricalCovariance(); // Calculate S_k
        if (!MyMixture.EmpCovOK()) break;

        MyMixture.UpdateCovarianceMatrices(); // Update Omega_k / Sigma_k (weighted penalized)
        if (!(MyMixture.OmegaOK() && MyMixture.SigmaOK())) break;

        // --- Calculate Penalized Log-Likelihood ---
        PenLogLik_1 = MyMixture.PenLogLik();
         if (!std::isfinite(PenLogLik_1)) {
             // Rprintf("Warning: PenLogLik non-finite (iter %d, lambda=%.4f, rho=%.4f). Stopping.\n", itr + 1, l, r);
             PenLogLik_1 = PenLogLik_0; // Revert for relative diff calc
             break;
         }

        // --- Check Convergence ---
        if (std::abs(PenLogLik_0) < tol || !std::isfinite(PenLogLik_0)) {
             relative_diff = (std::abs(PenLogLik_1 - PenLogLik_0) < tol) ? 0.0 : 1.0;
        } else {
             relative_diff = std::abs(PenLogLik_1 - PenLogLik_0) / (1.0 + std::abs(PenLogLik_0));
        }

        itr++;
    } // End EM loop

    // Return the variable role vector (0 if mean is zero across all clusters, 1 otherwise)
    return wrap(MyMixture.VarRole());
};




