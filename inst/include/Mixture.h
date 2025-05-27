#ifndef MIXTURE_H
#define MIXTURE_H
#include <RcppArmadillo.h>
#include <limits> 

class Mixture
{
private:
    /*  data and dimensions  */
    arma::mat        Xd;                // n × p data matrix
    int              n;                 // number of observations
    int              p;                 // number of variables
    int              nbClust;           // number of mixture components

    /*  model parameters  */
    arma::rowvec     prop;              // 1 × K   mixing proportions
    arma::mat        Mu;                // p × K   means
    arma::cube       CovarianceMatrix;  // p × p × K   Sigma_k
    arma::cube       PrecisionMatrix;   // p × p × K   Omega_k = Sigma_k^{-1}

    /*  helpers updated inside the EM  */
    arma::cube       EmpiricalCovariance; // p × p × K   weighted S_k
    arma::cube       Pk_cube;             // p × p × K   Casa weights
    arma::mat        ProbCond;            // n × K        posterior probs

    /*  penalties  */
    double           lambda;            // Lambda_mu  (mean penalty)
    double           rho;               // Lambda_Omega  (precision penalty)

    // numerically-robust log(det(M)) via Cholesky; returns −Inf on error
    double safe_log_det(const arma::mat& M, bool& ok) const;
public:
    /* ---------- constructors ---------- */

    // Common-shrinkage (Zhou-Pan-Shen) version
    explicit Mixture(Rcpp::List InputList);

    // Adaptive Casa weighting: pass P_k cube from R
    Mixture(Rcpp::List    InputList,
            double        lambda_mu,
            double        lambda_omega,
            arma::cube    Pk_in);

    ~Mixture(){}  // Destructor                                           

    /* ---------- EM core ---------- */
    double PenLogLik( void );
    void   GetProbCond( void );
    void   GetClassesSizes( void );
    void   UpdateMeans( void );
    void   GetEmpiricalCovariance( void );
    void   UpdateCovarianceMatrices( void );

    /* ---------- variable-ranking helper ---------- */
    arma::rowvec VarRole( void );

    /* ---------- getters for EM driver ---------- */
    inline int  getDimP()      const { return p; }
    inline bool ProbOK()       const { return ProbCond.is_finite(); }
    inline bool MuOK()         const { return Mu.is_finite(); }
    inline bool SigmaOK()      const { return CovarianceMatrix.is_finite(); }
    inline bool OmegaOK()      const { return PrecisionMatrix.is_finite(); }
    inline bool EmpCovOK()     const { return EmpiricalCovariance.is_finite(); }
    inline bool PropOK()       const { return arma::all(prop > 0.0); }

    /* ---------- misc ---------- */
    inline int  GetNbClust()   const { return nbClust; }
};

#endif /* MIXTURE_H */