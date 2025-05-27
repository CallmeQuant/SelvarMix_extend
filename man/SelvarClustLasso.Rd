\name{SelvarClustLasso}
\alias{SelvarClustLasso}
\title{
Enhanced Regularization for Variable Selection in Model-Based Clustering with Advanced Imputation
}
\description{
This function provides an extended implementation of variable selection in model-based clustering 
with enhanced penalization integrated from Casa et al (2022) and missing values handling from Sportisse et al (2024). 
The method extends the original approach of Celeux et al (2019) by incorporating
imputation, flexible penalty parameter selection, and support for 
missing not at random (MNAR) mechanisms. The variable ranking employs penalized EM algorithms 
with customizable regularization strategies.
}
\usage{
SelvarClustLasso(x, nbcluster, strategy = NULL, lambda = NULL, rho = NULL, 
                 num_vals_penalty = 5, type = "lasso", hsize = 3, 
                 criterion = "BIC", models = "gaussian_pk_sjk", 
                 rmodel = c("LI", "LB"), imodel = c("LI", "LB"), 
                 nbcores = min(2, detectCores(all.tests = FALSE, logical = FALSE)),
                 impute_missing = TRUE, use_copula = TRUE, scale_data = TRUE,
                 scale_check_method = "pairwise.complete.obs",
                 use_missing_pattern = FALSE, use_diag = TRUE,
                 true_labels = NULL, sd_ratio_threshold = 10,
                 cond_number_threshold = 30, rank = NULL,
                 rank_control = list(), mnarz_control = list())
}

\arguments{
  \item{x}{
    Matrix or data frame of continuous data. Rows correspond to observations 
    and columns correspond to variables. Missing values are supported.
  }
  \item{nbcluster}{
    Numeric vector specifying the number of clusters to consider (must be positive integers).
  }
  \item{strategy}{
    Optional MixAll strategy object for EM algorithm fitting. If NULL, a default strategy 
    with initialization is used.
  }
  \item{lambda}{
    Numeric vector of tuning parameters for L1 mean penalty. If NULL, automatically 
    determined using specified number of grid values.
  }
  \item{rho}{
    Numeric vector of tuning parameters for L1 precision matrix penalty. If NULL, 
    automatically determined using specified number of grid values.
  }
  \item{num_vals_penalty}{
    Integer specifying the number of penalty parameter values to consider in grid  
    when lambda or rho are not provided. Default is 5.
  }
  \item{type}{
    Character string defining the ranking procedure type. Must be "lasso" or "likelihood". 
    Default is "lasso". Currently, likelihood is not implemented.
  }
  \item{hsize}{
    Integer parameter controlling the strength of forward and backward algorithms 
    for selecting S and W variable sets. Default is 3. 
  }
  \item{criterion}{
    Character vector defining model selection criteria. Possible values: "BIC", "ICL", 
    or c("BIC", "ICL"). Default is "BIC".
  }
  \item{models}{
    Character string or Rmixmod Model object specifying the mixture models to fit. 
    Default is "gaussian_pk_sjk". Supports both Rmixmod, MixAll and mclust model specifications.
  }
  \item{rmodel}{
    Character vector defining covariance matrix forms for linear regression of U variables 
    on R variables. Options: "LI" (spherical), "LB" (diagonal), "LC" (general). 
    Default is c("LI", "LB").
  }
  \item{imodel}{
    Character vector defining covariance matrix forms for independent W variables. 
    Options: "LI" (spherical), "LB" (diagonal). Default is c("LI", "LB").
  }
  \item{nbcores}{
    Integer specifying number of CPU cores for parallel computation. Default uses 
    min(2, detectCores()).
  }
  \item{impute_missing}{
    Logical indicating whether to impute missing values. Default is TRUE.
  }
  \item{use_copula}{
    Logical indicating whether to use Gaussian copula for imputation at ranking stage (requires gcimputeR) 
    or missRanger. Default is TRUE.
  }
  \item{scale_data}{
    Logical indicating whether to automatically scale data based on condition number 
    and standard deviation ratios. Default is TRUE.
  }
  \item{scale_check_method}{
    Character string specifying method for checking scaling necessity. Options: 
    "pairwise.complete.obs" or "median". Default is "pairwise.complete.obs".
  }
  \item{use_missing_pattern}{
    Logical indicating whether to use missing pattern analysis with MNAR mechanisms. 
    Default is FALSE.
  }
  \item{use_diag}{
    Logical indicating whether to use diagonal elements in MNAR analysis. Default is TRUE.
  }
  \item{true_labels}{
    Optional vector of true cluster labels for validation purposes. Used to compare 
    clustering performance when available.
  }
  \item{sd_ratio_threshold}{
    Numeric threshold for standard deviation ratio to trigger automatic scaling. 
    Default is 10.
  }
  \item{cond_number_threshold}{
    Numeric threshold for condition number to trigger automatic scaling. Default is 30.
  }
  \item{rank}{
    Optional integer vector specifying variable ranks. If NULL, ranks are computed 
    automatically using the specified ranking method.
  }
  \item{rank_control}{
    List of control parameters for ranking procedure. See details for available options.
  }
  \item{mnarz_control}{
    List of control parameters for MNAR analysis. Includes options for mechanism 
    specification and convergence criteria.
  }
}

\details{
The \code{rank_control} list can contain the following elements:
\itemize{
\item \code{group_shrinkage_method}: Method for group shrinkage ("weighted_by_dist_to_diag_W0", etc.)
\item \code{distance_method}: Distance metric available from shapes::distcov ("Euclidean", "Riemannian", etc.)
\item \code{lambda_omega_0}: Lambda parameter for omega regularization
\item \code{penalize_diag}: Whether to penalize diagonal elements
\item \code{epsilon_weighted_by_W0}: Epsilon value for weighted distance calculations
\item \code{initialize}: Initialization method ("hc", "kmeans")
\item \code{laplacian_target_type}: Target laplacian matrix for laplacian_spectral only  ("diag_Omega_hat", "identity")
\item \code{adj_threshold}: Threshold for adjacency matrix in spectral distance
\item \code{laplacian_norm_type}: Type of laplacian normalization ("symmetric", "unsymmetric")
}

The \code{mnarz_control} list supports:
\itemize{
\item \code{mecha}: Missing mechanism ("MNARz", "mixed")
\item \code{diag}: Use diagonal elements
\item \code{rmax}: Maximum iterations
\item \code{tol}: Convergence tolerance
}
}

\value{
Returns a selvarmixext object (or list of objects if multiple criteria specified) containing:
\item{S}{Selected set of relevant clustering variables}
\item{R}{Selected subset of regressor variables}
\item{U}{Selected set of redundant variables}
\item{W}{Selected set of independent variables}
\item{criterionValue}{Criterion value for the selected model}
\item{nbcluster}{Selected number of clusters}
\item{model}{Selected mixture model specification}
\item{rmodel}{Selected covariance form for regression}
\item{imodel}{Selected covariance form for independent variables}
\item{parameters}{Object containing all mixture parameters}
\item{regparameters}{Matrix of regression coefficients for redundant variables}
\item{proba}{Matrix of conditional cluster membership probabilities}
\item{partition}{Vector of cluster assignments using MAP rule}
\item{imputedData}{Matrix of imputed data in original scale}
\item{parametersMNARz}{MNAR analysis results (if applicable)}
}

\author{
Binh H. Ho <binhho660@gmail.com>

Extended from original SelvarMix package by Mohammed Sedki, Gilles Celeux, and Cathy Maugis-Rabusseau.
}

\references{
  A. Sportisse, M. Marbac, F. Laporte, G. Celeux, C. Boyer, J. Josse, and C. Biernacki. Model-1617
  based clustering with missing not at random data. Statistics and Computing, 34(4):135, 2024. 

  A. Casa, A. Cappozzo, and M. Fop. Group-wise shrinkage estimation in penalized model-based1496
  clustering. Journal of Classification, 39(3):648-674, 2022. 
  
  Zhao, Y. and Udell, M. Missing value imputation for mixed data via Gaussian copula, KDD 2020.
  
  Sedki, M., Celeux, G., Maugis-Rabusseau, C., 2014. 
  "SelvarMix: A R package for variable selection in 
  model-based clustering and discriminant analysis with 
  a regularization approach". Inria Research Report.
  
  Zhou, H., Pan, W., and Shen, X., 2009. "Penalized model-based 
  clustering with unconstrained covariance matrices". 
  Electronic Journal of Statistics, vol. 3, pp.1473-1496.
  
  Maugis, C., Celeux, G., and Martin-Magniette, M. L., 2009. 
  "Variable selection in model-based clustering: 
  A general variable role modeling". Computational 
  Statistics and Data Analysis, vol. 53/11, pp. 3872-3882.
}

\seealso{
\code{\link{SelvarLearnLasso}}, \code{\link{SortvarClust}}, \code{\link{SortvarLearn}}

Required packages: \code{gcimputeR}, \code{missRanger}, \code{MixAll}, \code{mclust}
}

\examples{
\dontrun{
# Basic usage with automatic parameter selection
data(wine)
set.seed(123)
result <- SelvarClustLasso(x = wine[,1:27], 
                          nbcluster = 2:4, 
                          nbcores = 2)
summary(result)

# With Gaussian copula imputation and custom scaling
result_copula <- SelvarClustLasso(x = wine[,1:27], 
                                 nbcluster = 2:4,
                                 use_copula = TRUE,
                                 scale_data = TRUE,
                                 sd_ratio_threshold = 15)

# Using multiple criteria
result_multi <- SelvarClustLasso(x = wine[,1:27], 
                                nbcluster = 2:4,
                                criterion = c("BIC", "ICL"))

# With custom penalty parameters
result_custom <- SelvarClustLasso(x = wine[,1:27], 
                                 nbcluster = 3,
                                 lambda = seq(0.1, 1.0, 0.1),
                                 rho = seq(0.05, 0.5, 0.05))

# Advanced usage with MNAR analysis
result_mnar <- SelvarClustLasso(x = wine[,1:27], 
                               nbcluster = 3,
                               use_missing_pattern = TRUE,
                               mnarz_control = list(rmax = 200, tol = 1e-5))
}
}

\keyword{clustering}
\keyword{variable selection}
\keyword{regularization}
\keyword{missing values}