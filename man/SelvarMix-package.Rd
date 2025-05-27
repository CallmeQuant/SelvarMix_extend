\docType{package}
\name{SelvarMixExt-package}
\alias{SelvarMixExt-package}
\title{Enhanced Regularization for Variable Selection in Model-Based Clustering with Advanced Missing Value Handling}
\description{
  SelvarMixExt provides an extended framework for regularized variable selection 
  in model-based clustering and discriminant analysis with enhanced capabilities 
  for missing value imputation and MNAR (Missing Not At Random) mechanisms. 
  This package builds upon the foundational work of Sedki et al (2014) by 
  incorporating advanced imputation techniques including Gaussian copula methods,
   robust penalty parameter selection and MNAR handling.
  
  The enhanced procedure maintains the original two-step approach: first, variables are 
  ranked using a lasso-like procedure with flexible regularization 
  strategies; second, an adapted SRUW methodology is applied with support for 
  complex missing data patterns and modern clustering algorithms from Rmixmod, MixAll 
  and mclust frameworks.
}
\details{
  \tabular{ll}{
    Package: \tab SelvarMixExt\cr
    Type: \tab Package\cr
    Version: \tab 1.3.0\cr
    Date: \tab 2025-05-27\cr
    License: \tab GPL (>= 3)\cr
    LazyLoad: \tab yes\cr
  }
  
  The primary objective of SelvarMixExt is to perform robust variable selection 
  in model-based clustering frameworks, with particular 
  emphasis on handling complex missing data scenarios and high-dimensional datasets.
  
  The package uses the classical SRUW modeling framework of Maugis et al (2009) 
  which categorizes variables into four distinct roles:
  \itemize{
    \item \strong{S}: Relevant clustering variables that directly contribute to cluster structure
    \item \strong{R}: Subset of S variables used as regressors for redundant variables  
    \item \strong{U}: Redundant variables that can be explained by R variables
    \item \strong{W}: Independent variables that are uninformative for clustering
  }
  
  \strong{Key Enhancements:}
  \itemize{
    \item \strong{Advanced Pre-Imputation}: Gaussian copula-based imputation via gcimputeR 
          and robust alternatives using missRanger
    \item \strong{MNAR Support}: Handles Missing Not At Random mechanisms using 
          state-of-the-art methodologies
    \item \strong{Automatic Scaling}: Intelligent data conditioning based on 
          condition numbers and variance ratios
    \item \strong{Flexible Penalties}: Automatic grid search for optimal lambda 
          and rho parameters with customizable strategies
    \item \strong{Multi-Framework Support}: Integration with MixAll, mclust, and 
          Rmixmod clustering engines
    \item \strong{Robust Validation}: Enhanced model selection with multiple criteria 
          and performance metrics
  }
  
  The enhanced ranking procedure supports multiple regularization approaches including 
  group-wise shrinkage estimation with distance-based penalty weighting, making it 
  suitable for modern high-dimensional clustering challenges while maintaining 
  computational efficiency.
}

\examples{
\dontrun{
## Wine dataset analysis with enhanced features
## n = 178 observations, p = 27 variables
data(wine)

# Basic enhanced clustering with automatic parameter selection
set.seed(123)
result_basic <- SelvarClustLasso(x = wine[,1:27], 
                                nbcluster = 2:4, 
                                nbcores = 2)
summary(result_basic)
print(result_basic)

# Advanced usage with Gaussian copula imputation
result_copula <- SelvarClustLasso(x = wine[,1:27], 
                                 nbcluster = 2:4,
                                 use_copula = TRUE,
                                 scale_data = TRUE,
                                 impute_missing = TRUE)

# Multi-criteria selection with custom penalties  
result_multi <- SelvarClustLasso(x = wine[,1:27], 
                                nbcluster = 3,
                                criterion = c("BIC", "ICL"),
                                lambda = seq(0.1, 1.0, 0.2),
                                rho = seq(0.05, 0.5, 0.1))

# MNAR analysis for complex missing patterns
result_mnar <- SelvarClustLasso(x = wine[,1:27], 
                               nbcluster = 3,
                               use_missing_pattern = TRUE,
                               mnarz_control = list(mecha = "MNARz", 
                                                   rmax = 150))
}
}

\author{
  Binh H. Ho <binhho660@gmail.com>
  
  Extended from SelvarMix by Mohammed Sedki, Gilles Celeux, and Cathy Maugis-Rabusseau
}

\references{
  A. Sportisse, M. Marbac, F. Laporte, G. Celeux, C. Boyer, J. Josse, and C. Biernacki. 
  Model-based clustering with missing not at random data. Statistics and Computing, 
  34(4):135, 2024.
  
  A. Casa, A. Cappozzo, and M. Fop. Group-wise shrinkage estimation in penalized 
  model-based clustering. Journal of Classification, 39(3):648-674, 2022.
  
  Zhao, Y. and Udell, M. Missing value imputation for mixed data via Gaussian copula, 
  KDD 2020.
  
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

\keyword{package}
\keyword{clustering}
\keyword{variable-selection}
\keyword{missing-values}
\keyword{regularization}