# SelvarMixExt
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

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
library(devtools)
install_github("CallmeQuant/SelvarMix_extend")
library(SelvarMixExt)
```

## Reference (To be updated)
