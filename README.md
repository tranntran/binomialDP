binomialDP: Differentially Private Inference for Binomial Data
================

## Overview

The `binomialDP` package makes it easy to perform differentially private
inference for a population proportion in R. The package offers a quick
solution to find differentially private uniformly most powerful (UMP)
tests, private p-values, and private confidence intervals. Also included
in the package are the CDF function and the random deviate generator of
the Truncated-Uniform-Laplace (Tulap) distribution.

## Installation

To install the package, run one of the following:

``` r
## The latest official version from CRAN
install.packages("binomialDP")
library(binomialDP)

## The latest development version from Github
install.packages("devtools")
devtools::install_github("tranntran/binomialDP")
```

## Examples

Please see **Get Started** for the package demonstration and working
examples.

## Citation

If you use this package for academic purposes, please cite this package
along with the paper on differentially private inference for binomial
data:

Jordan Awan, Aleksandra Slavkovic and Tran Tran (2020). binomialDP:
Differentially Private Inference for Binomial Data. R package version
1.0.0.

Awan, Jordan Alexander, and Aleksandra Slavkovic. “Differentially
Private Inference for Binomial Data.” Journal of Privacy and
Confidentiality 10, no. 1 (2020).
