
# eDNAdetect

<!-- badges: start -->

<!-- badges: end -->

The goal of eDNAdetect is to predict the probability of detection given
a set of eDNA reads.

## Installation

You can install the development version of eDNAdetect from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("FreddieJH/eDNAdetect")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(eDNAdetect)
## basic example code

pred_detect(reads_vec = c(1000, 200, 400, 0, 0, 5000, 1500))
```
