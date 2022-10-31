
# edgeCounter

<!-- badges: start -->
[![R-CMD-check](https://github.com/yeyuan98/edgeCounter/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yeyuan98/edgeCounter/actions/workflows/R-CMD-check.yaml)
[![BioC status](http://www.bioconductor.org/shields/build/release/bioc/edgeCounter.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/edgeCounter)
<!-- badges: end -->

`edgeCounter` is a R-package for counting edges (5'/3') of BAM alignment reads. 

For certain sequencing applications such as ATAC-seq, fragment edges, but not 
full-length fragments themselves, are physically meaningful. In these cases, simple 
genome coverage functions offered by widely adopted tools are not well-suited.

This package provides a streamlined solution for counting alignment read edges.

## Installation

You can install the development version of edgeCounter from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("yeyuan98/edgeCounter")
```

## Examples

Please refer to package vignettes for examples of `edgeCounter`.

## Issues

Please submit issues and/or pull requests as you see fit.
