
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Build
Status](https://api.travis-ci.org/raff-k/RainSlide.svg?branch=master)](https://travis-ci.org/raff-k/RainSlide)

# RainSlide: Rainfall-induced landslide analysis tools

This package provides valuable tools to analyze rainfall-induced
landslides. In RainSlide the following analytic tools are implemented:

  - **survey** based on [Bornaetxea et
    al. 2018](https://doi.org/10.5194/nhess-18-2455-2018): Calculation
    of the effective surveyed area. The original Python script was
    translated to R and adapted to fit the raster-requirements. For
    computation GRASS GIS is required.
  - **event** based on [Melillo et
    al. 2014](https://doi.org/10.1007/s10346-014-0471-3) and
    [2018](https://doi.org/10.1016/j.envsoft.2018.03.024): Objective
    reconstruction of rainfall events responsible for landslides. The
    rainfall threshold is calculated using Rainslide::thresh.
  - **thresh** based on [Peruccacci et
    al. 2012](https://doi.org/10.1016/j.geomorph.2011.10.005) and
    [Rossi et al. 2017](https://doi.org/10.1016/j.geomorph.2017.02.001):
    Calculation of cumulated event rainfall–duration (ED) thresholds
    using the frequentist method.

The programming code is based on the references, but may vary slightly
due to R-specifications or computational performance reasons.

## Installation

You can install the latest version of RainSlide from
[github](https://github.com/raff-k/RainSlide) with:

``` r
# install.packages("devtools")
devtools::install_github("raff-k/RainSlide")
```
