
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GEE Whittaker

<!-- badges: start -->

[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/kongdd/gee_whittaker?branch=master&svg=true)](https://ci.appveyor.com/project/kongdd/gee_whittaker)
<!-- badges: end -->

Non-parametric weighted Whittaker smoothing in GEE

> ### 为何gee_whittaker核心代码转为闭源
> 
> - 无人遵守开源代码协议GPL3
> - 同类文章作者多不愿公开源代码

## Calibrate and Validate in R

The development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
install.packages("phenofit")
devtools::install_github("kongdd/gee_whittaker")
library(whittaker)
```

## Whittaker in Google Earth Engine (GEE)

The following is the main GEE script of the simpler version Whittaker
and an examples which smoothed 4-day MODIS LAI images in PML\_V2
model.

Please note that there are four necessary steps when using this method,
also shown in the above example:

1.  Pre-process, unmask NA values and initialize weights
    
    If skip this step, it will lead to matrix dimensions not equal,
    matrix can’t be inversed …

2.  apply Whittaker method (a 2d image array returned)

3.  convert 2d array into multi-bands (every single date corresponds to
    band)
    
    There is a small trick you should know:
    
    you should not convert 2d image into ImageCollection. If so, for
    each exporting task (one date one task), matrix operation in
    Whittaker (the most time-consuming part) will be executed
    repeatedly. Then, it will lead to n times slower (n is the number of
    images).
    
    Just export multi-bands image directly\!

4.  EXPORT the smoothed result

Please note that smoothing algorithm costs lots of computing resource.
You can’t smooth imagecollection and do further calculation or analysis
right now in the same script in GEE. The best option is exporting
smoothed images first.

At last, you should select a appropriate lambda parameter carefully when
using this method\! But if you are process MODIS, you can have a try
about wWHd in my previous blog.

For Chinese users, you might interested about my another blog in
[知乎](https://zhuanlan.zhihu.com/p/76278313), I have explained some
technique details in it.

Finally, if you not satisfied the smoothed result, you can have a try
about [other weights updating
methods](https://github.com/kongdd/phenofit/blob/master/R/wFUN.R).

# Calibration and Validation

## Calibrate lambda equation

1.  `test/whit_lambda/02_whit_lambda_main.R`: Prepare input data,
    calibrate optional lambda based on `V-curve` theory.

2.  `test/whit_lambda/03_whit_lambda_formula.R`: calibrate the empirical
    lambda equation based on Multiple Linear Regression

## Validate Whittaker performance

1.  `test/s1_reference_curve.R:` Get reference curve, and used as
    benchmark to evaluate smoothing methods’ performance.

2.  `test/s2_evaluate_performance.R`: Evaluate performance of different
    smoothing methods

3.  `test/s4_parameters_sensitivity.R`: Parameter sensitivity analysis

# **References**

> \[1\]. Kong, D., Zhang, Y., Gu, X., & Wang, D. (2019). A robust method
> for reconstructing global MODIS EVI time series on the Google Earth
> Engine. ISPRS Journal of Photogrammetry and Remote Sensing, *155*
> (May), 13–24. <https://doi.org/10.1016/j.isprsjprs.2019.06.014>

# Acknowledgements

Keep in mind that this repository is released under a GPL3 license.
