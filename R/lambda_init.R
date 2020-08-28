#' Initial lambda value of Whittaker smoother
#' 
#' This function is only suitable for 16-day EVI time-series.
#' 
#' @param y Numeric vector
#' @param param linear Regression coefficients of lambda simulating Equation.
#' 
#' @examples
#' library(phenofit)
#' data("MOD13A1")
#' 
#' dt <- tidy_MOD13.gee(MOD13A1$dt)
#' st <- MOD13A1$st
#' 
#' sitename <- dt$site[1]
#' d      <- dt[site == sitename, ] # get the first site data
#' lambda <- lambda_init(d$y)
#' @export
lambda_init  <- function(y, param = NULL){
    y        <- y[!is.na(y)] #rm NA values
    mean     <- mean(y)
    sd       <- sd(y)
    cv       <- sd/mean
    skewness <- skewness(y, type = 2)
    kurtosis <- kurtosis(y, type = 2)
    # lambda was transformed by log10
    # lambda   <- 0.555484 + 1.671514*mean - 3.434064*sd - 0.052609*skewness + 0.009057*kurtosis
    # lambda   <- 0.555465 + 1.501239*mean - 3.204295*sd - 0.031902*skewness # Just three year result
    # lambda <- 0.831120 + 1.599970*mean - 4.094027*sd - 0.035160*cv - 0.063533*skewness # all year
    ## update 2018-07-31
    # lambda <- 0.7835 +1.5959*mean -4.0371*sd +0.0048*cv -0.1032*skewness +0.0036*kurtosis # yearly
    # lambda <- 0.8209 +1.5008*mean -4.0286*sd -0.1017*skewness -0.0041*kurtosis            # 4-year

    if (is.null(param)) {
        lambda <- 0.8209 + 1.5008 * mean - 4.0286 * sd - 0.1017 * skewness - 0.0041 * kurtosis # 4-year
        # lambda <- 0.9809 + 0.7247*mean - 2.6752*sd - 0.3854*skewness - 0.0604*kurtosis
    } else {
        lambda <- param$`(Intercept)` + param$mean*mean + param$sd*sd +
            param$skewness*skewness + param$kurtosis*kurtosis + param$cv*cv
    }
    # lambda <- 0.817783 + 1.803588*mean - 4.263469*sd - 0.038240*cv - 0.066914*skewness - 0.011289*kurtosis  #3y
    return(10^lambda)
}


# #' Initial lambda value of whittaker taker
# #' @param y Numeric vector
# #' @export
# lambda_init <- function(y) {
#     y <- y[!is.na(y)] # rm NA values
#     mean <- mean(y)
#     sd <- sd(y)
#     cv <- sd / mean
#     skewness <- skewness(y, type = 2)
#     kurtosis <- kurtosis(y, type = 2)
#     # lambda was transformed by log10
#     # lambda   <- 0.555484 + 1.671514*mean - 3.434064*sd - 0.052609*skewness + 0.009057*kurtosis
#     # lambda   <- 0.555465 + 1.501239*mean - 3.204295*sd - 0.031902*skewness # Just three year result
#     # lambda <- 0.831120 + 1.599970*mean - 4.094027*sd - 0.035160*cv - 0.063533*skewness # all year

#     ## update 2018-07-31
#     # lambda <- 0.7835 +1.5959*mean -4.0371*sd +0.0048*cv -0.1032*skewness +0.0036*kurtosis # yearly
#     # lambda <- 0.8209 + 1.5008 * mean - 4.0286 * sd - 0.1017 * skewness - 0.0041 * kurtosis # 4-year
#     # lambda <- 0.817783 + 1.803588*mean - 4.263469*sd - 0.038240*cv - 0.066914*skewness - 0.011289*kurtosis  #3y
#     return(10^lambda)
# }


## 4 year group, 2018-07-31
# Call:
# lm_fun(formula = lambda ~ mean + sd + skewness + kurtosis, data = d,
#     na.action = na.exclude)
#
# Residuals:
#      Min       1Q   Median       3Q      Max
# -2.32247 -0.41889  0.05531  0.40772  3.06478
#
# Coefficients:
#              Estimate Std. Error t value Pr(>|t|)
# (Intercept)  0.820927   0.006844 119.950  < 2e-16 ***
# mean         1.500799   0.017574  85.399  < 2e-16 ***
# sd          -4.028618   0.045536 -88.471  < 2e-16 ***
# skewness    -0.101746   0.003102 -32.796  < 2e-16 ***
# kurtosis    -0.004148   0.001112  -3.729 0.000192 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.5689 on 75480 degrees of freedom
# Multiple R-squared:  0.2165,    Adjusted R-squared:  0.2164
# F-statistic:  5213 on 4 and 75480 DF,  p-value: < 2.2e-16
################################################################################
################################################################################
################################################################################
## All year togather
# Call:
# lm_fun(formula = lambda ~ mean + sd + cv + skewness, data = d,
#     na.action = na.exclude)

# Residuals:
#     Min      1Q  Median      3Q     Max
# -2.4662 -0.4267  0.1394  0.4144  2.5824

# Coefficients:
#              Estimate Std. Error t value Pr(>|t|)
# (Intercept)  0.831120   0.021897  37.956  < 2e-16 ***
# mean         1.599970   0.089914  17.794  < 2e-16 ***
# sd          -4.094027   0.168844 -24.247  < 2e-16 ***
# cv          -0.035160   0.008459  -4.157 3.25e-05 ***
# skewness    -0.063533   0.007966  -7.976 1.62e-15 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.5851 on 15135 degrees of freedom
# Multiple R-squared:  0.1572,    Adjusted R-squared:  0.1569
# F-statistic: 705.5 on 4 and 15135 DF,  p-value: < 2.2e-16
#
#          term    estimate   std.error  statistic       p.value
# 1 (Intercept)  0.83112003 0.021897138  37.955647 3.306152e-301
# 2        mean  1.59997023 0.089914112  17.794428  4.039120e-70
# 3          sd -4.09402663 0.168843860 -24.247412 1.876008e-127
# 4          cv -0.03516045 0.008458787  -4.156677  3.246863e-05
# 5    skewness -0.06353256 0.007965580  -7.975886  1.620571e-15

## Three year fitting result
# Call:
# lm_fun(formula = lambda ~ (mean + sd + cv + skewness + kurtosis),
#     data = d, na.action = na.exclude)

# Residuals:
#     Min      1Q  Median      3Q     Max
# -2.6601 -0.4564  0.0490  0.4418  2.7551

# Coefficients:
#              Estimate Std. Error t value Pr(>|t|)
# (Intercept)  0.817783   0.010569  77.379  < 2e-16 ***
# mean         1.803588   0.042016  42.926  < 2e-16 ***
# sd          -4.263469   0.081937 -52.033  < 2e-16 ***
# cv          -0.038240   0.004041  -9.462  < 2e-16 ***
# skewness    -0.066914   0.003762 -17.785  < 2e-16 ***
# kurtosis     0.011289   0.001506   7.496 6.62e-14 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.6639 on 90639 degrees of freedom
# Multiple R-squared:  0.1541,    Adjusted R-squared:  0.1541
# F-statistic:  3303 on 5 and 90639 DF,  p-value: < 2.2e-16

#          term    estimate   std.error  statistic      p.value
# 1 (Intercept)  0.81778299 0.010568590  77.378628 0.000000e+00
# 2        mean  1.80358830 0.042016401  42.925816 0.000000e+00
# 3          sd -4.26346883 0.081937136 -52.033413 0.000000e+00
# 4          cv -0.03823967 0.004041253  -9.462331 3.080333e-21
# 5    skewness -0.06691403 0.003762368 -17.785083 1.216689e-70
# 6    kurtosis  0.01128888 0.001505916   7.496350 6.621350e-14
