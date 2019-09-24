#' adjust_envelope
#'
#' @description Approach to upper envelope. This function works site by site.
#'
#' @details Smoothed time-series by Whittaker is often underestimated in the peak of
#' growing season (abbreviated as peak season). In R, we can further constrain
#' by yearly ylu. But, it's difficult to implement in GEE.
#'
#' 1. Not outlier.
#' 2. In middle growing season.
#' 3. only good values.
#'
#' @inheritParams phenofit::check_input
#' @param t Date object, corresponding date of y. If specified, yearly $ylu$ will
#' be calculated, which is used to filter peak season.
#' @param QC_flag Quality Control (QC) flag of original vegetation time-series y.
#' @param ylu Numeric, lower and upper boundary of y.
#' @param type Character, ylu get by "y" or "z".
#'
#' @export
adjust_envelope <- function(t, y, z, w, QC_flag, type = "y", trs = 0.5,
    alpha = 0.02, wmin = 0.2)
{
    # two solution:
    # (1) ylu by z
    # (2) ylu by y
    if (type == "y") {
        INPUT <- check_input(t, y, w, QC_flag = QC_flag,
            # nptperyear = nptperyear, south = south, maxgap = nptperyear/4,
            alpha = alpha, wmin = wmin)
        ylu <- INPUT$ylu
    } else if (type == "z") {
        ylu <- c(min(z), max(z))
        # INPUT <- check_input(t, z, QC_flag = QC_flag,
        #     # nptperyear = nptperyear, south = south, maxgap = nptperyear/4,
        #     alpha = alpha, wmin = wmin)
    }

    ## 1&2. good values and in middle growing season
    I_good <- QC_flag=='good'
    I_peak <- z >= (ylu[1] + (ylu[2] - ylu[1])*trs)

    ## 2.1 outlier constrain: value
    re <- y - z
    sd_re   <- sd(re, na.rm = T)
    mean_re <- mean(re, na.rm = T)
    I_not_outlier <- abs(re - mean_re) <= 3*sd_re
    # this step not work for US-KS2

    ## 2.2 outlier constrain: 1th order der
    z_temp <- z
    I_temp <- which(I_good & I_peak)
    I_fix  <- I_temp[z[I_temp] < y[I_temp]]
    z_temp[I_fix] <- y[I_fix]

    der_sd   <- diff(z) %>% abs() %>% sd(na.rm = T)
    der_mean <- diff(z) %>% abs() %>% mean(na.rm = T)
    # print(mean_der)
    A = range(z, na.rm = T) %>% diff()

    der   <- c(0, diff(z_temp)) %>% abs()
    der_z <- der - der_mean   # normalized by mean value
    I_der <- der < pmin(A/3, max(der, na.rm = T)) # abs(der_z) <= 1*der_sd #&

    ## 3. Finally, replace values
    I <- which(I_good & I_peak & I_not_outlier & I_der) #  & I_der

    # browser(), also check for outlier
    I_fix <- I[z[I] < y[I]]
    z[I_fix] <- y[I_fix]
    z

    # visual_z <- function(){
    #     I=1:length(y); plot(I, y[I]);
    #     lines(z, col = "red", lwd = 2);
    #     lines(z_temp, col = "green", lwd = 1.5); grid()
    # }
    # browser()
}

# d <- df[site == sitename]
# d[iters == "iter2", adjust_envelope(t, y, value, w, SummaryQA, type, trs)]

# d[iters == "iter1", ] %>% plot(, .)
# devtools::install_github("kongdd/plyr")

# df
#    site       date          t      y   w SummaryQA iters     value             meth GPP_NT GPP_DT
# 1: AR-SLu 2000-02-18 2000-02-25 0.2668 1.0      good iter1 0.2776958           whit_R     NA     NA
# 2: AR-SLu 2000-02-18 2000-02-25 0.2668 1.0      good iter2 0.2832229           whit_R     NA     NA
# 3: AR-SLu 2000-02-18 2000-02-25 0.2668 1.0      good iter1 0.2865657 whit_fluxcam_wWH     NA     NA
# 4: AR-SLu 2000-02-18 2000-02-25 0.2668 1.0      good iter2 0.3003137 whit_fluxcam_wWH     NA     NA
# 5: AR-SLu 2000-02-18 2000-02-25 0.2668 1.0      good iter1 0.2862013           wHANTS     NA     NA
