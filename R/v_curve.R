v_curve <- function (INPUT, lg_lambdas, d = 2, IsPlot = FALSE, wFUN = wTSM, 
    iters = 2) 
{
    y <- INPUT$y
    w <- INPUT$w
    nptperyear <- INPUT$nptperyear
    if (length(unique(y)) == 0) 
        return(NULL)
    param <- c(INPUT, nptperyear = nptperyear, wFUN = wFUN, iters = iters, 
        second = FALSE, lambda = NA)
    fits = pens = NULL
    for (lla in lg_lambdas) {
        z = whit2(y, 10^lla, w)
        fit = log(sum(w * (y - z)^2))
        pen = log(sum(diff(z, diff = d)^2))
        fits = c(fits, fit)
        pens = c(pens, pen)
    }
    dfits = diff(fits)
    dpens = diff(pens)
    llastep = lg_lambdas[2] - lg_lambdas[1]
    v = sqrt(dfits^2 + dpens^2)/(log(10) * llastep)
    nla = length(lg_lambdas)
    lamids = (lg_lambdas[-1] + lg_lambdas[-nla])/2
    k = which.min(v)
    lambda = 10^lamids[k]
    z <- whit2(y, lambda, w)
    d_sm <- data.table(t = INPUT$t, z)
    vc <- list(lambda = lambda, vmin = v[k], fit = d_sm, optim = data.table(lg_lambda = lamids, 
        v = v))
    cal_COEF <- function(y) {
        y <- y[!is.na(y)]
        list(mean = mean(y), sd = sd(y), kurtosis = kurtosis(y, 
            type = 2), skewness = skewness(y, type = 2))
    }
    vc$coef_all <- cal_COEF(INPUT$y)
    if (IsPlot) {
        par(mfrow = c(2, 1), mar = c(2.5, 2.5, 1, 0.2), mgp = c(1.3, 
            0.6, 0), oma = c(0, 0, 0.5, 0))
        ylim = c(0, max(v))
        plot(lamids, v, type = "l", col = "blue", ylim = ylim, 
            xlab = "log10(lambda)")
        points(lamids, v, pch = 16, cex = 0.5, col = "blue")
        abline(h = 0, lty = 2, col = "gray")
        abline(v = lamids[k], lty = 2, col = "gray", lwd = 2)
        title(sprintf("v-curve, lambda = %5.2f", lambda))
        grid()
        plot_input(INPUT, wmin = 0.2)
        colors <- c("blue", "red")
        lines(vc$fit$t, last(vc$fit), col = "blue", lwd = 1.2)
    }
    vc
}
