lambda_vcurve <- function (INPUT, lg_lambdas, d = 2, IsPlot = FALSE, wFUN = wTSM, 
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

# Functions for working with the V-curve
v_point = function(y, w = 0 * y + 1, lambda = 100, d = 2) {
    # Compute the value of the normalized V-curve for one value of lambda
    # Prepare for smoothing
    n = length(y)
    E = diag.spam(n)
    D = diff(E, diff = d)
    P = t(D) %*% D

    # Smooth for  log-lambdas to the left and to the right
    z     = whit2(y, lambda, w)
    pz    = P %*% z #D' * D * z
    zgrad = lambda * log(10) * whit2(- pz/ w, lambda, w) #whit2(- pz * lambda, lambda, 1)
    # zgrad1 = whit2(-lambda * pz, lambda, w)

    fit   = sum(w * (y - z) ^ 2)
    dlfit = 2 * sum(-zgrad * w * (y - z)) / fit
    pen   = sum(z * pz)
    dlpen = 2 * sum(pz * zgrad) / pen

    # Take distance
    v = sqrt(dlfit ^ 2 + dlpen ^ 2)
    return(v)
}

# sometimes not converge
v_opt = function(y, w = 0 * y + 1, d = 2, lambdas = c(0, 4), tol = 0.01) {
    # Locate the optimal value of log10(lambda) with optimizer
    # Specify bounds of search range for log10(lambda) in paramter 'lambdas'

    v_fun = function(lla, y, w, d) v_point(y, w, 10 ^ lla, d)
    op = optimize(v_fun, lambdas, y, w, d, tol = tol)
    return(op$minimum)
}
