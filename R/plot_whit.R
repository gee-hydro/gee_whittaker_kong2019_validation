# EBF的表现最差，r2只有0.01
plot_mls <- function(l, title = NULL, ...) {
    # ymax = 200
    # yobs, pred
    xyplot(l$model[[1]], l$fitted.values, title, ...)
}

xyplot <- function(x, y, title = NULL, ...) {
    lim = c(0, 3)
    smoothScatter(x, y,
        xlim = lim, ylim = lim,
        ...,
        nrpoints = 10, 
        xlab = "Observation", ylab = "Predicted"
    )
    abline(a = 0, b = 1, col = "red", lwd = 1)
    if (!is.null(title)) legend("topleft", title, border = "transparent", bty = "n")
}
