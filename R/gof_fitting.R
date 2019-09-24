#' Good-of-fitting (GOF) of fitted time-series
#'
#' @param mat_ref matrix of reference curvr `[ngrid, ntime]`
#' @param mat_fit matrix of fitted time-series `[ngrid, ntime]`
#' @param I_rem index of succeed fitted pixels.
#' @inheritParams simu_noise
#'
#' @seealso [phenofit::GOF]
#'
#' @examples
#' data("mat_ref")
#' data("lst_noise")
#' d_gof <- gof_fitting(mat_ref, lst_noise$real$y)
#' @export
gof_fitting <- function(mat_ref, mat_fit, I_rem, times = 100){
    ngrid <- nrow(mat_ref)
    # I_all <- 1:(ngrid*100)

    if (missing(I_rem)) {
        mat_fit2 <- mat_fit %>% as.matrix()
    } else {
        mat_fit2 <- matrix(NA, ngrid*100, 23)
        mat_fit2[I_rem, ] <- mat_fit %>% as.matrix()
    }

    # For each site, repeat 100 times.
    n_pixels <- ceiling(nrow(mat_fit2)/times)

    d_gof <- foreach(yobs = iter(mat_ref, "row"), i = icount(n_pixels),
        .combine = "rbind", .final = rm_rownames) %do% {
        # runningId(i, 100, ngrid)
        I <- ((i-1)*times+1):(i*times)

        mat_fit_I <- mat_fit2[I, ]
        foreach(ysim = iter(mat_fit_I, "row"), .combine = "rbind") %do% {
            GOF(yobs, ysim)
        }
    }

    I_row <- ceiling((1:nrow(d_gof))/times)

    Roughness <- coef_roughness(mat_fit2)
    d_gof %<>% cbind(I=I_row, Roughness, .)
    d_gof
}


#' Roughness of matrix
#'
#' @return `Roughness = diff(x)^2`
#'
#' @param x matrix, `[ngrid, ntime]`
#' @examples
#' data("lst_noise")
#' coef_roughness(lst_noise$`50%`$y)
#' @export
coef_roughness <- function(x){
    if (!is.matrix(x))
        x <- as.matrix(x)
    rowDiffs(x, differences = 2)^2 %>% rowSums2(na.rm = T)
}
