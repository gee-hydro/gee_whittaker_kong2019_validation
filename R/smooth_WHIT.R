#' Whittaker smoother
#'
#' @inheritParams check_input
#' @param ... other parameters to [check_input()], e.g. `perc_wc`, `wmin`.
#' `wmin` also is passed to `wFUN`.
#'
#' @seealso [check_input()]
#'
#' @importFrom data.table last
#' @export
smooth_WHIT <- function(y, w, lambda, iters = 3,
    IsPlot = FALSE, ...)
{
    l <- check_input(t, y, w, nptperyear = 23, south = F, ...)

    if (missing(lambda) || is.na(lambda)){
        # if missing, init according to v-curve
        lambda <- lambda_init(l$y)
    }

    r <- wWHIT(l$y, l$w, l$ylu, nptperyear, wFUN = wBisquare,
               iters = iters, lambda = lambda, ...)

    if (IsPlot) plot_fitting(data.frame(y, w), r)
    # GOF(d_ref$ref, dplyr::last(r$zs))
    last(r$zs)
}


#' @rdname smooth_WHIT
#' @export
smooth_WHIT_df <- function(d, lambda, iters = 3, IsPlot = FALSE, ...)
{
    d <- c(d, qc_summary(d$QC))
    l <- check_input(t, d$y, d$w, d$QC_flag, nptperyear = 23, south = F, ...)

    if (missing(lambda) || is.null(lambda)){
        # if missing, init according to v-curve
        lambda <- lambda_init(l$y)
    }

    r <- wWHIT(l$y, l$w, l$ylu, nptperyear, wFUN = wBisquare,
               iters = iters, lambda = lambda, ...)

    if (IsPlot) plot_fitting(d, r)
    GOF(d_ref$ref, last(r$zs))
}

#' Weigthed Whittaker Smoother
#'
#' @inheritParams check_input
#' @param lambda whittaker parameter (2-15 is suitable for 16-day VI). Multiple
#' lambda values also are accept, then a list object return.
#' @param second If true, in every iteration, Whittaker will be implemented
#' twice to make sure curve fitting is smooth. If curve has been smoothed
#' enough, it will not care about the second smooth. If no, the second one is
#' just prepared for this situation. If lambda value has been optimized, second
#' smoothing is unnecessary.
#' @param ... other parameters to `wFUN`
#'
#' @references
#' 1. Eilers, P.H.C., 2003. A perfect smoother. Anal. Chem. https://doi.org/10.1021/ac034173t \cr
#' 2. Frasso, G., Eilers, P.H.C., 2015. L- and V-curves for optimal smoothing. Stat.
#'      Modelling 15, 91-111. https://doi.org/10.1177/1471082X14549288
#'
#' @seealso [check_input]
#'
#' @examples
#' \dontrun{
#' library(phenofit)
#' data("MOD13A1")
#' dt <- tidy_MOD13.gee(MOD13A1$dt)
#' d <- dt[site == "AT-Neu", ]
#'
#' l <- check_input(d$t, d$y, d$w, nptperyear=23)
#' r_wWHIT <- wWHIT(l$y, l$w, l$ylu, nptperyear = 23, iters = 2)
#' }
#' @export
wWHIT <- function(y, w, ylu, nptperyear, wFUN = wBisquare, iters=1, lambda=15,
    second = FALSE, ...) #, df = NULL
{
    trs <- 0.5
    if (all(is.na(y))) return(y)
    n <- sum(w)

    OUT   <- list()
    yiter <- y
    for (j in seq_along(lambda)){
        lambda_j <- lambda[j]

        fits <- list()
        ws   <- list()
        for (i in 1:iters){
            ws[[i]] <- w
            z <- whit2(yiter, lambda_j, w)
            w <- wFUN(y, z, w, i, nptperyear, ...)

            # If curve has been smoothed enough, it will not care about the
            # second smooth. If no, the second one is just prepared for this
            # situation.
            if (second) z <- whit2(z, lambda_j, w) #genius move

            ## Based on our test, check_ylu and upper envelope will decrease
            # `wWWHIT`'s performance (20181029).
            z <- check_ylu(z, ylu) # check ylu

            ylu <- range(z)
            zc <- ylu[1] + (ylu[2] - ylu[1])*trs

            I_fix <- z > yiter & z > zc
            yiter[I_fix] <- z[I_fix] # upper envelope

            # browser()
            fits[[i]] <- z
            # wnew <- wFUN(y, z, w, i, nptperyear, ...)
            # yiter <- z# update y with smooth values
        }
        fits %<>% set_names(paste0('ziter', 1:iters))
        ws   %<>% set_names(paste0('witer', 1:iters))

        OUT[[j]] <- list(zs = fits, ws = ws)
    }
    if (length(lambda) == 1) OUT <- OUT[[1]]
    return(OUT)
}

# # CROSS validation
# if (validation){
#     h   <- fit$dhat
#     df  <- sum(h)
#     r   <- (y - z)/(1 - h)
#     cv  <- sqrt( sum( r^2*w ) /n )
#     gcv <- sqrt( sum( (r/(n-df))^2*w ))
#     LV  <- whit_V(y, z, w) #L curve, D is missing now
#     OUT[[j]] <- c(list(data = as_tibble(c(list(w = w), fits)),
#         df = df, cv = cv, gcv = gcv), LV)
# }else{
# }
