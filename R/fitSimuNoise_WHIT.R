#' fitSimuNoise_WHIT
#'
#' @param lst_noise list of simulated random and keypoint noise, `list(list(y, QC))`
#' @param lambdas lamdas for Whittaker
#' @inheritParams smooth_WHIT
#' @inheritParams phenofit::check_input
#' @param limits If not null, only process the first `limits` pixel severing
#' to debug.
#'
#' @seealso [check_input]
#' @examples
#' data('lst_noise')
#' r <- fitSimuNoise_WHIT(lst_noise)
#' @export
fitSimuNoise_WHIT <- function(
    lst_noise,
    lambdas = NULL,
    nptperyear = 23, wmin = 0.2, ...,
    outdir = "./OUTPUT", is_save = FALSE,
    limits = NULL)
{
    if (is.null(lambdas))
        lambdas <- c(NA, 2, 15) %>% set_names(c("wWHd", "wWH2", "wWH15"))

    pkgs <- c("foreach", "iterators", "phenofit", "magrittr")
    # 1. random
    NGRID <- lst_noise[[1]]$y %>% nrow()

    fit_wWHd <- foreach(lambda = lambdas, j = icount()) %do% {
        foreach(df_sim = lst_noise, i = icount()) %do% {
            x <- fitSimuNoise_WHIT_sub(df_sim, lambda, NGRID, nptperyear, wmin,
                ..., limits = limits)
        }
    }

    if (is_save) {
        outfile <- sprintf("%s/%s", outdir, "fitting_wWHd_noises3&keypoint.rda")
        save(fit_wWHd, file = outfile)
    }
    fit_wWHd
}

# ' @examples
# ' \dontrun{
# ' fitSimuNoise_WHIT_sub(df_sim, lambda, NGRID, nptperyear, ...)
# ' }
fitSimuNoise_WHIT_sub <- function(df_sim, lambda, NGRID,
    nptperyear, wmin, ...,
    limits = NULL)
{
    if (is.null(limits)) limits <- NGRID

    foreach(y = iter(df_sim$y, "row"),
        qc = iter(df_sim$QC, "row"),
        .combine = "rbind",
        # .final = . %>% set_rownames(NULL),
        j = icount(limits)) %do%
    {
        runningId(j, 10000, NGRID)
        # w <- rep(1, nptperyear)
        # w[qc == 3] <- wmin
        tryCatch({
            w <- qc_summary(qc, wmin = 0.2, wmid = 0.5, wmax = 1)$w
            smooth_WHIT(y[1, ], w, lambda, 3, FALSE, wmin = wmin, ...)
        }, error = function(e){
            cat(sprintf('[e] %d: %s\n', j, e$message))
            y[1, ] * NA_real_
        })
    }
    # do.call(rbind, res)
}
