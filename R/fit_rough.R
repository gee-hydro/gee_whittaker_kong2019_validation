#' rough_fitting
#' @param lambda Unless lambda is constant, lambda should be null.
#' @export
rough_fitting <- function(sitename, df, st, .FUN = wWHIT, lambda = NULL,
    IsPlot = FALSE, print = FALSE, ...)
{
    sp    <- st[site == sitename, ] # station point
    d     <- df[site == sitename, ] # get the first site data
    south <- sp$lat < 0
    titlestr <- with(sp, sprintf('[%03d,%s] %s, ', ID, as.character(site), IGBPname))
    cat(titlestr, "\n")

    tryCatch({
        # fill missing values
        # date : image date
        # t    : compositing date
        if (is_flux) {
            d <- merge(d, data.table(date), by = "date", all = T)
        } else {
            d <- merge(d, data.table(t = date), by = "t", all = T)
        }
        ############################################################################
        dnew  <- add_HeadTail(d, south, nptperyear)
        # 1. Check input data and initial parameters for phenofit
        INPUT <- check_input(dnew$t, dnew$y, dnew$w,
                             nptperyear, south = south,
                             maxgap = nptperyear/4, alpha = 0.02, wmin = 0.2)
        INPUT$y0 <- dnew$y

        ## 20180819 fixed lambda bug, lambda will overwrite new lambda in season_mov
        # if (is.null(lambda)) lambda <- init_lambda(INPUT$y)#*2w
        brks2 <- season_mov(INPUT,
                           rFUN = .FUN,
                           lambda = lambda, nf = nf, frame = frame,
                           adj.param = adj.param, # default is true
                           plotdat = d, IsPlot = IsPlot, print = print,
                           titlestr = titlestr, IsOnlyPlotbad = F, ...)
        brks2$GOF <- GOF_season3y(brks2)
        brks2
    # }, error = function(e){
    #     message(sprintf("[e]: %s, %s", titlestr, e$message))
    }, warning = function(w){
        message(sprintf("[w]: %s, %s", titlestr, w$message))
    })
}

#' GOF_season3y
#' GOF of season3y object
#' @export
GOF_season3y <- function(brks2){
    # browser()
    GOF_fun <- function(d){
        vars_iter <- d %>% contain("ziter")
        varnames  <- gsub("z", "", vars_iter)

        res <- map(vars_iter, function(varname){
            GOF(d$y, d[[varname]])
        }) %>% do.call(rbind, .) %>% 
            as.data.table() %>% cbind(iter = varnames, .)
        res
    }

    all  <- brks2$whit %>% GOF_fun()
    good <- brks2$whit[witer1 >= 1] %>% GOF_fun()

    list(all = all, good = good) %>% melt_list("type")
}
