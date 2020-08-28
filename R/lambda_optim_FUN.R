
lambda_optim_FUN <- function(sitename, wFUN = wSELF, deltaT = 1, extend = FALSE, 
    IsPlot = FALSE, is_normalize = FALSE){
    lambda_optim(sitename, df = df_org, deltaT, extend,
                 IsPlot = IsPlot, IsSave = F, file = "whit_formual_wBisquare.pdf",
                 wFUN = wFUN, is_normalize = is_normalize)
}

#' According to v-curve theory, get the optimal lambda value.
#'
#' Whittaker balanced the fidelity and smooth. The agreement index maybe poor
#' than others. But it'is much smoothing.
#'
#' @param deltaT int, nyears chunk
#' @export
lambda_optim <- function(sitename, df, deltaT, extend = T,
    IsPlot = F, IsSave = F, file = "test_whit_lambda.pdf",
    wFUN = wBisquare, iters = 2, is_normalize = FALSE){
    # sitename <- sites[i]#; grp = 1
    YEAR_END = 2019
    nperiod <- ceiling(length(2000:YEAR_END) / deltaT)

    d     <- df[site == sitename]

    # normalize y
    if (is_normalize) {
        alpha = 0.05
        range = quantile(d$y, probs = c(alpha/2, 1-alpha/2), na.rm = TRUE)
        ynorm = (d$y - range[1])/(range[2] - range[1])
        # ynorm <- scale(d$y)
        d$y   <- ynorm
    }
    
    dnew  <- add_HeadTail(d) #
    INPUT <- check_input(dnew$t, dnew$y, dnew$w, dnew$QC_flag,
                         maxgap = nptperyear/4, alpha = 0.02, wmin = 0.2)
    if (length(unique(INPUT$y)) <= 5) return(NULL)

    years <- year(ymd(dnew$t))
    # cat(sprintf('site: %s ...\n', sitename))

    # res <- numeric(nperiod)*NA_real_
    res <- list()
    for (i in 1:nperiod){
        year      <- (i - 1)*deltaT + 2000
        year_beg  <- year
        year_end <- min(year + deltaT - 1, YEAR_END)

        year_beg_ext <- ifelse(extend, year_beg-1, year_beg)
        year_end_ext <- ifelse(extend, year_end+1, year_end)

        I     <- which(years >= year_beg & years <= year_end)
        I_ext <- which(years >= year_beg_ext & years <= year_end_ext)

        INPUT_i <- lapply(INPUT[c("t", "y0", "y", "w", "QC_flag")], `[`, I_ext) %>%
            c(INPUT[c("ylu", "nptperyear")])

        res[[i]] <- tryCatch({
            if (IsPlot) par(mfrow = c(2, 1), mar = c(2.5, 2.5, 1, 0.2),
                mgp = c(1.3, 0.6, 0), oma = c(0, 0, 0.5, 0))
            
            # lambda_vcurve theory used here
            vc <- lambda_vcurve(INPUT_i, lg_lambdas = seq(-1, 3, by = 0.01), d = 2,
                          wFUN = wFUN, iters = iters,
                IsPlot = IsPlot)

            ind <- match(I, I_ext)
            # browser()
            # vc$fit <- vc$fit[ind, ]
            vc$fit <- cbind(vc$fit, y = INPUT_i$y0)[ind, ]

            # coefficient to construct Whittaker lambda formula
            y <- INPUT_i$y
            # y <- dnew$y[I_ext] # fixed 20180824
            y <- y[!is.na(y)] # should be NA values now
            vc$coef <- list(mean = mean(y),
                            sd = sd(y),
                            kurtosis = kurtosis(y, type = 2),
                            skewness = skewness(y, type = 2)) #%>% as.list()
            vc
            # listk(lambda = vc$lambda) #, vc
        }, error = function(e){
            message(sprintf("[e] %s, %d: %s", as.character(sitename), i, e$message))
            #return(NA)
        })
    }

    tryCatch({
        ## visualization
        res %<>% rm_empty()
        # For bare land, all y is equal will lead lambda is empty
        df_sm  <- res %>% map_df("fit") %>% data.table()
        coefs  <- res %>% map_df("coef")
        coefs$lambda <- map_dbl(res, ~first(.$lambda, default = NA_real_))

        info <- merge(df_sm[, .(t, z)], d, by = "t") %$% GOF(y, z) %>% as.list() %>% as.data.table()

        if (IsSave){
            titlestr <- info[c(1:3, 5, 7)] %>%
                {sprintf("%s = %.2f", names(.), .) %>% paste(collapse = ", ") }
            cairo_pdf(file, 10, 4)
            par(mar = c(2.5, 2.5, 1, 0.2),
                mgp = c(1.3, 0.6, 0), oma = c(0, 0, 0.5, 0))
            plot_input(d, nptperyear)
            lines(ziter1~t, df_sm, col = "blue", lwd = 1.2)
            lines(ziter2~t, df_sm, col = "red", lwd = 1.2)
            title(titlestr)
            dev.off()
            file.show(file)
        }
        listk(data = df_sm, coef = coefs, gof = info) # return
    }, error = function(e){
        message(sprintf("[e] %s, %d: %s", as.character(sitename), i, e$message))
        #return(NULL)
    })
}
