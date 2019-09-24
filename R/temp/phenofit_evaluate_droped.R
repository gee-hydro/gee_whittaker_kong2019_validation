GOF_fineFitting <- function(fit){
    d <- getFittings2(x)
    d[, as.list(GOF_extra2(y, value)), .(iters, meth)]
}

################################################################################
#############################   GOF functions  #################################
################################################################################
#' @export
tidy_rough_fitting <- function(x){
    # if (is.character(file)){
    #     x <- readRDS(file)
    # } else{
    #     x <- file
    # }
    x <- rm_empty(x)
    d_rough <- map(x, "whit") %>% melt_list("site") #%>% melt_list("meth")

    d_melt <- d_rough %>% .[, .(site, t, iter1 = ziter1, iter2 = ziter2)] %>%
        melt(id.vars = c("site", "t"),
             measure.vars = c("iter1", "iter2"), variable.name = "iters") %>%
        .[, .(site, t, iters, value)]

    # list(rough = d_rough, melt = d_melt)
    d_melt
}


# The neckness of speed is GOF
get_GOF3 <- function(d, iter = "iter1"){
    # d <- merge(df_org[, .(site, t, y0, w0, I_valid)], df_fit, by = c("site", "t"))
    # setkeyv(d, c("site", "t"))
    byname = c("site", "meth") %>% intersect(colnames(d)) # make sure by name exist

    # user  system elapsed
    # 18.87    0.03   19.09
    get_info <- function(is_valid){
        ddply_dt(d[I_valid == is_valid & iters == iter], .(GOF(y0, value)), byname)
    }

    list(cal = get_info(is_valid = 0), val = get_info(is_valid = 1)) %>%
        {.[sapply(., nrow) > 0]} %>% melt_list("type")
}

# get curve fitting results from phenofit object
getFittings2 <- function(fit){
    df_fit <- getFittings(fit) %>% data.table()

    whit <- melt(fit$seasons$whit[, .(t, y, w = witer2, iter1 = ziter1, iter2 = ziter2)],
                 measure.vars = c("iter1", "iter2"),
                 variable.name = "iters")
    whit$meth <- "whit_R"
    df_fit <- rbind(df_fit, whit)
    df_fit <- unique(df_fit) # remove duplicated value

    return(df_fit)
}



# Get GOF info for every
get_Fitting <- function(file){
    if (is(file, "list")){
        lst <- file
        is_fine <- FALSE
    } else if (is(file, "character")){
        lst  <- readRDS(file)
        is_fine <- grepl("phenofit", basename(dirname(file))) # Is fine curve fitting?
    }

    lst %<>% rm_empty()
    if (is_fine){
        # 1. read phenofit object
        df_fit <- llply(lst, getFittings2, .progress = "text") %>% melt_list("site")
    }else{
        # 2. read rough fitting
        df_fit <- tidy_rough_fitting(lst)
    }
    # info <- tryCatch(
    #     get_GOF3(df_org, df_fit), # GOF info
    #     error = function(e){ message(sprintf("%s : %s", file, e$message)) }
    # )
    df_fit
    # return(list(fit = df_fit, info = info)) # df_fit also input
}


get_GOF_fromFitting_I <- function(df_fit, df_org){
    d_perc <- df_org[, .(perc_good     = sum(w0 == 1)/.N,
                         perc_good_val = sum(w  == 1)/.N), .(site)]

    varnames <- c("site", "t", "iters", "value", "meth") %>%
        intersect(colnames(df_fit))
    df_fit <- df_fit[, ..varnames] %>% setkeyv(c("site", "t"))
    d      <- merge(df_org[, .(site, t, y0, w0, I_valid)], df_fit, by = c("site", "t"))

    setkeyv(d, c("site", "t"))

    iter1 <- suppressMessages(get_GOF3(d, "iter1"))
    iter2 <- suppressMessages(get_GOF3(d, "iter2"))

    ## 2. Roughness index and acf
    # if not include iters, info_rough will be error
    byname = c("site", "meth", "iters") %>% intersect(colnames(d)) # make sure `byname` exist

    info_rough <- ddply_dt(d, .(GOF_extra(y0, value)), byname)
    # vars <- c("Rg", "Rg_0")
    vars <- c("Rg", "Rg_norm_by_obs", "Rg_norm_by_pred" , "cv") #
    info_rough[, (vars) := lapply(.SD, function(x) map_dbl(x, first, default = NA)), .SDcols = vars]

    info <- listk(iter1, iter2) %>% melt_list("iters")
    info <- merge(info, d_perc, by = c("site")) # add percentage infomation
    list(info = info, rough = info_rough)
}

# get GOF from df_fitting RDS files
get_GOF_fromFitting <- function(file, df_org){
    lst <- readRDS(file)
    # a <- get_GOF_fromFitting_I(lst[[1]], df_org) # debug
    res <- llply(lst, get_GOF_fromFitting_I, df_org = df_org, .progress = "text")
    res
}


################################################################################
# merge phenofit INPUT, OUTPUT, GEE whittaker result and validation data
#
# Update 2018-08-09
# -----------------
# separate data preparing and tidying procedure

# get phenofit curve fitting result
# (input, fitting, validation data)
#' @export
get_phenofit_fitting <- function(infile, outfile, overwrite = F){
    if (file.exists(outfile) && !overwrite){
        load(outfile)
    } else {
        prefix  <- str_extract(infile, "\\w*(?=_MOD)")
        file_whit <- sprintf('data_test/gee_whit_%s.csv', prefix)
        ## 1. INPUT
        df_in   <- fread(infile, strip.white = T)                               # phenofit INPUT
        df_in$t    %<>% ymd()
        df_in$date %<>% ymd()

        ## 2.1 phenofit curve fitting
        lst     <- get_sbatch(paste0("Y:/github/phenofit_cluster/result/", prefix, "/"))
        df_temp <- llply(lst, getFittings2, .progress = "text")
        df_out  <- melt_list(df_temp, "site") %>% data.table() # phenofit OUTPUT

        ## 3. validation data
        valid_file <- sprintf("%svalid_%s_16day.csv", dir_data, prefix)
        df_valid   <- fread(valid_file)
        df_valid$date %<>% ymd()

        # remove NA values
        varnames <- colnames(df_valid)[4:5]
        eval(parse(text = sprintf("df_valid = df_valid[!(is.na(%s) & is.na(%s))]",
            varnames[1], varnames[2])))

        res <- list(data = df_in, fits = df_out, valid = df_valid)
        save(res, file = outfile)
    }
    res
}

get_phenofit_result <- function(infile, df_whit, df_out){
    ## 2.2 whit curve fitting
    # df_whit  <- readwhitMAT(dir_gdrive, prefix)            # GEE whit smoothing
    # add df_whit and merge
    df_list <- listk(df_in, df_out, df_valid)
    get_phenofit_update_whit(df_list, df_whit)

    ## 4. station info
    # st_file <- sprintf("%sst_%s.csv", dir_data, prefix)
    # st <- fread(st_file)   # station info
}

# merge gee_whit into phenofit
get_phenofit_update_whit <- function(df_list, df_whit){
    df_list$fits_merge <- with(df_list, {
        # df_whit      <- fread(file_whit); df_whit$date %<>% ymd()
        measure.vars <- colnames(df_whit) %>% .[grep("iter", .)]

        # unify curve fitting output format as phenofit
        df_whit <- merge(data[, .(site, t, date)], df_whit[, -1], by = c("site", "date")) #rm raw data in gee
        df_whit <- melt(df_whit, id.vars = c("site", "t", "meth"), measure.vars = measure.vars,
             variable.name = "iters")
        # df_whit$meth <- "whit_gee"

        rbind(fits[, .(site, t, iters, value, meth)],
             df_whit[, .(site, t, iters, value, meth)])
    })
    get_phenofit_result_merge(df_list)
}

#' merge phenofit list result
#' @param df_list List object returned from get_phenofit_result
get_phenofit_result_merge <- function(df_list){
    df_list$all <- with(df_list, {
        merge(data, fits_merge, c("site", "t"), all.x = T) %>%
            .[date < ymd(20180101) & date >= ymd(20000218)] %>%
            merge(valid[, c(1, 6, 4:5)], c("site", "date"), all.x = T)
    })
    df_list
}
