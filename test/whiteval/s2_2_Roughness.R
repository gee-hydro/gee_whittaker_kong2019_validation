# Dongdong Kong (2019-03-06)
#
source('R/main_whit_eval.R')
# HANTS -------------------------------------------------------------------
files <- dir(dir_whiteval, pattern = "*.csv$", full.names = T, recursive = T)

rough_HANTS <- llply(files, function(file){
    x <- data.table::fread(file)
    coef_roughness(x)
}, .progress = "text") %>% set_names(basename(files))

save(rough_HANTS, file = "rough_HANTS.rda")

# TSF ---------------------------------------------------------------------
files = dir(dir_whiteval, pattern = "fitting.*.txt$", full.names = T, recursive = T)
# meths = str_extract(files, ".{2}(?=_fit)")
# files_lst <- split(files, meths)

rough_TSF <- function(file){
    dat <- data.table::fread(file, skip = 1)
    res <- rep(NA, ngrid)
    I_rem <- dat$V1
    mat_sim <- as.matrix(dat[, 3:25])

    rough <- coef_roughness(mat_sim)
    res[I_rem] <- rough
    res
}

rough_TSF <- llply(files, rough_TSF, .progress = "text") %>% set_names(basename(files))
save(rough_TSF, file = "rough_TSF.rda")

# wWHd --------------------------------------------------------------------

load("F:/whit_eval/noise_random/fitting_wWHd_random.rda")
fits_1 <- fit_wWHd %>% purrr::transpose()
load("F:/whit_eval/fitting_wWHd_keypoint.rda")
fits_2 <- fit_wWHd %>% purrr::transpose()

fits <- c(fits_1, fits_2)

rough_wWHd <- foreach(fitI = fits, i = icount()) %do% {
    runningId(i)
    foreach(mat = fitI) %do% {
        coef_roughness(mat)
    } %>% set_names(names(fitI))
} %>% set_names(names(fits))

save(rough_wWHd, file = "rough_wWHd.rda")

# summary -----------------------------------------------------------------

files <- dir(".", "rough.*.rda", full.names = T)
for (file in files) {
    load(file)
}

# according to methods
meths <- names(rough_HANTS) %>% str_extract("(?<=_).{4,5}(?=_)")
rough_HANTS_fix <- rough_HANTS %>% split(meths)

meths <- names(rough_TSF) %>% str_extract("(?<=_).{2}(?=_)")
rough_TSF_fix <- rough_TSF %>% split(meths)

rough_wWHd_fix <- purrr::transpose(rough_wWHd) %>% map(~.[c(6:4, 1:3)])


ROWNAMES <- c(c("10%", "30%", "50%"), c("maxDer", "maxK", "real"))
roughs <- c(rough_wWHd_fix, rough_HANTS_fix, rough_TSF_fix) %>%
    map(~set_names(., ROWNAMES)) %>%
    purrr::transpose()


rough <- map_dfr(roughs, ~map_dfr(.x, mean, na.rm = TRUE)) %>%
    .[, c(1, 4, 8, 6, 7, 2:3, 5)] %>%
    cbind(type = ROWNAMES, .)
rough_perc <- map_dfr(roughs, ~map_dfr(.x, ~sum(.<0.01, na.rm=TRUE)/ngrid)) %>%
    .[, c(1, 4, 8, 6, 7, 2:3, 5)] %>%
    cbind(type = ROWNAMES, .)

list(rough = rough,
     rough_perc = rough_perc) %>%
    writelist_ToXlsx("gof_whit_rough.xlsx")

# 3.2 diff ----------------------------------------------------------------

diff_lst <- foreach(lst = roughs) %do% {
    wWHd <- lst[[1]]
    map(lst[-1], ~ (.- wWHd)) #%>% melt_list("meth")
} %>% set_names(names(roughs))


trans_diff <- function(d){
    d <- data.table(d)
    methods <- c("HANTS", "SG", "AG", "DL")

    # browser()
    d <- d[meth %in% methods, ]
    d$meth %<>% factor(methods)
    d[order(type, meth), ]
}

{
    delta <- 0.002
    rough_levels <- c(-Inf, -delta, delta, Inf)

    diff_rough <- map(diff_lst, function(lst){
        r <- llply(lst, function(x){
            table(cut(x, rough_levels))/ngrid
        })
        do.call(rbind, r) %>% data.table() %>% cbind(meth = names(r), .)
    }) %>% melt_list("type")


    writelist_ToXlsx(list(rough = trans_diff(diff_rough)), "diff_rough.xlsx")
}

# diff_NSE  <- map(diff_lst, ~.x[, as.list(table(cut(NSE , NSE_levels))/ngrid*100), .(meth)])
