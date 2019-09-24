# Dongdong Kong (2019-03-06)
source("test/main_pkgs.R")

types <- c("maxDer", "maxK", "real")
add_types <- function(d){
    cbind(type = c(names(percs), types), d)
}

split_TSF <- function(lst){
    meth <- names(lst) %>% str_extract(".{2}(?=_fit)")
    split(lst, meth) %>% map(~set_names(., names(percs)))
}

add_percs <- function(d){
    cbind(type = names(percs), d[, c(1, 4, 6:8, 2,3, 5)])
}
# MAIN SCRIPTS ------------------------------------------------------------

# 1. random
files_random   <- dir(paste0(dir_whiteval, "/noise_random/"), "gof_*", full.names = T)
for (i in seq_along(files_random)){
    load(files_random[i])
}

gof_TSF_fix <- split_TSF(gof_TSF)[c(3, 1, 2)]
gof_random  <- c(gof_wWHd, list(HANTS = gof_HANTS[1:3], MWHA = gof_HANTS[4:6]), gof_TSF_fix) %>%
    map(~set_names(., names(percs))) %>%
    .[c(1, 4, 6:8, 2:3, 5)] %>%
    purrr::transpose()

# 2. keypoint
files_keypoint <- dir(dir_whiteval, "gof_*", full.names = T)
for (i in seq_along(files_keypoint)){
    load(files_keypoint[i])
}

gof_wWHd_fix <- map(gof_wWHd, ~.[3:1]) # fix the order of wWHd, names is incorrect

gof_keypoint <- c(list(HANTS = gof_HANTS[1:3], MWHA = gof_HANTS[4:6]),
                  split(gof_TSF, rep(c("AG", "DL", "SG"), each = 3))[c(3, 1, 2)],
                  gof_wWHd_fix) %>%
    map(~set_names(., types)) %>%
    .[c(6, 1, 3:5, 7:8, 2)] %>%
    purrr::transpose()

gof_lst <- c(gof_random, gof_keypoint)

rm(gof_random, gof_keypoint, gof_TSF, gof_TSF_fix, gof_wWHd, gof_wWHd_fix)
# main scripts ------------------------------------------------------------

## 3.1 mean
# map_dfr(gof_keypoint, ~map_dfr(., ~mean(.x$NSE, na.rm = T)))
RMSE      <- map_dfr(gof_lst, ~map_dfr(., ~mean(.x$RMSE, na.rm = T))) %>% add_types()
RMSE_perc <- map_dfr(gof_lst, ~map_dfr(., ~nrow(.x[RMSE < 0.05])/1294700)) %>% add_types()

## 3.2 diff
diff_lst <- foreach(lst = gof_lst) %do% {
    wWHd <- lst[[1]]
    map(lst[-1], ~ (.- wWHd)) %>% melt_list("meth")
} %>% set_names(names(gof_lst))


delta <- 0.02
RMSE_levels <- c(-Inf, -delta, delta, Inf)

delta <- 0.1
NSE_levels <- c(-Inf, -delta, delta, Inf)

diff_RMSE <- map(diff_lst, ~.x[, as.list(table(cut(RMSE, RMSE_levels))/ngrid), .(meth)]) %>% melt_list("type")
diff_NSE  <- map(diff_lst, ~.x[, as.list(table(cut(NSE , NSE_levels))/ngrid), .(meth)]) %>% melt_list("type")

trans_diff <- function(d){
    d <- data.table(d)
    methods <- c("HANTS", "SG", "AG", "DL")

    # browser()
    d <- d[meth %in% methods, ]
    d$meth %<>% factor(methods)
    d[order(meth), ]
}

list(RMSE = diff_RMSE, NSE = diff_NSE) %>% lapply(trans_diff) %>%
    writelist_ToXlsx("diff_RMSE.xlsx")

# ggplot(r$maxDer, aes(RMSE)) + geom_histogram() + facet_wrap(~meth) +
#     scale_x_continuous(limits = c(-0.1, 0.1))

# 1. 只取RMSE即可
# 2. Global percentage of pixels with bias b 0.03 for the five schemes


# 2. random

# list(RMSE = rbind(RMSE.random, RMSE),
#      perc = rbind(RMSE_perc.random, RMSE_perc)) %>%
#     writelist_ToXlsx("eval_whit.xlsx")

## add figures

get_envelope <- function(
    alphas = c(.05, .1, .25, .5),
    group  = .(meth, index, iters, grp_perc))
{
    alphas <- c(.05, .1, .25, .5) %>% set_names(., .)
    d_envelope <- llply(alphas, function(alpha){
        # print(alpha)
        expr <- substitute(quote(
            res <- ddply_dt(d, .(quantile_envelope(RMSE, alpha)), group)
        ), list(alpha = alpha, group = group))
        # print(eval(expr))
        browser()
        eval(eval(expr))
    })
    d_envelope
}

d_enve <- d_envelope %>% melt_list("alpha")



info <- df_ref[, .(N = sum(SummaryQA == "good")/23), .(site)]
perc_good <- info$N

d  <- gof_keypoint$real$wWHd %>% cbind(site = rep(info$site, each = 100), perc_good = rep(info$N, each = 100), .)

ggplot(d, aes(perc_good, RMSE)) + geom_point() +
    geom_smooth()
# # The percentage of MAE < 0.03
# map_dfr(gof, ~map_dbl(., ~nrow(.x[MAE < 0.03])/1294700))

# x <- purrr::transpose(gof_wWHd)
# x$`50%` %>% melt_list("meth") %>% cbind(I = 1:1294700, .) %>%
# dcast(I~meth, value.var = "NSE") %>% plot(wWHd~wWH2, .)
