# Evaluate performance according to Reference curve
# Dongdong Kong (2019-03-06)
#
source('R/main_whit_eval.R')

add_types <- function(d){
    cbind(type = types, d[, c(6, 1, 5, 4, 3, 2, 7, 8)])
}

split_TSF <- function(lst){
    meth <- names(lst) %>% str_extract(".{2}(?=_fit)")
    split(lst, meth) %>% map(~set_names(., names(percs)))
}

add_percs <- function(d){
    cbind(type = names(percs), d[, c(1, 4, 6:8, 2,3, 5)])
}
# MAIN SCRIPTS ------------------------------------------------------------

# 1. keypoint
files_keypoint <- dir(dir_whiteval, "gof_*", full.names = T)
for (i in seq_along(files_keypoint)){
    load(files_keypoint[i])
}

gof_wWHd_fix <- map(gof_wWHd, ~.[3:1]) # fix the order of wWHd

types <- c("maxDer", "maxK", "real")
gof_keypoint <- c(list(HANTS = gof_HANTS[1:3], MWHA = gof_HANTS[4:6]),
                  split(gof_TSF, rep(c("AG", "DL", "SG"), each = 3)),
                  gof_wWHd_fix) %>% map(~set_names(., types)) %>%
    purrr::transpose()


# 1. 只取RMSE即可
# 2. Global percentage of pixels with bias b 0.03 for the five schemes

# map_dfr(gof_keypoint, ~map_dfr(., ~mean(.x$NSE, na.rm = T)))
RMSE      <- map_dfr(gof_keypoint, ~map_dfr(., ~mean(.x$RMSE, na.rm = T))) %>% add_types()
RMSE_perc <- map_dfr(gof_keypoint, ~map_dfr(., ~nrow(.x[RMSE < 0.05])/1294700)) %>% add_types()

# 2. random
files_random   <- dir(paste0(dir_whiteval, "/noise_random/"), "gof_*", full.names = T)
for (i in seq_along(files_random)){
    load(files_random[i])
}

gof_TSF_fix <- split_TSF(gof_TSF) %>% rev()
gof_random  <- c(gof_wWHd, list(HANTS = gof_HANTS[1:3], MWHA = gof_HANTS[4:6]), gof_TSF_fix) %>%
    map(~set_names(., names(percs))) %>%
    purrr::transpose()


RMSE.random      <- map_dfr(gof_random, ~map_dfr(., ~mean(.x$RMSE, na.rm = T))) %>% add_percs()
RMSE_perc.random <- map_dfr(gof_random, ~map_dfr(., ~nrow(.x[RMSE < 0.05])/1294700)) %>% add_percs()


list(RMSE = rbind(RMSE.random, RMSE),
     perc = rbind(RMSE_perc.random, RMSE_perc)) %>%
    writelist_ToXlsx("eval_whit.xlsx")

# # The percentage of MAE < 0.03
# map_dfr(gof, ~map_dbl(., ~nrow(.x[MAE < 0.03])/1294700))

# x <- purrr::transpose(gof_wWHd)
# x$`50%` %>% melt_list("meth") %>% cbind(I = 1:1294700, .) %>%
# dcast(I~meth, value.var = "NSE") %>% plot(wWHd~wWH2, .)
