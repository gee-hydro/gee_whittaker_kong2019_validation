# Evaluate performance according to Reference curve
# Dongdong Kong (2019-03-03)
#
# Input data was prepared by:
# 'test/07_whit/whit_lambda/02_whit_lambda_main.R'
source('R/main_whit_eval.R')

gof_fitting <- function(mat_ref.sm, mat_sim_raw, I_rem){
    ngrid <- nrow(mat_ref.sm)
    # I_all <- 1:(ngrid*100)

    if (missing(I_rem)) {
        mat_sim <- mat_sim_raw %>% as.matrix()
    } else {
        mat_sim <- matrix(NA, ngrid*100, 23)
        mat_sim[I_rem, ] <- mat_sim_raw %>% as.matrix()
    }

    foreach(yobs = iter(mat_ref.sm, "row"), i = icount(),
        .combine = "rbind", .final = rm_rownames) %do% {

        phenofit::runningId(i, 100, ngrid)
        I <- ((i-1)*100+1):(i*100)

        mat_sim_I <- mat_sim[I, ]
        foreach(ysim = iter(mat_sim_I, "row"), .combine = "rbind") %do% {
            phenofit::GOF(yobs, ysim)
        }
    }
}

## MAIN SCRIPTS ----------------------------------------------------------------
load(sprintf('%s/noise3_random.rda', dir_whiteval))

mat_ref.sm <- mat_ref[, 2:24] %>% as.matrix()
ngrid_org  <- nrow(mat_ref.sm)

## 1. wWHd ---------------------------------------------------------------------
s1_wWHd = FALSE
outfile_wWHd <- sprintf('%s/gof_wWHd_noises3.rda', dir_whiteval)

if (s1_wWHd) {
    load(sprintf('%s/fitting_wWHd_noises3.rda', dir_whiteval))

    gof_wWHd <- foreach(fit_meth = fit_wWHd) %do% {
        foreach(perc = percs, mat_sim = fit_meth, .final = set_names2,
                .packages = c("magrittr", "phenofit", "foreach", "iterators")) %dopar% {
            gof_fitting(mat_ref.sm, mat_sim)
        }
    }

    names(gof_wWHd) <- names(fit_wWHd)
    save(gof_wWHd, file = outfile_wWHd)
} else {
    load(outfile_wWHd)
}

map_dbl(gof_wWHd, ~mean(.$RMSE, na.rm = T))
# 10%        30%        50%
# 0.03241993 0.03294515 0.03603228
## 2. TSF ----------------------------------------------------------------------
s2_TSF = TRUE
if (s2_TSF) {

    split_TSF <- function(lst){
        meth <- names(lst) %>% str_extract(".{2}(?=_fit)")
        split(lst, meth) %>% map(~set_names(., names(percs)))
    }

    gof <- split_TSF(gof_TSF) %>% rev()
    gof <- c(gof_wWHd, list(HANTS = gof_HANTS), gof)

    # The percentage of MAE < 0.03
    map_dfr(gof, ~map_dbl(., ~nrow(.x[MAE < 0.03])/1294700))

    x <- purrr::transpose(gof_wWHd)
    x$`50%` %>% melt_list("meth") %>% cbind(I = 1:1294700, .) %>%
        dcast(I~meth, value.var = "NSE") %>% plot(wWHd~wWH2, .)

    # InitCluster(5)
    # split into different methods
    files = dir(dir_whiteval, pattern = "*.txt$", full.names = T)
    meths = str_extract(files, ".{2}(?=_fit)")
    files_lst <- split(files, meths)

    lst_TS <- llply(files, data.table::fread, skip = 1, .progress = "text")
    outfile_TSF <- sprintf('%s/gof_TSF_noises3_keypoint.rda', dir_whiteval)

    gof_TSF <- foreach(dat = lst_TS,
                       .packages = c("foreach", "iterators", "magrittr")) %dopar% {
        # foreach(file = files) %do% {
            # dat <- data.table::fread(file, skip = 1)
            I_rem <- dat$V1
            mat_sim <- as.matrix(dat[, 3:25])

            gof_fitting(mat_ref.sm, mat_sim, I_rem)
        # }
    } %>% set_names(basename(files) %>% gsub("fitting_|.txt", "", .))
    save(gof_TSF, file = outfile_TSF)
}
map_dbl(gof_TSF, ~mean(.$RMSE, na.rm = T))
# TSF_whit_eval_noise10_AG_fit.txt TSF_whit_eval_noise10_DL_fit.txt TSF_whit_eval_noise10_SG_fit.txt
# 0.05700367                       0.04537478                       0.03228564
# TSF_whit_eval_noise30_AG_fit.txt TSF_whit_eval_noise30_DL_fit.txt TSF_whit_eval_noise30_SG_fit.txt
# 0.05657278                       0.04659388                       0.03333651
# TSF_whit_eval_noise50_AG_fit.txt TSF_whit_eval_noise50_DL_fit.txt TSF_whit_eval_noise50_SG_fit.txt
# 0.05778334                       0.05086983                       0.03605201

## 3. TSF ----------------------------------------------------------------------
s3_HANTS = TRUE
if (s3_HANTS) {
    # InitCluster(5)
    # split into different methods
    files = dir(dir_whiteval, pattern = "*.csv$", full.names = T)
    # meths = str_extract(files, ".{2}(?=_fit)")
    # files_lst <- split(files, meths)

    lst <- llply(files, data.table::fread, .progress = "text")

    outfile_HANTS <- sprintf('%s/gof_HANTS_noises3_keypoints.rda', dir_whiteval)

    gof_HANTS <- foreach(mat_sim = lst,
                       .packages = c("foreach", "iterators", "magrittr")) %dopar% {
                           # foreach(file = files) %do% {
                           # dat <- data.table::fread(file, skip = 1)
                           # I_rem <- dat$V1
                           # mat_sim <- as.matrix(dat[, 3:25])

                           gof_fitting(mat_ref.sm, mat_sim)
                           # }
                       }
    save(gof_HANTS, file = outfile_HANTS)
}
map_dbl(gof_HANTS, ~mean(.$RMSE, na.rm = T))
