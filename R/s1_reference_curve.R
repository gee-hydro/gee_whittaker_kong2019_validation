# Reference curve
# Dongdong Kong (2019-03-03)
#
# Input data was prepared by:
# 'test/07_whit/whit_lambda/02_whit_lambda_main.R'
source('R/main_whit_eval.R')

# file ="data_test/whit_lambda/MOD13A1_st_1e3_20180725.rda"
file ="data_test/whit_lambda/MOD13A1_st_1e3_20180731.rda"
load(file)

df[, doy := yday(t)]
setkeyv(df, c("site", "doy"))

## remove sites has no complete DOYs
info <- df[, .N, .(site, doy)][, .N, .(site)]
df   <- merge(df, info[N == 23, .(site)])

# time of reference curve
t          <- seq(1, 366, 16)
nyear      <- 1
nptperyear <- 23
wmin       <- 0.2

# noise in 2010 is regarded as real noise.

## -----------------------------------------------------------------------------
Fig1 = FALSE
if (Fig1){
    ## illustration of reference EVI
    data('MOD13A1')
    dt <- MOD13A1$dt
    d <- dt[site == "CA-NS6"] %>% tidy_MOD13.gee()
    d[, doy := yday(date)]

    ## get ref
    d_ref <- d[, .(ref=get_reference(y, QC_flag == "good", 0.85),
                   ref9=get_reference(y, QC_flag == "good", 0.9)), .(doy)]
    plot_ref() # visual refer curve

    # test noise simulation
    perc <- 0.5
    s <- simu_noise(d_ref, perc, trans = F)
    x <- s[[1]]

    # (r <- smooth_WHIT(x, 2, 4, TRUE))
    r <- ldply(s, function(x){
        smooth_WHIT_raw(x, NULL, 3, FALSE)
    })
}

## 1. Get reference curve and simulate noise  ----------------------------------
file_noise_random   <- sprintf("%s/noise3_random.rda", dir_whiteval)
file_noise_keypoint <- sprintf("%s/noise3_keypoint.rda", dir_whiteval) # and real gap

is_simu_noise <- TRUE
if (is_simu_noise){
    # df[, is_good := SummaryQA == 'good']
    df_2010 <- df[year(t) == 2010, .(site, doy, SummaryQA)]
    df_ref  <- df[, .(ref=get_reference(y, SummaryQA == "good", 0.85)), .(site, doy)]
    info_notBare <- df_ref[, .(max = max(ref),
                               A   = max(ref) - min(ref)), .(site)][max >= 0.1 & A >= 0.05, ]
    df_ref <- merge(df_ref, info_notBare)
    # merge real gap
    df_ref <- merge(df_ref[, .(site, doy, ref)], df_2010, by = c("site", "doy"), all.x = T)
    df_ref[is.na(SummaryQA), SummaryQA := "cloud"] # missing QC is regarded as "cloud"

    # merge(st, info_notBare[A >= 0.05])[, .N, .(IGBPname)]
    mat_ref <- dcast(df_ref, site~doy, value.var = "ref")

    ## SIMULATE NOISE
    # 1. random, 加入perc噪点
    lst_noise_random <- llply(1:3, function(i){
        perc <- percs[i]
        res <- dlply(df_ref, .(site), simu_noise, perc = perc, .progress = "text")
        df_sim <- purrr::transpose(res) %>% map(~do.call(rbind, .))
    }) %>% set_names(names(percs))
    save(df_ref, mat_ref, lst_noise_random, file = file_noise_random)

    # 2. types 2-4
    types <- c("real", "maxK", "maxDer")
    lst_noise_keypoint <- foreach(type = types, i = icount()) %do% {
        res <- dlply(df_ref, .(site), simu_noise, type = type, .progress = "text")
        df_sim <- purrr::transpose(res) %>% map(~do.call(rbind, .))
    } %>% set_names(types)
    save(df_ref, mat_ref, lst_noise_keypoint, file = file_noise_keypoint)
}
load(file_noise_random)
load(file_noise_keypoint)

## 2. TSM curve fitting methods (SG, DL, AG) -----------------------------------

Fig2 = FALSE
if (Fig2){
    InitCluster("log.txt")

    # 1. random
    r <- foreach(perc = percs, df_sim = lst_sim, i = icount(),
                 .packages = c("rTIMESAT")) %dopar% {
        # perc <- percs[i]
        # print(perc)
        TSF_main_df(df_sim, perc, nyear, nptperyear, overwrite = FALSE, wait = TRUE, indir = dir_whiteval)
    }

    # 2. keypoint
    perc = 0.1
    r <- foreach(df_sim = lst_noise_keypoint,
                 type   = names(lst_noise_keypoint),
                 i = icount(),
                 .packages = c("rTIMESAT")) %dopar% {
        # perc <- percs[i]
        # print(perc)
        indir <- sprintf("%s/%s", dir_whiteval, type)
        TSF_main_df(df_sim, perc, nyear, nptperyear, overwrite = FALSE, wait = TRUE, indir = indir)
    }

    ## tts to txt
    files <- dir(dir_whiteval, pattern = "*.tts", recursive = T, full.names = T)
    # keypoint files
    files <- files[-(7:15)]
    outfiles <- files %>% gsub("F:/whit_eval/|(/perc_10/TSF_whit_eval_noise10)|_fit.tts", "", .) %>%
        strsplit("_") %>%
        {sprintf("%s/fitting_%s_%s.txt", dir_whiteval, map_chr(., ~.[2]), map_chr(., ~.[1]))}


    r <- foreach(file = files, outfile = outfiles, i = icount()) %do% {
        # print(file)
        TSF_fit2time(file, 1, 1e8, 1, 1, wait = F,
                     outdir = dir_whiteval, outfile = outfile)
    }
}

# files_fix <- basename(files) %>% gsub(".txt", "", .) %>% strsplit("_") %>%
#     {sprintf("%s/fitting_%s_%s.txt", dirname(files), map_chr(., ~.[2]), map_chr(., ~.[1]))}

## 3. wWHIT curve fitting methods -----------------------------------
rm_rownames <- . %>% set_rownames(NULL) %>% data.table::data.table()

Fig3 = TRUE
if (Fig3){
    # 1. random
    lst_sim <- lst_noise_random
    NGRID <- lst_sim[[1]]$y %>% nrow()
    lambdas <- c(NA, 2, 15) %>% set_names(c("wWHd", "wWH2", "wWH15"))
    fit_wWHd <- foreach(lambda = lambdas, j = icount()) %do% {
        foreach(perc = percs, df_sim = lst_sim, i = icount(),
                 .packages = c("foreach", "iterators", "phenofit", "magrittr")) %dopar% {
            # perc <- percs[i]
            # print(perc)
            # TSF_main_df(df_sim, perc, nyear, nptperyear, overwrite = FALSE, wait = TRUE)
            foreach(y = iter(df_sim$y, "row"), qc = iter(df_sim$QC, "row"),
                    j = icount(), .combine = "rbind", .final = rm_rownames) %do% {
                phenofit::runningId(j, 10000, NGRID)
                w <- rep(1, nptperyear)
                w[qc == 3] <- wmin

                tryCatch(
                    # stop(simpleError(1)),
                    smooth_WHIT(y, w, lambda, 3, FALSE),
                    error = function(e){
                        cat(sprintf('[e] %d: %s', j, e$message))
                        return(y*NA_real_)
                    })
            }
        } %>% set_names(names(percs))
    } %>% set_names(names(lambdas))

    ## 2. key point
    lst_sim <- lst_noise_keypoint
    NGRID <- lst_sim[[1]]$y %>% nrow()
    lambdas <- c(NA, 2, 15) %>% set_names(c("wWHd", "wWH2", "wWH15"))
    fit_wWHd <- foreach(lambda = lambdas, j = icount()) %do% {
        foreach(types = types, df_sim = lst_sim, i = icount(),
                .packages = c("foreach", "iterators", "phenofit", "magrittr")) %dopar% {
                    # perc <- percs[i]
                    # print(perc)
                    # TSF_main_df(df_sim, perc, nyear, nptperyear, overwrite = FALSE, wait = TRUE)
                    foreach(y = iter(df_sim$y, "row"), qc = iter(df_sim$QC, "row"),
                            j = icount(), .combine = "rbind", .final = rm_rownames) %do% {
                                phenofit::runningId(j, 10000, NGRID)
                                w <- rep(1, nptperyear)
                                w[qc == 3] <- wmin

                                tryCatch(
                                    # stop(simpleError(1)),
                                    smooth_WHIT(y, w, lambda, 3, FALSE),
                                    error = function(e){
                                        cat(sprintf('[e] %d: %s', j, e$message))
                                        return(y*NA_real_)
                                    })
                            }
                } %>% set_names(types)
    } %>% set_names(names(lambdas))

    # 1208502 failed
    save(fit_wWHd, file = "fitting_wWHd_noises3.rda")
    save(fit_wWHd, file = "fitting_wWHd_keypoints.rda")
}

# d <- dt[site == "CA-NS6", .(y = EVI/1e4, date, doy = yday(date),
#                             SummaryQA, is_good = SummaryQA == 0)]
# d[, c("QC_flag", "w") := qc_summary(SummaryQA)]


fix_whit <- function(){
    I_missing <- c(1208502, 1212102, 1217502, 1219302, 1225502, 1247823, 1251202) # perc_50

    d_new <- matrix(NA, nrow = ngrid, ncol = 23)
    d_new[setdiff(I_all, I_missing), ] <- r_whit[[3]] %>% as.matrix()
    r_whit[[3]] <- data.table(d_new)
    names(r_whit) <- names(percs)
    fit_wWHd <- r_whit
}
