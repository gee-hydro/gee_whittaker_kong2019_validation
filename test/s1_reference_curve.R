# Reference curve
# Dongdong Kong (2019-03-03)
#
# 'test/07_whit/whit_lambda/02_whit_lambda_main.R'
source("test/main_pkgs.R")

# time of reference curve
t          <- seq(1, 366, 16)
nyear      <- 1
nptperyear <- 23
wmin       <- 0.2

## -----------------------------------------------------------------------------

## 1. Get reference curve and simulate noise  ----------------------------------
#  Noise in 2010 is regarded as real noise.
s1_simu_noise <- FALSE
if (s1_simu_noise){
    ## 1.1 load data
    # file ="data_test/whit_lambda/MOD13A1_st_1e3_20180725.rda"
    file <- "data_test/whit_lambda/MOD13A1_st_1e3_20180731.rda"
    load(file)

    df[, doy := yday(t)]
    setkeyv(df, c("site", "doy"))

    ## remove sites has no complete DOYs
    info <- df[, .N, .(site, doy)][, .N, .(site)]
    df   <- merge(df, info[N == 23, .(site)])

    ## 1.2 check reference curve
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
        r <- ldply(s, function(x){ smooth_WHIT_df(x, NULL, 3, FALSE) })
    }

    ## 1.3 simulate noise
    # df[, is_good := SummaryQA == 'good']
    df_2010 <- df[year(t) == 2010, .(site, doy, SummaryQA)]
    df_ref0  <- df[, .(ref=get_reference(y, SummaryQA == "good", 0.85)), .(site, doy)]
    info_notBare <- df_ref0[, .(max = max(ref),
                               A   = max(ref) - min(ref)), .(site)][max >= 0.1 & A >= 0.05, ]
    df_ref0 <- merge(df_ref0, info_notBare)
    # merge real gap
    df_ref0 <- merge(df_ref0[, .(site, doy, ref)], df_2010, by = c("site", "doy"), all.x = T)
    df_ref0[is.na(SummaryQA), SummaryQA := "cloud"] # missing QC is regarded as "cloud"

    # merge(st, info_notBare[A >= 0.05])[, .N, .(IGBPname)]
    df_ref  <- dcast(df_ref0, site~doy, value.var = "ref")
    mat_ref <- df_ref[, -1] %>% as.matrix()

    ## SIMULATE NOISE
    # 1. random, 加入perc噪点
    lst_noise_random <- llply(1:3, function(i){
        perc <- percs[i]
        res <- dlply(df_ref0, .(site), simu_noise, perc = perc, .progress = "text")
        df_sim <- purrr::transpose(res) %>% map(~do.call(rbind, .))
    }) %>% set_names(names(percs))
    save(mat_ref, df_ref, lst_noise_random, file = file_noise_random)

    # 2. types 2-4
    types <- c("real", "maxK", "maxDer")
    lst_noise_keypoint <- foreach(type = types, i = icount()) %do% {
        res <- dlply(df_ref0, .(site), simu_noise, type = type, .progress = "text")
        df_sim <- purrr::transpose(res) %>% map(~do.call(rbind, .))
    } %>% set_names(types)
    save(mat_ref, df_ref, lst_noise_keypoint, file = file_noise_keypoint)
}

load(file_noise_random)
load(file_noise_keypoint)

mat_ref <- df_ref[1:100, -1] %>% as.matrix()
use_data(mat_ref)
lst_noise <- c(map(lst_noise_random, function(x) {map(x, ~.[1:100, ])}),
       map(lst_noise_keypoint, function(x) {map(x, ~.[1:100, ])}))

#  2 curve fitting
## 2.1 TSM curve fitting methods (SG, DL, AG) -----------------------------------
# fitSimuNoise_TSF
is_fitSimuNoise_TSF = FALSE
if (is_fitSimuNoise_TSF) {
    InitCluster(3, "log.txt")
    fitSimuNoise_TSF(lst_noise_random, lst_noise_keypoint, percs,
        nyear, nptperyear, indir = dir_whiteval)
}

# files_fix <- basename(files) %>% gsub(".txt", "", .) %>% strsplit("_") %>%
#     {sprintf("%s/fitting_%s_%s.txt", dirname(files), map_chr(., ~.[2]), map_chr(., ~.[1]))}

## 3. wWHIT curve fitting methods -----------------------------------

Fig3 = TRUE
if (Fig3){
    # 1208502 failed
    res <- fitSimuNoise_WHIT(lst_noise_random, lst_noise_keypoint,
        nptperyear, wmin,
        outdir = "./OUTPUT", limits = NULL)
}

# d <- dt[site == "CA-NS6", .(y = EVI/1e4, date, doy = yday(date),
#                             SummaryQA, is_good = SummaryQA == 0)]
# d[, c("QC_flag", "w") := qc_summary(SummaryQA)]
