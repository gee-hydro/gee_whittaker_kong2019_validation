# Evaluate performance according to Reference curve
# Dongdong Kong (2019-03-03)
#
# Input data was prepared by:
source("test/main_pkgs.R")

## MAIN SCRIPTS ----------------------------------------------------------------
load(sprintf('%s/noise3_random.rda', dir_whiteval))

mat_ref.sm <- mat_ref[, 2:24] %>% as.matrix()
ngrid_org  <- nrow(mat_ref.sm)

## 1. wWHd ---------------------------------------------------------------------
s1_wWHd = FALSE
outfile_wWHd <- sprintf('%s/gof_wWHd_noises3.rda', dir_whiteval)
outfile_wWHd <- sprintf('%s/gof_wWHd_keypoint.rda', dir_whiteval)

if (s1_wWHd) {
    # load(sprintf('%s/fitting_wWHd_noises3.rda', dir_whiteval))
    load(sprintf('%s/fitting_wWHd_keypoint.rda', dir_whiteval))

    gof_wWHd <- foreach(fit_meth = fit_wWHd) %do% {
        foreach(perc = percs, mat_sim = fit_meth, .final = set_names2,
                .packages = c("magrittr", "phenofit", "foreach", "iterators")) %dopar% {
            gof_fitting(mat_ref.sm, mat_sim)
        }
    }

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
# map_dbl(gof_TSF, ~mean(.$RMSE, na.rm = T))

## 3. HANTS ----------------------------------------------------------------------
eval_gof_HANTS <- function(indir, postfix = "keypoint"){
    files <- dir(indir, pattern = "*.csv$", full.names = T)
    names <- gsub('.csv$', '', basename(files))

    lst <- llply(files, data.table::fread, .progress = "text")
    # lst_rough <- llply(lst, coef_roughness) %>% set_names(names)

    outfile_HANTS <- sprintf('%s/gof_HANTS_%s.rda', indir, postfix)

    gof_HANTS <- foreach(mat_sim = lst,
        .export = c("gof_fitting", 'rm_rownames', "mat_ref.sm"),
        .packages = c("foreach", "iterators", "magrittr")) %dopar% {
           # foreach(file = files) %do% {
           # dat <- data.table::fread(file, skip = 1)
           # I_rem <- dat$V1
           # mat_sim <- as.matrix(dat[, 3:25])

           gof_fitting(mat_ref.sm, mat_sim)
           # }
        } %>% set_names(names)
    save(gof_HANTS, file = outfile_HANTS)
}

s3_HANTS = TRUE
if (s3_HANTS) {
    # InitCluster(5)
    # split into different methods
    # meths = str_extract(files, ".{2}(?=_fit)")
    # files_lst <- split(files, meths)
    # indir <- dir_whiteval

    eval_gof_HANTS(paste0(dir_whiteval, '/noise_random'), postfix = 'random')
    # eval_gof_HANTS(dir_whiteval, postfix = 'keypoint')
}
# map_dbl(gof_HANTS, ~mean(.$RMSE, na.rm = T))
