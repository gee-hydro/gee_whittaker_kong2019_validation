#! Rscript
# 本程序用于：
# 1. 读取MODIS Terra-LAI在16000测试点上的数据，测试optimal lambda公式
# 2. 随后调用`s2_calibrate_lambda_equation`获得optimal lambda公式
# @author
# Dongdong Kong, 2020-08-29

# setwd("/media/kong/Various Data/Research/GEE_repos/gee_whittaker/")
library(phenofit)
library(rgdal)
devtools::load_all()
# args <- commandArgs(TRUE)
# print(args)
source("scripts/main_pkgs.R")
# nptperyear <- 23

# file_poly = "~/ArcGIS/huabei.shp"
file_poly = "O:/ChinaBasins/dem_raw/NorthChina/bou2_NorthChina_4p.shp"
file_sp = "data-raw/whit_lambda/shp/st_1e3_mask.shp"
poly <- read_sf(file_poly) %>% st_dissolve()
sp <- read_sf(file_sp)

inds = sapply(st_intersects(sp, poly),
       function(z) if (length(z)==0) NA_integer_ else z[1]) %>%
    which.notna()
st_test = st_data[inds, ] %>% data.table() %>% select(-index) %>%
    cbind(I = 1:nrow(.), .)
fwrite(st_test, "data-raw/st_test-NorthChina&GuangDong.csv")
sp[inds, ] %>% select(-index) %>%
    write_shp("data-raw/st_test-NorthChina&GuangDong.shp")

write_fig({
    par(mar = c(3, 3, 1, 1))
    plot(poly %>% st_geometry(), axes = TRUE)
    plot(sp[inds, 1], col = "red", axes = TRUE, add = TRUE)
}, "Figs1_Distribution of NorthChina&GuangDong tested points (sp111).pdf", 5, 5)

## retrieve MOD13A2 data from GEE

## 1. tidy GEE exported data ---------------------------------------------------
# file ="data_test/whit_lambda/MOD13A1_st_1e3_20180725.rda"
file = "data-raw/whit_lambda/terra_LAI_st-1e3_v20200610.RDS"
if (!file.exists(file)) {
    library(sf)
    # 原始数据较大，读取之后，放在了如下压缩文件中：
    # "data-raw/whit_lambda/raw-csv-v20200608/raw-csv-v20200608.zip"
    indir <- "data-raw/whit_lambda/raw-csv-v20200608/combined_LAI/"
    files <- dir(indir, "EVI.*.csv", full.names = T)
    files <- dir(indir, "LAI.*.csv", full.names = T)

    lst <- llply(files, fread, nrows = 5, .progress = "text")
    ncols = sapply(lst, ncol)
    I_bad = which(ncols < 11)
    files[I_bad] %>% basename()
    # "EVI_terra_MOD13A1_st_1e3_11_0m_buffer.csv"
    # [1] "LAI_terra_MOD15A2H_006_st_1e3_11_0m_buffer.csv"
    # [2] "LAI_terra_MOD15A2H_006_st_1e3_13_0m_buffer.csv"
    # [3] "LAI_terra_MOD15A2H_006_st_1e3_15_0m_buffer.csv"
    # [4] "LAI_terra_MOD15A2H_006_st_1e3_16_0m_buffer.csv"
    dt  <- do.call(rbind, lst) %>% unique()

    # d <- dt[1:1e3, ]
    # d[, c("QC_flag", "w") := qc_summary(SummaryQA)]
    dt[, `:=`(
        # t = str_sub(`system:index`, 1, 10) %>% ymd,
        index = str_sub(`system:index`, 12, 31),
        # w = qc_summary(SummaryQA, wmin = 0.2),
        # SummaryQA = factor(SummaryQA, levels = qc_values, labels = qc_levels),
        `system:index` = NULL,
        .geo = NULL
    )]
    setkeyv(dt, c("site", "date"))

    # site info
    st <- read_sf("data-raw/whit_lambda/shp/st_1e3_mask.shp") %>% data.table() %>% .[, 1:2]
    df <- merge(st, dt, by = c("index", "site"))[, -1]
    setkeyv(df, c("site", "date"))
    # df <- df[order(site), .(site, y = EVI/1e4, t, w, SummaryQA)]

    # remove sites less then 3y valid data
    info <- df[, .N, site][order(N)]
    # site_rm <- info[N < nptperyear*3, site]
    # I_rm <- which(df$site %in% site_rm)
    # df <- df[-I_rm, ]
    # saveRDS(df, "terra_EVI_st-1e3_v20200610.RDS")
    saveRDS(df, file)
    # save(df, st, file = file)
}

{
    df = readRDS(file)
    df = df[date <= "2019-12-31", .(site, date, y = Lai_500m / 10, FparExtra_QC)] # , FparLai_QC
    df[, c("QC_flag", "w") := qc_FparLai(FparExtra_QC)]

    # rm sites with valid VIobs <= 30%
    info = df[, .(perc = sum(!is.na(y))/.N), .(site)]
    sites = info[perc >= 0.3, site] # 2048 sites removed
    # df$date %<>% ymd()

    ################################################################################
    # fill missing template
    date = fill_missdate(8)
    temp = data.table(date, site = rep(sites, rep(length(date), length(sites))))
    # t as here is image date, other than pixel data.
    df_org = merge(df, temp, by = c("date", "site"), all.y = T) # fill missing values
    colnames(df_org)[1] = "t"
    # 13751
    # info = df_org[, .(perc = sum(!is.na(y)) / .N), .(site)]
    ################################################################################
    # 1. test whether lambda values are significant different among
    # different \code{dt}.
    ## parameters
    # deltaT    = 1
    # is_extent = F
    # subfix <- sprintf("_grp%d", ifelse(is_extent, deltaT, 0))
    deltaTs = c(1, 2, 3, 4, 6, 20) %>% set_names(., paste0("chunk", .))
    is_extends = c(TRUE, FALSE) %>% set_names(c("extend", "non-extend"))
    pars = expand.grid(deltaT = deltaTs, is_extend = is_extends)
    # i <- 1:nrow(pars) %>% set_names(pars$deltaT, )
    # set extent = false, it will not enclude previous and subsequent year' data.
    # sites    <- unique(df$site) %>% set_names(., .)
    sitename = sites[2]
    extend = FALSE
    # lambda_optim_FUN(sitename)
}

InitCluster(6, kill = FALSE)
is_normalize = FALSE
## the parameters' senetivity need to be validated -----------------------------
temp = foreach(deltaT = deltaTs, i = icount()) %do% {
    # runningId(i, prefix = glue("chunksize={deltaT}"))
    foreach(is_extend = is_extends, j = icount()) %do% {
        # runningId(j, prefix = names(is_extends)[j])
        extend_str <- ifelse(is_extend, "Extend", "nonExtend")
        subfix <- sprintf("chunk%02d_%s", deltaTs[i], extend_str)
        # print(subfix)
        if (is_normalize) {
            outfile = glue("OUTPUT/whit_lambda/normalized/v020_{subfix}.RDS")
        } else {
            outfile = glue("OUTPUT/whit_lambda/original/v020_{subfix}.RDS")
        }

        if (file.exists(outfile)) return()
        ans <- foreach(sitename = sites %>% set_names(., .), k = icount()) %dopar% {
            runningId(k, 10, prefix = subfix)
            tryCatch({
                d <- lambda_optim_FUN(sitename, deltaT = deltaT, extend = is_extend)
            }, error = function(e) {
                message(sprintf('%s', e))
            })
        }
        saveRDS(ans, outfile)
    }
}

# info <- map_depth(temp, 3, "gof") %>% melt_tree(c("chunk", "extend", "site"))
# info[, chunksize := str_extract(chunk, "\\d{1,2}") %>% as.numeric() %>% as.factor()]

# ggplot(info, aes(extend, NSE)) + geom_boxplot()
# ggplot(info, aes(chunksize, NSE)) + geom_boxplot()

# res <- par_sbatch(sites, lambda_optim_FUN, wFUN = wTSM, deltaT, is_extend,
#                   return.res = TRUE, Save = TRUE,
#                   outdir = paste0("OUTDIR/whit_lambda/v020", subfix))

# a <- lambda_optim_FUN(sitename)
# res <- lambda_optim_FUN(102)
# deltaT <- 1 # current is 4 at GEE
# res.bisquare <- lambda_optim(sitename, df, deltaT = 1, extent = T, IsPlot = F, IsSave = F,
#                     wFUN = wBisquare, file = "whit_formual_wBisquare.pdf")
# res.TSM <- lambda_optim(sitename, df, deltaT = 1, extent = T, IsPlot = F, IsSave = T,
#                     wFUN = wTSM, file = "whit_formual_wTSM.pdf")
# res.self <- lambda_optim(sitename, df, deltaT = 1, extent = T, IsPlot = F, IsSave = T,
#                     wFUN = wSELF, file = "whit_formual_wSELF.pdf")

# group = F # three year group
# outdir <- ifelse(group, '_grp', '')
# outdir <- paste0("result", outdir)

# cpus_per_node <- as.numeric(Sys.getenv('SLURM_CPUS_ON_NODE'))
#par_sbatch(sites, lambda_optim_FUN, wFUN = wBisquare, Save = T,
 #          outdir = paste0("result/whit_lambda/wBisquare", subfix) )
# par_sbatch(sites, lambda_optim_FUN, wFUN = wTSM, Save = T,
#            outdir = paste0("result/whit_lambda/wTSM", subfix) )

# res <- par_sbatch(sites, lambda_optim_FUN, wFUN = wBisquare,
#                   return.res = F, Save = T,
#             outdir = paste0("result/whit_lambda/wBisquare", subfix))
# res <- par_sbatch(sites, lambda_optim_FUN, wFUN = wTSM,
#                   return.res = F, Save = T,
#             outdir = paste0("result/whit_lambda/wTSM", subfix))
# res <- par_sbatch(sites, lambda_optim_FUN, wFUN = wSELF,
#                   return.res = F, Save = T,
#             outdir = paste0("result/whit_lambda/wSELF", subfix))

# a <- readRDS("result/whit_lambda/wSELF_grp1/results_0.RDS")
# a <- get_sbatch("result/whit_lambda/wSELF_grp1/")
#
# sapply(a, length) %>% {which(. == 1)}
#l
# res <- list()
# for (i in seq_along(sites)[101:1000]){
#     sitename <- sites[i]
#     res[[i]] <- lambda_optim_FUN(sitename, wSELF)
# }
# system.time({ res <- pbmclapply(sites, lambda_optim, df = df, mc.cores = cpus_per_node) })
# x$`system:index` %<>% str_sub( 1, 31)
