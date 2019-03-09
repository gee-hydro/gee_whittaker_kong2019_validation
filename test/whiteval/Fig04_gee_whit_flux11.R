library(maptools)
# install_github('kongdd/ggplot2')
# install_github('kongdd/plyr')
# source("test/load_pkgs.R")
source('G:/Github/phenology/phenology/phenofit/test/load_pkgs.R')
source("R/main_phenofit_test.R")
# source("test/stable/ggplot/geom_boxplot_no_outlier.R")
# source('R/plot_phenofit.R')

qc_values <- c("0", "1", "2", "3")
qc_levels <- c("good", "margin", "snow/ice", "cloud")
qc_colors <- c("grey60", "#00BFC4", "#F8766D", "#C77CFF") %>% set_names(qc_levels)
qc_shapes <- c(19, 15, 4, 17) %>% set_names(qc_levels)

infile <- file_cam
dir_gdrive <- "D:/Documents/Google 云端硬盘/whit"

# gee_Whittaker version:
# ------------------------------------------------------------------------------
# whit_fluxcam_wWH: latest

file <- "data_test/phenofit.rda"
if (!file.exists(file)) {
    # 1. gee whit
    # ---------------------------------------------------------
    # Modified 20180807, multiple Whittakers version were tested.
    # 'gcesapelo' and 'coaloilpoint' are located at the seaside.
    # It's normal sometimes nodata.

    ## rename meth
    sub_levels <- . %>% {
        level <- levels(.)
        mapvalues(., level, gsub("phenoflux166_|phenocam133_|_v105",
            "", level))
    }

    # dir_gdrive <- "D:/Document/GoogleDrive/whit"  #data/gee_phenofit/v2/

    files <- dir(dir_gdrive, "*.geojson", full.names = T)
    patterns <- str_extract(basename(files), ".*(?=_\\d{4}_)") %>%
        unique()

    df <- llply(patterns, function(pattern) readwhitMAT(dir_gdrive,
        pattern), .progress = "text") %>% set_names(patterns) %>%
        melt_list("meth")
    df$meth %<>% as.factor()

    whit_flux <- df[grep("phenoflux", meth), ] %>% data.frame() %>%
        data.table()
    whit_cam <- df[grep("phenocam", meth), ] %>% data.frame() %>%
        data.table()

    whit_flux$meth %<>% sub_levels
    whit_cam$meth %<>% sub_levels
    # fwrite(df_flux, 'data_test/gee_whit_phenoflux166.csv')
    # fwrite(df_cam , 'data_test/gee_whit_phenocam133.csv')

    # 2. merge others
    # ---------------------------------------------------------
    lst_cam  <- get_phenofit_fitting(file_cam, "data_test/phenofit_fitting_phenocam133.rda")
    lst_flux <- get_phenofit_fitting(file_flux, "data_test/phenofit_fitting_phenoflux166.rda")

    df_cam  <- get_phenofit_update_whit(lst_cam, whit_cam)
    df_flux <- get_phenofit_update_whit(lst_flux, whit_flux)

    # load
    load("data_test/phenofit_rough.rda")
    df_cam$fits_merge %<>% rbind(d_rough_cam$melt)
    df_cam %<>% get_phenofit_result_merge()

    df_flux$fits_merge %<>% rbind(d_rough_flux$melt)
    df_flux %<>% get_phenofit_result_merge()

    st_cam <- fread(file_st_cam)
    st_cam <- st_cam[!(site %in% c("cedarcreek", "rosemount"))]
    st_flux <- fread(file_st_flux)
    # save(df_cam, df_flux, st_cam, st_flux, file = file)
} else {
    load(file)
}

# FR-Pue, AU-Dry, US-KS2
sites_sel <- c("RU-Fyo", "GF-Guy", "US-UMB", "CN-Cha",
    "US-Whs", "AU-How", "ZA-Kru", "CN-HaM", "US-Los", "DE-Geb") # rm
{
    # st_flux <- st[site %in% sites]
    # st_flux[match(sites_sel, site), `:=`(selected = 1, label = sprintf("(%s) %s",
    #     letters[1:.N], IGBPname))]
    # df2sp(st_flux) %>%
    # rgdal::writeOGR('.', "st_flux101", drive = "ESRI Shapefile")
    # df2sp(st_cam) %>% writePointsShape('st_cam133.shp')
}

# load('D:/Documents/GoogleDrive/phenofit.rda')
methods <- c("AG", "BECK", "ELMORE", "ZHANG", "whit_R", "whit_gee")[-5]

##
i <- 1
prefix <- c("phenoflux", "phenocam")[i]
df <- if (i == 1) df_flux else df_cam
st <- if (i == 1) st_flux else st_cam

df <- df$all
df$SummaryQA %<>% factor(qc_levels)

st$IGBPname %<>% factor(IGBPnames_006)
st <- st[order(IGBPname, site), ]  # reorder according to IGBP
st[site %in% sites_sel, `:=`(titlestr, sprintf("(%s) %s, %s",
    letters[1:.N], site, IGBPname))]

# make sure different curve fitting methods have the same
# length fitting
formula <- if (i == 1) formula(site + date + t + y + GPP_NT +
    GPP_DT + SummaryQA + iters ~ meth) else formula(site + date +
    t + y + gcc + vci + SummaryQA + iters ~ meth)
IGBP.all = T
ylim2 <- if (i == 1) c(35, 100) else c(64, 100)

xlab <- st[, .N, IGBPname] %>% rbind(data.table(IGBPname = "all", N = nrow(st)))
xlab[, `:=`(label, sprintf("%s\n(n=%2d)", IGBPname, N))]
outfile <- sprintf("Fig8_valid_%s.pdf", prefix)

source("R/main_phenofit_test.R")
# over_perform(df[iters == 'iter2'], st, 'gpp', formula, prefix, xlab, ylim2, IGBP.all = IGBP.all, outfile)

# site figure data input
# df_trim <- melt(df_trim, measure.vars = methods, variable.name = 'meth')

# sites  <- unique(df$site)
sites    <- st$site
sitename <- sites[100]

# vars <- c('get_range', 'save_pdf', 'lgd_vci', 'lgd_gpp', 'qc_levels', 'qc_colors', 'qc_shapes', 'methods')
# cl <- cluster_Init(pkgs = c('data.table', 'ggplot2', 'magrittr'))
# clusterExport(cl, vars) res <- parLapplyLB(cl, sites,
# plot_whit, df_trim, st, prefix_fig = paste0('whit_', prefix))

st[, `:=`(ID, 1:.N)]
# check_Whittaker(sites, df_trim, 'gee_whit_flux166_v2.pdf') # all
# sites

## 2. select representative points for fluxnet

adj_envelope <- function(){
    df2 <- as.data.frame(df) %>% data.table()
    type = "y"; trs = 0.5
    df_adj <- ddply(df2, .(site), function(d){
        d[iters == "iter2", value := adjust_envelope(t, y, value, w, SummaryQA, type, trs)]
        return(d)
    })
    # all.equal(df_adj, df)
    return(df_adj)
}

# land cover homogeneous
st2 <- fread("file:///G:/Github/phenology/phenology/data/flux/station/st_flux97.csv")
sites2 <- st2[order(factor(IGBP, IGBPnames_006))]$site
sites <- c(sites2, sites_sel) %>% unique()
sites <- st[site %in% sites, ]$site

# adjust version
source('R/improve_Whittaker.R')
df_adj <- adj_envelope()
df_trim_adj <- dcast(df_adj, formula, value.var = "value", fun.aggregate = mean)  # %>% na.omit()

check_Whittaker(sites_sel, df_trim_adj, "gee_whit_flux11_adj2.pdf", "whit_fluxcam_wWH")  # all sites
# check_Whittaker(sites, df_trim_adj, "gee_whit_flux166_adj2.pdf", "whit_fluxcam_wWH")  # all sites

# Original version
Fig4_old = FALSE
if (Fig4_old) {
    df_trim     <- dcast(df    , formula, value.var = "value", fun.aggregate = mean)  # %>% na.omit()
    check_Whittaker(sites_sel, df_trim    , "gee_whit_flux11.pdf", "whit_fluxcam_wWH")  # all sites
}

Fig5 = TRUE
if (Fig5) {
    library(rTIMESAT)
    source('R/TSF_main.R')

    lim_year <- c(2001, 2017)
    lim_date <- lim_year %>% paste(c("-01-01", "-12-31"), sep = "") %>% as.Date()
    type = "y"; trs = 0.5 # adjust_envelope

    d_whit <- df[site %in% sites_sel & meth == "whit_fluxcam_wWH" & iters == "iter2" &
                     (date >= lim_date[1] & date <= lim_date[2])] # site == sitename &
    d_whit[, value := adjust_envelope(t, y, value, w, SummaryQA, type, trs)]

    d <- d_whit[, .(site, date, t, y, w, SummaryQA)] %>% unique()
    d$SummaryQA %<>% as.numeric()

    mat_y  <- dcast(d, date~site, value.var = "y") %>% .[, -1] %>% as.matrix() %>% t()
    mat_qc <- dcast(d, date~site, value.var = "SummaryQA") %>% .[, -1] %>% as.matrix() %>% t()


    sitename <- "US-KS2"
    d <- df[iters == "iter2" & site == sitename & meth == "whit_fluxcam_wWH", ]
    d <- d[date >= lim_date[1] & date <= lim_date[2] & !duplicated(date), ]
    d$SummaryQA %<>% as.numeric() # 1:4, good margin snow/ice cloud

    ## convert to QC_flag
    d[, c("QC_flag", "w") := qc_summary(SummaryQA-1, wmin = 0.2, wmid = 0.5, wmax = 0.8)]

    r_SG <- TSF_main(d, nptperyear = 23, FUN = 1, half_win = 10, cache = F)
    r_AG <- TSF_main(d, nptperyear = 23, FUN = 2, cache = F)
    r_DL <- TSF_main(d, nptperyear = 23, FUN = 3, cache = F)
}
## wWHIT
l <- check_input(d$t, d$y, d$w, nptperyear=23)
r_wWHIT <- wWHIT(l$y, l$w, l$ylu, nptperyear = 23,
                 wFUN = wBisquare, iters = 3, lambda = 10)
df_whit <- r_wWHIT$zs %>% as.data.frame()

pdat <- d[, wWHd := value]
pdat$SG <- r_SG$fit$v1
pdat$AG <- r_AG$fit$v1
pdat <- cbind(pdat, df_whit)
p <- show_fitting(pdat, methods = c("SG", "ziter1", "ziter2", "ziter3", "AG"),
                  colors = c("red", "blue", "green", "purple",  "black"), show.legend = T)
write_fig(p, "lambda_10.pdf", 10, 4)
# check_Whittaker(sitename, df_trim, "gee_whit_flux.pdf", "whit_fluxcam_wWH")
# merge_pdf('../whit_phenoflux166.pdf', indir = './')
# merge_pdf('whit_phenoflux166.pdf', indir = 'Figure/')
# merge_pdf('whit_phenocam133.pdf', indir = 'Figure',
# 'whit_phenocam.*.pdf', del = F)
# merge_pdf('whit_phenoflux166.pdf', indir = 'Figure',
# 'whit_phenoflux.*.pdf', del = F)

################################################################################ colnames(df)[4] <- 'raw'

## check impact of lambda and weigth updating
## names(table(df$meth)) %>% .[grep('*WH*', .)] %>%
## paste(collapse = '', '') %>% paste0(''', ., ''') %>% cat
## methods <- c('wWH', 'WH_p15', 'WH_p2', 'wWH_p15', 'wWH_p2')

# df[, `:=`( re = value - y, doy = yday(date) )]


## 1. autocorrelation shows that weighted function works
## better, and wWH_p2 is the best in all methods.  Initial
## lambda is useless here.  df_acf <- ddply(df[, ], .(meth,
## site), function(d){ with(d, { acf(-re, lag.max = 10, plot =
## F, na.action = na.pass)$acf[,,1][-1] }) }) %>% data.table()
## %>% melt(c('meth', 'site'), variable.name = 'lag')
## df_acf$lag %<>% factor() %>% as.numeric() %>% as.factor()
## ggplot(df_acf, aes(lag, value, color = meth)) +
## geom_boxplot()

# grid.newpage(); grid.draw(lgd) arrangeGrob(p, lgd, nrow =
# 2, heights = c(15, 1), padding = unit(1, 'line')) #return,

# p0 <- ggplot(d[meth == 'raw'], aes(date, value, shape =
# SummaryQA, color = SummaryQA)) + geom_point() +
# theme(legend.position = 'none') + scale_color_manual(values
# = c('good' = 'grey60', 'margin' = '#00BFC4', 'snow/ice' =
# '#F8766D', 'cloud' = '#C77CFF'), drop = F) +
# scale_shape_manual(values = c(19, 15, 4, 17), drop = F) +
# scale_y_continuous(lim = lim_raw) p3 <-
# ggplot_dual_axis(p1, p2) #%>% as.ggplot() p <-
# ggplot_dual_axis(p3, p0, add_yaxis_r = F)


## check marginal point
d <- df[meth == "AG" & iters == "iter1"]
d[, doy := yday(date)]

d_sel <- d[site %in% sites_sel, ]
ggplot(d_sel, aes(doy, y, color = SummaryQA, fill = SummaryQA)) +
    geom_point(aes(shape = SummaryQA)) +
    geom_smooth(data = d_sel[SummaryQA %in% c("good", "margin")]) +
    facet_wrap(~site, scales = "free_y") +
    scale_shape_manual(values = qc_shapes) +
    scale_color_manual(values = qc_colors) +
    scale_fill_manual(values = qc_colors)


