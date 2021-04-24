source("test/main_pkgs.R")

indir = "/mnt/i/Research/gee_Whittaker version2/ChinaCropPhen1km/"
files = dir(indir, "*.tif", full.names = TRUE)

library(stringr)
{
    kinds_crop <- str_extract(basename(files), "(?<=CHN_).{2,12}(?=_)")
    phenophase <- str_extract(basename(files), "(?<=_).{2,6}(?=_\\d{4})")
    table(phenophase)
}

lst_files <- split(files, kinds_crop) %>%
    {.[grep("Maize|Wheat", names(.))]}

InitCluster(6, kill = FALSE)
res = foreach(files = lst_files, i = icount()) %dopar% {
    # runningId(i)
    l <- foreach(infile = files, j = icount()) %do% {
        runningId(j, prefix = i)
        raster(infile) %>% as.array()
    }
    arr <- abind::abind(l, along = 3)
    arr_mean = apply_3d(arr, FUN = rowMeans2)
    arr_nvalid = length(l) - apply_3d(arr %>% is.na, FUN = rowSums2) # n_valid

    list(doy = arr_mean, n_valid = arr_nvalid)
}
# arr重新转为raster格式

# GR&EM: Greenup or Emergence
# TR: Transplanting
# ET: early rice
# LR: long rice
# SR: single rice

# 在这里挑选春小麦和夏玉米的区域, double cropping rice识别难度较大
inds <- files %>% .[grep("Maize|Wheat", basename(.))]
# save(res, file = "mean_phenodate_Maize&wheat-(2000-2015).rda", compress = FALSE)

loc = coordinates(r) %>% as.data.table() %>% set_colnames(c("lon", "lat"))
Maize_LOS <- res$Maize_MA$doy - res$Maize_V3$doy
Wheat_LOS <- res$Wheat_MA$doy - res$`Wheat_GR&EM`$doy

d_nvalid = map(res, "n_valid") %>% map(~ t(.x) %>% as.numeric) %>% as.data.table()
d_doy = map(res, "doy") %>% c(listk(Maize_LOS, Wheat_LOS)) %>%
    map(~ t(.x) %>% as.numeric) %>% as.data.table()

inds_good = {rowSums2(is.na(as.matrix(d_doy))) != ncol(d_doy) } %>% which()
df = cbind(loc[inds_good], d_doy[inds_good], n = d_nvalid[inds_good])

# 挑选双生长季的作物
d_dGS = df[!is.na(Maize_LOS + Wheat_LOS) & n.Maize_MA >= 10 & n.Wheat_MA >= 10]

prj = CRS("+proj=aea +lat_0=30 +lon_0=95 +lat_1=15 +lat_2=65 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
sp_dGS <- df2sp(d_dGS, prj = crs(r))
sp_84 = spTransform(sp_dGS, prj84)

{
    set.seed(1)
    inds_sm <- sample(1:nrow(sp_84), 1e2)
    sp_sm <- sp_84[inds_sm, ] %>%
        {.[which(.@coords[,1] > 100), ]}
}
# plot()
## 小麦和玉米种植都超过10年的pixels
write_shp(sp_84, "data-raw/ChinaPheno1km_10y_DoubleSeasonPixels_Maize&Wheat.shp")
write_shp(sp_sm, "data-raw/ChinaPheno1km_10y_DoubleSeasonPixels_Maize&Wheat_st1e2.shp")



# write_fig(spplot(sp_dGS), "a.pdf")
r <- raster(files[1])
# sp1 = as_SpatialPixelsDataFrame(r)
sp = as_SpatialGridDataFrame(r)# %>% as_SpatialPixelsDataFrame()
# sp2 = sp[inds_good, ]
{
    sp2 = sp
    sp2@data <- d_nvalid
    r2 <- brick(sp2)
    writeRaster(r2, "ChinaPheno1km_nvalid.tif")

    sp2@data <- d_doy
    writeRaster(brick(sp2), "ChinaPheno1km_doy.tif")
    # sp2 <- as_SpatialPixelsDataFrame(sp2)
    # write_fig(spplot(sp2), "sp_nvalid.png")
}
