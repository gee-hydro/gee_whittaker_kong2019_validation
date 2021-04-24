# source("test/main_pkgs.R")
library(sf)
library(rgee)
library(rgee2)
ee_Initialize(drive = TRUE)

# InitCluster(1)
sf = read_sf("data-raw/ChinaPheno1km_10y_DoubleSeasonPixels_Maize&Wheat_st1e2.shp")
# sf %>% as.data.table()
sf2 = sf %>% {cbind(Id = 1:nrow(.), .)} %>% .[, "Id"]

# load_all("/mnt/i/Research/rpkgs/rgee2")
## LAI
imgcol = ee$ImageCollection$Dataset$MODIS_006_MCD15A3H$
    filterDate("2010-01-01", "2016-12-31")
df = ee_extract2(imgcol, sf2, prefix = "ChinaPheno_sp1e2_", via = "drive")
## EVI
imgcol = ee$ImageCollection$Dataset$MODIS_006_MOD13A1$
    filterDate("2010-01-01", "2016-12-31")
df = ee_extract2(imgcol, sf2, prefix = "ChinaPheno_sp1e2_", via = "drive")

# sp2 = sf2
drive_csv_clean("data-raw/temp/ChinaPheno_sp1e2_MODIS_006_MCD15A3H.csv", sf2)
drive_csv_clean("data-raw/temp/ChinaPheno_sp1e2_MODIS_006_MOD13A1.csv", sf2)
