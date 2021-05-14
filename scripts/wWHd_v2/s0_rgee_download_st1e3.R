library(rgee2)
library(rgee)

ee_Initialize(email = "cuijian426@gmail.com")
# data-raw/whit_lambda/shp/st_1e3_mask.shp, mv to
st <- read_sf("shp/st_1e3_mask.shp") %>% select(c(1, 3))

lst_inds = chunk(1:nrow(st), 16)
for(i in 1:length(lst_inds)) {
    runningId(i)
    ind = lst_inds[[i]]
    sp = st[ind, ]
    prefix = sprintf("st1e3_2000-2020_ith%02d_", i)
    # ee_extract2(ee$ImageCollection$Dataset$MODIS_006_MOD15A2H,
    #             sp, via = "drive", lazy = TRUE,
    #             prefix = prefix)
    ee_extract2(ee$ImageCollection$Dataset$MODIS_006_MCD15A3H,
                sp, via = "drive", lazy = TRUE,
                prefix = prefix)
}

# ee_Initialize(email = "")

