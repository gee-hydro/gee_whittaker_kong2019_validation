library(raster)
library(purrr)

# files = dir("n:/MODIS/Terra_LAI_nc/", "_0_1", full.names = "TRUE")
# ndate <- map_int(files, ~brick(.x) %>% {dim(.)[3]})

d = read.table("n:/MODIS/MODIS.jl/whittaker2/examples/LAI_dates.csv", header = TRUE) %>% data.table()
d$date %<>% ymd()

dates_all = get_modis_date(ymd("2000-02-26"), ymd("2019-12-31"))
# seq(1, 365, 16)
## Three date missing
c("2001-06-10", "2001-07-04", "2016-02-10") %>% match(format(d$date))
# [1:60, 60:61, 61:733, 733:end]
# "2001-06-18" "2001-06-26" "2016-02-18"
