df_lai = fread("data-raw/wWHd2/ChinaPheno_sp1e2_MODIS_006_MCD15A3H.csv")
df_evi = fread("data-raw/wWHd2/ChinaPheno_sp1e2_MODIS_006_MOD13A1.csv")

df_lai[, c("QC_flag", "w") := qc_FparLai(FparExtra_QC)]
library(tidyverse)
d = df_lai[, .(y = mean(Lai, na.rm = TRUE)), .(Id, yday(date))]
p <- ggplot(d, aes(yday, y)) +
    geom_point() +
    geom_line() +
    facet_wrap(~Id)

d = df_evi[, .(y = mean(EVI, na.rm = TRUE)), .(Id, yday(date))]
p <- ggplot(d, aes(yday, y)) +
    geom_point() +
    geom_line() +
    facet_wrap(~Id)

write_fig(p, "a.pdf", 20, 20)
