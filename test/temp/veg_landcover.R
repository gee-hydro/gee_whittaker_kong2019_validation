library(rgdal)

# rgdal::readGDAL
d <- read_xlsx("N:/DATA/China/植被类型/veg/植被类型代码表.xlsx")

# 1-km resolution
d_multi <- d[grep("二|三|两", `植被亚类`)]
d_single <- d[-grep("二|三|两", `植被亚类`)]


# rm: (30)温带一年生草本荒漠
names_single <- table(d_single$植被亚类) %>% {.[grep("年", names(.))]} %>% .[-1] %>% names()
names_single_new <- rep("一年一熟", length(names_single))

names_multi <- d_multi$植被亚类 %>% table() %>% names()
names_multi_new <- c("两年三熟或一年两熟", "一年两熟", "一年两熟", "一年两熟或三熟水旱轮作", "一年三熟")

d$植被亚类 %<>% factor(c(names_single, names_multi),
                   c(names_single_new, names_multi_new))
r <- rgdal::readGDAL("N:/DATA/China/植被类型/veg/veg-100w")
r$vegcode <- d$植被亚类[match(r$band1, d$Value)]

write_fig(spplot(r[2]), "crop_land.jpg", 10, 6)

