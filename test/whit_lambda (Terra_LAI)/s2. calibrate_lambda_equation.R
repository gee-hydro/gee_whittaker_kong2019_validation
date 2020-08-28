# !/usr/bin/Rscript
source('test/main_pkgs.R')
library(sf)
library(broom)
# setwd("/media/kong/Various Data/Research/GEE_repos/gee_whittaker/")
# print(getwd())
{
    ## 1. station info
    st   <- sf::read_sf("data-raw/whit_lambda/shp/st_1e3_mask.shp")
    coor <- st_geometry(st) %>% do.call(rbind, .) %>% data.table() %>%  set_colnames(c("lon", "lat"))
    st <- as.data.table(st)[, c(1, 3)] %>% cbind(coor)
    st$site %<>% as.character()
    st$IGBP <- factor(st$IGBPcode, levels = IGBPcodes_006, labels = IGBPnames_006)
    # colnames(st)[2] <- "IGBP"
}

# args <- commandArgs(TRUE)
# print(args)
files = dir("OUTPUT/whit_lambda/", full.names = TRUE, recursive = TRUE)

type = dirname(files) %>% basename()
chunksize = basename(files) %>% str_extract("(?<=chunk)\\d{1,}_.*(?=\\.)")
prefix = paste0(type, "_", chunksize)

#
lst <- llply(files, readRDS, .progress = "text")
names(lst) <- prefix
# l <- lst$`normalized-01`
gof <- map(lst, ~map(.x, "gof") %>% do.call(rbind, .)) %>% melt_list("id_str")
gof[, `:=`(type = str_extract(id_str, ".*(?=_\\d)"),
           chunksize = str_extract(id_str, "(?<=_).*(?=_)"),
           is_extend = str_extract(id_str, "(?<=\\d{2}_).*"))]

ggplot(gof, aes(chunksize, R2, color = type)) +
    geom_boxplot2()

# InitCluster(8, kill = FALSE)
lst_full = foreach(l = lst, i = icount()) %do% {
    runningId(i)
    d = map(l, ~ .x$coef) %>%
        melt_list("site") %>%
        data.table()
}
#  melt_list("id_str")
# st <- read_sf("data-raw/whit_lambda/shp/st_1e3_mask.shp") %>% data.table() %>% .[, -c(2, 4)]

# 结果显示前后扩展一年非常有必要 -----------------------------------------------
# 无扩展时R2 = 0.0905, 扩展时R2 = 0.167
nyears = str_extract(prefix, "(?<=_).*(?=_)") %>% as.numeric()
ngrp = length(nyears)/2
# 不考虑植被类型发生突变

df = foreach(d = lst_full,
             i = icount(),
             nyear = nyears,
             is_extend = rep(c(1, 0), ngrp),
             kind = rep(c("normalized", "original"), each = 12)) %do% {
    # d$site %<>% as.numeric()
    cbind(d, nyear, is_extend, kind)
    # lm(lambda ~ mean + sd + kurtosis +skewness, d) %>% glance()
} %>% do.call(rbind, .) %>% merge(st)


## 不同植被类型的模型表现-------------------------------------------------------
{
    par(mfrow = c(1, 2))
    IGBPs_bad = c("EBF", "CNV", "URB")#[1]
    d = df[kind == c("normalized", "original")[2] & is_extend == 0 & nyear == 4]
    d[, lambda2 := mark_outlier(lambda), .(site, nyear, is_extend, kind)]
    # d[lambda < 1, lambda := 1]
    l <- lm(log10(lambda2) ~ mean + sd + skewness + kurtosis, d[!(IGBP %in% IGBPs_bad), ]) %T>% {print(glance(.))}
    plot_mls(l)
    l$coefficients
    # lm(log10(lambda) ~ sd  + kurtosis + skewness, d)
    # lm(log10(lambda) ~ sd + mean + kurtosis + skewness, d)
    ## 检查采用该模型对c("EBF", "CNV", "URB")的预测情况
    d$ypred = predict(l, newdata = d)
    with(d[IGBP %in% IGBPs_bad[1]], xyplot(log10(lambda), ypred))
}

# # A tibble: 1 x 11
# r.squared adj.r.squared sigma statistic p.value    df   logLik     AIC     BIC deviance df.residual
# <dbl>         <dbl> <dbl>     <dbl>   <dbl> <int>    <dbl>   <dbl>   <dbl>    <dbl>       <int>
#     1     0.500         0.500 0.596    27944.       0     5 -100956. 201924. 201982.   39796.      111944

# 以此编写GEE版本的函数
# >     l$coefficients
# (Intercept)        mean          sd    skewness    kurtosis
# 1.77365505  0.43062881 -0.34192178 -0.30107590  0.03221195

{
    r <- NULL
    write_fig({
        par(mfrow = c(4, 4), mar = c(2, 2, 1, 1), mgp = c(3, 0.6, 0))
        igbp_codes = IGBPcodes_006[2:17]
        r <<- foreach(igbp_code = igbp_codes, igbp_name = names(igbp_codes), i = icount()) %do% {
            d = df[kind == c("normalized", "original")[2] & is_extend == 1 & IGBPcode == igbp_code]
            l <- lm(log10(lambda) ~ mean + sd + kurtosis + skewness + lat, d)

            title = sprintf("(%s) %s", letters[i], igbp_name)
            plot_mls(l, title = title)
            glance(l)
        } %>% melt_list("IGBPname") %>% data.table
    }, "whit_lambda_equation.pdf", 10, 10)
    # l$coefficients
}
r[order(adj.r.squared)]


# m = aov(lambda ~ factor(nyear), df[is_extend == 1, ])
# TukeyHSD(m)

# lm(lambda ~ mean + sd + kurtosis + skewness + nyear, df[is_extend == 0]) %>% glance()
# lm(lambda ~ mean + sd + kurtosis + skewness + nyear, df[is_extend == 1]) %>% glance()

# # 是否考虑年影响不大
# lm(lambda ~ mean + sd + kurtosis + skewness, df[is_extend == 0]) %>% glance()
# lm(lambda ~ mean + sd + kurtosis + skewness, df[is_extend == 1]) %>% glance()

# {
#     l <- lm(lambda ~ mean + sd + kurtosis +skewness, df[kind == c("normalized", "original")[1] & is_extend == 1])
#     info = data.table(ypred = l$fitted.values, yobs = l$model$lambda)

# }

# info = foreach(d = res, i = icount()) %do% {
#     lm(lambda ~ mean + sd + kurtosis +skewness, d) %>% glance()
# } %>% melt_list("id_str")

