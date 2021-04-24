# source("test/main_pkgs.R")
suppressMessages({
    library(tidyverse)
    library(data.table)
    library(plyr)
    library(magrittr)
    library(maptools)

    library(glue)
    library(lubridate)

    library(foreach)
    library(iterators)

    library(Ipaper)
    # library(phenofit)
    # library(whittaker)
})

# install_github('kongdd/ggplot2')
# install_github('kongdd/plyr')
# source("test/load_pkgs.R")
# source('G:/Github/phenology/phenology/phenofit/test/load_pkgs.R')
# source("R/main_phenofit_test.R")
# source("test/stable/ggplot/geom_boxplot_no_outlier.R")
# source('R/plot_phenofit.R')
dir_whiteval <- "data_test/whit_eval"

file_noise_random   <- sprintf("%s/noise3_random.rda", dir_whiteval)
file_noise_keypoint <- sprintf("%s/noise3_keypoint.rda", dir_whiteval) # and real gap

{
    theme_kong <- theme_grey(base_size = 14) +
        theme(panel.grid.minor = element_blank(),
            panel.grid.major = element_line(size = 0.3))
    # theme(legend.position = "none", panel.grid = element_blank()))
    theme_set(theme_kong)
}

# ------------------------------------------------------------------------------
fix_whit <- function(){
    I_missing <- c(1208502, 1212102, 1217502, 1219302, 1225502, 1247823, 1251202) # perc_50

    d_new <- matrix(NA, nrow = ngrid, ncol = 23)
    d_new[setdiff(I_all, I_missing), ] <- r_whit[[3]] %>% as.matrix()
    r_whit[[3]] <- data.table(d_new)
    names(r_whit) <- names(percs)
    fit_wWHd <- r_whit
}

# tidy gof of parameter sensitivity result
tidy_gof <- function(res){
    df.RMSE <- map(res, "RMSE") %>% as.data.table() %>% cbind(I = 1:nrow(.), .) %>%
        melt("I", variable.name = "var")

    df.Roughness <- map(res, "Roughness") %>% as.data.table() %>% cbind(I = 1:nrow(.), .) %>%
        melt("I", variable.name = "var")

    df  <- list(RMSE = df.RMSE, Roughness = df.Roughness) %>% melt_list("index")
    df2 <- df[, .(value = median(value)),.(I2 = ceiling(I/100), var, index)]
    df2
}
