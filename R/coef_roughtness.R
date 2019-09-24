# source('test/07_whit/whit_eval/main_whit_eval.R')
# library(rTIMESAT)
# library(matrixStats)
# source('G:/Github/phenology/phenology/phenofit/test/load_pkgs.R')
# dir_whiteval <- "G:/Github/phenology/phenology/phenofit/data_test/whit_eval"

ngrid <- 1294700
I_all <- 1:1294700
percs <- c(0.1, 0.3, 0.5) %>% set_names(paste0(.*100, "%"))

rm_rownames <- . %>% set_rownames(NULL) %>% data.table::data.table()
set_names2  <- . %>% set_names(names(percs))


#' @export
plot_fitting <- function(d, fit){
    par(mar = c(3.5, 3, 1, 1), mgp = c(1.2, 0.6, 0))
    t <- d$t
    if (is.null(t)) t <- 1:length(d$y)

    plot_input(d)
    iters <- length(fit$zs)

    colors <- c("blue", "red", "green")
    if (iters < 3) colors <- c("blue", "red")

    for (i in 1:iters){
        lines(t, fit$zs[[i]], col = colors[i], lwd = 2)
    }
}
