#' check the performance of Whittaer, and compared with other methods
#' 
#' @export
check_Whittaker <- function(sites, df_trim, file, methods){
    # methods <- c('WH_p2',  "WH_p15",'wWH_p2', "wWH_p15",'wWH') #'WH_p15', 'wWH_p15', 'wWH_v13'
    if (missing(methods)) methods <- c("AG", "ZHANG", "wHANTS", "wSG", "whit_fluxcam_wWH") #"wWH2", "wWH",

    # methods <- "wWH_v13"
    nmeth   <- length(methods)
    show.legend <- ifelse(nmeth == 1, FALSE, TRUE)
    # Cairo::CairoPDF("gee_whit_flux166_all_meth_v2.pdf", 9, nmeth*1.8)

    ps <- list()
    for (i in seq_along(sites)){
        runningId(i)
        sitename <- sites[i]
        # for single method, ggplot obj return
        p <- plot_methods(sitename, df_trim, st, prefix_fig = paste0("whit_", prefix),
                          methods = methods, show.legend = show.legend)
        # print(p)
        ps[[i]] <- p
        # for all methods, grob obj return
        # p <- plot_methods(sitename, df_trim, st, prefix_fig = paste0("whit_", prefix), methods)
    }
    # dev.off()
    ylab_r  <- expression("GPP ( gC "*mm^-1*d^-1*" )")
    # ylab_r  <- "VCI"
    # if (nmeth > 1){
        write_fig2ps(ps, lgd_gpp, ylab_r, file, width = 10, nrow=5)
    # }
}
