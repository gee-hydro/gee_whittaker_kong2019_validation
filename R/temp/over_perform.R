methods <- c('AG', 'BECK', 'ELMORE', 'ZHANG', 'whit_R', 'whit_gee')[-5]

#' Overall performance of different curve fitting methods
#' 
#' @param st data.frame, station info
#' @param var_VI variable name of vegetation time series, can be one of "gpp" 
#' and "vci".
#' 
#' @examples
#' \dontrun{
#' over_perform(df, formula, prefix)
#' }
#' @export
over_perform <- function(d, st, var_VI, formula, prefix, xlab, ylim2 = NULL, IGBP.all = F, outfile){
    # only period when all curve fitting methods have result is kept.
    # browser()

    df_trim <- dcast(d, formula, value.var = "value", fun.aggregate = mean) %>% na.omit()
    df_trim[, raw := y]

    I_var   <- 1:8

    id_vars <- colnames(df_trim)[I_var]
    methods <- colnames(df_trim)[-I_var] %>% sort() %>% .[c(4, 1:3, 12, 5:11)] %>% rm_empty()
    methods %<>% setdiff(c("whit_R", "wWH_p2", "wWH_p15", "WH_p2","WH_p15")) #,"AG", "BECK", "ELMORE", "ZHANG"
    methods <- c("raw", "AG", "ZHANG", "wHANTS", "wSG", "wWH", "wWH2","whit_fluxcam_wWH")
    nmeth   <- length(methods)

    df_trim <- melt(df_trim, id.vars = id_vars, measure.vars = methods, variable.name = "meth")

    # new_levels <- c("raw ", "AG ", "Beck ", "Elmore ", "Zhang ", "WH")
    # methods2 <- methods %>% map_chr(~paste0(.x, " "))
    cols  <- scales::hue_pal()(nmeth-1) %>% c("black", .) %>% set_names(methods)
    cols2 <- cols; cols2[1] <- "transparent"

    # df_trim$meth %<>% mapvalues(methods, methods2)

    # visualization,
    info_ai <- df_trim[SummaryQA == "good", # & meth != "raw "
                       .(ai = agree_index(y, value)), .(site, meth)] %>% merge(st)
    if (var_VI == "gpp") {
        info_r  <- df_trim[, .(R = stat_fun(value, GPP_DT)), .(site, meth)] %>% merge(st)
    } else {
        info_r  <- df_trim[, .(R = stat_fun(value, vci)), .(site, meth)] %>% merge(st)
    }

    # add a column 'all', independent of IGBP
    add_IGBPall <- . %>% {
        temp <- .; temp$IGBPname <- "all"
        rbind(., temp)
    }
    if (IGBP.all){
        info_ai %<>% add_IGBPall()
        info_r  %<>% add_IGBPall()
    }

    th <- theme(panel.grid.major.x = element_blank(),
                panel.grid.major.y = element_line(size = 0.2))
    grid_x <- geom_vline(data = xlab, xintercept = (1:(nrow(xlab) - 1)) + 0.5,
        linetype = 2, size = 0.2, color = "grey70")

    # s <- ggplot(x, aes(x = meth, y = value, fill = kind)) +
    #     geom_bar(width = 1, stat = "identity",
    #              position = position_stack(reverse = TRUE)) +
    #     geom_text(aes(label = value))
        # coord_flip()
    # p + geom_subview(s)
        # ggrepel::geom_text_repel(data = d[raw - wWH > 0.05], # | whit_fluxcam_wWH < 0
        #                          aes(label = site))

    # 1. show correlation, , alpha = 0.6, fill = meth
    p1 <- { ggplot(info_r, aes(IGBPname, R, color = meth), position = "dodge") +
                scale_colour_manual(values = cols,
                                    guide = guide_legend(direction = "horizontal", nrow = 1, keywidth = 1)) +
                # scale_fill_manual(values = cols) +
                labs(x = "IGBP", y = "Correlation (r)")
          } %>% boxplot()
    p1 <- p1 + grid_x +
        # geom_text(data = info_r[1, ], aes(fontface = "bold"), x = -Inf, y = -Inf,
        #           label = c("(a)"), hjust = -1, vjust = -2.6,
        #           show.legend = F, size = 4) +
        th +
        theme(axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              plot.margin = margin(3, 3, 2, 3))

    p2 <- {ggplot(info_ai, aes(IGBPname, ai, colour = meth), position = "dodge") +
        scale_colour_manual(values = cols2) } %>%
        boxplot() %>% `+`(labs(x = "IGBP", y = "Agreement Index (AI)"))
    p2 <- p2 + grid_x + th +
        theme(legend.position = 'none',
                    plot.margin = margin(3, 3, 2, 3)) +
        scale_x_discrete(breaks = xlab$IGBPname, labels = xlab$label)

    if (!is.null(ylim2)) p2 <- p2 + coord_cartesian(ylim = ylim2)
    p <- gridExtra::arrangeGrob(p1, p2, nrow = 2, heights = c(0.9, 1), padding = unit(1, "line"))
    # grid.draw(p)
    if (missing(outfile)) outfile <- sprintf("valid_%s.pdf", prefix)
    write_fig(p, outfile, 10, 8, show = T)

    return(info_r)
    # save_pdf(sprintf("valid_%s_AI.pdf", prefix), 12, 5, p2)
}
