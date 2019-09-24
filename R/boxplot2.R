#' @export
boxplot2 <- function(p, width = 0.95, size = 0.7){
    # width  <- 0.95
    width2 <- width - 0.15
    dodge <- position_dodge(width = width)

    p + stat_summary(fun.data = box_qtl,
                     position = dodge, size = size,
                     geom = "errorbar", width = width2) +
        geom_boxplot2(coef = 0,
                      width = width2,
                      lwd = size - 0.2,
                      notch = F, outlier.shape = NA, position=dodge) +
        grid_x +
        geom_text(data = d_lab, aes(x = "ENF",
                                    y = Inf, color = NULL, label = label),
                  vjust = 1.5, hjust = 1.1, fontface = "bold", size =5, show.legend = F)
}
