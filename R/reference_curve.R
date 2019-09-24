#' get reference of site-doy
#' @export
get_reference <- function(y, is_good, prob = 0.85){
    ngood <- sum(is_good, na.rm = T)
    if (ngood >= 4){
        ref <- median(y[is_good], na.rm = TRUE)
    } else {
        ref <- quantile(y, prob, na.rm = T)
    }
}

#' plot reference EVI curve of one site.
#' @export
plot_ref <- function(){
    qc_colors <- phenofit:::qc_colors
    qc_levels <- phenofit:::qc_levels
    qc_shapes <- phenofit:::qc_shapes

    font.size <- 16
    p <- ggplot(d, aes_string("doy", "y")) +
            geom_point(size = 2, alpha = 0.75,
                aes_string(shape="QC_flag", color = "QC_flag", fill = "QC_flag")) +
            scale_color_manual(values = c(qc_colors, "iter1" = "blue", "iter2" = "red"), drop = F) +
            scale_fill_manual(values = qc_colors, drop = F) +
            scale_shape_manual(values = qc_shapes, drop = F) +
            # scale_x_date(date_labels = "%Y/%m", breaks = seq(xlim_date[1], xlim_date[2], 'year')) +
            theme_gray(base_size = font.size) +
                theme(legend.position="none",
                    axis.title = element_text(size = font.size),
                    axis.text = element_text(size = font.size - 2)
                    # axis.text.x = element_text(angle = 10, hjust = 1, vjust = 1)
                ) +
            labs(x = 'DOY', y = 'EVI') +
        geom_line(data = d_ref, aes(doy, ref), color = "black", size = 1)+
        # geom_line(data = d_ref, aes(doy, ref9)) +
        theme(plot.margin = margin(2, 2, -2, 2, unit = "pt"))

            # geom_vline(xintercept = ymd(0101 + seq(xlim_year[1], xlim_year[2], 1)*1e4),
                       # color = "grey", linetype = 1, size = 0.4)
    # print(p)
    lgd <- phenofit:::make_legend(linename = "Reference",
        linecolor = "black", 4, cex = 1.2)
    g <- arrangeGrob(p, lgd, nrow = 2, heights = c(8, 1),
                padding = unit(0, "line"))
    # grid.newpage()
    # grid.draw(g)

    write_fig(g, "Fig2_ref.pdf", 7, 4)
    # write_fig(g, "Fig2_ref.tif", 7, 4)
}

