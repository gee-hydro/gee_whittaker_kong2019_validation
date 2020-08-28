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


## MODIS date functions --------------------------------------------------------

get_modis_date <- function(date_begin, date_end, dn = 8) {
    years <- seq(year(date_begin), year(date_end))

    dates <- foreach(year = years, i = icount(), .combine = c) %do% {
        str <- sprintf("%d%03d", year, seq(1, 365, dn))
        as.Date(str, "%Y%j")
    }
    dates[dates >= date_begin & dates <= date_end]
}

#' fill_missdate
#' fill missing date
#' @export
fill_missdate <- function(dn = 8) {
    years <- 2000:2019
    # dn = 8
    doy <- seq(1, 366, dn)
    nptperyear <- length(doy)
    date <- sprintf("%4d%03d", rep(years, each = nptperyear), doy) %>%
        parse_date_time("%Y%j") %>%
        lubridate::date()
    if (years[1] == 2000) date <- date[date >= "2000-02-26"]
    # date  <- date[1:(length(date)-12)] # for 2018
    date
}
