get_modis_date <- function(date_begin, date_end, dn = 8) {
    years = seq(year(date_begin), year(date_end))

    dates <- foreach(year = years, i = icount(), .combine = c) %do% {
        str = sprintf("%d%03d", year, seq(1, 365, dn))
        as.Date(str, "%Y%j")
    }
    dates[dates >= date_begin & dates <= date_end]
}
