
#' fill_missdate
#' fill missing date
#' @export
fill_missdate <- function(dn = 8){
    years <- 2000:2019
    # dn = 8
    doy   <- seq(1, 366, dn)
    nptperyear = length(doy)
    date  <- sprintf("%4d%03d", rep(years, each = nptperyear), doy) %>% parse_date_time("%Y%j") %>% lubridate::date()
    if (years[1] == 2000) date <- date[date >= "2000-02-26"]
    # date  <- date[1:(length(date)-12)] # for 2018
    date
}
