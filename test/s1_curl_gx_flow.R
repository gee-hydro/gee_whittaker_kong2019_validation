library(magrittr)
library(httr)
library(xml2)
library(rvest)
library(purrr)
library(tidyverse)
library(lubridate)
library(data.table)
library(plyr)

header <- "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/71.0.3578.98 Safari/537.36"
set_config(c(
    # verbose(),
    timeout(60),
    add_headers(
        Connection =  "keep-alive",
        `Accept-Encoding` = "gzip, deflate, br",
        `Accept-Language` = "zh-CN,zh;q=0.8,en;q=0.6",
        # Host = "modis.ornl.gov",
        # Origin = "https://modis.ornl.gov",
        # Referer = url,
        # `Upgrade-Insecure-Requests` = 1,
        `User-Agent` = header
    )
))
# handle_reset("https://modis.ornl.gov") #quite important
# Sys.setlocale("LC_TIME", "english")#

xml_check <- function(x){
    if(class(x)[1] %in% c("xml_document", "xml_node")) x else read_html(x)
}

html_inputs <- function(p, xpath = "//form/input"){
    xml_check(p) %>% xml_find_all(xpath) %>%
    {setNames(as.list(xml_attr(., "value")), xml_attr(., "name"))}
}


url_root <- "http://slt.gxzf.gov.cn:9000/page/index.html?act=3"
p <- GET(url_root) %>% content()

url <- "http://222.216.6.171:3379/arcgis/rest/services/gxsh/mapCommon/MapServer/12/query?where=1%3D1&text=&objectIds=&time=&geometry=&geometryType=esriGeometryEnvelope&inSR=&spatialRel=esriSpatialRelIntersects&relationParam=&outFields=*&returnGeometry=true&maxAllowableOffset=&geometryPrecision=&outSR=&returnIdsOnly=false&returnCountOnly=false&orderByFields=&groupByFieldsForStatistics=&outStatistics=&returnZ=false&returnM=false&gdbVersion=&returnDistinctValues=false&f=pjson"

params <- getFormParams(url) %>% as.list()

p <- POST("http://222.216.6.171:3379/arcgis/rest/services/gxsh/mapCommon/MapServer/12/query",
          body = params, encode = "form") %>% content()

p1 <- GET("http://slt.gxzf.gov.cn:9000/gxsl/japi/api/sl323/realtime/river/") %>% content()

p1 <- GET("http://slt.gxzf.gov.cn:9000/gxsl/japi/api/sl323/realtime/river?time=2019-03-05") %>% content()

write_xml(p, "a.html")

p1 <- GET("http://slt.gxzf.gov.cn:9000/gxsl/japi/api/sl323/realtime/river?time=2018-03-06") %>% content()
p1$result[[1]] %>% as.data.frame()

x1 <- GET("http://slt.gxzf.gov.cn:9000/gxsl/japi/api/sl323/realtime/river/gcx/80805600") %>% content()
a <- x1$result %>% purrr::transpose() %>% as.data.table()
print(a)
