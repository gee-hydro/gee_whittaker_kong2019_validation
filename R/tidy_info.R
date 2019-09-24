
fix_name <- function(x) {
    names(x) %<>% str_extract(".*(?=_)")
    # names(x)
    x
}

#' @export
tidy_info <- function(file){
    x <- readRDS(file)
    a <- llply(x, fix_name) %>% purrr::transpose() %>%
        llply(function(l) {
            info <- map(l, "info") %>% do.call(rbind, .)
            rough <- map(l, "rough") %>% do.call(rbind, .)
            list(info = info, rough = rough)
        })
    names <- names(a)

    for (i in 1:length(a)){
        name <- names[i]
        d_info  <- a[[i]]$info
        d_rough <- a[[i]]$rough

        if (name != "phenofit"){
            d_info$meth <- name
            d_rough$meth <- name
        }
        d_info  %<>% reorder_name(c("site", "meth"))
        d_rough %<>% reorder_name(c("site", "meth"))

        a[[i]] <- merge(d_info, d_rough) %>%
            reorder_name(c("site", "meth", "type", "iters"))#list(info = d_info, rough = d_rough)
    }
    a %<>% do.call(rbind, .) #transpose() %>% map(~
    # a <- fix_name(a)
    # d <- .[NSE > 0, ] %>%
    #     melt(id.vars = c("site", "meth", "type"), variable.name = "index") %>%
    #     .[index %in% c("R2", "NSE", "RMSE")]

    # d$meth %<>% factor(methods)
    # a$Rg %<>% unlist()
    return(a)
}

