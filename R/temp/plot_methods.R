# source("test/07_whit/main_phenofit_test.R")
# windowsFonts(Times = windowsFont("Times New Roman"),
#              Arial = windowsFont("Arial"))
fontsize = 14

# re-select colors
cols <- colors()[c(74, 134, 50)]
cols <- c("blue", "red", cols[3]) #"#B2B302"
color_valid <- cols[3] #"green4"


lgd_null  <- phenofit:::make_legend(nmax_points = 4, linename = c())
linecolor <- c("blue", "red", "black")
lgd_name.gpp <- c(expression(1^{th}*"iteration of wWHd smoothed EVI"),
               expression(2^{th}*"iteration of wWHd smoothed EVI"), "flux site GPP")
lgd_lab.gpp  <- legendGrob(lgd_name.gpp, nrow = 1,
           gp = grid::gpar(cex = 1.2, lty = 1, lwd = 3, col = linecolor))
lgd_gpp <- arrangeGrob(lgd_null, lgd_lab.gpp)
# lgd_gpp <- phenofit:::make_legend(linename = c("iter1", "iter2", "GPP"), nmax_points = 4, 
#                    linecolor = cols)

lgd_vci <- phenofit:::make_legend(linename = c("iter1", "iter2", "VCI"), nmax_points = 4,
                   linecolor = cols)

stat_fun <- function(Y_obs, Y_sim){
    R      <- NA_real_
    pvalue <- NA_real_
    tryCatch({
        cor.obj <- cor.test(Y_obs, Y_sim, use = "complete.obs")
        R       <- cor.obj$estimate[[1]] # statistic
        pvalue  <- cor.obj$p.value
    }, error = function(e){
        message(e$message)
    })
    # c(R = R, pvalue = pvalue)[1]
    # pvalue
    R^2
}

table_count <- function(x, levels){
    t <- table(x)
    I <- match(names(t), levels)

    res <- rep(NA, length(levels))
    res[I] <- t
    set_names(res, levels)
}

# agreement index
agree_index <- function(Y_obs, Y_sim){
    I <- which(!(is.na(Y_sim) | is.na(Y_obs))) # | is.na(w)))
    # n_obs <- length(Y_obs)
    n_sim <- length(I)

    Y_sim <- Y_sim[I]
    Y_obs <- Y_obs[I]

    u <- mean(Y_sim)

    d = 100 - sum((Y_sim - Y_obs)^2) / sum( (abs(Y_sim - u) + abs(Y_obs - u))^2 )*100
    d # agreement index
}

box_qtl <- function(x){
    x <- stats::na.omit(x)
    quantile(x, c(0.1, 0.9)) %>% set_names(c("ymin", "ymax"))
}

# boxplot for over all correlation and agreement index
boxplot <- function(p, width = 0.95){
    # width  <- 0.95
    width2 <- width - 0.15
    dodge <- position_dodge(width = width)

    p + stat_summary(fun.data = box_qtl,
                     position = dodge,
                     geom = "errorbar", width = width2) +
        geom_boxplot2(coef = 0,
                  width = width2,
                 lwd = 0.3,
                 notch = F, outlier.shape = NA, position=dodge) +
    theme_light(base_size = fontsize, base_family = "Arial") +
    theme(legend.position = c(1-0.01, 0.01), legend.justification = c(1, 0),
          panel.grid.major = element_line(linetype = 2),
          panel.grid.minor = element_blank(),
          legend.title=element_blank(),
          axis.text = element_text(color = "black"))
}


# geom_point(aes(fill = meth), pch = 21,
#            position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9))
# ggplot(info_df, aes(meth, R), position = "dodge") +
#     stat_summary(fun.data = box_qtl,
#                  position = position_dodge(width = 0.9),
#                  geom = "errorbar", width = 0.9) +
#     geom_boxplot2(coef = 0, width = 0.9, notch = F, outlier.shape = NA)

# get dominant method occurred times
#' @export
stat_dominant <- function(){
    info <- dcast(info_r, site+lat+lon+IGBPname~meth, value.var = "R")
    cols_del <- c(1:4, 9) # 9:whit_R
    methods <- colnames(info)[-cols_del]

    mat <- info[, -cols_del, with = F] %>% as.matrix()
    I <- rowSums(is.na(mat)) == 0

    best <- mat[I, ] %>% {
        data.table(min = methods[apply(., 1, which.min)],
                   max = methods[apply(., 1, which.max)])
    } %>% cbind(info[I, 1:4], .)
    best
    info %<>% cbind(best)

    a_min <- ddply(info, .(IGBP), function(d) table_count(d$min, methods))
    a_max <- ddply(info, .(IGBP), function(d) table_count(d$max, methods))

    t_min <- table(best$min)
    t_max <- table(best$max)
    listk(info, a_min, a_max, t_min, t_max)
    # writelist_ToXlsx(listk(c_min, c_max), "gee_info_count.xlsx")
}

#' @export
get_range <- function(d, alpha = c(0, 1)){
    if (length(alpha) == 1) alpha %<>% rep(2)
    res <- d[, .(min = quantile(value, alpha[1], na.rm = T),
                 max = quantile(value, alpha[2], na.rm = T))]
    unlist(res)
}

#' plot_methods
#'
#' plot curve fitting series and validation data (i.e.GPP or VCI ) to check the
#' smoothing performance
#'
#' @param df_trim A data.table of two iters curve fitting result.
#' For MODIS product, the current year last value maybe same as
#' the coming year first value. Need to remove the duplicated data. Besides,
#' We don't constrain the equal length of different curve fitting series as
#' performance index part.
#' @param st A dataframe of station information, ID, site, IGBPname, lat, lon.
#' @param methods one of 'AG', 'BECK', 'ELMORE', 'ZHANG', 'whit_R' and 'whit_gee'.
#'
#' @examples
#' \dontrun{
#' plot_whit(sitename, df_trim, st, prefix_fig = "whit")
#' }
#' @export
plot_methods <- function(sitename, df_trim, st, prefix_fig = "whit", methods, 
    show.legend = T, base_size = 16)
{
    ## figure title and filename
    sp    <- st[site == sitename, ] # station point
    # titlestr <- with(sp, sprintf('[%03d,%s] %s, lat = %5.2f, lon = %6.2f',
    #                                  ID, site, IGBPname, lat, lon))
    titlestr <- sp$titlestr
    if ( length(titlestr) == 0 || is.na(titlestr) ){
        titlestr <- with(sp, sprintf('[%03d,%s] %s, lat = %5.2f, lon = %6.2f',
                                     ID, site, IGBPname, lat, lon))
    }
    file_pdf <- sprintf('Figure/%s_[%03d]_%s.pdf', prefix_fig, sp$ID[1], sp$site[1])

    ##
    x <- df_trim[site == sitename , ]
    #if (all(is.na(x$whit_gee))) return()
    # d <- melt(x, measure.vars = methods, variable.name = "meth")
    d <- melt(x,
              measure.vars = c(contain(x, "^y$|GPP_NT|GPP_DT|vci|gcc"), methods),
              variable.name = "meth")
    pdat  <- d[meth %in% methods]

    ## scale validation variable (e.g. GPP or VCI)
    # The range of validation variable and \code{whit_gee} should be equal.
    var_wh <- "whit_fluxcam_wWH" #"wWH" #,"whit_fluxcam_wWH"
    d_valid_scale <- d[meth %in% c(var_wh, "GPP_DT", "vci") & iters == "iter2"] %>%
        dcast(site+date~meth, value.var = "value") %>% na.omit() %>%
        melt(id.vars = c("site", "date"), variable.name = "meth")
    lim_fit   <- get_range(d_valid_scale[ grep(var_wh, meth)], alpha = c(0.01, 0.99))
    lim_valid <- get_range(d_valid_scale[ grep("GPP|vci", meth)], alpha = c(0.01, 0.99))
    lim_raw   <- get_range(d[ grep("y", meth)]) # ylim

    # r <- phenofit:::.normalize(y_fit, y_raw)
    # lim_valid_adj <- lm(lim_valid~r) %>% predict(data.frame(r = c(0, 1)))
    coef <- lm(lim_valid~lim_fit) %>% coef() # coefs used to rescale validation data

    ## only keep iter2
    x <- df_trim[site == sitename & iters == "iter2", ]
    if ("vci" %in% colnames(x)){
        d_valid <- x[, .(valid = (vci - coef[1])/coef[2]), .(site, date, t)]
        ylab_r  <- "VCI"
        lgd     <- lgd_vci
    }else{
        d_valid <- x[, .(valid = (GPP_DT - coef[1])/coef[2]), .(site, date, t)]
        ylab_r  <- expression("GPP ( gC "*mm^-1*d^-1*" )")
        lgd     <- lgd_gpp
    }

    d_raw   <- x[, .(date, y, SummaryQA)]

    ## ggplot, not only whit_gee, I also need to know all curve fitting methods
    #  performance
    # y|GPP|vci|gcc|

    IsSingle <- length(methods) == 1
    if (IsSingle){
        d_lab <- data.frame(meth = methods, lab = titlestr)
    }else{
        d_lab <- data.frame(meth = methods,
                         lab = sprintf("(%s) %-8s", letters[1:length(methods)], methods))
    }
    lwd <- 0.65
    # color_valid <- "green4" # set as global variable
    p1 <- ggplot(pdat, aes(date, value)) +
        geom_vline(xintercept = ymd(0101 + (2001:2017)*1e4),
                   color = "white", linetype = 3, size = 0.4) +
        geom_point(data = d_raw, aes(date, y, shape = SummaryQA, color = SummaryQA), size = 1.4) +
        geom_line(data = pdat[iters == "iter1"], color = "blue", size = lwd) +
        geom_line(data = pdat[iters == "iter2"], color = "red", size = lwd) +
        geom_line(data = d_valid, aes(date, valid), size = lwd, color = color_valid) +
        labs(y = "EVI") + 
        theme_grey(base_size = base_size)+
        theme(legend.position = "none",
              # axis.text.y.left = element_text(color = cols[2]), #iter2 color
              # axis.title.y.left = element_text(color = cols[2]),
              axis.text.y.right = element_text(color = color_valid),
              axis.title.y.right = element_text(color = color_valid),
              plot.margin = margin(t = 4, r = 2, b = 0, l = 2, unit = "pt"),
              legend.margin = margin(),
              strip.text = element_blank(),
              panel.grid.minor = element_blank(),
              # panel.grid.major.x = element_blank(),
              # panel.grid.major.y = element_line(size = 0.2),
              panel.grid = element_line(size = 0.4) #linetype = 2
              # axis.ticks.y.right = element_text(color = "blue"),
              ) +
        scale_y_continuous(lim = lim_raw,
                           sec.axis = sec_axis(~.*coef[2]+coef[1], name = ylab_r)) +
        scale_x_date(breaks = ymd(0101 + seq(2001, 2017, 4)*1e4)) +
        facet_wrap(~meth, ncol = 1) +
        scale_color_manual(values = qc_colors, drop = F) +
        scale_shape_manual(values = qc_shapes, drop = F) +
        geom_text(data = d_lab, aes(label = lab), fontface = "bold",
            x = -Inf, y =Inf, vjust = 1.5, hjust = -0.08, size = 4)

    if (!IsSingle) p1 <- p1 + ggtitle(titlestr)

    if (show.legend)
        p1 <- gridExtra::arrangeGrob(p1, lgd_short, nrow = 2,
            heights = c(20, 1), padding = unit(0.5, "line")) #return,

    if (!IsSingle) write_fig(p1, file_pdf, width = 11, height = 7, show = F)
    p1
}
