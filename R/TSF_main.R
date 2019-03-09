# source('test/07_whit/main_TSF.R')

#' @param FUN 
TSF_main <- function(d, nptperyear = 23, job_name,
  FUN = 2, iters = 3, half_win = floor(nptperyear/4), cache = T)
{
    if (missing(job_name)) {
        job_name <- d$site[1]
    }
    nyear <- floor(nrow(d)/nptperyear)
    npt   <- nyear * nptperyear
    d <- d[1:npt, ]
    # pls make sure it's complete year in input

    file_y   <- sprintf("TSM_%s_y.txt", job_name)
    file_w   <- sprintf("TSM_%s_w.txt", job_name)
    file_set <- sprintf("TSM_%s.set", job_name)

    write_input(d$y  , file_y, nptperyear)
    write_input(d$SummaryQA, file_w, nptperyear)

    ## 2. Update options
    options <- list(
       job_name            = job_name,
       file_y              = file_y,             # Data file list/name
       file_w              = file_w,             # Mask file list/name
       nyear_and_nptperear = c(nyear, nptperyear),      # No. years and no. points per year
       ylu                 = c(0, 9999),     # Valid data range (lower upper)
       qc_1                = c(1, 1, 1),     # Quality range 1 and weight
       qc_2                = c(2, 2, 0.5),   # Quality range 2 and weight
       qc_3                = c(3, 4, 0.2),   # Quality range 3 and weight
       A                   = 0.1,            # Amplitude cutoff value
       output_type         = c(1, 1, 0),     # Output files (1/0 1/0 1/0)
       seasonpar           = 1.0,            # Seasonality parameter (0-1)
       iters               = 2,              # No. of envelope iterations (3/2/1)
       FUN                 = 2,              # Fitting method (3/2/1)
       half_win            = 6,              # half Window size for Sav-Gol.
       meth_pheno          = 1,              # Season start / end method (4/3/2/1)
       trs                 = c(0.5, 0.5)     # Season start / end values
    )

    options$job_name <- job_name
    options$FUN <- FUN
    options$half_win <- half_win

    # update setting
    opt <- update_setting(options)
    write_setting(opt, file_set)

    TSF_process(file_set) # call TSF_process.exe

    file_tts <- sprintf("%s_fit.tts", opt$job_name)
    file_tpa <- sprintf("%s_TS.tpa", opt$job_name)

    # note: only suit for ascii
    tidy_tts <- function(d_tts){
        sites <- d_tts$row %>% paste0("v", .)
        npt   <- ncol(d_tts) - 2
        d <- d_tts %>% {.[, 3:ncol(.)]} %>% as.matrix() %>% t() %>% data.frame() %>%
            set_colnames(sites) %>% cbind(t = 1:npt, .) %>% 
            set_rownames(NULL)
        d
    }

    d_tts <- read_tts(file_tts) %>% tidy_tts()
    d_tpa <- read_tpa(file_tpa)

    if (!cache){
        status1 <- file.remove(c(file_tts, file_tpa, file_y, file_w, file_set))
        status2 <- dir(".", "*.ndx", full.names = T) %>% file.remove()
    }
    list(fit = d_tts, pheno = d_tpa)
}

## pdat
#             t      y QC_flag        SG        AG      wWHd
# 1: 2004-01-06 0.3807    good 0.3760034 0.3713598 0.3729616
# 2: 2004-01-22 0.3254    good 0.3730653 0.3615378 0.3625175
# 3: 2004-02-16 0.3801    good 0.3776879 0.3588165 0.3557184
# 4: 2004-03-01 0.3556    good 0.3792459 0.3589701 0.3524135
# 5: 2004-03-19 0.2936    good 0.3830829 0.3617738 0.3541164
show_fitting <- function(pdat, methods = c("SG", "AG", "value"), colors, show.legend = T){
    qc_colors <- phenofit:::qc_colors
    qc_shapes <- phenofit:::qc_shapes

    xlim_year <- pdat$date %>% {c(first(.), last(.))} %>% year()
    xlim_date <- paste0(xlim_year, "-01-01") %>% as.Date()

    nmeth   <- length(methods)

    font.size <- 14
    if (missing(colors)){
        colors <- hue_pal()(nmeth) # c("red", "blue", "black")  
    }

    # self defined legend
    nmax_points <- ifelse(max(pdat$QC_flag %>% as.numeric()) <= 4, 4, 6)    
    lgd <- phenofit:::make_legend(linename = methods,
        linecolor = colors, nmax_points, cex = 1.2)

    p <- ggplot(pdat, aes_string("t", "y")) +
        geom_point(size = 3, alpha = 0.75,
            aes_string(shape="QC_flag", color = "QC_flag", fill = "QC_flag")) +
        scale_color_manual(values = c(qc_colors, "iter1" = "blue", "iter2" = "red"), drop = F) +
        scale_fill_manual(values = qc_colors, drop = F) +
        scale_shape_manual(values = qc_shapes, drop = F) +
        scale_x_date(date_labels = "%Y/%m", breaks = seq(xlim_date[1], xlim_date[2], 'year')) +
        theme_gray(base_size = font.size) +
            theme(legend.position="none",
                axis.title = element_text(size = font.size),
                axis.text = element_text(size = font.size - 2)
                # axis.text.x = element_text(angle = 10, hjust = 1, vjust = 1)
            ) +
        labs(x = 'Time', y = 'Vegetation Index') +
        geom_vline(xintercept = ymd(0101 + seq(xlim_year[1], xlim_year[2], 1)*1e4),
                   color = "grey", linetype = 1, size = 0.4)

    for (i in seq_along(methods)){
        method <- methods[i]    
        p <- p + geom_line(aes_string(y = method), size = 0.8, alpha = 0.7, color = colors[i])
    }
    
    if (show.legend){
        arrangeGrob(p, lgd, nrow = 2, heights = c(10, 1),
                padding = unit(1, "line"))
    } else {
        p
    }
}
