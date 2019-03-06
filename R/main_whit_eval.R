# source('test/07_whit/whit_eval/main_whit_eval.R')
library(rTIMESAT)
library(matrixStats)

source('G:/Github/phenology/phenology/load_pkgs.R')

# dir_whiteval <- "G:/Github/phenology/phenology/phenofit/data_test/whit_eval"
dir_whiteval <- "F:/whit_eval"

ngrid <- 1294700
I_all <- 1:1294700
percs <- c(0.1, 0.3, 0.5) %>% set_names(paste0(.*100, "%"))


rm_rownames <- . %>% set_rownames(NULL) %>% data.table::data.table()
set_names2  <- . %>% set_names(names(percs))

#' coef_roughness
#' D2
#'
#' @param x matrix, [ngrid, ntime]
coef_roughness <- function(x){
    x <- as.matrix(x)
    rowDiffs(x, differences = 2)^2 %>% rowSums2(na.rm = T)
}

#' get reference of site-doy
get_reference <- function(y, is_good, prob = 0.85){
    ngood <- sum(is_good, na.rm = T)
    if (ngood >= 4){
        ref <- median(y[is_good], na.rm = TRUE)
    } else {
        ref <- quantile(y, prob, na.rm = T)
    }
}


#' plot reference EVI curve of one site.
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


#' Qiang Zhang, 2018, AFM, Eq.2
der_k <- function(y) {
    f1 <- c(0, diff(y))
    f2 <- c(0, 0, diff(y, differences = 2))
    k  <- f2/((1+f1^2)^1.5)
    k
}


#' add_noise
#'
#' QC variable also returned.
#' Same as MOD13A1 SummaryQA, 0: good value, 3: cloud contaminated
add_noise <- function(y, perc = 0.1, seed, I){
    n <- length(y)

    if (!missing(seed)) set.seed(seed)
    if (missing(I) || is.null(I)) {
        I <- sample(1:n, ceiling(n*perc))
    }

    n_noise <- length(I)
    # was random reduced by 5% - 40% with an interval of 5%

    QC <- rep(0, n)
    if (length(I) > 0){
        A_noise <- ceiling(runif(n_noise, 0.05, 0.4)/0.05)*0.05

        y[I] <- y[I] * (1 - A_noise)
        QC[I] <- 3
    }
    list(y = y, QC = QC)
}

#' simu_noise
#'
#' simulate 100 times about \code{add_noise}
#' @param type One of c("random", "real", "maxK", "maxDer")
#' @param ntime Simulate how many times
simu_noise <- function(d, perc = 0.1, trans = TRUE, type, ntime = 100){
    y <- d$ref
    # 生成两个噪点
    if (type == "random") {
        I <- NULL
    } else if (type == "real") {
        I <- which(d$SummaryQA != "good")
    } else if (type == "maxK") {
        k <- der_k(y)
        I <- c(which.max(k), which.min(k))
    } else if (type == "maxDer") {
        f1 <- c(0, diff(y))
        I  <- c(which.max(f1), which.min(f1))
    }

    r <- llply(1:ntime, function(i){ add_noise(y, perc, i, I) })

    if (trans) {
        purrr::transpose(r) %>% map(~do.call(rbind, .))
    } else {
        r
    }
}


#' TSF_process of TIMESAT
#' Three curve fitting function (i.e. SG, AG and DL) in TIMESAT are used.
TSF_main_df <- function(df_sim, perc = 0.1, nyear = 1, nptperyear = 23,
    overwrite = TRUE, wait = TRUE, indir)
{
    dir_root <- sprintf("%s/perc_%d", indir, perc*100)
    if (!dir.exists(dir_root)) dir.create(dir_root, recursive = TRUE)
    setwd(dir_root)

    file_y <- sprintf("TSF_whit_eval_y_noise%d.txt", perc*100)
    file_w <- sprintf("TSF_whit_eval_w_noise%d.txt", perc*100)

    if (!file.exists(file_y) || overwrite) {
        system.time(write_input(df_sim$y , file = file_y, nptperyear = nptperyear))
    }

    if (!file.exists(file_w) || overwrite) {
        system.time(write_input(df_sim$QC, file = file_w, nptperyear = nptperyear))
    }

    # nyear <- floor(ncol(df_sim$y)/nptperyear)
    TSF_meths <- c("SG", "AG", "DL")

    ## FOR TIMESAT
    ## 2. Update options
    options <- list(
       job_name            = "whit_eval",
       file_y              = file_y,             # Data file list/name
       file_w              = file_w,             # Mask file list/name
       nyear_and_nptperear = c(nyear, nptperyear),      # No. years and no. points per year
       ylu                 = c(0, 1),        # Valid data range (lower upper)
       qc_1                = c(0, 0, 1),     # Quality range 1 and weight
       qc_2                = c(1, 1, 0.5),   # Quality range 2 and weight
       qc_3                = c(2, 3, 0.2),   # Quality range 3 and weight
       A                   = 0.05,            # Amplitude cutoff value
       output_type         = c(0, 1, 0),     # Output files (1/0 1/0 1/0), 1: seasonality data; 2: smoothed time-series; 3: original time-series
       seasonpar           = 0.5,            # Seasonality parameter (0-1)
       iters               = 2,              # No. of envelope iterations (3/2/1)
       FUN                 = 2,              # Fitting method (1/2/3): (SG/AG/DL)
       half_win            = 6,              # half Window size for Sav-Gol.
       meth_pheno          = 1,              # Season start / end method (4/3/2/1)
       trs                 = c(0.5, 0.5)     # Season start / end values
    )

    for (i in 1:3){
        options$FUN <- i
        meth <- TSF_meths[options$FUN]

        options$job_name <- sprintf('TSF_whit_eval_noise%d_%s', perc*100, meth)

        file_set <- sprintf("TSF_whit_eval_noise%d_%s.set", perc*100, meth)
        file_tts <- sprintf("%s_fit.tts", options$job_name)
        file_tpa <- sprintf("%s_TS.tpa", options$job_name)

        if (!file.exists(file_tts) || overwrite){
            opt <- update_setting(options)
            write_setting(opt, file_set)
            TSF_process(file_set, 4, wait=wait) # call TSF_process.exe
        }
    }
}

plot_fitting <- function(d, fit){
    par(mar = c(3.5, 3, 1, 1), mgp = c(1.2, 0.6, 0))
    t <- d$t
    if (is.null(t)) t <- 1:length(d$y)

    plot_input(d)
    iters <- length(fit$zs)

    colors <- c("blue", "red", "green")
    if (iters < 3) colors <- c("blue", "red")

    for (i in 1:iters){
        lines(t, fit$zs[[i]], col = colors[i], lwd = 2)
    }
}

smooth_WHIT_raw <- function(d, lambda, iters = 3, IsPlot = FALSE){
    d <- c(d, qc_summary(d$QC))
    l <- check_input(t, d$y, d$w, d$QC_flag, nptperyear = 23, south = F)

    if (missing(lambda) || is.null(lambda)){
        # if missing, init according to v-curve
        lambda <- init_lambda(l$y)
    }

    r <- wWHIT(l$y, l$w, l$ylu, nptperyear, wFUN = wBisquare,
               iters = iters, lambda = lambda)

    if (IsPlot) plot_fitting(d, r)
    GOF(d_ref$ref, last(r$zs))
}

smooth_WHIT <- function(y, w, lambda, iters = 3, IsPlot = FALSE){
    l <- check_input(t, y, w, nptperyear = 23, south = F)

    if (missing(lambda) || is.na(lambda)){
        # if missing, init according to v-curve
        lambda <- init_lambda(l$y)
    }

    r <- wWHIT(l$y, l$w, l$ylu, nptperyear, wFUN = wBisquare,
               iters = iters, lambda = lambda)

    if (IsPlot) plot_fitting(data.frame(y, w), r)
    # GOF(d_ref$ref, dplyr::last(r$zs))
    dplyr::last(r$zs)
}
