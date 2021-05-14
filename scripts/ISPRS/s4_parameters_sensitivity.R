source("scripts/main_pkgs.R")
# data("mat_ref")
# data('lst_noise')

## parameters
{
    limits = NULL

    lambdas  <- c(NA, 2, 15) %>% set_names(c("wWHd", "wWH2", "wWH15")) %>% .[1]
    perc_wcs <- seq(0.1, 0.9, 0.05) %>% set_names(., paste0(.*100, "%"))
    ws <- seq(0.1, 0.5, 0.05) %>% set_names(., paste0(.*100, "%"))
    # mat_ref1 <- mat_ref[limits/100, , drop = FALSE]
}

## 1. Examples -------------------------------------------------------------
# only dynamic lambda was tested here.
s1_para = FALSE
if (s1_para) {
    load(file_noise_random)
    # load(file_noise_keypoint)

    InitCluster(17)
    # lst_noise <- lst_noise_keypoint[1]
    lst_noise <- lst_noise_random %>% map(function(x){
        list(y = x$y[seq(10, 12947*100, 10), ],
             QC = x$QC[seq(10, 12947*100, 10), ])
    })
    times <- 10 # noise simulated times

    # 1. perc_wc
    res.perc_wc <- foreach(perc_wc = perc_wcs,
                       .packages = c("whittaker", "magrittr"),
                       i = icount()) %dopar% {
        Ipaper::runningId(i, prefix = "perc_wc")
        r <- fitSimuNoise_WHIT(lst_noise, lambdas, limits = NULL, wmin = 0.2, perc_wc = perc_wc)$wWHd
        d_gof <- purrr::map(r, ~gof_fitting(mat_ref, .x, times = times)) %>%
            Ipaper::melt_list("type")
    }
    killCluster()
    saveRDS(res.perc_wc, "OUTPUT/para_perc_wc.RDS")

    # 2. wmin
    res.wmin <- foreach(wmin = ws,
                       .packages = c("whittaker", "magrittr"),
                       i = icount()) %dopar% {
        Ipaper::runningId(i, prefix = "wmin")
        r <- fitSimuNoise_WHIT(lst_noise, lambdas, limits = NULL, wmin = wmin, perc_wc = 0.4)$wWHd
        d_gof <- purrr::map(r, ~gof_fitting(mat_ref, .x, times = times)) %>%
           Ipaper::melt_list("type")
    }
    saveRDS(res.wmin, "OUTPUT/para_wmin.RDS")
}

## 2. Figures -------------------------------------------------------------
Fig9_para_sensitivity <- function(){
    # 1. prepare data
    res.perc_wc <- readRDS("OUTPUT/para_perc_wc.RDS")
    res.wmin    <- readRDS("OUTPUT/para_wmin.RDS")

    tidy_gof2 <- function(res){
        map(res, ~.x[, .(I, Roughness, RMSE, type)]) %>%
            melt_list("param") %>%
            melt(measure.vars = c("Roughness", "RMSE"))
    }
    df.perc_wc <- tidy_gof2(res.perc_wc)
    df.wmin    <- tidy_gof2(res.wmin)

    level.wmin <- df.wmin$param %>% unique() %>% {as.numeric(gsub("\\%", "", .) )/100} %>% as.character()
    level.perc_wc <- df.perc_wc$param %>% unique() %>% gsub("\\%", "", .)
    level.varname <- c("bold(perc[wc] * ' (%)')", "bold(w[min])")
    level.index   <- c("bold(RMSE)", "bold(roughness)")

    df.wmin[, param := factor(param, labels = level.wmin)]
    df.perc_wc[, param := factor(param, labels = level.perc_wc)]
    df <- list(wmin = df.wmin, perc_wc = df.perc_wc) %>% melt_list("varname")
    df[, varname := factor(varname, labels = level.varname)]
    df[, variable   := factor(variable, c("RMSE", "Roughness"), labels = level.index)]

    d_lab <- expand.grid(variable = level.index, varname = level.varname) %>%
        cbind(lab = sprintf("(%s)", letters[1:4]), .)

    # 2. produce figure
    p <- ggplot(df, aes(param, value, fill = type)) +
        stat_summary(fun.data = box_qtl, geom = "errorbar", width = 0.5, color = "black") +
        geom_boxplot2(show.errorbar = FALSE) +
        facet_grid(variable~varname, scale = "free", labeller = label_parsed) +
        geom_text(data = d_lab, aes(label = lab, fill = NULL), x = -Inf, y = Inf,
                  hjust = -0.2, vjust = 1.5, size = 5.5) +
        labs(y = NULL, x = NULL) +
        theme(
            strip.text.x = element_text(margin = margin(1,0,1,0, "pt")*5, size = 14),
            strip.text.y = element_text(face = "bold"),
            strip.text = element_text(size = 13, face = "bold"))

    # curve -------------------------------------------------------------------
    alpha = 0.25
    df2 <- df[, .(value = median(value),
                  # sd = sd(value),
                  ymin = quantile(value, alpha),
                  ymax = quantile(value, 1 - alpha)),
              .(type, param, variable, varname)]
    df2[, param := as.numeric(as.character(param))]
    # df2 <- df2 %>% mutate(ymin = pmax(0, value - sd), ymax = value + sd)
    name_legend <- "Percentage of random gaps"
    p <- ggplot(df2, aes(param, value, color = type, shape = type)) +
        geom_line(size = 0.5) +
        geom_point(size = 1.7) +
        geom_text(data = d_lab, aes(color = NULL, shape = NULL, label = lab),
                  x = -Inf, y = Inf, hjust = -0.3, vjust = 2.5, size = 5, show.legend = FALSE) +
        # geom_errorbar(aes(ymin = ymin, ymax = ymax)) +
        # stat_summary(fun.data = box_qtl, geom = "errorbar", width = 0.5, color = "black") +
        facet_grid(variable~varname, scale = "free", labeller = label_parsed) +
        theme(legend.position = "bottom",
              legend.margin = margin(-7, 0, -5, 0, "pt"),
              axis.text = element_text(color = "black"),
              strip.text.x = element_text(size = 14)) +
        labs(x = NULL, y = NULL, color = name_legend, shape = name_legend)

    gt = ggplot_gtable(ggplot_build(p))
    gt$widths[7] <- unit(15/17, "null")
    write_fig(gt, "Fig9_parameter_sensitivity2.tif", 10, 5)

    write_fig(gt, "Fig9_parameter_sensitivity2.pdf", 10, 5)
    write_fig(gt, "Fig9_parameter_sensitivity2.emf", 10, 5)

}

## real-noise ------------------------------------------------------------------
# 1. prepare data
{
    res.perc_wc <- readRDS("OUTPUT/real-gaps/para_perc_wc.RDS")
    res.wmin    <- readRDS("OUTPUT/real-gaps/para_wmin.RDS")

    tidy_gof2 <- function(res){
        data.table(
            x = names(res),
            RMSE = res %>% map_dbl(~mean(.$RMSE, na.rm = TRUE)),
            Roughness = res %>% map_dbl(~mean(.$Roughness, na.rm = TRUE))
        )
    }
    df.perc_wc <- tidy_gof2(res.perc_wc)
    df.wmin    <- tidy_gof2(res.wmin)

    level.wmin <- df.wmin$x %>% unique() %>% {as.numeric(gsub("\\%", "", .) )/100} %>% as.character()
    level.perc_wc <- df.perc_wc$x %>% unique() %>% gsub("\\%", "", .)
    level.varname <- c("bold(perc[wc] * ' (%)')", "bold(w[min])")
    level.index   <- c("bold(RMSE)", "bold(Roughness)")

    df.wmin[, x := factor(x, labels = level.wmin)]
    df.perc_wc[, x := factor(x, labels = level.perc_wc)]
    df <- list(wmin = df.wmin, perc_wc = df.perc_wc) %>% melt_list("varname") %>%
        melt(measure.vars = c("RMSE", "Roughness"))
    df[, varname := factor(varname, labels = level.varname)]
    df[, variable   := factor(variable, c("RMSE", "Roughness"), labels = level.index)]
    df[, x := as.numeric(as.character(x))]

    d_lab <- expand.grid(variable = level.index, varname = level.varname) %>%
        cbind(lab = sprintf("(%s)", letters[1:4]), .)

    # 2. produce figure
    p <- ggplot(df, aes(x, value)) +
        # stat_summary(fun.data = box_qtl, geom = "errorbar", width = 0.5, color = "black") +
        # geom_boxplot2(show.errorbar = FALSE) +
        geom_line(size = 0.5) +
        geom_point(size = 1.7) +
        facet_grid(variable~varname, scale = "free", labeller = label_parsed) +
        geom_text(data = d_lab, aes(label = lab, fill = NULL), x = -Inf, y = Inf,
                  hjust = -0.2, vjust = 1.5, size = 5.5) +
        labs(y = NULL, x = NULL) +
        theme(
            strip.text.x = element_text(margin = margin(1,0,1,0, "pt")*5, size = 14),
            strip.text.y = element_text(face = "bold"),
            strip.text = element_text(size = 13, face = "bold"))

}

Fig9_para_sensitivity_try2 <- function(){
    # Figure_s1, parameter sensitivity envelope
    alphas <- c(.05, .1, .25, .5) %>% set_names(., .)
    d_enve <- llply(alphas, function(alpha){
        df[, as.list(quantile_envelope(value, alpha)), .(index, var)]
    }) %>% melt_list("alpha")
    d_enve$xmid <- as.character(d_enve$var) %>% str_extract("\\d{1,}") %>% as.numeric() %>% {./100}

    ggplot(d_enve, aes(xmid, ymin)) +
        # geom_point() + geom_density2d() +
        geom_ribbon(aes(x = xmid, ymin = ymin, ymax = ymax, fill = alpha)) +
        facet_wrap(~index, scales = "free", labeller=label_parsed) +
        geom_line(data = d_enve[alpha == 0.5, ], color ="white")
}

# MAIN -------------------------------------------------------------------------
Fig9_para_sensitivity()
