#' TSF_process of TIMESAT
#' Three curve fitting function (i.e. SG, AG and DL) in TIMESAT are used.
#' @export
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
