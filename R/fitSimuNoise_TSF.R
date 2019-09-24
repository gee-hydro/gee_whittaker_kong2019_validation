#' fitSimuNoise_TSF
#' 
#' @param ... ignored
#' 
#' @export
fitSimuNoise_TSF <- function(
    lst_noise_random, lst_noise_keypoint, 
    nptperyear, nyear, 
    indir, ...)
{
    # 1. random
    r <- foreach(df_sim = lst_noise_random, i = icount(),
                 .packages = c("rTIMESAT")) %dopar% 
    {
         # print(perc)
         TSF_main_df(df_sim, perc, nyear, nptperyear, overwrite = FALSE, wait = TRUE, indir = indir)
    }

    # 2. keypoint
    perc = 0.1
    r <- foreach(df_sim = lst_noise_keypoint,
                 type   = names(lst_noise_keypoint), i = icount(),
                 .packages = c("rTIMESAT")) %dopar% 
    {
         # perc <- percs[i]
         # print(perc)
         indir <- sprintf("%s/%s", indir, type)
         TSF_main_df(df_sim, perc, nyear, nptperyear, overwrite = FALSE, wait = TRUE, indir = indir)
     }

    ## tts to txt
    files <- dir(indir, pattern = "*.tts", recursive = T, full.names = T)
    # keypoint files
    files <- files[-(7:15)]
    outfiles <- files %>% gsub("(/perc_10/TSF_whit_eval_noise10)|_fit.tts", "", .) %>%
        strsplit("_") %>%
        {sprintf("%s/fitting_%s_%s.txt", indir, map_chr(., ~.[2]), map_chr(., ~.[1]))}

    r <- foreach(file = files, outfile = outfiles, i = icount()) %do% {
        # print(file)
        TSF_fit2time(file, 1, 1e8, 1, 1, wait = F,
                     outdir = indir, outfile = outfile)
    }
}
