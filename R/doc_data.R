#' Simulated random noise and keypoint noise.
#' 
#' 10\%, 30\% and 50\% simulated random noise and 3 type keypoint noise.
#' 
#' Each noise contains two variables:
#' * __`y`:__ vegetation enhanced index (EVI) time-series
#' * __`QC`:__ in the range of `[0-3]`, see `qc_summary`
#'  
#' @examples
#' data("lst_noise")
#' str(lst_noise, 2)
"lst_noise"

#' data.frame of Reference curve used to simulate noises.
#' 
#' @seealso [simu_noise]
#' 
#' @examples
#' data("d_ref")
#' print(head(d_ref))
"d_ref"

#' matirx of Reference curve used to evaluate GOF.
#' 
#' @seealso [gof_fitting]
#' 
#' @examples
#' data("mat_ref")
#' print(head(mat_ref))
"mat_ref"
