#' add_noise
#'
#' QC variable also returned.
#' Same as MOD13A1 SummaryQA, 0: good value, 3: cloud contaminated.
#' 
#' For each site, repeat 100 times.
#' @param I If provide, noise position is fixed at `I`.
#' 
#' @export
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
#' simulate specific type noise 100 times through `add_noise`.
#' 
#' @param d data.frame, at least with the column of `ref` and `SummaryQA`.
#' @param type One of c("random", "real", "maxK", "maxDer")
#' @param times For each, how many times noise simulated?
#' 
#' @note Each time only process one site
#' 
#' @seealso [add_noise]
#' 
#' @examples
#' data("d_ref")
#' r <- simu_noise(d_ref, type = "real")
#' @export
simu_noise <- function(d, perc = 0.1, type = "random", times = 100){
    y         <- d$ref
    QC_factor <- d$SummaryQA
    # 生成两个噪点
    if (type == "random") {
        I <- NULL
    } else if (type == "real") {
        I <- which(QC_factor != "good")
    } else if (type == "maxK") {
        k <- der_k(y)
        I <- c(which.max(k), which.min(k))
    } else if (type == "maxDer") {
        f1 <- c(0, diff(y))
        I  <- c(which.max(f1), which.min(f1))
    }

    r <- llply(1:times, function(i){ add_noise(y, perc, i, I) }) %>% 
        purrr::transpose() %>% map(~do.call(rbind, .))
    if (type == "real") {
        QC <- as.numeric(QC_factor) - 1
        r$QC <- matrix(QC, ncol = length(QC), nrow = times, byrow = TRUE)
    }
    r
}

#' @references
#' Qiang Zhang, 2018, AFM, Eq.2
#' @export
der_k <- function(y) {
    f1 <- c(0, diff(y))
    f2 <- c(0, 0, diff(y, differences = 2))
    k  <- f2/((1+f1^2)^1.5)
    k
}
