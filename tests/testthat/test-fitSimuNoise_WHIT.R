test_that("fitSimuNoise_WHIT and GOF works", {
    data("mat_ref")
    data('lst_noise')
    r <- fitSimuNoise_WHIT(lst_noise, limits = 100, wmin = 0.1, perc_wc = 0.5)

    d_gof <- gof_fitting(mat_ref[1, , drop = FALSE], r$wWHd$`10%`)
    expect_equal(nrow(r$wWHd$real), 100)
})
