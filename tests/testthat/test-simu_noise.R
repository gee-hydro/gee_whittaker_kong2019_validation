test_that("multiplication works", {
    data("d_ref")

    expect_silent(r_real <- simu_noise(d_ref, type = "real"))
    expect_silent(r_rand <- simu_noise(d_ref, type = "random"))

    expect_equal(names(r_real), c("y", "QC"))
})
