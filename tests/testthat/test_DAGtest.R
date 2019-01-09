test_that("test_testDAU_group", {
    bg <- buildBackgroundModel(seq.example, proteome = proteome.example, 
        numSubsamples = 10L)
    t <- testDAU(seq.example, bg)
    expect_true(t@statistics["K", "AA0"] > 100)
    expect_equal(0, t@pvalue["K", "AA0"])
})

