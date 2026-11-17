test_that("test_bg_p_bg", {
    set.seed(1)
    bg <- buildBackgroundModel(seq.example, background = "inputSet", 
                               numSubsamples = 2L)
    set.seed(1)
    bg <- buildBackgroundModel(seq.example, background = "nonInputSet", 
        proteome = proteome.example, numSubsamples = 2L)
})

test_that("test_bg_p_combinations", {
  set.seed(1)
    bg <- buildBackgroundModel(seq.example, background = "inputSet", 
        model = "anchored", numSubsamples = 2L)
})

test_that("test_bg_p_model", {
  set.seed(1)
    bg <- buildBackgroundModel(seq.example, targetPosition = "Nterminus", 
        proteome = proteome.example, numSubsamples = 2L)
    set.seed(1)
    bg <- buildBackgroundModel(seq.example, targetPosition = "Cterminus", 
        proteome = proteome.example, numSubsamples = 2L)
})

test_that("test_bg_p_others", {
  set.seed(1)
    bg <- buildBackgroundModel(seq.example, uniqueSeq = FALSE, 
        proteome = proteome.example, numSubsamples = 2L)
    set.seed(1)
    bg <- buildBackgroundModel(seq.example, replacement = TRUE, 
        proteome = proteome.example, numSubsamples = 2L)
})

