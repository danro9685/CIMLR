cimlr = CIMLR(X = GliomasReduced$in_X, c = 3, cores.ratio = 0)
input_data = rbind(GliomasReduced$in_X$point_mutations,GliomasReduced$in_X$copy_numbers,GliomasReduced$in_X$methylations,GliomasReduced$in_X$expression_values)
ranks = CIMLR_Feature_Ranking(A = cimlr$S, X = input_data)

context("CIMLR")
test_that("structure of output is compliant", {
    expect_equal(names(cimlr), c("y", "S", "F", "ydata",
        "alphaK", "execution.time", "converge", "LF"))
})

context("CIMLR ranking")
test_that("structure of output is compliant", {
    expect_equal(names(ranks), c("pval", "aggR"))
})
