# Unit test for combine_estimators() (in sim_functions.R)

# Checks precision-weighted pooling of two independent estimates against
# hand-computed values: the pooled point estimate is the inverse-variance
# weighted average, and each pooled SE is sqrt() of the weighted VARIANCES
# (not the weighted SEs) -- see the "+" units bug this guards against.

library( testthat )

source( "scripts/sim_functions.R" )

est1 <- list( est = 1, SEhat = 2, SEhat.full = 1.5, SEhat.samp = 1.2,
             S2.0 = 3, S2.1 = 4, SE = 2, method = "A" )
est2 <- list( est = 2, SEhat = 4, SEhat.full = 3.0, SEhat.samp = 2.6,
             S2.0 = 5, S2.1 = 6, SE = 4, method = "B" )

out <- combine_estimators( est1, est2 )

# prec = var2 / (var1 + var2) = 16 / (4 + 16) = 0.8
prec <- 0.8

test_that( "point estimate is the inverse-variance weighted average", {
    expect_equal( out$est, prec * 1 + (1 - prec) * 2 )
} )

test_that( "SEhat combines VARIANCES, not raw SEs, under the sqrt", {
    expect_equal( out$SEhat, sqrt( prec^2 * 2^2 + (1 - prec)^2 * 4^2 ) )
    # this is what the pre-fix bug would have given -- confirm we do NOT match it
    expect_false( isTRUE( all.equal(
        out$SEhat, sqrt( prec^2 * 2 + (1 - prec)^2 * 4 ) ) ) )
} )

test_that( "SEhat.full and SEhat.samp use the SAME weight as SEhat (additive decomposition)", {
    expect_equal( out$SEhat.full, sqrt( prec^2 * 1.5^2 + (1 - prec)^2 * 3.0^2 ) )
    expect_equal( out$SEhat.samp, sqrt( prec^2 * 1.2^2 + (1 - prec)^2 * 2.6^2 ) )
} )

test_that( "SE (oracle) uses its own weight, and combines variances", {
    # var1T = var2T = SE^2, both estimates have SE = 2 and SE = 4 here too,
    # so this happens to reuse the same prec, but is computed independently
    expect_equal( out$SE, sqrt( prec^2 * 2^2 + (1 - prec)^2 * 4^2 ) )
} )

test_that( "S2.0 / S2.1 are plain precision-weighted averages (not sqrt'd)", {
    expect_equal( out$S2.0, prec * 3 + (1 - prec) * 5 )
    expect_equal( out$S2.1, prec * 4 + (1 - prec) * 6 )
} )

test_that( "method labels are concatenated", {
    expect_equal( out$method, "A+B" )
} )

test_that( "equal-precision inputs average down to the simple mean", {
    e1 <- list( est = 1, SEhat = 2, SEhat.full = 1, SEhat.samp = 1,
               S2.0 = 3, S2.1 = 3, SE = 2, method = "A" )
    e2 <- list( est = 3, SEhat = 2, SEhat.full = 1, SEhat.samp = 1,
               S2.0 = 5, S2.1 = 5, SE = 2, method = "B" )
    o <- combine_estimators( e1, e2 )
    expect_equal( o$est, 2 )
    expect_equal( o$SEhat, sqrt( 0.5^2 * 4 + 0.5^2 * 4 ) )
    expect_equal( o$SEhat, sqrt(2) )
} )

cat( "All combine_estimators() checks passed.\n" )
