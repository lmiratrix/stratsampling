



library( tidyverse )

dta = generate_population(N=1000, n_strat=4 )

sdta = generate_observed_sample( dta, 
                                 Z="Z", strata="strata",
                                 n_code=200 )


test_that( "data generated correctly", {
  expect_true( all( c("strata","ID", "Yobs", "Yhat", "Z" ) %in% colnames(dta) ) )
  
  tt = table( dta$strata )
  expect_true( all(tt == 250 ) )
  
  head( dta )
  
  
  expect_true( is.data.frame(sdta) )
  expect_true( nrow(sdta) == 1000 )
  expect_true( all( c("strata","ID", "Yobs", "Yhat", "Z", "S" ) %in% colnames(sdta) ) )
  
  head( sdta )
  
  tt = table( sdta$Z, sdta$S )
  expect_true( all( tt > 0 ) )
} )


test_that("core srs estimators call correctly", {
  #undebug( est_ma_strat )
  e1 = est_ma_srs( data=sdta, 
                   Yhat="Yhat",
                   S = "S" )
  
  e1
  expect_true( is.na( e1$SE ) )
  expect_true( !is.na( e1$SEhat ) )
  expect_true( e1$SEhat > 0 )
  
  
  # Do we get oracle if we don't zap all the outcomes?
  e1b = est_ma_srs( data=dta, 
                    Yhat="Yhat",
                    S = sdta$S )
  e1b
  expect_equal( e1b$SEhat, e1$SEhat )
  expect_true( !is.na( e1$est ) )
  
  # Passing a new variable not in the dataframe
  e4 = est_ma_srs( data = sdta,
                   Yobs = ifelse( sdta$Z == 1, sdta$Yobs + 10, sdta$Yobs ),
                   Yhat = "Yhat",
                   S = sdta$S,
                   Z = sdta$Z )
  
  e4
  e1
  expect_equal( e4$est - e1$est, 10 )
  
  
  # Check we get actual variances when we give oracle fully observed outcomes
  e7 = est_ma_srs( data = dta,
                   Yhat = "Yhat",
                   S = sdta$S )
  e7
  expect_true( !is.na( e7$SE ) )
  expect_equal( e1$est, e7$est )
  
} )



test_that( "core strat estimators call correctly", {
  
  e2 = est_ma_strat( data = sdta,
                     Yhat = "Yhat" )
  
  e2b = est_ma_strat( data = dta,
                      Yhat = "Yhat",
                      S = sdta$S )
  
  expect_true( is.na( e2$SE ) )
  expect_true( e2$SEhat > 0 )
  expect_equal( e2b$SEhat, e2$SEhat )
  
  
  # Check different names
  d2 <- sdta %>%
    rename( ZZ = Z,
            SS = S,
            strata2 = strata,
            Yhat2 = Yhat )
  
  e3 = est_ma_strat( data = d2,
                     Yhat = "Yhat2",
                     S = "SS",
                     strata = "strata2",
                     Z = "ZZ" )
  expect_equal( e2, e3 )
  
  
  
  
  e4 = est_ma_strat( data = sdta,
                     Yobs = ifelse( sdta$Z == 1, sdta$Yobs + 10, sdta$Yobs ),
                     Yhat = "Yhat",
                     S = sdta$S,
                     strata = sdta$strata,
                     Z = sdta$Z )
  e2
  e4
  expect_equal( e2$est + 10, e4$est )
  
  
 
  e4
  
  
  # Shifting all predictions should not shift final estimate due to bias correction
  e5 = est_ma_strat( data = sdta,
                     Yhat = ifelse( sdta$Z == 1, sdta$Yhat + 10, sdta$Yhat ),
                     S = sdta$S,
                     strata = sdta$strata,
                     Z = sdta$Z )
  e2
  e5
  expect_equal( e2$est, e5$est )
  
  
  
  # Check we get actual variances when we give oracle fully observed outcomes
  e6 = est_ma_strat( data = dta,
                     Yhat = "Yhat",
                     S = sdta$S )
  e6
  expect_true( !is.na( e6$SE ) )
  expect_equal( e2$est, e6$est )
  expect_equal( e2$SEhat, e6$SEhat )
  
  
  
  
  
})



test_that( "power calculator works", {
  
  set.seed( 433443 )
  
  sdta = generate_observed_sample( dta, 
                                   Z="Z", strata="strata",
                                   n_code=200 )
  
  head( sdta )
  sdta$Yhat = sdta$Yhat + 3*sdta$strata - 0.5*sdta$Z
  

  sdta$Yobs_full = ifelse( sdta$Z == 1, sdta$Y1, sdta$Y0 )
  
  # Check power calcs if they align with estimation
  es = est_ma_strat( data = sdta,
                     Yobs = "Yobs_full",
                     Yhat = "Yhat", return_pieces = TRUE )
  
  es_srs = est_ma_strat( data = sdta,
                     Yobs = "Yobs_full", strata = NULL,
                     Yhat = "Yhat", return_pieces = TRUE )
  
  est_subset(data = sdta,
             Yobs = "Yobs_full" )
  
  es_srs$res
  
  es$res
  
  es$est.part2
  
  # Check full coding values
  hnd <- sdta %>%
    group_by( Z ) %>%
    summarize( varY = var( Yobs, na.rm=TRUE ),
               ney.hat = varY / n() )
  hnd
  es$est.part2$ney.hat
  expect_equal( sqrt( sum( es$est.part2$ney.hat ) ), sqrt( sum( hnd$ney.hat ) ), tolerance=0.001 )
  
  expect_equal( es$res$SEhat.full, sqrt( sum( hnd$ney.hat ) ), tolerance=0.001 )
  
  
  head( sdta )
  cor( sdta$Yobs_full, sdta$Yhat )^2
  
  R2s <- estimate_R2( data = sdta,
                      Yobs = "Yobs_full",
                      Yhat = "Yhat",
                      strata = "strata" )
  
  R2s
  
  R2 = 1 - var( sdta$Yobs_full - sdta$Yhat ) / var( sdta$Yobs_full )
  R2
  
  expect_equal( R2, R2s$R2_ML )
  # It's negative!
  
  var( sdta$Y0 )
  var( sdta$Y0[sdta$Z==0] )
  
  cc <- calc_nominal_var( N=nrow(sdta), 
                          p = mean(sdta$Z), 
                          h = 1,
                          var_Y0 = var( sdta$Yobs[sdta$Z==0], na.rm=TRUE ), 
                          R2_ML = R2, R2_strat = 0 )
  sqrt( cc )
  es$res$SEhat.full
  expect_equal( cc, es$res$SEhat.full^2, tolerance=0.05 )
  
  
  
  cc <- calc_nominal_var( N=nrow(sdta), 
                          p = mean(sdta$Z), 
                          h = 0.2,
                          var_Y0 = var( sdta$Yobs[sdta$Z==0], na.rm=TRUE ), 
                          R2_ML = R2, R2_strat = 0 )
  sqrt( cc )
  expect_equal( sqrt(cc), es_srs$res$SEhat, tolerance=0.1 )
  
  
  # Aside: how do our variance estimates align?
  var( sdta$Yobs[sdta$Z==0], na.rm=TRUE ) / var( sdta$Y0 )
  var( sdta$Yobs[sdta$Z==0], na.rm=TRUE ) / var( sdta$Y0[ sdta$Z==0 ] )
  var( sdta$Y0 ) / var( sdta$Y0[ sdta$Z==0 ] )
  
  
  pp <- plot_power_curve( N=nrow(sdta), 
                          p = mean(sdta$Z), 
                          var_Y0 = var( sdta$Yobs[sdta$Z==0], na.rm=TRUE ), 
                          R2_ML = R2, R2_strat = 0, h_min = 0.05 )
  
  tbl <- attr( pp, "table" )
  tbl
  min( tbl$SE )
  expect_equal( min( tbl$SE ), es$res$SEhat.full, tolerance=0.05 )

  R2_str <- sdta %>%
    group_by( strata, Z ) %>%
    summarise( var_e = var(Yobs_full - Yhat ),
               .groups = "drop")
  
  R2_str
  
  # How much within variance shrinks vs. overall (due to strata bias)
  R2_strat <- 1 - mean( R2_str$var_e ) / var( sdta$Yobs_full - sdta$Yhat )
  
  R2_strat
  expect_equal( R2_strat, R2s$R2_strat )
  (1-R2_strat) * (1-R2)
  
  pp <- plot_power_curve( N=nrow(sdta), 
                          p = mean(sdta$Z), 
                          var_Y0 = var( sdta$Yobs[sdta$Z==0], na.rm=TRUE ), 
                          R2_ML = R2, R2_strat = R2_strat, h_min = 0.05 )
  
  tbl <- attr( pp, "table" )
  tbl
  min( tbl$SE )
  expect_equal( min( tbl$SE ), es$res$SEhat.full, tolerance=0.02 )
  
  tt <- filter( tbl, h == 0.2 )
  tt
  
  es$res
  expect_equal( tt$SE[[1]], es$res$SE, tolerance = 0.1 )
  
  
  
  # Check very different h
  sdta = generate_observed_sample( dta, 
                                   Z="Z", strata="strata",
                                   n_code=500 )
  hh = mean( sdta$S )
  hh
  es = est_ma_strat( data = sdta,
                     Yobs = "Yobs",
                     Yhat = "Yhat")
  es
  ff = filter( tbl, abs(h - hh) < 0.005 )
  ff
  expect_equal( ff$SE[[1]], es$SEhat, tolerance = 0.10 )  
})




test_that( "different call structures work", {
  
  # All vector, no DF
  e2a = est_ma_strat( Yhat = sdta$Yhat,
                      Yobs = sdta$Yobs,
                      S = sdta$S,
                      Z = sdta$Z,
                      strata = sdta$strata )
  
  # All character names
  e2b = est_ma_strat( data = sdta, 
                      Yhat = "Yhat",
                      Yobs = "Yobs",
                      S = "S",
                      Z = "Z",
                      strata = "strata" )
  
  expect_equal( e2a, e2b )
  
  
  # Mix of both
  e2c = est_ma_strat( data = sdta, 
                      Yhat = "Yhat",
                      Yobs = sdta$Yobs,
                      S = "S",
                      Z = sdta$Z,
                      strata = "strata" )
  expect_equal( e2a, e2c )  
  
})


test_that( "no treatment works", {
  
  mn = mean( dta$Yobs, na.rm=TRUE )
  mn
  e2a = est_ma_strat( Yhat = sdta$Yhat,
                      Yobs = sdta$Yobs,
                      S = sdta$S,
                      Z = NULL,
                      strata = sdta$strata )
  
  expect_true( !is.na( e2a$SEhat ) )
  e2a$est
  mn
  expect_equal( e2a$est, mn, tolerance = 0.2 )  
  
})





test_that( "Small sample srs = 1 strat", {

  dta = generate_population(N=100, n_strat=2 )

  sdta = generate_observed_sample( dta,
                                   Z="Z", strata=NULL,
                                   n_code=10 )


  e1b = est_ma_srs( data=sdta,
                    Yhat="Yhat" )

  expect_true( !is.na( e1b$SEhat ) )
  expect_true( e1b$SEhat > 0 )

})


