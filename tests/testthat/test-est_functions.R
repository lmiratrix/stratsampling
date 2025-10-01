



library( tidyverse )

dta = generate_population(N=1000, n_strat=4, oracle=FALSE )

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
  
  e1c = est_ma_srs_old( data=sdta, 
                        Yhat="Yhat",
                        S = "S" )
  
  e1
  e1c
  
  # Should they be exact or just close?
  expect_equal( e1, e1c, tolerance = 0.0001 )
  
  
  # Passing a new variable not in the dataframe
  e4 = est_ma_srs( data = sdta,
                   Yobs = ifelse( sdta$Z == 1, sdta$Yobs + 10, sdta$Yobs ),
                   Yhat = "Yhat",
                   S = sdta$S,
                   Z = sdta$Z )
  
  e4
  e1
  expect_equal( e4$est - e1$est, 10 )
  

  if ( FALSE ) {
    e1 = est_ma_strat( data=sdta, 
                       Yhat="Yhat",
                       Yobs =dta$Yobs,
                       strata = NULL,
                       S = "S", return_pieces = TRUE  )
    
    
    e1c = est_ma_srs_old( data=sdta, 
                          Yobs = dta$Yobs,
                          Yhat="Yhat",
                          S = "S", return_pieces = TRUE )
    
    e1
    e1c
    e1$est.part2
    e1c$est.part
    
    
  }
  
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
  
  dta = generate_population(N=100, n_strat=2, oracle=FALSE )
  
  sdta = generate_observed_sample( dta, 
                                   Z="Z", strata=NULL,
                                   n_code=10 )
  
  
  e1b = est_ma_srs( data=sdta, 
                    Yhat="Yhat" )
  
  e1c = est_ma_srs_old( data=sdta, 
                        Yhat="Yhat" )
  
  e1b
  e1c
  
  # Should this be more exact?   The "full" is differnet---the correction?
  expect_equal( e1b$SEhat, e1c$SEhat, tolerance = 0.001 )
  
  # Should they be exact or just close?
  expect_equal( e1b, e1c, tolerance = 0.01 )
  
})


