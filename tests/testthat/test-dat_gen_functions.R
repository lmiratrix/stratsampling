
library( tidyverse )

test_that("DGP works", {
  
  pop <- generate_population( 24, 4, 
                              strata_config="balanced-exact", 
                              tau=2, 
                              bias="moderate", 
                              residual_var="heterogeneous",
                              seed=123 )
  expect_true( is.data.frame(pop) )
  expect_true( all( c("strata","ID", "Yobs", "Yhat", "Z" ) %in% colnames(pop) ) )
  
  S = generate_sample_index( pop, n_code=6, strata="strata", seed=123 )
  expect_true( is.numeric(S) )
  expect_equal( length(S), nrow(pop) )
  sum(S)
  tt <- table( S, pop$strata )
  tt
  expect_true( all( tt[2,] == 2 ) )
  
  set.seed( 404404 )
  pop <- generate_population( 24, 3, 
                              tau=2, 
                              bias="moderate", 
                              residual_var="heterogeneous" )
  table( pop$strata, pop$Z )
  S = generate_sample_index( pop, n_code=6, Z = "Z", strata="strata", seed=123, p_tx = 0.5 )
  sum(S)
  expect_equal( sum(S), 6 )
  
  tt <- table( pop$strata[S==1], pop$Z[S==1] )
  tt
  expect_true( all( tt == 1 ) )
  
  dd <- generate_observed_sample( pop, n_code=6, strata="strata", by_tx=FALSE, seed=123 )  
  dd
  expect_true( is.data.frame(dd) )
  
})



test_that( "DGP controls variance correctly", {
  
  set.seed( 53434 )
  
  library( tidyverse )
  dd <- generate_population(N=1000,
                            n_strat=4,
                            strata_config="unbalanced",
                            tau=2,
                            bias="large",
                            residual_var="heterogeneous",
                            sd_y=3,
                            R2_target=0.8,
                            no_tx=FALSE )
  tt <- table( dd$strata, dd$Z )
  tt
  
  expect_true( sd( tt ) > 20 )
  
  ds <- dd %>% group_by( Z ) %>%
    summarise( sd = sd( Yobs ) )
  ds
  
  expect_equal( ds$sd, c( 3, 3 ), tolerance = 0.05 )  

  M = lm( Yobs ~ Z + Yhat, data=dd )
  s <- summary( M )  
  s
  #expect_equal( s$r.squared, 0.8, tolerance = 0.1 )  
  var( dd$Yhat - dd$Yobs ) / var( dd$Yobs )
  expect_equal( 1 - var( dd$Yhat - dd$Yobs ) / var( dd$Yobs ), 0.8, tolerance = 0.025 )
  
  
  set.seed( 53434 )
  
  dd <- generate_population(N=10000,
                            n_strat=4,
                            strata_config="balanced-exact",
                            tau=2,
                            bias="extreme-contrast",
                            residual_var="extreme-contrast",
                            sd_y=10,
                            R2_target=0.30,
                            no_tx=TRUE )
  tt <- table( dd$strata, dd$Z )
  expect_equal( sd( tt ), 0 )
  
  ds <- dd %>% group_by( Z ) %>%
    summarise( sd = sd( Yobs ) )
  expect_equal( ds$sd, 10, tolerance = 0.01 )  
  
  M = lm( Yobs ~ Yhat, data=dd )
  s <- summary( M )  
  s
  
  e = dd$Yobs - dd$Yhat
  var( e )
  var( e ) / var( dd$Yobs )
  #expect_equal( s$r.squared, 0.8, tolerance = 0.1 )  
  var( dd$Yhat - dd$Yobs ) / var( dd$Yobs )
  expect_equal( 1 - var( dd$Yhat - dd$Yobs ) / var( dd$Yobs ), 0.3, tolerance = 0.05 )
  
  
  
  
  # Third, most extreme test
  set.seed( 534344 )
  
  dd <- generate_population(N=5000,
                            n_strat=4,
                            strata_config="unbalanced",
                            tau=2,
                            bias="extreme-contrast",
                            residual_var="extreme-contrast",
                            sd_y=10,
                            R2_target=0.55,
                            no_tx=TRUE )
  tt <- table( dd$strata, dd$Z )
  tt
  
  expect_true( sd( tt ) > 40 )
  
  ds <- dd %>% group_by( Z ) %>%
    summarise( sd = sd( Yobs ) )
  ds
  
  expect_equal( ds$sd, 10, tolerance = 0.01 )  
  
  M = lm( Yobs ~ Yhat, data=dd )
  s <- summary( M )  
  s
  
  #expect_equal( s$r.squared, 0.8, tolerance = 0.1 )  
  var( dd$Yhat - dd$Yobs ) / var( dd$Yobs )
  expect_equal( 1 - var( dd$Yhat - dd$Yobs ) / var( dd$Yobs ), 0.55, tolerance = 0.025 )
  
  
})