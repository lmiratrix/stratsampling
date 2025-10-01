
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



