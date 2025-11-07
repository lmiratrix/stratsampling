
library( tidyverse )

test_that("the sampling function executes", {
  
  dta = generate_population(N=200, n_strat=4, oracle=FALSE )
  head( dta )
  
  #debug( generate_sample_index )
  S = generate_sample_index( dta, n_code = 30,
                             strata="strata", seed=123 )
  expect_true( is.numeric(S) )
  expect_equal( length(S), nrow(dta) )
  
  sdta = generate_observed_sample( dta, 
                                   Z="Z", strata="strata",
                                   n_code=30 )
  expect_true( is.data.frame(sdta ) )  
  
})





test_that("Sample index by strata works", {
  
  set.seed( 404404)
  pop <- generate_population( 30, 3, 
                              strata_config="unbalanced", 
                              tau=2, 
                              bias="moderate", 
                              residual_var="heterogeneous",
                              seed=123 )
  pop$strata = rep( c(1,2,3), c(15,10,5) )
  table( pop$strata )
  expect_true( is.data.frame(pop) )
  expect_true( all( c("strata","ID", "Yobs", "Yhat", "Z" ) %in% colnames(pop) ) )
  table( pop$strata, pop$Z )
  
  
  S = generate_sample_index( pop, n_code=c(2,2,2), strata="strata", seed=123 )
  tt <- table( S, pop$strata )
  tt[2,]
  expect_true( all( tt[2,] == 2 ) )
  
  S = generate_sample_index( pop, n_code=c(2,1,1,1,2,1), 
                             strata="strata", Z = "Z",
                             seed=123 )
  tt <- table( S, pop$strata )
  tt
  expect_equal( as.numeric( tt[2,] ), c( 3, 3, 2 ) )
  ss = paste( pop$Z, pop$strata, sep="-" )
  tt = table( ss, S )
  tt
  expect_equal( as.numeric( tt[,2] ), c( 2, 1, 1, 1, 2, 1 ) )
  
  
  # SRS by treatment
  S = generate_sample_index( pop, n_code=c(3, 6), strata=NULL,
                             Z = "Z", seed=123 )
  
  tt <- table( S, pop$Z  )
  tt
  expect_equal( as.numeric( tt[2,] ), c( 3, 6 ) )
  
  
  # SRS by treatment with total sample size
  pop$Z = rep( c(0,1), c(20,10) )
  #debug( generate_sample_index )
  S = generate_sample_index( pop, n_code=12, strata=NULL,
                             p_tx = 0.5,
                             Z = "Z", seed=123 )
  
  tt <- table( S=S, Z=pop$Z  )
  tt
  expect_equal( as.numeric( tt[2,] ), c( 6, 6 ) )
  
  # propotional to treatment group size
  S = generate_sample_index( pop, n_code=12, strata=NULL,
                             p_tx = NULL,
                             Z = "Z", seed=123 )
  
  tt <- table( S=S, Z=pop$Z  )
  tt
  expect_equal( as.numeric( tt[2,] ), c( 8, 4 ) )

    
})




test_that("the optimal allocation sampling functions work", {
  
  dta = generate_population(N=1000, n_strat=4, oracle=FALSE ) %>%
    dplyr::rename( ZZ = Z,
            SS = strata )
  
  head( dta )
  hits = dta$SS == 1 & dta$ZZ == 1
  dta$Yobs[ hits ] = dta$Yobs[ hits ] + rnorm( sum(hits), sd=10 )
  
 # debug(calc_optimal_allocation_oracle)
  calc <- calc_optimal_allocation_oracle( dta, n_code = 100, 
                                          strata="SS", 
                                          Z = "ZZ",
                                          return_counts = TRUE )
  calc
  expect_true( is.data.frame(calc) )
  expect_equal( nrow(calc), 8 )
  expect_true( which.max(calc$n.zk) == 5 )
  
  dd <- calc_optimal_allocation_oracle( dta, n_code = 100, strata="SS", Z = "ZZ",
                                        return_counts = FALSE )
  
  table( dd$.strata, dd$pi )
  head( dd )
  expect_equal( sum(dd$S), 100, tolerance = 0.02 )
  
  dd$.s2 = interaction( dd$SS, dd$ZZ )
  tt <- table( dd$.s2, dd$S )
  tt
  
  expect_equal( as.numeric(tt[,2]),
                calc$n.zk )
  

  
  # No strata testing
  dd <- calc_optimal_allocation_oracle( dta, n_code = 100, strata="SS", Z = NULL,
                                        return_counts = FALSE )
  head( dd )
  table( dd$.strata, dd$pi )
  d3 <- dd %>% group_by( .strata ) %>%
    summarize( n_un = n_distinct(pi),
               n = n(), 
               n_samp = sum(S) )
  d3
  expect_equal( d3$n_un, c(1,1,1,1) )
  expect_equal( sum(dd$S), 100, tolerance = 0.02 )
})





test_that("weight calculations work", {
  
  
  dta = generate_population(N=30, n_strat=3, oracle=FALSE ) %>%
    dplyr::rename( ZZ = Z,
            SS = strata )
  dta
  
  S = generate_sample_index( dta, n_code = 10 ,
                             strata="SS", Z = "ZZ", seed=1234 )
  
  dta$pi = calc_sampling_weights(dta, Z="ZZ", S=S, strata="SS" )
  dta$S = S
  
  sums <- dta %>% 
    group_by( SS, ZZ, pi ) %>%
    summarise( N = n(),
               n = sum(S),
               rt = n / N, 
               .groups = "drop")
  
  expect_equal( sums$pi, sums$rt )
  
} )





