

source( here::here( "scripts/sim_functions.R") )

test_that("sim functions call", {
  
  pop <- generate_population( 300, 5, 
                              strata_config="balanced-exact", 
                              tau=2, 
                              bias="moderate", 
                              residual_var="heterogeneous",
                              seed=12344 )
  
  one_iter_srs(pop, n_code=60)
  
  #undebug( one_iter_optimal_strat )
  one_iter_optimal_strat(pop, n_code=60)
  
  alls <- one_iter_all(pop, n_code=60)  
  alls  

  expect_true( is.data.frame(alls) )
  expect_true( nrow(alls) == 4 )
  
})
