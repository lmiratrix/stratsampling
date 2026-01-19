

# Calibration simulation

# This simulation just runs the code to check our variance estimators


library( tidyverse )

library( stratsampling )


one_run <- function( ) {

  dta = generate_population(N=1000, n_strat=4, oracle=FALSE,
                            tau = 1,
                            bias = "large",
                            residual_var = "heterogeneous" )

  d_srs = generate_observed_sample( dta,
                                    Z="Z", strata=NULL,
                                    n_code=200 )

  e_srs = est_ma_strat( data=d_srs,
                        Yhat="Yhat",
                        strata = NULL,
                        S = "S"  )

  e_srs_old = est_ma_srs_old( data=d_srs,
                              Yhat="Yhat",
                              S = "S"  )


  d_strat = generate_observed_sample( dta,
                                      Z="Z", strata="strata",
                                      n_code=200 )


  e_strat = est_ma_strat( data=d_strat,
                          Yhat="Yhat",
                          strata = "strata",
                          S = "S"  )

  # Neyman
  SEfull = sqrt( var( dta$Yobs[dta$Z==1] )/ sum(dta$Z==1) +
                   var( dta$Yobs[dta$Z==0] )/ sum(dta$Z==0) )

  oracle = tibble( est = mean( dta$Yobs[dta$Z==1] ) - mean( dta$Yobs[dta$Z==0] ),
                   SE = NA,
                   SEhat = SEfull,
                   SEfull = SEfull,
                   S2.1 = var( dta$Yobs[dta$Z==1] ),
                   S2.0 = var( dta$Yobs[dta$Z==0] ) )

  bind_rows( oracle = oracle,
             SRS = e_srs,
             SRS_old = e_srs_old,
             Strat = e_strat, .id = "method" ) %>%
    as_tibble()
}


R = 100
rps = map( 1:R, ~ one_run(), .progress = TRUE)

rps <- bind_rows( rps, .id="RunID" )
rps

rps %>%
  group_by( method ) %>%
  summarise( SE = sd( est ),
             SEhat = sqrt( mean( SEhat^2 )),
             SEfull = sqrt( mean( SEfull^2 )),
             mean.est = mean( est ),
             S2.0 = mean( S2.0 ),
             S2.1 = mean( S2.1 ) )


