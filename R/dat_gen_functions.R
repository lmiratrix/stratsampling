# Code to generate fake data for simulation




generate_surrogate_noise <- function( N,
                                      n_strat,
                                      strata_config = c("balanced-exact", "balanced-approx", "unbalanced"),
                                      bias = c( "none", "moderate", "large", "extreme-contrast"),
                                      residual_var = c("homogeneous", "heterogeneous", "extreme-contrast"),
                                      sd_y = 3,
                                      R2_target = 0.5,
                                      rescale_v = TRUE ) {
    
    strata_config = match.arg(strata_config)
    bias = match.arg(bias)
    residual_var = match.arg(residual_var)
    
    # Generate strata (separately within treatment arms)
    strat = generate_strata_within_arm(N, n_strat, strata_config)
    
    # Calculate weights for strata
    wts_tab = as.numeric(table(factor(strat, levels = 1:n_strat)))
    wts =  wts_tab / sum(wts_tab)
    
    # Specify residual variance of predictions per scenario
    if ( residual_var == "extreme-contrast" & n_strat < 3 ) {
        stop("n_strat must be at least 3 for extreme-contrast residual variance configuration")
    }
    
    err_var_raw = switch(residual_var,  # variance multipliers across strata
                         "homogeneous" = rep(1, n_strat),
                         "heterogeneous" = seq(0.25, 4, length=n_strat),
                         "extreme-contrast"=c(0.1, rep(1, n_strat-2), 10))
    
    # Normalize to mean 1 (so we can preserve the target R^2)
    if ( rescale_v ) {
        err_var = err_var_raw / sum(err_var_raw * wts) 
    } else {
        err_var = err_var_raw
    }
    
    
    # Specify bias magnitude (units are relative to outcome standard deviation)
    bias_raw = switch(bias,
                      "none" = rep(0, n_strat),
                      "moderate" = seq(-1, 1, length.out = n_strat),
                      "large" = seq(-1, 1, length.out = n_strat), # up to 100% of outcome SD
                      "extreme-contrast" = c(-1, rep(0, n_strat - 2), 1) # up to 100% of outcome SD focused in only 2 strata
    )
    
    # Weight-centering so overall mean bias = 0 under unbalanced strata
    bias_centered  = bias_raw - sum(bias_raw * wts)
    max_abs_bias  = max(abs(bias_centered))
    b_scaled = if (max_abs_bias > 0) bias_centered / max_abs_bias else bias_centered
    
    if (bias=="none"){
        b = rep(0, n_strat)
    } else if (bias=="moderate"){
        b = b_scaled * 0.5 * sd_y  # max bias = 50% of outcome SD
    }
    else if (bias %in% c("large", "extreme-contrast")){
        b = b_scaled * sd_y # max bias up to 100% of outcome SD
    } else {
        stop( "unrecognized bias option" )
    }
    
    
    # Calibrate target R2
    target_err_var = sd_y^2*(1-R2_target)
    
    # Weighted variance of bias and weighted mean of noise variance
    b_bar_w = sum(wts * b)
    Var_b_w = sum(wts * (b - b_bar_w)^2)
    E_sig2_w = sum(wts * err_var)
    curr_err_var = Var_b_w + E_sig2_w
    
    c = if (curr_err_var > 0) sqrt(target_err_var / curr_err_var) else 0
    
    b = c * b
    sigma_eps = c * sqrt(err_var)
    
    # stratum-specific residual variance
    eps = rnorm(N, 0, sd=sigma_eps[as.integer(strat)])
    
    # ML predictions with bias
    strat_noise = b[as.integer(strat)] + eps
    
    dat = data.frame( strata=strat, 
                      e = strat_noise ) %>%
        dplyr::arrange( strata ) 
    
    return(dat)
}


#' Generate a population with strata and outcomes and
#' surrogate/predicted outcomes that depend on strata.
#'
#' These strata are generated post-treatment, in that we generate
#' n_strat strata within the treated and control units separately, and
#' then generate outcomes.
#'
#' @param no_tx If TRUE, do not generate treatment assignment (Z) and
#'   set all outcomes to control potential outcome.
#'
#' @export
generate_population <- function(N,
                                n_strat,
                                strata_config = c("balanced-exact", "balanced-approx", "unbalanced"),
                                tau = 0,
                                bias = c( "none", "moderate", "large", "extreme-contrast"),
                                residual_var = c("homogeneous", "heterogeneous", "extreme-contrast"),
                                sd_y = 3,
                                R2_target = 0.5,
                                no_tx = FALSE,
                                rescale_v = TRUE,
                                seed=NULL,
                                oracle = TRUE # what does this parameter do?
) {
    
    strata_config = match.arg(strata_config)
    bias = match.arg(bias)
    residual_var = match.arg(residual_var)
    
    if (!is.null(seed)){set.seed(seed)}
    
    # Generate the human-coded data
    Y0 = rnorm(N, 0, sd=sd_y) 
    Y1 = Y0 + tau
    
    if ( no_tx ) {
        Z = rep( 0, N ) 
    } else {
        # Randomize 50% to treatment
        Z = as.numeric( sample( N ) <= floor(N/2 ))
    }
    
    Yobs = ifelse(Z==1, Y1, Y0)
    
    dat = data.frame( Y0=Y0, Y1=Y1, Yobs=Yobs, Z=Z ) %>%
        arrange( Z )
    
    N1 = sum(dat$Z == 1)
    N0 = sum(dat$Z == 0)
    
    dat$Yhat = rep( NA, N )
    dat$strata = rep( NA, N )
    
    strat = generate_surrogate_noise( N = N0,
                                      n_strat = n_strat,
                                      strata_config = strata_config,
                                      bias = bias,
                                      residual_var = residual_var,
                                      sd_y = sd_y,
                                      R2_target = R2_target,
                                      rescale_v = rescale_v )
    stopifnot( n_distinct(strat$strata) == n_strat )
    
    # ML predictions with bias
    dat$Yhat[ dat$Z == 0 ] = dat$Yobs[ dat$Z == 0 ] + strat$e
    dat$strata[ dat$Z == 0 ] = strat$strata
    
    if ( N1 > 0 ) {
        strat1 = generate_surrogate_noise( N = N1,
                                           n_strat = n_strat,
                                           strata_config = strata_config,
                                           bias = bias,
                                           residual_var = residual_var,
                                           sd_y = sd_y,
                                           R2_target = R2_target,
                                           rescale_v = rescale_v  )
        
        stopifnot( n_distinct(strat1$strata) == n_strat )
        # ML predictions with bias
        dat$Yhat[ dat$Z == 1 ] = dat$Yobs[ dat$Z == 1 ] + strat1$e
        dat$strata[ dat$Z == 1 ] = strat1$strata
    }
    
    dat = dat %>%
        dplyr::arrange( Z, strata ) %>%
        mutate( ID = 1:n() ) %>%
        dplyr::relocate( Z, strata, ID ) 
    
    return(dat)
}



# Helper function to generate strata within a single treatment arm
generate_strata_within_arm <- function(N_arm, n_strat, strata_config) {
    
    if (N_arm %% n_strat != 0 & strata_config == "balanced-exact") {
        stop("N_arm must be divisible by n_strat for balanced-exact configuration")
    }
    
    strat_arm = switch(strata_config,
                       "balanced-exact" = sample(rep(1:n_strat, N_arm/n_strat)),
                       "balanced-approx" = sample(1:n_strat, size=N_arm, replace=TRUE,
                                                  prob=rep(1/n_strat, n_strat)),
                       "unbalanced" = sample(1:n_strat, size=N_arm, replace=TRUE,
                                             prob=runif(n_strat, min=0.2, max=0.8))
    )
    
    return(strat_arm)
}



#' Describe the DGP population defined by set parameters
#' 
#' Utility function to explore the different simulation factors
#' 
#' @export
describe_population <- function(n_strat,
                                strata_config = c("balanced-exact", "balanced-approx", "unbalanced"),
                                tau = 0,
                                bias = c( "none", "moderate", "large", "extreme-contrast"),
                                residual_var = c("homogeneous", "heterogeneous", "extreme-contrast"),
                                sd_y = 3,
                                R2_target = 0.5,
                                no_tx = FALSE,
                                rescale_v = TRUE,
                                seed=NULL,
                                verbose = FALSE ) {
    
    if ( is.data.frame(n_strat) ) {
        dat = n_strat 
    } else {
        dat = generate_population(N=100000,
                                  n_strat=n_strat,
                                  strata_config=strata_config,
                                  tau=tau,
                                  bias=bias,
                                  residual_var=residual_var,
                                  sd_y=sd_y,
                                  R2_target=R2_target,
                                  no_tx=no_tx,
                                  rescale_v = rescale_v,
                                  seed=seed )
    }
    
    sum <- dat %>%
        group_by( Z, strata ) %>%
        summarize( b = mean( Yhat - Yobs ),
                   sd_e = sd( Yhat - Yobs ),
                   n = n(), .groups = "drop" ) %>%
        group_by( Z ) %>%
        mutate( w = n / sum(n) )
    
    if ( verbose ) {
        print( knitr::kable( sum, digits = 2) )
    }
    
    agg <- sum %>%
        group_by( Z ) %>%
        summarise( Var_b = sum( w * b^2 ),
                   E_sig2 = sum( w * sd_e^2 ), .groups = "drop" ) %>%
        mutate( prop_bias = Var_b / (Var_b + E_sig2),
                R2 = 1 - (Var_b+E_sig2) / var( dat$Yobs ) ) %>%
        summarize( across( everything(), mean ) )
    
    if ( verbose ) {
        print( knitr::kable( agg, digits = 3) )
    }
    agg        
}





# Demo code ---- 

if ( FALSE ) {
    
    library( tidyverse )
    library(stratsampling)
    
    
    # Generate data
    dd <- generate_population(N=1000,
                              n_strat=4,
                              strata_config="unbalanced",
                              tau=2,
                              bias="large",
                              residual_var="heterogeneous",
                              sd_y=3,
                              R2_target=0.8,
                              no_tx=FALSE,
                              seed=1234)
    table( dd$strata, dd$Z )
    
    
    # Illustrate describe_population
    
    describe_population(4, 
                        bias = "extreme-contrast",
                        residual_var = "extreme-contrast",
                        rescale_v = TRUE )
    describe_population(4, 
                        bias = "extreme-contrast",
                        residual_var = "extreme-contrast",
                        rescale_v = FALSE )
    
    
    # Call describe_population on all the values for the arguments
    
    grid <- expand_grid(
        n_strat = 4,
        strata_config = c("balanced-exact", "balanced-approx", "unbalanced"),
        #strata_config = "balanced-exact",
        tau = c(0),
        bias = c( "none", "moderate", "large", "extreme-contrast"),
        residual_var = c("homogeneous", "heterogeneous", "extreme-contrast"),
        sd_y = c(3),
        R2_target = c(0.8),
        rescale_v = TRUE
    )
    
    
  
    
    results <- grid %>%
        rowwise() %>%
        mutate( desc = list( describe_population( n_strat = n_strat,
                                                  strata_config = strata_config,
                                                  tau = tau,
                                                  bias = bias,
                                                  residual_var = residual_var,
                                                  sd_y = sd_y,
                                                  R2_target = R2_target,
                                                  no_tx = FALSE,
                                                  rescale_v = rescale_v,
                                                  seed=1234,
                                                  verbose = FALSE ) ) ) %>%
        unnest( desc )
    
    results <- results %>%
        dplyr::select( -Z, -n_strat, -tau, - sd_y, -R2_target )
    
    results %>%
        knitr::kable( digits=2 )
    
    ggplot( results, aes( bias, prop_bias, color = strata_config, 
                          lty = rescale_v, group=interaction( strata_config, rescale_v) )) +
        facet_wrap( ~ residual_var ) +
        geom_point() + geom_line() +
        theme_minimal() + 
        coord_flip() +
        scale_y_continuous( limits=c(0,1) )
    
    
}
