# Code to generate fake data for simulation


#' Generate a population with strata and outcomes and
#' surrogate/predicted outcomes that depend on strata.
#'
#' These strata are generated post-treatment, in that we generate
#' n_strat strata within the treated and control units separately, and
#' then generate outcomes.   That said, DGP makes strata parameters
#' the same for corresponding tx and co strata.
#'
#' TODO: Review the DGP -- is what we want?
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
                                 seed=NULL,
                                 oracle = TRUE # what does this parameter do?
) {
    
    strata_config = match.arg(strata_config)
    bias = match.arg(bias)
    residual_var = match.arg(residual_var)
    
    if (!is.null(seed)){set.seed(seed)}
    
    if ( no_tx ) {
        Z = rep( 0, N ) 
    } else { 
        Z = as.numeric( sample( N ) <= floor(N/2 ))
    }
    N1 = sum(Z == 1)
    N0 = sum(Z == 0)
    
    
    # Generate strata (separately within treatment arms)
    strat = numeric(N)
    
    control_ids = which(Z == 0)
    strat[control_ids] = generate_strata_within_arm(N0, n_strat, strata_config)
    
    if ( N1 > 0 ) {
        treatment_ids = which(Z == 1)
        strat[treatment_ids] = generate_strata_within_arm(N1, n_strat, strata_config)
    }
    
    # Calculate weights for when strata are unbalanced
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
    
    # Normalize to mean 1 (to preserve the target R^2)
    err_var = err_var_raw / sum(err_var_raw * wts) 
    
    
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
    
    # human-coded outcomes
    Y0 = rnorm(N, 0, sd=sd_y) 
    Y1 = Y0 + tau
    
    Yobs = ifelse(Z==1, Y1, Y0)
    
    # stratum-specific residual variance
    eps = rnorm(N, 0, sd=sigma_eps[as.integer(strat)])
    
    # ML predictions with bias
    Yhat = Yobs + b[as.integer(strat)] + eps
    
    
    #dat = data.frame(ID=1:N, Y0=Y0, Y1=Y1, Yhat0=Yhat0, Yhat1=Yhat1, strata=strat)
    dat = data.frame( Yobs=Yobs, Yhat=Yhat, Z=Z, strata=strat) %>%
        dplyr::arrange( strata ) %>%
        mutate( ID = 1:n() ) %>%
        dplyr::relocate( strata, ID ) 
    
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


