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
#' @export
generate_population <- function(N, n_strat,
                                strata_config = c("balanced-exact", "balanced-approx", "unbalanced"),
                                tau = 0,
                                bias = c( "none", "moderate", "large", "extreme-contrast"),
                                residual_var = c("homogeneous", "heterogeneous", "extreme-contrast"),
                                sd_y = 3,
                                R2_target = 0.5,
                                seed=NULL,
                                oracle = TRUE # what does this parameter do?
) {
    
    strata_config = match.arg(strata_config)
    bias = match.arg(bias)
    residual_var = match.arg(residual_var)
    
    if (!is.null(seed)){set.seed(seed)}
    
    
    Z = as.numeric( sample( N ) <= floor(N/2 ))
    N1 = sum(Z == 1)
    N0 = sum(Z == 0)
    
    # Calibrate baseline residual variance of predictions to achieve target R^2 between Yobs and Yhat (overall)
    s2_sigma = R2_target * sd_y^2
    s2_eps = (1 - R2_target)*sd_y^2
    
    
    # Generate strata (separately within treatment arms)
    strat = numeric(N)
    
    control_ids = which(Z == 0)
    strat[control_ids] = generate_strata_within_arm(N0, n_strat, strata_config)
    
    treatment_ids = which(Z == 1)
    strat[treatment_ids] = generate_strata_within_arm(N1, n_strat, strata_config)
    
    
    # Calculate weights for when strata are unbalanced
    wts_tab = as.numeric(table(factor(strat, levels = 1:n_strat)))
    wts =  wts_tab / sum(wts_tab)
    
    # Specify residual variance of predictions per scenario
    if ( residual_var == "extreme-contrast" & n_strat < 3 ) {
        stop("n_strat must be at least 3 for extreme-contrast residual variance configuration")
    }
    
    var_mult_raw = switch(residual_var,  # variance multipliers across strata
                          "homogeneous" = rep(1, n_strat),
                          "heterogeneous" = seq(0.25, 4, length=n_strat),
                          "extreme-contrast"=c(0.1, rep(1, n_strat-2), 10))
    
    # Normalize to mean 1 (to preserve the target R^2)
    var_mult = var_mult_raw / sum(var_mult_raw * wts) 
    
    
    # Specify bias magnitude (units are relative to outcome standard deviation)
    bias_raw = switch(bias,
                      "none" = rep(0, n_strat),
                      "moderate" = seq(-1, 1, length.out = n_strat),
                      "large" = seq(-1, 1, length.out = n_strat), # up to 100% of outcome SD
                      "extreme-contrast" = c(-1, rep(0, n_strat - 2), 1) # up to 100% of outcome SD focused in only 2 strata
    )
    
    # Weight-centering so overall mean bias = 0 under unbalanced strata
    bias_centered  = bias_raw - sum(bias_raw * wts)
    b_scaled = bias_centered/max(abs(bias_centered))
    
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
    
    
    # Outcome generation
    # Yobs = Y.true + eps
    # Yhat = Y.true + bias
    Y.true = rnorm(N, 0, sd=sqrt(s2_sigma)) # latent signal
    
    # stratum-specific residual variance
    eps = rnorm(N, 0, sqrt(s2_eps * var_mult[strat]))
    
    # human-coded outcomes
    Y0 = Y.true + eps
    Y1 = Y0 + tau
    
    # ML predictions with bias
    Yhat = Y.true + b[as.integer(strat)]
    
    Yobs = ifelse(Z==1, Y1, Y0)
    
    #dat = data.frame(ID=1:N, Y0=Y0, Y1=Y1, Yhat0=Yhat0, Yhat1=Yhat1, strata=strat)
    dat = data.frame(ID=1:N, Yobs=Yobs, Yhat=Yhat, Z=Z, strata=strat) %>%
        dplyr::relocate( strata, ID ) %>%
        dplyr::arrange( strata )
    
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


