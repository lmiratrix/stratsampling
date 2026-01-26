#' Compute oracle estimate (full coding) using shared functions
#' @return Data frame with oracle estimate and SE
#' 
#' @export
est_oracle <- function(data, Yobs = "Yobs", Z="Z") {
    data$Yobs = get_col_or_vector(data, Yobs)
    
    if ( !is.null(Z) ) {
        data$Z = get_col_or_vector(data, Z)
    } else {
        data$Z = 1
    }
    
    oracle_part <- data %>%
        group_by(Z) %>%
        summarise(
            N_z = n(),
            Ybar = mean(Yobs),
            var_Y = var(Yobs),
            .groups = "drop"
        )
    
    if (!is.null(Z)){ 
        est <- diff(oracle_part$Ybar)  # Treatment - Control
    }
    else{ est=oracle_part$Ybar}
    
    se <- sqrt(sum(oracle_part$var_Y / oracle_part$N_z))
    
    data.frame(
        method = "Oracle",
        est = est,
        SE = se
    )
}


#' Simple subset coding estimator
#'
#' This is the (Baseline) function to calculate tau_hat_subset (the
#' simple subset treatment effect estimator) assuming we have
#' human-coded scores, Yobs, for a subsample of n<N chosen using SRS
#'
#' @return A data.frame with columns:
#' \item{est}{the point estimate}
#' \item{SE}{The SE of the estimator if everything were coded}
#' \item{SEhat}{the estimated standard error of the estimator}
#' \item{SEhat.full}{NA (not applicable for this estimator)}
#' \item{SEhat.samp}{NA (not applicable for this estimator)}
#' \item{S2.0}{the sample variance of the coded outcomes in the control group}
#' \item{S2.1}{the sample variance of the coded outcomes in the treatment group}
#' 
#' @export
est_subset = function( data, Yobs = "Yobs",
                       S = "S", Z = "Z", finite_pop=FALSE ) {
  
  data$Yobs = get_col_or_vector(data, Yobs)
  data$S    = get_col_or_vector(data, S)
  data$Z    = get_col_or_vector(data, Z)
  
  
  if ( !is.null(Z) ) {
    data$Z = get_col_or_vector(data, Z)
  } else {
    data$Z = 1
  }
  
  n.code=sum(data$S)
  N=nrow(data)
  
  est.part = data %>%
    group_by( Z ) %>%
    dplyr::summarise(
      N=n(),
      var.y=var(Yobs),
      n.z = sum(S),
      Ybar.coded=mean(Yobs[S==1]),
      var.y.coded=var(Yobs[S==1]),
      .groups="drop"
    )
  
  
  if (finite_pop) { 

    # add finite population correction when Z is fixed
    est.part <- est.part %>% mutate( fpc = (N-n.z)/N )
    var.hat = sum( est.part$fpc * est.part$var.y.coded/est.part$n.z)
    
  } else {
    
    var.hat = sum(est.part$var.y.coded/est.part$n.z)

  } 
  
  # If we see all the outcomes we can use the oracle est of the 
  # variance of y terms.
  var = sum(est.part$var.y/est.part$n.z)
  
  if (nrow(est.part) == 2) {
    est = diff(est.part$Ybar.coded)
    
    res <- data.frame( est=est,
                       SE = sqrt(var),
                       SEhat = sqrt(var.hat),
                       SEhat.full=NA,
                       SEhat.samp=NA,
                       S2.0 = est.part$var.y.coded[[1]],
                       S2.1 = est.part$var.y.coded[[2]])
    
  } else if (nrow(est.part)==1) {
    est = est.part$Ybar.coded
    
    res <- data.frame( est=est,
                       SE = sqrt(var),
                       SEhat = sqrt(var.hat),
                       SEhat.full=NA,
                       SEhat.samp=NA,
                       S2 = est.part$var.y.coded[[1]] )
  } else {
    stop( "No rows or too many rows in est.part" )
  }
  
  return( res )
  
}




#' Calculate treatment impacts with stratified sample
#'
#' Function to calculate the model-assisted treatment effect estimator
#' assuming we have K fixed strata and hand-coded scores for a
#' subsample of n<N documents chosen by stratified random sampling
#' (with fixed and known sampling probabilities).
#' 
#' If strata set to null, then do the SRS version.
#' 
#' @return A data.frame with columns:
#' \item{est}{the point estimate}
#' \item{SE}{Not sure what this is.}
#' \item{SEhat}{the estimated standard error of the estimator}
#' \item{SEhat.full}{the estimated standard error component due to full coding variance
#' (i.e., if we had coded everything)}
#' \item{SEhat.samp}{the estimated standard error component due to sampling variance
#' (i.e., due to only coding a subset of documents) HOW DIFFERENT FROM SEhat?}
#' \item{S2.0}{the estimated population variance of the outcomes in the
#' control group}
#' \item{S2.1}{the estimated population variance of the outcomes in the
#' treatment group}
#' 
#' @export
est_ma_strat = function(data = NULL,
                        Yhat = "Yhat", Yobs ="Yobs",
                        S = "S", Z = "Z",
                        strata = "strata",
                        return_pieces = FALSE, finite_pop=FALSE ) {
  
  if ( is.null( data ) ) {
    data = data.frame( Yhat = Yhat,
                       Yobs = Yobs,
                       S = S )
  } else {
    data$Yhat = get_col_or_vector(data, Yhat)
    data$Yobs = get_col_or_vector(data, Yobs)
    data$S    = get_col_or_vector(data, S)
    #data$Z    = get_col_or_vector(data, Z)
  }
  
  if ( !is.null(Z) ) {
    data$Z = get_col_or_vector(data, Z)
  } else {
    data$Z = 1
  }
  
  if ( !is.null( strata ) ) {
    data$strata = get_col_or_vector(data, strata)
  } else {
    data$strata = "one"
  }
  
  # Calculate within-strata quantities
  est.part = data %>%
    group_by( Z, strata ) %>%
    dplyr::summarise(
      # First calculate effect estimate using model-assisted estimator
      N.zk = n(),
      n.zk = sum(S),
      Ybar = mean( Yhat ),
      adj = mean( Yobs[S==1] - Yhat[S==1] ),
      Ybar.adj = Ybar + adj,
      
      # variance components
      frac = (N.zk-n.zk)/n.zk,
      
      var.y = var(Yobs),
      var.resid = var(Yobs-Yhat),
      
      var.resid.coded = var(Yobs[S==1] - Yhat[S==1]),
      Ybar2.coded = mean(Yobs[S==1]^2),
      
      Ybar.obs.zk = mean(Yobs[S==1]),
      s2.obs.zk = var(Yobs[S==1]),
      
      .groups="drop" )
  
  # Aggregate over strata within each treatment arm
  est.part2 = est.part %>%
    group_by(Z) %>%
    dplyr::summarise(
      N.z=sum(N.zk),
      
      # Point estimates
      Ybar = weighted.mean( Ybar, w = N.zk ),
      adj = weighted.mean( adj, w = N.zk ),
      Ybar.adj = weighted.mean( Ybar.adj, w = N.zk ),
      
      # variance components
      inner.z = weighted.mean( frac * var.resid, w = N.zk ) / N.z,
      inner.z.coded = weighted.mean( frac * var.resid.coded, w = N.zk ) / N.z,
      
      
      # For calculating S2.hat
      Y2bar.obs = weighted.mean( Ybar2.coded, w = N.zk ),
      Ybar.obs = weighted.mean( Ybar.obs.zk, w = N.zk ),
      pt.3 = weighted.mean( (N.zk/N.z) * ((N.zk - n.zk)/(N.zk-1)) * (s2.obs.zk / n.zk),
                            w = N.zk ),
      
      .groups="drop" )
  
  # Variance estimator
  est.part2 <- est.part2 %>%
    mutate( S2.hat = (N.z/(N.z-1)) * ( Y2bar.obs - Ybar.obs^2 + pt.3 ),
            ney.hat = S2.hat / N.z )
  
  
  if (finite_pop) {
    # Only sampling variance component (fixed treatment assignments)
    var.hat = sum(est.part2$inner.z.coded)
    var = sum(est.part2$inner.z)  # Oracle: only sampling component
    
  } else {
    var.hat.A = sum( est.part2$ney.hat )
    var.hat.B = sum( est.part2$inner.z.coded )
    var.hat = var.hat.A + var.hat.B
    
    # Oracle variance estimate
    # For oracle, calculate full (not stratified) statistics
    oracle.var = data %>%
      group_by(Z) %>%
      dplyr::summarise( N=n(),
                        ney.y = var(Yobs) / N,
                        .groups = "drop" )
    var = sum( est.part2$inner.z ) + sum( oracle.var$ney.y )
  }
  
  # Point estimate
  if ( !is.null( Z ) & length(unique(data$Z))>1) {
    stopifnot( nrow(est.part2) == 2 )
    est = diff(est.part2$Ybar.adj)
  } else {
    stopifnot( nrow(est.part2) == 1 )
    est = est.part2$Ybar.adj
  }
  
  res <-  data.frame( est=est,
                      SE=sqrt(var),
                      SEhat = sqrt(var.hat),
                      SEhat.full = ifelse(finite_pop,0,sqrt(var.hat.A)),
                      SEhat.samp = ifelse(finite_pop, sqrt(var.hat), sqrt(var.hat.B) ))
  
  if ( !is.null( Z ) & length(unique(data$Z))>1) {
    res$S2.0 = est.part2$S2.hat[[1]]
    res$S2.1 = est.part2$S2.hat[[2]] # Est population variance
  } else {
    res$S2 = est.part2$S2.hat[[1]] # Est population variance
  }
  
  if ( return_pieces ) {
    return( list( res = res,
                  est.part = est.part,
                  est.part2 = est.part2,
                  oracle.var = oracle.var ) )
  }
  
  return( res )
  
  
}








#' Calculate treatment impacts with SRS sample
#' 
#' (Baseline) function to calculate tau_hat_ma (the model-assisted
#' treatment effect estimator) assuming we have human-coded scores,
#' Yobs, for a subsample of n<N documents plus predicted scores for all
#' N documents, with the subsample chosen using SRS
#'
#' @export
est_ma_srs = function(data,
                      Yhat = "Yhat", Yobs = "Yobs",
                      S = "S", Z = "Z",finite_pop=FALSE) {
  
  est_ma_strat( data, Yhat=Yhat, Yobs=Yobs, S=S, Z=Z, strata=NULL, finite_pop=finite_pop )
}







# For power calculations ----


#' Estimate what the SE would be given a stratification design
#'
#' This method calculates the SE given some possible population
#' summary statistics for the proportional stratified estimator. This
#' method is used for making power calculations when planning to
#' subsample documents.
#'
#' @param N The total population size
#' @param p The proportion of units assigned to treatment
#' @param h The proportion of units that will be human coded (i.e.,
#'   the sampling fraction)
#' @param var_Y0 The variance of the outcome in the control group
#' @param R2_ML The pseudo R-squared of the model's surrogate outcomes
#'   for predicting actual outcomes.
#' @param R2_strat The proportion of the residual variance explained by the
#'  stratification variable.
#'
#' @return The estimated variance of the model-assisted proportional
#'   stratified estimator.
#'
#' @export
calc_nominal_var = function( N, p, h, var_Y0, R2_ML, R2_strat ) {
  
  N1 = N * p
  N0 = N * (1-p)
  
  n1 = N1 * h
  n0 = N0 * h
  
  var_e = (1-R2_ML) * (1-R2_strat) * var_Y0
  
  SE2_full = (var_Y0 / N1) + (var_Y0 / N0)
  
  SE2_samp = ( (1/N1^2)*( N1*(N1-n1)/n1 ) + (1/N0^2)*( N0*(N0-n0)/n0 ) ) * var_e
  
  SE2_full + SE2_samp

}


#' Power curve table and plot for stratified vs SRS sampling
#'
#' Build the power-curve table (via `make_power_table()`) and/or a
#' ggplot (via `plot_power_curve()`) comparing stratified sampling to
#' SRS, given population parameters.
#'
#' You can estimate parameters from pilot data with [estimate_R2()].
#'
#' @name make_power_table
#'
#' @inheritParams calc_nominal_var
#' @seealso [calc_nominal_var()] for calculating for a single scenario
#'   and fixed h, and [estimate_R2()] for estimating design
#'   parameterd.
#' @param h_min Minimum proportion of documents coded to include in
#'   the plot/table.
#' @return `make_power_table()` returns a tibble with columns `h`,
#' `SE`, and `MDES`.
#'
#' `plot_power_curve()` returns a ggplot object with power curves. The
#' underlying data table is stored as an attribute `"table"`.
#' @export
#' 
make_power_table <- function( N, p, var_Y0, R2_ML, R2_strat, h_min = 0.05 ) {
  
  h_seq = seq( h_min, 1, by = 0.01 )
  
  SE_seq = sapply( h_seq, function(h) {
    sqrt( calc_nominal_var( N, p, h, var_Y0, R2_ML, R2_strat ) )
  } )
  
  MDES = 2.8*SE_seq
  
  pow_table = tibble( h = h_seq,
                      SE = SE_seq,
                      MDES = MDES )
  
  pow_table
}


#' @rdname make_power_table
#' @export
plot_power_curve <- function( N, p, var_Y0, R2_ML, R2_strat, h_min = 0.05 ) {
  
  strat_pow = make_power_table( N, p, var_Y0, R2_ML, R2_strat, h_min=h_min )
  srs_pow = make_power_table( N, p, var_Y0, R2_ML=R2_ML, R2_strat=0, h_min = h_min )
  
  pow_table = bind_rows( strat = strat_pow,
                    srs = srs_pow,
                    .id = "design" )
  
  full_MDES = min( pow_table$MDES )
    
  plt <- ggplot( pow_table, aes( h, MDES, color=design ) ) +
    geom_line() +
    labs( x = "Proportion of documents coded (h)",
          y = "MDES" ) +
    geom_hline( yintercept = full_MDES, lty = 2 ) +
    # expand y-axis to include 0
    expand_limits( y = 0 ) +
    theme_minimal()
  
  attr( plt, "table" ) <- pow_table
  
  plt
}



#' Given data, estimate the R2 of the stratification and ML
#'
#' This can be used on pilot data to get a rough sense of values for a
#' power calculation.
#'
#' This is a back of the envelope calculation that does not take
#' varying proportions of human coded outcomes within strata into
#' account.
#'
#' @seealso [calc_nominal_var()], [plot_power_curve()],
#'   [make_power_table()]
#' @inheritParams est_ma_strat
#' @return A tibble containing the estimated R2 of the stratification
#'   (how much machine predictions vary across strata) and of the
#'   surrogates for the outcome.
#'
#' @export
estimate_R2 <- function(data = NULL,
                                        Yhat = "Yhat", Yobs ="Yobs",
                                        S = NULL, Z = "Z",
                                        strata = "strata" ) { 
  
  if ( is.null( data ) ) {
    data = data.frame( Yhat = Yhat,
                       Yobs = Yobs,
                       S = S )
  } else {
    data$Yhat = get_col_or_vector(data, Yhat)
    data$Yobs = get_col_or_vector(data, Yobs)
    #data$Z    = get_col_or_vector(data, Z)
  }

  if ( !is.null(S) ) {
    data$S    = get_col_or_vector(data, S)
  } else {
    data$S = 1
  }
  
  if ( !is.null(Z) ) {
    data$Z = get_col_or_vector(data, Z)
  } else {
    data$Z = 1
  }
  
  if ( !is.null( strata ) ) {
    data$strata = get_col_or_vector(data, strata)
  } else {
    data$strata = "one"
  }
  
  var_e_full <- var( data$Yobs - data$Yhat, na.rm=TRUE )
  
  stats <- data %>%
    group_by( strata, Z ) %>%
    summarise( var_e = var( Yobs - Yhat, na.rm=TRUE ),
               var_Y = var( Yobs, na.rm=TRUE ),
               n = n(),
               .groups="drop" )
  
  R2s <- stats %>%
    summarise( R2_strat = 1 - weighted.mean( var_e, n ) / var_e_full,
               n = sum( n ),
               K = n() )
  
  R2s$R2_ML = 1 - var_e_full / var( data$Yobs, na.rm=TRUE )
  
  R2s
  
}




if ( FALSE ) {
  plot_power_curve( 1000, p=0.5, var_Y0=1, R2_ML=0.7, R2_strat=0.2 )
  
  plot_power_curve( 1000, p=0.5, var_Y0=1, R2_ML=0.7, R2_strat=0.5 )
  plot_power_curve( 1000, p=0.5, var_Y0=1, R2_ML=0.7, R2_strat=0 )
  
}




# Old code ----


#' (Baseline) function to calculate tau_hat_ma (the model-assisted
#' treatment effect estimator) assuming we have human-coded scores,
#' Yobs, for a subsample of n<N documents plus predicted scores for all
#' N documents, with the subsample chosen using SRS
#' 
#' This function was replaced.
#'
#' @noRd
est_ma_srs_old = function(data,
                          Yhat = "Yhat", Yobs = "Yobs",
                          S = "S", Z = "Z",
                          return_pieces = FALSE) {
  
  data$Yhat = get_col_or_vector(data, Yhat)
  data$Yobs = get_col_or_vector(data, Yobs)
  data$S    = get_col_or_vector(data, S)
  data$Z    = get_col_or_vector(data, Z)
  
  n.code=sum(data$S)
  N=nrow(data)
  data$pi = rep(n.code/N, N)
  
  est.part = data %>%
    group_by( Z ) %>%
    dplyr::summarise( # First calculate effect estimate using model-assisted estimator
      Ybar = mean( Yhat ),
      adj = weighted.mean( Yobs[S==1] - Yhat[S==1],  w=1/pi[S==1] ),
      Ybar.adj = Ybar + adj,
      
      # Next calculate "true" variance of tau.hat
      N = n(),
      var.y = var(Yobs),
      var.resid = var(Yobs-Yhat),
      
      # Finally, calculate the estimated variance of tau.hat based on partial coding
      var.y.coded = var(Yobs[S==1]),
      n.z = sum(S),
      var.resid.coded = var(Yobs[S==1] - Yhat[S==1]),
      .groups="drop"
    )
  
  est = diff(est.part$Ybar.adj)
  var = sum(est.part$var.y/est.part$N) + sum( (est.part$N - est.part$n.z)/est.part$N * (est.part$var.resid/est.part$n.z))
  
  var.hat.A = sum(est.part$var.y.coded/est.part$N)
  var.hat.B = sum( (est.part$N -est.part$n.z)/est.part$N * (est.part$var.resid.coded/est.part$n.z))
  
  
  
  
  res <- data.frame( est=est,
                     SE = sqrt(var),
                     SEhat = sqrt(var.hat.A + var.hat.B),
                     SEhat.full = sqrt(var.hat.A),
                     SEhat.samp = sqrt(var.hat.B),
                     S2.0 = est.part$var.y.coded[[1]] ,
                     S2.1 = est.part$var.y.coded[[2]] )
  
  if ( return_pieces ) {
    return( list( res = res,
                  est.part = est.part ) )
  }
  
  return( res )
  
}


