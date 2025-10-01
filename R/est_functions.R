
#' (Baseline) function to calculate tau_hat_subset (the simple subset
#' treatment effect estimator) assuming we have human-coded scores,
#' Yobs, for a subsample of n<N chosen using SRS
#'
#' @export
est_subset = function(data, Yobs = "Yobs",
                          S = "S", Z = "Z") {

    data$Yobs = get_col_or_vector(data, Yobs)
    data$S    = get_col_or_vector(data, S)
    data$Z    = get_col_or_vector(data, Z)

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

    est = diff(est.part$Ybar.coded)
    var = sum(est.part$var.y/est.part$N)

    var.hat = sum(est.part$var.y.coded/est.part$n.z)

    res <- data.frame( est=est,
                       SE = sqrt(var),
                       SEhat = sqrt(var.hat),
                       SEhat.full=NA,
                       SEhat.samp=NA,
                       S2.0 = est.part$var.y.coded[[1]] ,
                       S2.1 = est.part$var.y.coded[[2]] )


    return( res )

}

#' (Baseline) function to calculate tau_hat_ma (the model-assisted
#' treatment effect estimator) assuming we have human-coded scores,
#' Yobs, for a subsample of n<N documents plus predicted scores for all
#' N documents, with the subsample chosen using SRS
#'
#' @export
est_ma_srs_old = function(data,
                          Yhat = "Yhat", Yobs = "Yobs",
                          S = "S", Z = "Z",
                          return_pieces = FALSE ) {

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






#' Calculate treatment impacts with stratified sample
#'
#' Function to calculate the model-assisted treatment effect estimator
#' assuming we have K fixed strata and hand-coded scores for a subsample of n<N documents
#' chosen by stratified random sampling (with fixed and known sampling probabilities)
#'
#' @export
est_ma_strat = function(data = NULL,
                        Yhat = "Yhat", Yobs ="Yobs",
                        S = "S", Z = "Z",
                        strata = "strata",
                        return_pieces = FALSE ) {

  if ( is.null( data ) ) {
    data = data.frame( Yhat = Yhat,
                       Yobs = Yobs,
                       S = S )
  } else {
    data$Yhat = get_col_or_vector(data, Yhat)
    data$Yobs = get_col_or_vector(data, Yobs)
    data$S    = get_col_or_vector(data, S)
    data$Z    = get_col_or_vector(data, Z)
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

  var.hat.A = sum( est.part2$ney.hat )
  var.hat.B = sum( est.part2$inner.z.coded )


  # Oracle variance estimate
  # For oracle, calculate full (not stratified) statistics
  oracle.var = data %>%
    group_by(Z) %>%
    dplyr::summarise( N=n(),
                      ney.y = var(Yobs) / N,
                      .groups = "drop" )
  var = sum( est.part2$inner.z ) + sum( oracle.var$ney.y )


  # Point estimate
  if ( !is.null( Z ) ) {
    stopifnot( nrow(est.part2) == 2 )
    est = diff(est.part2$Ybar.adj)
  } else {
    stopifnot( nrow(est.part2) == 1 )
    est = est.part2$Ybar.adj
  }

  res <-  data.frame( est=est,
                      SE=sqrt(var),
                      SEhat = sqrt(var.hat.A + var.hat.B),
                      SEhat.full = sqrt(var.hat.A),
                      SEhat.samp = sqrt(var.hat.B) )

  if ( !is.null( Z ) ) {
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





#'
#' @export
est_ma_srs = function(data,
                      Yhat = "Yhat", Yobs = "Yobs",
                      S = "S", Z = "Z") {

  est_ma_strat( data, Yhat=Yhat, Yobs=Yobs, S=S, Z=Z, strata=NULL )
}
