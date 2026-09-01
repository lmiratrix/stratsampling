# Tests for pilot_functions.R: generalized optimal allocation with
# user-supplied variance parameters (var_alloc), pilot-aware budgeting
# (S_pilot), pilot variance estimation, and two-mode pre-coding
# diagnostics.

library( tidyverse )

# Shared fixture: heterogeneous residual variance so the optimal
# allocation is non-trivial.
set.seed( 20260810 )
dat <- generate_population( N = 1000, n_strat = 4,
                            strata_config = "balanced-exact",
                            tau = 2,
                            bias = "large",
                            residual_var = "heterogeneous",
                            R2_target = 0.5 )


test_that( "var_alloc equal to oracle variances reproduces oracle allocation", {
    
    counts_oracle <- calc_optimal_allocation_oracle( dat, n_code = 300,
                                                     Z = "Z", strata = "strata",
                                                     return_counts = TRUE )
    
    vt <- counts_oracle %>%
        dplyr::select( .Z, .strata, var.alloc = var.resid )
    
    counts_supplied <- calc_optimal_allocation( dat, n_code = 300,
                                                Z = "Z", strata = "strata",
                                                var_alloc = vt,
                                                return_counts = TRUE )
    
    expect_equal( counts_supplied$.strata, counts_oracle$.strata )
    expect_equal( counts_supplied$var.alloc, counts_oracle$var.alloc )
    expect_equal( counts_supplied$n.zk, counts_oracle$n.zk )
    
    # Wrapper agrees with the generalized function in oracle mode
    counts_direct <- calc_optimal_allocation( dat, n_code = 300,
                                              Z = "Z", strata = "strata",
                                              return_counts = TRUE )
    expect_equal( counts_direct$n.zk, counts_oracle$n.zk )
    
    # Same check with no treatment arm
    counts_oracle_noZ <- calc_optimal_allocation_oracle( dat, n_code = 300,
                                                         Z = NULL, strata = "strata",
                                                         return_counts = TRUE )
    vt_noZ <- counts_oracle_noZ %>%
        dplyr::select( .Z, .strata, var.alloc = var.resid )
    counts_supplied_noZ <- calc_optimal_allocation( dat, n_code = 300,
                                                    Z = NULL, strata = "strata",
                                                    var_alloc = vt_noZ,
                                                    return_counts = TRUE )
    expect_equal( counts_supplied_noZ$n.zk, counts_oracle_noZ$n.zk )
    
    # Missing cell in var_alloc is an error
    expect_error( calc_optimal_allocation( dat, n_code = 300,
                                           Z = "Z", strata = "strata",
                                           var_alloc = vt[-1, ],
                                           return_counts = TRUE ) )
} )


test_that( "pilot mode: budget accounting, floors, and pilot inclusion", {
    
    set.seed( 4044 )
    dat$S_pilot = generate_sample_index( dat, n_code = 50,
                                         Z = "Z", strata = "strata" )
    
    vt = estimate_pilot_variances( dat, S = "S_pilot",
                                   Z = "Z", strata = "strata" )
    expect_true( all( is.finite( vt$var.alloc ) ) )
    # generate_sample_index() rounds per cell and does not reconcile
    # totals, so compare against the realized pilot count
    expect_equal( sum( vt$n_pilot.zk ), sum( dat$S_pilot ) )
    
    # h = 0.30 total budget, pilot counted inside it
    n_code = round( 0.30 * nrow(dat) )
    
    set.seed( 505 )
    dd <- calc_optimal_allocation( dat, n_code = n_code,
                                   Z = "Z", strata = "strata",
                                   var_alloc = vt, S_pilot = "S_pilot" )
    
    # Total coded = round(h * N), and every pilot unit is in the sample
    expect_equal( sum(dd$S), n_code )
    expect_true( all( dd$S[ dd$S_pilot == 1 ] == 1 ) )
    
    # Realized per-cell counts match the (deterministic) allocation table,
    # respect the min_n floor, and never fall below the pilot counts
    counts <- calc_optimal_allocation( dat, n_code = n_code,
                                       Z = "Z", strata = "strata",
                                       var_alloc = vt, S_pilot = "S_pilot",
                                       return_counts = TRUE )
    realized <- dd %>%
        group_by( .Z, .strata ) %>%
        summarise( n = sum(S), n_pilot = sum(S_pilot), .groups = "drop" ) %>%
        arrange( .Z, .strata )
    
    expect_equal( realized$n, counts$n.zk )
    expect_equal( realized$n_pilot, counts$n_pilot.zk )
    expect_true( all( counts$n.zk >= pmin( 2, counts$N.zk ) ) )
    expect_true( all( counts$n.zk >= counts$n_pilot.zk ) )
    expect_equal( counts$n_new.zk, counts$n.zk - counts$n_pilot.zk )
    
    # Downstream estimation runs and returns finite results
    e <- est_ma_strat( dd, Yhat = "Yhat", S = "S", Z = "Z", strata = "strata" )
    expect_true( is.finite( e$est ) )
    expect_true( is.finite( e$SEhat ) & e$SEhat > 0 )
    
    # Over-piloted cell: pilot = ALL of stratum 1; allocation must keep
    # those units and still hit the total budget
    S_conc = as.integer( dat$strata == 1 )
    expect_warning(
        vt2 <- estimate_pilot_variances( dat, S = S_conc,
                                         Z = "Z", strata = "strata" )
    )
    dd2 <- calc_optimal_allocation( dat, n_code = 600,
                                    Z = "Z", strata = "strata",
                                    var_alloc = vt2, S_pilot = S_conc )
    expect_equal( sum(dd2$S), 600 )
    cc2 <- calc_optimal_allocation( dat, n_code = 600,
                                    Z = "Z", strata = "strata",
                                    var_alloc = vt2, S_pilot = S_conc,
                                    return_counts = TRUE )
    expect_true( all( cc2$n.zk >= cc2$n_pilot.zk ) )
    
    # Budget below the pilot floors triggers a warning
    expect_warning( calc_optimal_allocation( dat, n_code = 100,
                                             Z = "Z", strata = "strata",
                                             var_alloc = vt2, S_pilot = S_conc,
                                             return_counts = TRUE ) )
} )


test_that( "estimate_pilot_variances matches oracle when pilot = everyone", {
    
    vt_full <- estimate_pilot_variances( dat, S = rep( 1, nrow(dat) ),
                                         Z = "Z", strata = "strata" )
    
    oracle_counts <- calc_optimal_allocation_oracle( dat, n_code = 300,
                                                     Z = "Z", strata = "strata",
                                                     return_counts = TRUE )
    
    expect_equal( vt_full$var.alloc, oracle_counts$var.resid )
    expect_equal( vt_full$N.zk, oracle_counts$N.zk )
} )


test_that( "full pilot reproduces oracle stratification rankings", {
    
    dat2 <- dat %>%
        mutate( resid = Yobs - Yhat,
                W = round( rnorm( n() ), 1 ) )
    
    candidates <- create_stratification_candidates( dat2,
                                                    vars = c("Yhat", "W"),
                                                    Z = "Z" )
    
    cmp_oracle <- compare_stratifications( dat2, candidates, Z = "Z" )
    
    cmp_pilot <- compare_stratifications_precoding( dat2, candidates, Z = "Z",
                                                    Yobs = "Yobs",
                                                    S = rep( 1, nrow(dat2) ) )
    
    # With a 100% pilot, the pilot-based between-strata spread is the
    # oracle quantity, so values (and hence rankings) must agree exactly
    # on the candidates both methods evaluate.
    joined <- inner_join(
        cmp_oracle$summary %>% dplyr::select( stratification, oracle_spread = mean_spread ),
        cmp_pilot %>% dplyr::select( stratification, pilot_spread = mean_spread_pilot ),
        by = "stratification" )
    
    expect_true( nrow(joined) >= 5 )
    expect_equal( joined$pilot_spread, joined$oracle_spread, tolerance = 1e-10 )
    expect_equal( cor( joined$oracle_spread, joined$pilot_spread,
                       method = "spearman" ), 1 )
    
    # Proxy mode is unchanged: same columns as before, no pilot columns,
    # sorted by the Yhat-spread proxy
    cmp_proxy <- compare_stratifications_precoding( dat2, candidates, Z = "Z" )
    expect_false( "mean_spread_pilot" %in% names(cmp_proxy) )
    expect_true( all( c("stratification", "mean_spread",
                        "avg_within_var_Yhat") %in% names(cmp_proxy) ) )
    expect_equal( cmp_proxy$mean_spread,
                  sort( cmp_proxy$mean_spread, decreasing = TRUE ) )
} )