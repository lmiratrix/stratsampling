
test_that( "collapse ordering works", {

    x = sample( 1:10, replace=TRUE, 100 )
    x_collapsed = stratsampling:::collapse_ordered_adjacent( x, 4, "X" )
    table( x_collapsed )
    expect_equal( n_distinct(x_collapsed), 4 )

    x = sample( 1:3, replace=TRUE, 100 )
    x_collapsed = stratsampling:::collapse_ordered_adjacent( x, 4, "X" )
    table( x_collapsed )
    expect_equal( n_distinct(x_collapsed), 3 )

    x = sample( 1:4, replace=TRUE, 100 )
    x_collapsed = stratsampling:::collapse_ordered_adjacent( x, 4, "X" )
    table( x_collapsed )
    expect_equal( n_distinct(x_collapsed), 4 )  
    
    x_collapsed = stratsampling:::collapse_ordered_adjacent( x, 1, "X" )
    table( x_collapsed )
    expect_equal( n_distinct(x_collapsed), 1 )  
    
})


test_that("stratification works", {
    
    
    X1 = rnorm( 100 )
    X2 = rnorm( 100 )
    Z = sample( 0:1, replace=TRUE, 100 )
    
    str = stratsampling:::make_quantile_strata( X1, 2 )
    expect_equal( n_distinct(str), 2 )
    
    str = stratsampling:::make_quantile_strata( X1, 3 )
    expect_equal( n_distinct(str), 3 )
    
    str = stratsampling:::make_quantile_strata( X1, 4 )
    expect_equal( n_distinct(str), 4 )
    
    
    
    str = make_interaction_strata( X1, 3, "X1",
                                   X2, 4, "X2" )

    expect_equal( n_distinct(str), 12 )
        
    
    str = make_interaction_strata( X1, 3, "X1" )
    
    expect_equal( n_distinct(str), 3 )
    
    X3 = sample( 1:3, replace=TRUE, 100 )
    
    str = stratsampling:::make_quantile_strata( X3, 2 )
    expect_equal( n_distinct(str), 2 )
    
    str = stratsampling:::make_quantile_strata( X3, 2 )
    expect_equal( n_distinct(str), 2 )
    
    str = stratsampling:::make_quantile_strata( X3, 3 )
    expect_equal( n_distinct(str), 3 )
    
    str = stratsampling:::make_quantile_strata( X3, 3 )
    expect_equal( n_distinct(str), 3 )
    
    
    str = make_interaction_strata( X1, 3, "X1",
                                   X3 )
    expect_equal( n_distinct(str), 9 )
    
    str = make_interaction_strata( X1, 3, "X1",
                                   X3, 4 )
    expect_equal( n_distinct(str), 9 )

    str = make_interaction_strata( X1, 6, "X1",
                                   X3, 4 )
    expect_equal( n_distinct(str), 18 )
    
})


test_that( "calculating multiple stratificaton works", {
    
    dat = generate_population(N=1000, n_strat=4, 
                              bias = "large",
                              strata_config ="unbalanced",
                              tau = 2, R2_target = 0.8,
                              residual_var = "heterogeneous" )
    
    dat <- dat %>%
        rename( X2 = strata,
                Tx = Z ) %>%
        mutate( X1 = round( rnorm(nrow(dat)), digits = 2 ) )
    dat$X3 = cut( abs( dat$Yobs - dat$Yhat ),
                  breaks=quantile( abs( dat$Yobs - dat$Yhat ),
                                   probs=seq(0,1,length.out=4) ), include.lowest = TRUE,
                  labels = c("easy", "medium", "hard" ) )
    dat$X3 = as.numeric( dat$X3 )
    
    dat <- as_tibble( dat )

    str = make_interaction_strata( dat$X2, 3, "X2",
                                   dat$X1, 2, "X1" )
    
    expect_equal( length( str ), nrow( dat ) )
    
    # Create candidate stratifications
    candidates <- create_stratification_candidates(dat, 
                                                   vars = c("Yhat", "X2", "X1", "X3"),
                                                   Z = "Tx" )
    sapply( candidates, length )
        
    expect_true( all( sapply( candidates, length ) == nrow( dat ) ) )
    
})
