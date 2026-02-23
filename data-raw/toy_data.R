## code to prepare `toy_data` dataset

library( stratsampling )

set.seed( 4044 )
dat = generate_population(N=1000, n_strat=4, 
                          bias = "large",
                          strata_config ="unbalanced",
                          tau = 2, R2_target = 0.8,
                          residual_var = "heterogeneous" )

head( dat )
ggplot( dat, aes(x=interaction(strata,Z), y=Yhat, color=Z)) + geom_boxplot()


dat <- dat %>%
    rename( X2 = strata,
            Tx = Z ) %>%
    mutate( X1 = round( rnorm(nrow(dat)), digits = 2 ) ) # add X1 variable to test that it doesn't get used in stratification


dat$X3 = cut( abs( dat$Yobs - dat$Yhat ),
              breaks=quantile( abs( dat$Yobs - dat$Yhat ), probs=seq(0,1,length.out=4) ), include.lowest = TRUE,
              labels = c("easy", "medium", "hard" ) )
dat$X3 = as.numeric( dat$X3 )
table( dat$X3 )
dat <- as_tibble( dat )
dat
ggplot( dat, aes(x=interaction(X3,Tx), y=Yhat, color=Tx)) +
    geom_boxplot()

str( dat )

# Create candidate stratifications
candidates <- create_stratification_candidates(dat, 
                                               vars = c("Yhat", "X2", "X1", "X3"),
                                               Z = "Tx" )
length( candidates )
sapply(candidates, length)

# Evaluate all candidates
comparison_precoding <- compare_stratifications_precoding(dat, 
                                                          candidates,
                                                          Z="Tx") 
comparison_precoding



### Make strata, take the sample and record it to make fake coded data

# NOTE: This has to align with the vignette

toy_data = dat

toy_data$strata <- make_interaction_strata( x1 = toy_data$Yhat, n1 = 3, "Yhat",
                                            x2 = toy_data$X3, n2 = 3, "X3",
                                            z_vec = toy_data$Tx)

toy_data$sample = generate_sample_index( toy_data, 100,
                                         Z = "Tx", strata = "strata",
                                         seed = 4044 )

table( toy_data$strata, toy_data$sample )


coded_dat <- filter( toy_data, sample == 1 ) %>%
    dplyr::select( ID, Yobs ) %>%
    arrange( ID ) %>%
    as_tibble()

coded_dat



# Save all the data ----


toy_data = as_tibble( dat ) %>%
    dplyr::select( ID, X1, X2, X3, Tx, Yhat )

toy_data
usethis::use_data(toy_data, overwrite = TRUE)

toy_data_coded = coded_dat
usethis::use_data(toy_data_coded, overwrite = TRUE)
