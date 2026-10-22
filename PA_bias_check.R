


N = 100000
t1 = rexp( N ) - 1
t2A = rexp( N ) - 1
t2B = rnorm( N, 0, 0.1 )

big_samp = t1 < 2
t2 = ifelse( big_samp, t2B, t2A )

varB = ifelse( big_samp, 0.1^2, 1^2 )
t_fin = (1/(1+varB)) * t2 + (varB/(1+varB)) * t1
se_fin = sqrt( varB / (1+varB) )

library( tidyverse )

res <- tibble( t_fin = t_fin, se_fin, big_samp = big_samp ) %>%
    mutate( sig = abs( t_fin / se_fin ) > 2 )

summary( res$se_fin )

ggplot( res[1:1000,], aes( t_fin, se_fin ) ) +
    geom_point()

mean( res$sig )

res |>
    group_by( big_samp ) |>
    summarize( sig = mean( sig ),
               est = mean( t_fin ) )

ggplot( res, aes( t_fin, col = big_samp, group=big_samp ) ) +
    geom_density()

