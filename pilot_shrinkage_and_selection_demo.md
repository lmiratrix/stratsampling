Demo: pilot shrinkage, pilot-selected stratification, and estimate
pooling
================
2026-09-11

This is a working demo of three new pieces added to `stratsampling`,
adapted from the `if (FALSE) {}` sandbox blocks in
`R/pilot_functions.R`, `R/strata_eval_functions.R`, and
`scripts/sim_functions.R`. For the full technical derivation of the
shrinkage estimators, see `scripts/pilot_shrinkage_memo.md` – this
document only summarizes how the pieces fit together and shows them
running.

``` r
devtools::load_all( quiet = TRUE )   # picks up the in-development package code
library( tidyverse )
```

## A toy population with one real stratification and one useless covariate

Throughout, `strata` is the DGP’s real stratification (real
between-stratum bias/variance heterogeneity), and `X1` is pure noise
with no relationship to `Yobs`/`Yhat` – included so the demos below have
both a “good” and a “bad” candidate stratification to contrast.

``` r
set.seed( 42 )
dat <- generate_population( N = 1000, n_strat = 4,
                            strata_config = "unbalanced",
                            tau = 0.5, bias = "moderate",
                            residual_var = "heterogeneous",
                            R2_target = 0.6 )
dat$X1 <- rnorm( nrow(dat) )
```

We can then make a bunch of possible stratifications:

``` r
candidates <- create_stratification_candidates( dat, vars = c("Yhat", "X1", "strata"), Z = "Z" )
length( candidates )
```

    ## [1] 24

A small pilot (100 of the 1000 units) is coded; everyone else’s `Yobs`
is treated as not-yet-observed:

``` r
dat_full <- dat   # saving full data for sim stuff later on
dat$S <- as.numeric( sample( nrow(dat) ) <= 100 )
dat$Yobs[ dat$S == 0 ] <- NA
table( dat$S, dat$Z )
```

    ##    
    ##       0   1
    ##   0 447 453
    ##   1  53  47

# Part 1: Empirical-Bayes shrinkage of pilot means and variances

`estimate_pilot_variances(..., shrink = TRUE)` shrinks each cell’s raw
pilot mean and variance toward a pooled value, with the amount of
shrinkage set by how little pilot data that cell has (see the memo for
the full derivation: variances are shrunk first, fully pooled across
treatment arms; means are shrunk “hybrid” style – the amount of
shrinkage is based on all the strata across all the arms, but each cell
is still pulled toward its *own* arm’s grand mean).

How much shrinkage happens depends entirely on whether the candidate
stratification captures real signal. We compare three very different
candidates from the same pilot:

- **`strata_levels`** – the true stratification: real between-stratum
  heterogeneity, so shrinkage should be mild.
- **`X1_quintiles`** – cuts on the useless `X1`: any apparent
  between-cell spread here is pure noise, so shrinkage should be
  aggressive.
- **`Yhat_X1_2x2`** – crosses the meaningful `Yhat` with the useless
  `X1`: some real signal, diluted by noise and split into smaller cells,
  so shrinkage should land in between.

``` r
three_candidates <- c( "strata_levels", "X1_quintiles", "Yhat_X1_2x2" )

rr_three <- map_dfr( three_candidates, function( cname ) {
    estimate_pilot_variances( dat = NULL,
                              S = dat$S, Z = dat$Z,
                              strata = candidates[[ cname ]],
                              Yobs = dat$Yobs, Yhat = dat$Yhat,
                              shrink = TRUE ) %>%
        mutate( candidate = factor( cname, levels = three_candidates ) )
} )
```

Point size is scaled to each cell’s full-population stratum size
(`N.zk`), so it’s easy to see when a heavily-shrunk point belongs to a
small stratum:

``` r
ggplot( rr_three, aes( mean_resid_pilot, mean.alloc, col = as.factor(.Z), size = N.zk ) ) +
    facet_wrap( ~ candidate, nrow = 1 ) +
    geom_point() +
    scale_size_continuous( range = c(1, 6), name = "Stratum size (N.zk)" ) +
    geom_abline( intercept = 0, slope = 1, linetype = 2 ) +
    coord_fixed() + theme_minimal() +
    labs( title = "Mean shrinkage: raw vs shrunk cell means",
          x = "raw mean_resid_pilot", y = "shrunk mean.alloc", col = "Z" )
```

![](pilot_shrinkage_and_selection_demo_files/figure-gfm/shrink-mean-plot-1.png)<!-- -->

``` r
ggplot( rr_three, aes( var_resid_pilot, var.alloc, col = as.factor(.Z), size = N.zk ) ) +
    facet_wrap( ~ candidate, nrow = 1 ) +
    geom_point() +
    scale_size_continuous( range = c(1, 6), name = "Stratum size (N.zk)" ) +
    geom_abline( intercept = 0, slope = 1, linetype = 2 ) +
    coord_fixed() + theme_minimal() +
    labs( title = "Variance shrinkage: raw vs shrunk cell variances",
          x = "raw var_resid_pilot", y = "shrunk var.alloc", col = "Z" )
```

![](pilot_shrinkage_and_selection_demo_files/figure-gfm/shrink-var-plot-1.png)<!-- -->

Points on the dashed 45-degree line are unshrunk; points pulled off it
(toward a flat line, since all cells within a candidate then share close
to the same pooled value) show heavy shrinkage. `strata_levels` should
sit closest to the line, `X1_quintiles` should show the cells nearly
collapsing onto a single pooled value (both in the mean and, especially,
the variance plot), and `Yhat_X1_2x2` should fall in between – the
shrinkage machinery is reacting to how much real between-cell signal is
actually present, not just applying a fixed discount.

# Part 2: Using pilot data to choose a stratification

`compare_stratifications_precoding()` can now rank candidates two ways:

- **Proxy mode** (no `Yobs`/`S`): ranks by the spread of stratum-mean
  `Yhat` alone – no coding needed, but only a proxy for what actually
  matters (spread of the *residual*, `Yobs - Yhat`).
- **Pilot mode** (`Yobs` + `S` supplied): ranks by
  `mean_spread_pilot - avg_within_var_pilot`, a direct pilot estimate of
  the between-stratum residual spread net of within-stratum noise –
  computed via `estimate_pilot_variances()` internally, so it inherits
  the `shrink` argument.

``` r
cc_shrink   <- compare_stratifications_precoding( dat, candidates, Z = "Z", Yhat = "Yhat",
                                                  Yobs = dat$Yobs, S = dat$S, shrink = TRUE )
cc_noshrink <- compare_stratifications_precoding( dat, candidates, Z = "Z", Yhat = "Yhat",
                                                  Yobs = dat$Yobs, S = dat$S, shrink = FALSE )
cc_proxy    <- compare_stratifications_precoding( dat, candidates, Z = "Z", Yhat = "Yhat" )

tibble( rank = 1:5,
       pilot_shrink   = cc_shrink$stratification[1:5],
       pilot_noshrink = cc_noshrink$stratification[1:5],
       proxy_only     = cc_proxy$stratification[1:5] )
```

    ## # A tibble: 5 × 4
    ##    rank pilot_shrink         pilot_noshrink  proxy_only          
    ##   <int> <chr>                <chr>           <chr>               
    ## 1     1 Yhat_strata_3x3      Yhat_strata_3x3 Yhat_quintiles      
    ## 2     2 Yhat_strata_3x2      Yhat_strata_3x2 Yhat_quartiles      
    ## 3     3 Yhat_strata_2x3      Yhat_strata_2x3 Yhat_strata_3xlevels
    ## 4     4 Yhat_strata_2x2      Yhat_strata_2x2 Yhat_strata_3x3     
    ## 5     5 Yhat_strata_2xlevels strata_levels   Yhat_strata_3x2

With only 100 pilot units split across up to ~10 cells, the shrunk and
un-shrunk pilot rankings mostly agree here (both put the `Yhat_strata`
interactions on top), but they need not: with a small pilot, a candidate
with a couple of thin cells can get an artificially small or large raw
variance by chance and move up or down the un-shrunk ranking in a way
shrinkage would have damped out. The proxy-only ranking, having no
access to `Yobs` at all, is a different (weaker) signal again – notice
it never gets to see that the `Yhat_strata` interactions outperform
plain `Yhat` cuts, since that requires the pilot’s actual residuals.
Re-knit this document with a different `set.seed()` above, or a smaller
pilot, to see the shrunk/un-shrunk rankings diverge more.

# Part 3: Pilot-selected stratification, with precision-weighted pooling

`one_iter_pilot_select()` puts this together as a two-wave design:

1.  Draw a small SRS pilot; use `compare_stratifications_precoding()`
    (pilot mode, shrinkage on by default) to pick the best candidate
    stratification.
2.  Draw a wave-2 sample of the remaining budget from the *full*
    population (not the un-piloted remainder – see the comment above the
    function for why), stratified on the chosen candidate, both
    proportionally and with pilot-informed optimal (Neyman) allocation.
3.  Combine the pilot-wave estimate and each wave-2 estimate with
    `combine_estimators()` – an inverse-variance (“precision”) weighted
    pooling of two estimates of the same quantity (see its doc string
    for the exact formulas, and the note there on the independence
    assumption).

**Note:** `one_iter_pilot_select()` draws its own pilot and wave-2
samples internally, so (unlike Parts 1-2 above, which worked from a
single small pre-drawn pilot with `Yobs` masked outside it) it needs
everyone’s real `Yobs` available – we pass it `dat_full`, the unmasked
snapshot saved above, rather than the `dat` object from Part 1.

``` r
source( "scripts/sim_functions.R" )
set.seed( 2026 )
result <- one_iter_pilot_select( dat_full, n_code = 200, candidates = candidates,
                                 pilot_select_frac = 0.05,
                                 Z = "Z", Yhat = "Yhat", Yobs = "Yobs" )
result %>% select( method, est, SEhat, SE, strat_choice, n_pilot )
```

    ##                         method       est     SEhat        SE  strat_choice
    ## 1       Stratified-PilotSelect 0.3654152 0.2530734 0.2716972 strata_levels
    ## 2          Optimal-PilotSelect 0.6817672 0.2658350 0.2666485 strata_levels
    ## 3 Pilot+Stratified-PilotSelect 0.5453291 0.2272252 0.2440313 strata_levels
    ## 4    Pilot+Optimal-PilotSelect 0.8100495 0.2363269 0.2403534 strata_levels
    ##   n_pilot
    ## 1      50
    ## 2      50
    ## 3      50
    ## 4      50

Reading the four rows:

- **`Stratified-PilotSelect`** / **`Optimal-PilotSelect`** – the
  wave-2-only estimates (proportional vs. Neyman-optimal allocation) on
  the pilot-chosen stratum (`strat_choice`).
- **`Pilot+Stratified-PilotSelect`** / **`Pilot+Optimal-PilotSelect`** –
  each of those, precision-pooled with the pilot-wave (SRS) estimate via
  `combine_estimators()`. Since the pilot here is used only to *choose*
  the stratification and is otherwise excluded from the wave-2 draw, the
  two waves being combined are close to independent draws, which is what
  `combine_estimators()`’s pooling formula assumes.

The pooled rows should generally have smaller `SEhat` than either wave
alone – that’s the whole point of pooling two (roughly) independent
estimates of the same effect – though see the comment block above
`one_iter_pilot_select()` in `scripts/sim_functions.R` for the
acknowledged caveat: if a wave-2 draw happens to recode a unit already
in the pilot, the two waves are not *exactly* independent, so the pooled
SE here is a mild approximation.

# Further Information

- `scripts/pilot_shrinkage_memo.md` – full derivation of the two-stage
  empirical-Bayes shrinkage (variances first, fully pooled across arms;
  means “hybrid”, shared amount-of-shrinkage but arm-specific target),
  with citations.
- `R/pilot_functions.R` – `estimate_pilot_variances()` and the shrinkage
  internals; its own `if (FALSE) {}` sandbox has the source for Part 1
  above.
- `R/strata_eval_functions.R` – `evaluate_stratification_precoding()` /
  `compare_stratifications_precoding()`; its sandbox has the source for
  Part 2 above.
- `scripts/sim_functions.R` – `one_iter_pilot_select()` and
  `combine_estimators()` (see its doc string for the pooling formulas);
  its sandbox has the source for Part 3 above.
