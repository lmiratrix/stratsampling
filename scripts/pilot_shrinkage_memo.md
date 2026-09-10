# Empirical-Bayes shrinkage in `estimate_pilot_variances()`

This memo derives the `shrink = TRUE` option of `estimate_pilot_variances()` (`R/pilot_functions.R`): a two-stage empirical-Bayes procedure that shrinks the per-cell (treatment arm × stratum) residual variances first, then uses those shrunk variances — not the raw ones — to shrink the per-cell mean residuals. Each stage is a standard estimator from a different literature; citations are given inline and flagged `(TO VERIFY)` since they are from memory and should be checked against the original sources before this memo is cited anywhere.

Throughout, a "cell" is one treatment-arm × stratum combination, indexed by `k` within a fixed arm (shrinkage is done separately within each arm, matching the existing pooled-variance calculation in the code). Let `e_ik = Y_ik - Yhat_ik` be the pilot residual for coded unit `i` in cell `k`, `n_k` the pilot count in cell `k`, `ēbar_k = mean(e_ik)` and `S²_k = var(e_ik)` the raw per-cell mean and variance.

## Why shrink at all

With few pilot units per cell, both `ēbar_k` and `S²_k` are noisy. Pooling completely (using one variance/mean for everything) throws away real between-cell heterogeneity; using each cell's raw estimate on its own overreacts to sampling noise, especially in small cells. Empirical-Bayes shrinkage — pulling each cell's estimate partway toward a pooled value, with the amount of pulling set by how little information that cell has — is the standard middle ground.

## Stage 2 first, conceptually: shrinking the means

Model the cell means as a two-level (random-effects / small-area) structure:

```
ēbar_k | mu_k ~ (ψ_k)           [sampling variance ψ_k = S²_k / n_k]
mu_k         ~ (sigma²_u)       [between-cell variance of the true cell means]
```

The posterior mean of `mu_k` (equivalently, the best linear unbiased predictor under this model) is the classic shrinkage estimator:

```
A_k    = ψ_k / (ψ_k + sigma²_u)
e*_k   = (1 - A_k)·ēbar_k + A_k·ēbar_w
```

where `ēbar_w` is the precision-weighted grand mean. This is the Fay-Herriot small-area estimator (Fay & Herriot, 1979) (TO VERIFY), and is algebraically the same shrinkage rule underlying the James-Stein estimator (James & Stein, 1961) (TO VERIFY) and Morris's parametric empirical-Bayes intervals (Morris, 1983) (TO VERIFY). Cells with little precision of their own (`ψ_k` large relative to `sigma²_u`) get pulled hard toward the pooled mean; well-estimated cells barely move.

### Estimating `sigma²_u`: DerSimonian-Laird method of moments

`sigma²_u` (the between-cell variance) is estimated by a closed-form, non-iterative method-of-moments formula from the random-effects meta-analysis literature (DerSimonian & Laird, 1986) (TO VERIFY):

```
w_k     = 1 / ψ_k
ēbar_w  = sum(w_k · ēbar_k) / sum(w_k)
Q       = sum(w_k · (ēbar_k - ēbar_w)^2)
sigma²_u = max(0, (Q - (K - 1)) / (sum(w_k) - sum(w_k²)/sum(w_k)))
```

with `K` the number of cells (within the arm). This is a one-shot computation — no numerical optimization — and clamps at 0 when the cells look no more variable than sampling noise alone predicts (the "fully pool" case).

`A_k` and `e*_k` are then computed as above; implemented in `.shrink_cell_means()`.

### Why not just `lmer(resid ~ 1 + (1|strata))`?

A random-intercept model is a valid alternative but assumes a single, pooled residual variance across all cells, so the implicit sampling variance it uses for shrinkage is `sigma²/n_k`, not the per-cell `S²_k/n_k`. That conflicts with the setting this package is built for, where residual variance is deliberately heterogeneous across strata (`generate_population(residual_var = "heterogeneous")`). Getting a mixed model to respect per-cell variances requires `nlme::lme(weights = varIdent(...))` — a new dependency for something the closed-form estimator above already gives directly, and in a way that's naturally heteroskedasticity-aware since `ψ_k` already varies by cell.

## Stage 1: shrinking the variances first

The concern that motivates doing this stage *before* stage 2: with a small `n_k`, the raw `S²_k` can come out too small by chance. If it does, `ψ_k = S²_k/n_k` looks small, `A_k` comes out small, and the mean-shrinkage step concludes (wrongly) that cell `k`'s own mean should be trusted and left unshrunk — exactly backwards, since a small `n_k` is precisely the situation where the cell mean is least trustworthy, regardless of what its noisy variance estimate says. This is the identical failure mode that motivates variance moderation in limma (Smyth, 2004) (TO VERIFY): genes with by-chance tiny variance estimates from few replicates produce artificially inflated t-statistics; the fix is to moderate the variance before it enters any downstream statistic.

### Model

```
S²_k | sigma²_k ~ sigma²_k · chisq(d_k) / d_k         [d_k = n_k - 1 pilot degrees of freedom]
sigma²_k        ~ Scaled-Inv-chisq(d_0, s_0²)          [conjugate prior]
```

`s_0²` is the pooled within-arm variance already computed by the existing code (`var_pooled`). The posterior mean is a degrees-of-freedom-weighted average of the prior and the cell's own estimate:

```
sigma²_k,shrunk = (d_0·s_0² + d_k·S²_k) / (d_0 + d_k)
```

Low-df (small-pilot) cells get pulled hard toward `s_0²`; well-estimated cells barely move — this is exactly the posterior-mean form used for moderated variances in limma (Smyth, 2004) (TO VERIFY).

### Estimating the prior degrees of freedom `d_0`

`d_0` is estimated by matching the sample variance of `log(S²_k)` across cells to its value under the model above. Using the standard chi-squared result `Var(log(chisq(d)/d)) ≈ trigamma(d/2)`:

```
Var(log(S²_k)) ≈ trigamma(d_k/2) + trigamma(d_0/2)
```

so, using the average `d_k` across cells as a quick approximation for the (generally slightly different) exact unequal-df formula (Smyth, 2004) (TO VERIFY):

```
v_log   = var(log(S²_k))                      # empirical, across cells
avg_tg  = mean(trigamma(d_k / 2))
target  = v_log - avg_tg
d_0     = if (target <= 0) Inf                # cells look no more variable than sampling
                                               # noise predicts -> pool completely
          else solve trigamma(d_0/2) = target for d_0   [trigamma is monotone decreasing, so
                                                          this has a unique solution, found here
                                                          by uniroot()]
```

`trigamma()` is a base-R function, and `uniroot()` is base `stats`, so this needs no new package dependency (in particular, no dependency on `limma` itself, which is a Bioconductor package (TO VERIFY) and not appropriate to add here).

Implemented in `.shrink_cell_variances()`.

## Putting the two stages together

`estimate_pilot_variances(..., shrink = TRUE)`:

1. Computes raw `mean_resid_pilot` and `var_resid_pilot` per cell, and `var_pooled` per arm, exactly as before (unchanged, and always reported raw in the output).
2. **Stage 1**: shrinks `var_resid_pilot` toward `var_pooled` via the df-weighted posterior mean above, using a data-driven `d_0` &rarr; `var.alloc`.
3. **Stage 2**: uses `var.alloc` (the *shrunk* variance from stage 1, not the raw one) to form `ψ_k = var.alloc/n_k`, and shrinks `mean_resid_pilot` toward the precision-weighted grand mean via the Fay-Herriot/DerSimonian-Laird estimator above &rarr; `mean.alloc`.

When `shrink = FALSE` (the default, and the prior behavior of this function), both steps are skipped and `var.alloc`/`mean.alloc` are just the raw per-cell values, with the pre-existing `min_pilot` fallback to `var_pooled` still applied to `var.alloc` for cells with too little pilot data.

## Caveat

This is a sequential plug-in procedure, not a single fully joint hierarchical model: it estimates the variance-shrinkage hyperparameters (`d_0`, `s_0²`) and treats the resulting shrunk variances as known when shrinking the means, rather than jointly estimating everything (and propagating the extra uncertainty from stage 1 into stage 2). This is the standard practical simplification used throughout applied empirical-Bayes work (it's exactly what limma's moderated-variance-then-moderated-t-statistic pipeline does), and it directly fixes the "artificially confident small cell" pathology motivating this whole exercise — it just means the resulting standard errors, if ever read off this pipeline directly, would be mildly optimistic. For this package's use (feeding `var.alloc` into `calc_optimal_allocation()` for sample-allocation design), that is an acceptable simplification; a design-stage variance input does not need to carry its own uncertainty quantification the way a final treatment-effect standard error would.

## References

- Fay, R. E., & Herriot, R. A. (1979). Estimates of income for small places: An application of James-Stein procedures to census data. *Journal of the American Statistical Association*. (TO VERIFY)
- James, W., & Stein, C. (1961). Estimation with quadratic loss. *Proceedings of the Fourth Berkeley Symposium on Mathematical Statistics and Probability*. (TO VERIFY)
- Morris, C. N. (1983). Parametric empirical Bayes inference: Theory and applications. *Journal of the American Statistical Association*. (TO VERIFY)
- DerSimonian, R., & Laird, N. (1986). Meta-analysis in clinical trials. *Controlled Clinical Trials*. (TO VERIFY)
- Smyth, G. K. (2004). Linear models and empirical Bayes methods for assessing differential expression in microarray experiments. *Statistical Applications in Genetics and Molecular Biology*. (TO VERIFY)
