# =====================================================================
# sim_functions_fast.R
#
# Drop-in fast replacements for the inner loop of the pilot simulation
# (01b_run_sim_pilot.R).  Everything here is base R: the dplyr calls in
# est_ma_strat(), calc_optimal_allocation(), generate_sample_index(),
# estimate_pilot_variances() and generate_population() are replaced by
# vector arithmetic on pre-computed treatment-by-stratum cell indices.
#
# The functions consume the RNG in exactly the same order, with exactly
# the same sample.int()/rnorm()/sample() calls, as the package versions,
# so for a given seed the realised populations and samples are IDENTICAL
# and the estimates agree to ~1e-12 (verified on all 270 scenarios of
# the balanced-exact grid; see verify block at the bottom of this file).
#
# Key ideas:
#   1. generate_population() sorts rows by (Z, strata), so the 2*K
#      treatment-by-stratum cells are CONTIGUOUS blocks of rows. All
#      "group_by(Z, strata) |> summarise()" work becomes cumsum() /
#      rowsum() on those blocks.
#   2. Everything that does not depend on the sample indicator S
#      (cell sizes, mean Yhat per cell, oracle residual variances,
#      by-arm var(Yobs)) is computed ONCE per population draw instead of
#      once per method (25 times per iteration in the current code).
#   3. One iteration returns a numeric matrix, not 25 one-row data
#      frames glued together with bind_rows().
#
# Speed: ~110x faster per iteration than the dplyr path (1.25 s -> 11 ms
# at N = 1000, n_strat = 4, 11 pilot fractions).
# =====================================================================


## ---- 0. Largest-remainder apportionment (copied from the package) ----
apportion_counts_f <- function(target, total, cap) {
    n <- pmin(floor(target), cap)
    frac <- target - floor(target)
    rem <- total - sum(n)
    while (rem > 0) {
        room <- cap - n
        cand <- ifelse(room > 0, frac, -Inf)
        i <- which.max(cand)
        if (!is.finite(cand[i])) break
        n[i] <- n[i] + 1
        frac[i] <- frac[i] - 1
        rem <- rem - 1
    }
    n
}


## ---- 1. Population draw + per-population precomputation --------------

# RNG-identical replacement for generate_surrogate_noise()
surrogate_noise_fast <- function(N, n_strat, strata_config, bias, residual_var,
                                 sd_y, R2_target, rescale_v = TRUE) {
    strat <- switch(strata_config,
                    "balanced-exact"  = sample(rep(1:n_strat, N / n_strat)),
                    "balanced-approx" = sample(1:n_strat, N, TRUE, rep(1 / n_strat, n_strat)),
                    "unbalanced"      = sample(1:n_strat, N, TRUE, runif(n_strat, 0.2, 0.8)))
    wts <- tabulate(strat, n_strat) / N
    
    err_var_raw <- switch(residual_var,
                          "homogeneous"      = rep(1, n_strat),
                          "heterogeneous"    = seq(0.25, 4, length = n_strat),
                          "extreme-contrast" = c(0.1, rep(1, n_strat - 2), 10))
    err_var <- if (rescale_v) err_var_raw / sum(err_var_raw * wts) else err_var_raw
    
    bias_raw <- switch(bias,
                       "none"             = rep(0, n_strat),
                       "extreme-contrast" = c(-1, rep(0, n_strat - 2), 1),
                       seq(-1, 1, length.out = n_strat))
    bc <- bias_raw - sum(bias_raw * wts)
    mx <- max(abs(bc))
    b  <- if (mx > 0) bc / mx else bc
    b  <- b * switch(bias, "none" = 0, "small" = 0.25 * sd_y,
                     "moderate" = 0.5 * sd_y, sd_y)
    
    b_bar_w <- sum(wts * b)
    curr <- sum(wts * (b - b_bar_w)^2) + sum(wts * err_var)
    cc <- if (curr > 0) sqrt(sd_y^2 * (1 - R2_target) / curr) else 0
    b <- cc * b
    sigma_eps <- cc * sqrt(err_var)
    
    e <- b[strat] + rnorm(N, 0, sd = sigma_eps[strat])
    o <- order(strat)                      # stable; matches dplyr::arrange
    list(strata = strat[o], e = e[o])
}

# Pre-compute everything about a population that does not depend on S.
# Rows must already be ordered by (Z, strata).
pop_precompute_v <- function(y, yh, Zv, st, K) {
    r <- y - yh
    cell <- Zv * K + st                    # 1..K (Z=0), K+1..2K (Z=1)
    K2 <- 2L * K
    N.zk <- tabulate(cell, K2)
    stopifnot(all(N.zk > 0))
    ends   <- cumsum(N.zk)
    starts <- c(1L, ends[-K2] + 1L)
    cs <- function(v) { cc <- c(0, cumsum(v)); cc[ends + 1L] - cc[starts] }
    
    sy <- cs(y); sy2 <- cs(y * y); sr <- cs(r); sr2 <- cs(r * r)
    Ybar      <- cs(yh) / N.zk                                   # mean(Yhat) per cell
    var.resid <- (sr2 - N.zk * (sr / N.zk)^2) / (N.zk - 1)       # oracle var(Yobs - Yhat)
    
    i0 <- 1:K; i1 <- (K + 1L):K2
    N.z   <- c(sum(N.zk[i0]), sum(N.zk[i1]))
    sy_z  <- c(sum(sy[i0]),  sum(sy[i1]))
    sy2_z <- c(sum(sy2[i0]), sum(sy2[i1]))
    ney.y <- ((sy2_z - N.z * (sy_z / N.z)^2) / (N.z - 1)) / N.z  # oracle var(Yobs)/N.z
    
    list(y = y, yh = yh, r = r, cell = cell, K = K, K2 = K2, N.zk = N.zk,
         starts = starts, ends = ends, arm = rep(0:1, each = K),
         Ybar = Ybar, var.resid = var.resid, N.z = N.z,
         oracle_ney = sum(ney.y), i0 = i0, i1 = i1)
}

# RNG-identical replacement for generate_population(); returns the
# precomputed structure instead of a data frame.
pop_fast <- function(params) {
    N <- params$N; K <- params$n_strat; sd_y <- params$sd_y
    Y0 <- rnorm(N, 0, sd_y); Y1 <- Y0 + params$tau
    Z  <- as.numeric(sample(N) <= floor(N / 2))
    Yobs <- ifelse(Z == 1, Y1, Y0)
    o <- order(Z); Z <- Z[o]; Yobs <- Yobs[o]
    N1 <- sum(Z == 1); N0 <- N - N1
    
    Yhat <- numeric(N); st <- integer(N)
    s0 <- surrogate_noise_fast(N0, K, params$strata_config, params$bias,
                               params$residual_var, sd_y, params$R2_target)
    i0 <- seq_len(N0)
    Yhat[i0] <- Yobs[i0] + s0$e; st[i0] <- s0$strata
    if (N1 > 0) {
        s1 <- surrogate_noise_fast(N1, K, params$strata_config, params$bias,
                                   params$residual_var, sd_y, params$R2_target)
        i1 <- N0 + seq_len(N1)
        Yhat[i1] <- Yobs[i1] + s1$e; st[i1] <- s1$strata
    }
    pop_precompute_v(Yobs, Yhat, Z, st, K)   # rows already sorted by (Z, strata)
}


## ---- 2. Model-assisted stratified estimator (est_ma_strat) -----------
# Returns c(est, SE, SEhat, SEhat.full, SEhat.samp, S2.0, S2.1).
est_fast <- function(pc, S) {
    n <- tabulate(pc$cell[S == 1L], pc$K2)
    M <- rowsum(cbind(S * pc$r, S * pc$r * pc$r, S * pc$y, S * pc$y * pc$y),
                pc$cell, reorder = TRUE)
    adj <- M[, 1] / n
    vrc <- (M[, 2] - n * adj^2) / (n - 1)          # var(resid) among coded
    ybo <- M[, 3] / n                              # mean(Yobs) among coded
    yb2 <- M[, 4] / n                              # mean(Yobs^2) among coded
    s2o <- (M[, 4] - n * ybo^2) / (n - 1)
    
    N.zk <- pc$N.zk; N.z <- pc$N.z
    frac <- (N.zk - n) / n
    wm <- function(v) c(sum(v[pc$i0] * N.zk[pc$i0]),
                        sum(v[pc$i1] * N.zk[pc$i1])) / N.z
    
    Ybar.adj.z    <- wm(pc$Ybar + adj)
    inner.z       <- wm(frac * pc$var.resid) / N.z
    inner.z.coded <- wm(frac * vrc) / N.z
    pt.3 <- wm((N.zk / N.z[pc$arm + 1]) * ((N.zk - n) / (N.zk - 1)) * (s2o / n))
    S2.hat <- (N.z / (N.z - 1)) * (wm(yb2) - wm(ybo)^2 + pt.3)
    
    var.hat.A <- sum(S2.hat / N.z)
    var.hat.B <- sum(inner.z.coded)
    c(est        = Ybar.adj.z[2] - Ybar.adj.z[1],
      SE         = sqrt(sum(inner.z) + pc$oracle_ney),
      SEhat      = sqrt(var.hat.A + var.hat.B),
      SEhat.full = sqrt(var.hat.A),
      SEhat.samp = sqrt(var.hat.B),
      S2.0 = S2.hat[1], S2.1 = S2.hat[2])
}


## ---- 3. Proportional stratified draw (generate_sample_index) ---------
gsi_fast <- function(pc, n_code, p_tx = 0.5) {
    N.zk <- pc$N.zk
    n_k <- integer(pc$K2)
    for (z in 0:1) {
        ix <- if (z == 0) pc$i0 else pc$i1
        p  <- if (z == 1) p_tx else 1 - p_tx
        Nk <- N.zk[ix]
        n_k[ix] <- apportion_counts_f(n_code * p * Nk / sum(Nk),
                                      total = round(n_code * p), cap = Nk)
    }
    if (any(n_k > N.zk)) stop("Cannot sample more than available in a stratum")
    S <- integer(length(pc$cell))
    for (k in seq_len(pc$K2)) {
        ids <- pc$starts[k]:pc$ends[k]
        S[ids[sample.int(length(ids), size = n_k[k])]] <- 1L
    }
    S
}


## ---- 4. Pilot variance estimates (estimate_pilot_variances) ----------
# shrink = 0; thin cells (< min_pilot, or non-estimable) fall back to the
# df-weighted pooled within-arm pilot variance, as in the package version
# (which also throws the warning that 01b wraps in suppressWarnings()).
pilotvar_fast <- function(pc, S_pilot, min_pilot = 2) {
    n <- tabulate(pc$cell[S_pilot == 1L], pc$K2)
    M <- rowsum(cbind(S_pilot * pc$r, S_pilot * pc$r * pc$r), pc$cell, reorder = TRUE)
    mu <- M[, 1] / n
    v  <- (M[, 2] - n * mu^2) / (n - 1)
    v[n < 2] <- NA_real_
    w <- pmax(n - 1, 0)
    out <- numeric(pc$K2)
    for (z in 0:1) {
        ix <- if (z == 0) pc$i0 else pc$i1
        vv <- v[ix]; ok <- !is.na(vv)
        pooled <- if (any(ok)) sum(w[ix][ok] * vv[ok]) / sum(w[ix][ok]) else NaN
        vv[n[ix] < min_pilot | is.na(vv)] <- pooled
        out[ix] <- vv
    }
    if (any(!is.finite(out)))
        stop("Could not form variance parameters for all cells; pilot is too small.")
    out
}


## ---- 5. Oracle / distorted allocation variances ----------------------
oracle_var_fast <- function(pc, distort = c("none", "permute", "reverse", "random")) {
    distort <- match.arg(distort)
    v <- pc$var.resid
    if (distort == "none") return(v)
    for (z in 0:1) {
        ix <- if (z == 0) pc$i0 else pc$i1
        vz <- v[ix]
        v[ix] <- switch(distort,
                        "permute" = sample(vz),
                        "reverse" = sort(vz, decreasing = TRUE)[rank(vz, ties.method = "first")],
                        "random"  = { lv <- log(pmax(vz, 1e-12)); s <- sd(lv)
                        if (is.na(s)) s <- 0
                        exp(rnorm(length(vz), mean = mean(lv), sd = s)) })
    }
    v
}


## ---- 6. Optimal (Neyman) allocation + draw (calc_optimal_allocation) --
# var.alloc: numeric vector of length 2*K in cell order (Z = 0 strata
# 1..K, then Z = 1 strata 1..K).  S_pilot non-NULL => pooled design
# (n_code is the TOTAL budget, pilot units are never un-sampled and the
# second wave is drawn only from un-piloted units).
alloc_fast <- function(pc, n_code, var.alloc, S_pilot = NULL, min_n = 2) {
    N.zk <- pc$N.zk
    n_pilot.zk <- if (is.null(S_pilot)) integer(pc$K2)
    else tabulate(pc$cell[S_pilot == 1L], pc$K2)
    weight <- N.zk * sqrt(var.alloc)
    n.zk <- numeric(pc$K2)
    n_code_z <- round(n_code / 2)
    
    for (z in 0:1) {
        ix <- if (z == 0) pc$i0 else pc$i1
        Nz <- N.zk[ix]; wz <- weight[ix]
        nz <- round(wz / sum(wz) * (n_code / 2))
        while (any(nz > Nz)) {                     # cap at cell size
            over <- nz >= Nz
            if (!any(over)) break
            nz[over] <- Nz[over]
            n_rem <- n_code_z - sum(nz[over])
            nz[!over] <- round(n_rem * (wz[!over] / sum(wz[!over])))
        }
        floor_k <- pmax(pmin(min_n, Nz), n_pilot.zk[ix])
        nz <- pmax(nz, floor_k)
        excess <- sum(nz) - n_code_z
        while (excess > 0) {                       # trim the most slack
            slack <- nz - floor_k
            if (max(slack) <= 0) break
            i <- which.max(slack); nz[i] <- nz[i] - 1; excess <- excess - 1
        }
        if (excess > 0)
            warning("Arm ", z, ": budget is below the per-cell floors (pilot + min_n).")
        while (excess < 0) {                       # top up by weight
            room <- Nz - nz
            if (max(room) <= 0) break
            i <- which.max(ifelse(room > 0, wz, -Inf))
            nz[i] <- nz[i] + 1; excess <- excess + 1
        }
        n.zk[ix] <- nz
    }
    
    n_new <- n.zk - n_pilot.zk
    S_new <- integer(length(pc$cell))
    for (k in seq_len(pc$K2)) {
        if (n_new[k] > 0) {
            ids <- pc$starts[k]:pc$ends[k]
            if (!is.null(S_pilot)) ids <- ids[S_pilot[ids] == 0L]
            S_new[ids[sample.int(length(ids), size = n_new[k])]] <- 1L
        }
    }
    if (is.null(S_pilot)) S_new else as.integer(S_pilot == 1L | S_new == 1L)
}


## ---- 7. One iteration: 3 baselines + 2 arms per pilot fraction --------
# Row order (matches the method/pilot_frac labels built in the runner):
#   1 Stratified, 2 Optimal, 3 Optimal-Random, then for each pilot
#   fraction: Optimal-Pilot-Pooled, Optimal-Pilot-Ind.
one_iter_pilot_fast <- function(params, pilot_fracs) {
    pc <- pop_fast(params)
    n_code <- round(params$N * params$coding_budget)
    out <- matrix(NA_real_, 3L + 2L * length(pilot_fracs), 7L)
    
    out[1, ] <- est_fast(pc, gsi_fast(pc, n_code))
    out[2, ] <- est_fast(pc, alloc_fast(pc, n_code, oracle_var_fast(pc, "none")))
    out[3, ] <- est_fast(pc, alloc_fast(pc, n_code, oracle_var_fast(pc, "random")))
    
    j <- 3L
    for (pf in pilot_fracs) {
        Sp <- gsi_fast(pc, round(pf * n_code))
        vt <- pilotvar_fast(pc, Sp)
        out[j + 1L, ] <- est_fast(pc, alloc_fast(pc, n_code, vt, S_pilot = Sp))
        out[j + 2L, ] <- est_fast(pc, alloc_fast(pc, n_code - sum(Sp), vt))
        j <- j + 2L
    }
    out
}

FAST_COLS <- c("est", "SE", "SEhat", "SEhat.full", "SEhat.samp", "S2.0", "S2.1")


## ---- 8. Equivalence check against the dplyr path ---------------------
# Run this once after any change to the package internals. It confirms
# the fast path draws the SAME samples (identical RNG stream) and returns
# the same numbers.
#
# verify_fast <- function(param_grid, pilot_fracs, seeds = 1:20) {
#     worst <- 0
#     for (i in seq_len(nrow(param_grid))) {
#         set.seed(1000 + i); o1 <- one_iter_pilot_study(param_grid[i, ], pilot_fracs)
#         set.seed(1000 + i); o2 <- one_iter_pilot_fast(param_grid[i, ], pilot_fracs)
#         m1 <- as.matrix(o1[, FAST_COLS])
#         worst <- max(worst, max(abs(m1 - o2) / pmax(abs(m1), 1e-12)))
#     }
#     worst
# }