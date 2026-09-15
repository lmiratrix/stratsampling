# =====================================================================
# sim_functions_fast.R
#
# Drop-in fast replacements for the inner loops of the pilot simulation
# (01b_run_sim_pilot.R -> one_iter_pilot_fast) and of the main factorial
# simulation (01_run_full_sim.R -> one_iter_all_fast, oracle_fast).
# Everything here is base R: the dplyr calls in est_ma_strat(),
# est_subset(), calc_optimal_allocation(), generate_sample_index(),
# estimate_pilot_variances(), create_stratification_candidates(),
# compare_stratifications_precoding() and generate_population() are
# replaced by vector arithmetic on pre-computed treatment-by-stratum
# cell indices.
#
# The functions consume the RNG in exactly the same order, with exactly
# the same sample.int()/rnorm()/sample() calls, as the package versions,
# so for a given seed the realised populations and samples are IDENTICAL
# and the estimates agree to ~1e-12 (see verify_fast() and
# verify_all_fast() at the bottom of this file; rerun them after ANY
# change to the package internals or to sim_functions.R).  Verified
# 2026-09-14/15 on all 270 scenarios x 11 shares (01b) and on all 270
# scenarios x 3 strata configurations with pilot_frac = 1/3 (01).
#
# Key ideas:
#   1. generate_population() sorts rows by (Z, strata), so the 2*K
#      treatment-by-stratum cells of the TRUE strata are contiguous
#      blocks of rows.  Pilot-SELECTED strata (quantile cuts of Yhat)
#      are not contiguous, so every cell structure carries an explicit
#      list of row indices per cell (`idx`); the contiguous case is just
#      the special case idx[[k]] = starts[k]:ends[k].
#   2. Everything that does not depend on a sample indicator (cell
#      sizes, mean Yhat per cell, oracle residual variances, by-arm
#      var(Yobs), the candidate stratifications and their proxy
#      diagnostic) is computed ONCE per population draw instead of once
#      per method.
#   3. One iteration returns a numeric matrix, not ~60 one-row data
#      frames glued together with bind_rows().
#
# Method registry (row order within one iteration):
#   one_iter_pilot_fast (01b): 1 Stratified, 2 Optimal, 3 Optimal-Random,
#     4 Optimal-Permuted, 5 Optimal-Reversed, then for each pilot fraction:
#     Optimal-Pilot-Pooled, Optimal-Pilot-Ind,
#     Stratified-PilotSelect, Optimal-PilotSelect,
#     Pilot+Stratified-PilotSelect, Pilot+Optimal-PilotSelect
#   one_iter_all_fast (01): Subset, SRS, Stratified, Optimal,
#     Stratified-Random, Optimal-Permuted, Optimal-Reversed,
#     [Optimal-Random], then the same six pilot arms at one share.
#   (sim_functions.R of 2026-09-15: one_iter_pilot_select() returns four
#   rows -- proportional and Neyman wave 2 on the pilot-selected strata,
#   Neyman now always on EB-SHRUNK pilot variances, plus each of them
#   precision-combined with the SRS pilot estimate via
#   combine_estimators(); the raw-variance arm is gone.)
#
# Columns: est, SE, SEhat, SEhat.full, SEhat.samp, S2.0, S2.1,
#          strat_choice (integer code into CAND_NAMES; NA for the
#          fixed-strata arms), n_pilot (realised pilot size; NA for the
#          four baseline arms).
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

# Build the cell structure used by est_fast / gsi_fast / alloc_fast /
# pilotvar_fast for ANY stratification with K levels, given per-row
# integer codes `st` in 1..K (level order = the order dplyr sorts the
# factor, i.e. factor-level order) and arm indicator Zv in {0, 1}.
# Rows are assumed sorted by Z (arm 0 first), as generate_population()
# guarantees.  Cell k = Z*K + st, so cells are ordered (Z, strata)
# exactly like the package's group_by(.Z, .strata) / arrange(.Z, .strata).
# Every (arm, level) cell must be non-empty (true for the true strata
# and for within-arm quantile cuts).
cells_fast <- function(y, yh, Zv, st, K, contiguous = FALSE) {
    r <- y - yh
    cell <- Zv * K + st
    K2 <- 2L * K
    N.zk <- tabulate(cell, K2)
    stopifnot(all(N.zk > 0))
    
    if (contiguous) {
        ends   <- cumsum(N.zk)
        starts <- c(1L, ends[-K2] + 1L)
        idx    <- lapply(seq_len(K2), function(k) starts[k]:ends[k])
    } else {
        idx <- split(seq_along(cell), factor(cell, levels = seq_len(K2)))
        names(idx) <- NULL
    }
    
    M  <- rowsum(cbind(y, y * y, r, r * r, yh), cell, reorder = TRUE)
    sy <- M[, 1]; sy2 <- M[, 2]; sr <- M[, 3]; sr2 <- M[, 4]
    Ybar      <- M[, 5] / N.zk                                   # mean(Yhat) per cell
    var.resid <- (sr2 - N.zk * (sr / N.zk)^2) / (N.zk - 1)       # oracle var(Yobs - Yhat)
    
    i0 <- 1:K; i1 <- (K + 1L):K2
    N.z   <- c(sum(N.zk[i0]), sum(N.zk[i1]))
    sy_z  <- c(sum(sy[i0]),  sum(sy[i1]))
    sy2_z <- c(sum(sy2[i0]), sum(sy2[i1]))
    ney.y <- ((sy2_z - N.z * (sy_z / N.z)^2) / (N.z - 1)) / N.z  # oracle var(Yobs)/N.z
    
    list(y = y, yh = yh, r = r, cell = cell, st = st, Zv = Zv, K = K, K2 = K2, N.zk = N.zk,
         idx = idx, arm = rep(0:1, each = K),
         arm_rows = list(seq_len(N.z[1]), N.z[1] + seq_len(N.z[2])),
         Ybar = Ybar, var.resid = var.resid, N.z = N.z,
         var.y.z = ney.y * N.z,                 # var(Yobs) within arm
         oracle_ney = sum(ney.y), i0 = i0, i1 = i1)
}

# Pre-compute everything about a population that does not depend on S.
# Rows must already be ordered by (Z, strata).
pop_precompute_v <- function(y, yh, Zv, st, K) {
    cells_fast(y, yh, Zv, st, K, contiguous = TRUE)
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
        ids <- pc$idx[[k]]
        S[ids[sample.int(length(ids), size = n_k[k])]] <- 1L
    }
    S
}

# SRS within arm (generate_sample_index with strata = NULL): one cell
# per arm, arm 0 drawn first.  Used for the selection pilot.
srs_fast <- function(pc, n_code, p_tx = 0.5) {
    S <- integer(length(pc$cell))
    for (z in 0:1) {
        ids <- pc$arm_rows[[z + 1L]]
        p   <- if (z == 1) p_tx else 1 - p_tx
        n_z <- apportion_counts_f(n_code * p, total = round(n_code * p),
                                  cap = length(ids))
        if (n_z > length(ids)) stop("Cannot sample more than available in a stratum")
        S[ids[sample.int(length(ids), size = n_z)]] <- 1L
    }
    S
}


## ---- 4. Pilot variance estimates (estimate_pilot_variances) ----------
# Raw per-cell pilot residual mean / variance / count, plus the
# df-weighted pooled within-arm pilot variance (var_pooled).  Shared by
# the no-shrink and shrink paths below.
pilot_cells_fast <- function(pc, S_pilot) {
    n <- tabulate(pc$cell[S_pilot == 1L], pc$K2)
    M <- rowsum(cbind(S_pilot * pc$r, S_pilot * pc$r * pc$r), pc$cell, reorder = TRUE)
    mu <- M[, 1] / n
    v  <- (M[, 2] - n * mu^2) / (n - 1)
    v[n < 2] <- NA_real_
    w <- pmax(n - 1, 0)
    pooled <- numeric(2)
    for (z in 0:1) {
        ix <- if (z == 0) pc$i0 else pc$i1
        vv <- v[ix]; ok <- !is.na(vv)
        pooled[z + 1L] <- if (any(ok)) sum(w[ix][ok] * vv[ok]) / sum(w[ix][ok]) else NaN
    }
    list(n = n, mu = mu, v = v, pooled = pooled)
}

# shrink = FALSE; thin cells (< min_pilot, or non-estimable) fall back to
# the df-weighted pooled within-arm pilot variance, as in the package
# version (which also throws the warning that 01b wraps in
# suppressWarnings()).
pilotvar_fast <- function(pc, S_pilot, min_pilot = 2) {
    p <- pilot_cells_fast(pc, S_pilot)
    out <- numeric(pc$K2)
    for (z in 0:1) {
        ix <- if (z == 0) pc$i0 else pc$i1
        vv <- p$v[ix]
        vv[p$n[ix] < min_pilot | is.na(vv)] <- p$pooled[z + 1L]
        out[ix] <- vv
    }
    if (any(!is.finite(out)))
        stop("Could not form variance parameters for all cells; pilot is too small.")
    out
}

# Verbatim copy of stratsampling:::.shrink_cell_variances (R/pilot_functions.R).
# Re-copy if the package version changes.  NOTE: as written, the root of
# trigamma(x/2) = target is doubled (d0 = 2 * root), so the prior df is
# twice the value the memo derives; whatever is decided about that line
# must be changed HERE and in the package together.
shrink_cell_variances_f <- function( s2, n, s0 ) {
    
    d <- pmax( n - 1, 0 )
    ok <- d >= 1 & is.finite( s2 ) & s2 > 0
    
    if ( sum(ok) < 2 ) {
        # Not enough cells to estimate a prior df; pool completely.
        return( s0 )
    }
    
    v_log  <- var( log( s2[ok] ) )
    avg_tg <- mean( trigamma( d[ok] / 2 ) )
    target <- v_log - avg_tg
    
    d0 <- if ( !is.finite(target) || target <= 0 ) {
        Inf
    } else {
        root <- tryCatch(
            uniroot( function(x) trigamma(x / 2) - target, c(1e-4, 1e6) )$root,
            error = function(e) NA_real_ )
        if ( is.na(root) ) Inf else 2 * root
    }
    
    if ( !is.finite(d0) ) {
        return( s0 )
    }
    
    d_use  <- ifelse( d <= 0 | !is.finite(s2), 0, d )
    s2_use <- ifelse( d <= 0 | !is.finite(s2), 0, s2 )
    
    ( d0 * s0 + d_use * s2_use ) / ( d0 + d_use )
}

# shrink = TRUE: stage-1 empirical-Bayes shrinkage of the raw per-cell
# variances toward the pooled within-arm variance, applied within arm.
# (The package also shrinks the cell MEANS in stage 2, but
# calc_optimal_allocation() only consumes var.alloc, so mean.alloc is
# not needed here; it involves no RNG.)  The min_pilot fallback is NOT
# applied on this path, matching the package.
pilotvar_shrink_fast <- function(pc, S_pilot) {
    p <- pilot_cells_fast(pc, S_pilot)
    out <- numeric(pc$K2)
    for (z in 0:1) {
        ix <- if (z == 0) pc$i0 else pc$i1
        out[ix] <- shrink_cell_variances_f(p$v[ix], p$n[ix],
                                           rep(p$pooled[z + 1L], length(ix)))
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
            ids <- pc$idx[[k]]
            if (!is.null(S_pilot)) ids <- ids[S_pilot[ids] == 0L]
            S_new[ids[sample.int(length(ids), size = n_new[k])]] <- 1L
        }
    }
    if (is.null(S_pilot)) S_new else as.integer(S_pilot == 1L | S_new == 1L)
}


## ---- 7. Candidate stratifications of Yhat (create_stratification_candidates)
# Mirrors create_stratification_candidates(dat, vars = "Yhat", Z = "Z",
# verbose = FALSE) with its DEFAULTS (cut_by_arm = TRUE,
# n_groups = c(3, 4, 5), min_cell_size = 2, no breaks/custom/jitter).
# With a single continuous variable there are no interactions, so the
# candidate list is Yhat_tertiles, Yhat_quartiles, Yhat_quintiles, each a
# within-arm quantile cut (quantile() type 7, cut(include.lowest = TRUE)).
# No RNG is consumed by the package version either (jitter_vars = NULL,
# seed = NULL), so this is computed ONCE per population.
#
# If the reference call ever changes (other vars, breaks, interaction
# configs, cut_by_arm = FALSE), this function must be extended to match;
# it stops if Yhat looks categorical to the package (<= 30 distinct values).
CAND_N_GROUPS <- c(3L, 4L, 5L)
CAND_NAMES    <- paste0("Yhat_", c("tertiles", "quartiles", "quintiles"))

quantile_codes_by_arm_f <- function(pc, k) {
    x <- pc$yh
    st <- integer(length(x))
    probs <- seq(0, 1, length.out = k + 1)
    for (z in 0:1) {
        idx <- pc$arm_rows[[z + 1L]]
        breaks <- quantile(x[idx], probs = probs, na.rm = TRUE)
        if (length(unique(breaks)) < k + 1) return(NULL)
        st[idx] <- cut(x[idx], breaks = breaks, labels = FALSE, include.lowest = TRUE)
    }
    st
}

candidates_fast <- function(pc, n_groups = CAND_N_GROUPS, min_cell_size = 2) {
    if (length(unique(pc$yh)) <= 30)
        stop("Yhat has <= 30 distinct values; the package would treat it as ",
             "categorical and candidates_fast() does not mirror that path.")
    keep <- logical(length(n_groups))
    cells <- vector("list", length(n_groups))
    spread <- rep(NA_real_, length(n_groups))
    for (j in seq_along(n_groups)) {
        k  <- n_groups[j]
        st <- quantile_codes_by_arm_f(pc, k)
        if (is.null(st)) next                                   # non-unique breaks
        pcs <- cells_fast(pc$y, pc$yh, pc$arm[pc$cell], st, k)
        if (min(pcs$N.zk) < min_cell_size) next                 # screen_candidates()
        keep[j] <- TRUE
        cells[[j]] <- pcs
        # Proxy diagnostic (evaluate_stratification_precoding, proxy mode):
        # mean over arms of sum_k w_k (mean_Yhat_k - Yhat_bar_z)^2 with
        # w_k = N_zk / N_z and Yhat_bar_z = weighted.mean(mean_Yhat, w_k).
        vm <- numeric(2)
        for (z in 0:1) {
            ix <- if (z == 0) pcs$i0 else pcs$i1
            w_k <- pcs$N.zk[ix] / pcs$N.z[z + 1L]
            m_k <- pcs$Ybar[ix]
            vm[z + 1L] <- sum(w_k * (m_k - weighted.mean(m_k, w_k))^2)
        }
        spread[j] <- mean(vm)
    }
    list(keep = keep, cells = cells, mean_spread = spread,
         codes = seq_along(n_groups), n_groups = n_groups)
}


## ---- 8. Pilot-mode selection (compare_stratifications_precoding) -----
# Ranks the kept candidates by mean_spread_pilot = mean over arms of
# sum_k w_k (mean_resid_k - resid_bar_z)^2, where mean_resid_k is the
# mean pilot residual in cell k (NaN if the cell has no pilot units,
# which makes the whole candidate NA), w_k = N_zk / N_z and resid_bar_z
# = weighted.mean(mean_resid, w_k).  NAs sort last (dplyr::arrange).  If
# the top candidate is NA (no candidate is pilot-evaluable) the package
# falls back to which.max(mean_spread) over the (stably sorted, i.e.
# original-order) table.  Optional min_cell screen on min_stratum_size
# as in one_iter_pilot_select(); it is skipped if it would empty the table.
# Returns the integer code (position in CAND_NAMES) of the winner.
select_fast <- function(cands, S_pilot, min_cell = 0) {
    codes <- which(cands$keep)
    if (length(codes) == 0) stop("No valid stratifications found")
    ps <- rep(NA_real_, length(codes))
    for (j in seq_along(codes)) {
        pcs <- cands$cells[[codes[j]]]
        n <- tabulate(pcs$cell[S_pilot == 1L], pcs$K2)
        if (any(n == 0)) next                                   # NaN in the package
        mu <- rowsum(S_pilot * pcs$r, pcs$cell, reorder = TRUE)[, 1] / n
        vm <- numeric(2)
        for (z in 0:1) {
            ix <- if (z == 0) pcs$i0 else pcs$i1
            w_k <- pcs$N.zk[ix] / pcs$N.z[z + 1L]
            vm[z + 1L] <- sum(w_k * (mu[ix] - weighted.mean(mu[ix], w_k))^2)
        }
        ps[j] <- mean(vm)
    }
    if (min_cell > 0) {
        mss <- vapply(codes, function(cd) min(cands$cells[[cd]]$N.zk), numeric(1))
        if (any(mss >= min_cell)) { codes <- codes[mss >= min_cell]; ps <- ps[mss >= min_cell] }
    }
    o <- order(-ps, na.last = TRUE)                             # stable, like arrange(desc())
    best <- codes[o[1]]
    if (is.na(ps[o[1]])) {
        best <- codes[which.max(cands$mean_spread[codes])]     # proxy fallback
    }
    best
}


## ---- 9. Pilot-SELECTED stratification arm (one_iter_pilot_select) ----
# Mirrors sim_functions.R (2026-09-15): SRS pilot (within arm) -> choose
# the stratification -> wave-1 estimate from the pilot alone (est_ma_srs,
# i.e. one stratum per arm) -> two fresh wave-2 draws of
# n_main = n_code - n_pilot from the FULL population on the selected
# strata: proportional, then Neyman on EB-SHRUNK pilot variances -> each
# wave-2 estimate precision-combined with the wave-1 estimate.
# RNG order: srs_fast, gsi_fast, alloc_fast (identical to the package
# path; estimate_pilot_variances / est_ma_srs / combine consume no RNG).
#
# Rows (matches bind_rows(out.prop, out.opt, out.prop.pool, out.opt.pool)):
#   1 Stratified-PilotSelect        2 Optimal-PilotSelect
#   3 Pilot+Stratified-PilotSelect  4 Pilot+Optimal-PilotSelect
#
# If the pilot cannot yield shrunk variance parameters for every cell
# (an arm whose cells all hold <= 1 pilot unit; only possible at the very
# smallest pilot, e.g. n_code = 100 and share 0.10) the package function
# stop()s AFTER the proportional wave-2 draw.  Here rows 1 and 3 are
# still returned, rows 2 and 4 are NA and no allocation draw is made,
# so the RNG stream stays aligned with a reference that wraps the
# estimate_pilot_variances() call in tryCatch() at the same point.
#
# combine_estimators(): inverse-variance weights from SEhat^2 for est,
# SEhat, SEhat.full, SEhat.samp, S2.0, S2.1; separate weights from SE^2
# for the oracle SE.  Vectors are (est, SE, SEhat, SEhat.full,
# SEhat.samp, S2.0, S2.1).
combine_fast <- function(e1, e2) {
    var1 <- e1[3]^2; var2 <- e2[3]^2
    prec <- var2 / (var1 + var2)
    out <- e1
    out[1] <- prec * e1[1] + (1 - prec) * e2[1]
    out[3] <- sqrt(prec^2 * var1 + (1 - prec)^2 * var2)
    out[4] <- sqrt(prec^2 * e1[4]^2 + (1 - prec)^2 * e2[4]^2)
    out[5] <- sqrt(prec^2 * e1[5]^2 + (1 - prec)^2 * e2[5]^2)
    out[6] <- prec * e1[6] + (1 - prec) * e2[6]
    out[7] <- prec * e1[7] + (1 - prec) * e2[7]
    var1T <- e1[2]^2; var2T <- e2[2]^2
    precT <- var2T / (var1T + var2T)
    out[2] <- sqrt(precT^2 * var1T + (1 - precT)^2 * var2T)
    out
}

# pc1: the one-stratum-per-arm cell structure of the same population
# (one_stratum_fast(pc)); pass it in to avoid rebuilding it per share.
pilot_select_fast <- function(pc, cands, n_code, pilot_select_frac,
                              min_cell = 0, pc1 = one_stratum_fast(pc)) {
    n_pilot <- round(pilot_select_frac * n_code / 2) * 2
    n_main  <- n_code - n_pilot
    stopifnot(n_pilot >= 1, n_main >= 1)
    
    Sp   <- srs_fast(pc, n_pilot)                               # wave 1
    best <- select_fast(cands, Sp, min_cell = min_cell)
    pcs  <- cands$cells[[best]]
    
    out <- matrix(NA_real_, 4L, 9L)
    out[, 8L] <- best
    out[, 9L] <- n_pilot
    
    wave1 <- est_fast(pc1, Sp)                                  # est_ma_srs on the pilot
    
    prop <- est_fast(pcs, gsi_fast(pcs, n_main))                # (a) proportional wave 2
    out[1L, 1:7] <- prop
    out[3L, 1:7] <- combine_fast(wave1, prop)
    
    vt <- tryCatch(pilotvar_shrink_fast(pcs, Sp), error = function(e) NULL)
    if (is.null(vt)) return(out)                                # pilot too thin
    opt <- est_fast(pcs, alloc_fast(pcs, n_main, vt))           # (b) Neyman wave 2, shrunk variances
    out[2L, 1:7] <- opt
    out[4L, 1:7] <- combine_fast(wave1, opt)
    out
}


## ---- 10. One iteration: 5 baselines + 6 arms per pilot fraction -------
# Row order (matches the method/pilot_frac labels built in the runner):
#   1 Stratified, 2 Optimal, 3 Optimal-Random, 4 Optimal-Permuted,
#   5 Optimal-Reversed, then for each pilot fraction: Optimal-Pilot-Pooled,
#   Optimal-Pilot-Ind, Stratified-PilotSelect, Optimal-PilotSelect,
#   Pilot+Stratified-PilotSelect, Pilot+Optimal-PilotSelect.
# The reference (one_iter_pilot_study in 01b_run_sim_pilot.R) must
# compute the baselines in this order; "reverse" consumes no RNG in the
# distortion itself but its allocation draw does.
#
# paired_pilot: FALSE calls one_iter_optimal_pilot() twice per fraction and so draws
#   a FRESH proportional pilot for the Ind arm.  TRUE  mirrors the current 01b_run_sim_pilot.R,
#   (one pilot shared by Pooled and Ind).  Must match the
#   reference or verify_fast() fails from the Ind arm onward.
one_iter_pilot_fast <- function(params, pilot_fracs, paired_pilot = FALSE) {
    pc <- pop_fast(params)
    cands <- candidates_fast(pc)                                # no RNG
    pc1 <- one_stratum_fast(pc)                                 # for the wave-1 (SRS pilot) estimate
    n_code <- round(params$N * params$coding_budget)
    out <- matrix(NA_real_, 5L + 6L * length(pilot_fracs), 9L)
    
    out[1, 1:7] <- est_fast(pc, gsi_fast(pc, n_code))
    out[2, 1:7] <- est_fast(pc, alloc_fast(pc, n_code, oracle_var_fast(pc, "none")))
    out[3, 1:7] <- est_fast(pc, alloc_fast(pc, n_code, oracle_var_fast(pc, "random")))
    out[4, 1:7] <- est_fast(pc, alloc_fast(pc, n_code, oracle_var_fast(pc, "permute")))
    out[5, 1:7] <- est_fast(pc, alloc_fast(pc, n_code, oracle_var_fast(pc, "reverse")))
    
    j <- 5L
    for (pf in pilot_fracs) {
        Sp <- gsi_fast(pc, round(pf * n_code))
        vt <- pilotvar_fast(pc, Sp)
        out[j + 1L, 1:7] <- est_fast(pc, alloc_fast(pc, n_code, vt, S_pilot = Sp))
        out[j + 1L, 9L]  <- sum(Sp)
        if (!paired_pilot) {
            Sp <- gsi_fast(pc, round(pf * n_code))
            vt <- pilotvar_fast(pc, Sp)
        }
        out[j + 2L, 1:7] <- est_fast(pc, alloc_fast(pc, n_code - sum(Sp), vt))
        out[j + 2L, 9L]  <- sum(Sp)
        out[j + 3:6, ]   <- pilot_select_fast(pc, cands, n_code, pf, pc1 = pc1)
        j <- j + 6L
    }
    out
}

## ---- 10b. Main-grid arms (01_run_full_sim.R / one_iter_all) ----------
# The four arms one_iter_all() has that the pilot study does not:
#   Subset            SRS within arm, difference of coded means (est_subset)
#   SRS               SRS within arm, model-assisted estimator with one
#                     stratum per arm (est_ma_srs = est_ma_strat(strata = NULL))
#   Stratified-Random proportional allocation on strata labels permuted
#                     within arm (one_iter_random_strat)
#   Optimal-Reversed  Neyman allocation with rank-reversed variances

# dplyr::sample_n(n_code / 2) on group_by(Z): one sample.int() per arm,
# arm 0 first (slice() evaluates per group in key order).
srs_draw_fast <- function(pc, n_code) {
    S <- integer(length(pc$cell))
    for (z in 0:1) {
        ids <- pc$arm_rows[[z + 1L]]
        S[ids[sample.int(length(ids), size = n_code / 2)]] <- 1L
    }
    S
}

# est_subset(): difference of coded means; SE uses the oracle var(Yobs)
# per arm, SEhat the coded-sample variance; SEhat.full / SEhat.samp are NA.
subset_est_fast <- function(pc, S) {
    ybar <- vy <- vyc <- nz <- numeric(2)
    for (z in 0:1) {
        ids <- pc$arm_rows[[z + 1L]]
        yc  <- pc$y[ids][S[ids] == 1L]
        nz[z + 1L]   <- length(yc)
        ybar[z + 1L] <- mean(yc)
        vyc[z + 1L]  <- var(yc)
    }
    c(est = ybar[2] - ybar[1],
      SE = sqrt(sum(pc$var.y.z / nz)),
      SEhat = sqrt(sum(vyc / nz)),
      SEhat.full = NA_real_, SEhat.samp = NA_real_,
      S2.0 = vyc[1], S2.1 = vyc[2])
}

# One stratum per arm: the cell structure behind est_ma_srs()
one_stratum_fast <- function(pc) {
    cells_fast(pc$y, pc$yh, pc$Zv, rep(1L, length(pc$cell)), 1L, contiguous = TRUE)
}

# group_by(Z) |> mutate(strata_random = sample(strata)): one permutation
# per arm, arm 0 first; returns the cell structure on the permuted labels
# (non-contiguous cells, same cell sizes as the true strata).
random_strata_fast <- function(pc) {
    st <- integer(length(pc$cell))
    for (z in 0:1) {
        ids <- pc$arm_rows[[z + 1L]]
        st[ids] <- pc$st[ids][sample.int(length(ids))]
    }
    cells_fast(pc$y, pc$yh, pc$Zv, st, pc$K)
}

# RNG-identical replacement for one_iter_all(dat, n_code, include_rand,
# pilot_frac) applied to a population drawn by generate_population().
# Row order (matches one_iter_all):
#   1 Subset, 2 SRS, 3 Stratified, 4 Optimal, 5 Stratified-Random,
#   6 Optimal-Permuted, 7 Optimal-Reversed, [8 Optimal-Random if include_rand],
#   then if pilot_frac is set: Optimal-Pilot-Pooled, Optimal-Pilot-Ind,
#   Stratified-PilotSelect, Optimal-PilotSelect,
#   Pilot+Stratified-PilotSelect, Pilot+Optimal-PilotSelect.
# paired_pilot = FALSE mirrors one_iter_all() as written (two
# one_iter_optimal_pilot() calls, fresh pilot for Ind).
one_iter_all_fast <- function(params, pilot_frac = NULL, include_rand = TRUE,
                              paired_pilot = FALSE) {
    pc <- pop_fast(params)
    pc1 <- one_stratum_fast(pc)
    n_code <- round(params$N * params$coding_budget)
    nr <- 7L + as.integer(include_rand) + if (!is.null(pilot_frac)) 6L else 0L
    out <- matrix(NA_real_, nr, 9L)
    
    out[1L, 1:7] <- subset_est_fast(pc, srs_draw_fast(pc, n_code))
    out[2L, 1:7] <- est_fast(pc1, srs_draw_fast(pc, n_code))
    out[3L, 1:7] <- est_fast(pc, gsi_fast(pc, n_code))
    out[4L, 1:7] <- est_fast(pc, alloc_fast(pc, n_code, oracle_var_fast(pc, "none")))
    pcr <- random_strata_fast(pc)
    out[5L, 1:7] <- est_fast(pcr, gsi_fast(pcr, n_code))
    out[6L, 1:7] <- est_fast(pc, alloc_fast(pc, n_code, oracle_var_fast(pc, "permute")))
    out[7L, 1:7] <- est_fast(pc, alloc_fast(pc, n_code, oracle_var_fast(pc, "reverse")))
    j <- 7L
    if (include_rand) {
        j <- 8L
        out[8L, 1:7] <- est_fast(pc, alloc_fast(pc, n_code, oracle_var_fast(pc, "random")))
    }
    if (!is.null(pilot_frac)) {
        cands <- candidates_fast(pc)                            # no RNG
        Sp <- gsi_fast(pc, round(pilot_frac * n_code))
        vt <- pilotvar_fast(pc, Sp)
        out[j + 1L, 1:7] <- est_fast(pc, alloc_fast(pc, n_code, vt, S_pilot = Sp))
        out[j + 1L, 9L]  <- sum(Sp)
        if (!paired_pilot) {
            Sp <- gsi_fast(pc, round(pilot_frac * n_code))
            vt <- pilotvar_fast(pc, Sp)
        }
        out[j + 2L, 1:7] <- est_fast(pc, alloc_fast(pc, n_code - sum(Sp), vt))
        out[j + 2L, 9L]  <- sum(Sp)
        out[j + 3:6, ]   <- pilot_select_fast(pc, cands, n_code, pilot_frac, pc1 = pc1)
    }
    out
}

all_method_labels <- function(include_rand = TRUE, pilot_frac = NULL) {
    c("Subset", "SRS", "Stratified", "Optimal", "Stratified-Random",
      "Optimal-Permuted", "Optimal-Reversed",
      if (include_rand) "Optimal-Random",
      if (!is.null(pilot_frac)) c("Optimal-Pilot-Pooled", "Optimal-Pilot-Ind",
                                  "Stratified-PilotSelect", "Optimal-PilotSelect",
                                  "Pilot+Stratified-PilotSelect", "Pilot+Optimal-PilotSelect"))
}

# RNG-identical replacement for sim_oracle() in 01_run_full_sim.R (full
# coding: difference in arm means, SE from the arm variances).
oracle_fast <- function(params) {
    pc <- pop_fast(params)
    m <- c(mean(pc$y[pc$arm_rows[[1]]]), mean(pc$y[pc$arm_rows[[2]]]))
    c(est = m[2] - m[1], SE = sqrt(sum(pc$var.y.z / pc$N.z)))
}


FAST_COLS <- c("est", "SE", "SEhat", "SEhat.full", "SEhat.samp", "S2.0", "S2.1",
               "strat_choice", "n_pilot")

fast_method_labels <- function(pilot_fracs) {
    list(method = c("Stratified", "Optimal", "Optimal-Random", "Optimal-Permuted",
                    "Optimal-Reversed",
                    rep(c("Optimal-Pilot-Pooled", "Optimal-Pilot-Ind",
                          "Stratified-PilotSelect", "Optimal-PilotSelect",
                          "Pilot+Stratified-PilotSelect", "Pilot+Optimal-PilotSelect"),
                        times = length(pilot_fracs))),
         pilot_frac = c(rep(NA_real_, 5L), rep(pilot_fracs, each = 6L)))
}


## ---- 11. Equivalence check against the dplyr path --------------------
# Run this once after any change to the package internals, to
# sim_functions.R, or to this file.  Requires library(stratsampling),
# source("scripts/sim_functions.R") and one_iter_pilot_study() from
# 01b_run_sim_pilot.R in the calling environment.  Confirms the fast
# path draws the SAME samples (identical RNG stream) and returns the
# same numbers, the same chosen stratification and the same pilot size.
verify_fast <- function(param_grid, pilot_fracs, seeds = NULL, paired_pilot = FALSE,
                        tol = 1e-8, verbose = TRUE) {
    if (is.null(seeds)) seeds <- 1000 + seq_len(nrow(param_grid))
    lab <- fast_method_labels(pilot_fracs)
    worst <- 0; bad <- character(0)
    for (i in seq_len(nrow(param_grid))) {
        set.seed(seeds[i])
        o1 <- suppressWarnings(one_iter_pilot_study(param_grid[i, ], pilot_fracs))
        set.seed(seeds[i])
        o2 <- suppressWarnings(one_iter_pilot_fast(param_grid[i, ], pilot_fracs,
                                                   paired_pilot = paired_pilot))
        stopifnot(nrow(o1) == nrow(o2), identical(o1$method, lab$method))
        m1 <- as.matrix(o1[, FAST_COLS[1:7]])
        rel <- abs(m1 - o2[, 1:7]) / pmax(abs(m1), 1e-12)
        w_i <- max(rel, na.rm = TRUE)
        worst <- max(worst, w_i)
        sel <- !is.na(o2[, 8])
        same_choice <- identical(CAND_NAMES[o2[sel, 8]], o1$strat_choice[sel])
        same_npilot <- isTRUE(all.equal(o2[sel, 9], o1$n_pilot[sel]))
        if (w_i > tol || !same_choice || !same_npilot)
            bad <- c(bad, sprintf("scenario %d: max rel diff %.2e, choice %s, n_pilot %s",
                                  i, w_i, same_choice, same_npilot))
        if (verbose && i %% 10 == 0) cat(sprintf("  %d / %d scenarios, worst so far %.2e\n",
                                                 i, nrow(param_grid), worst))
    }
    if (length(bad)) { cat("MISMATCHES:\n"); cat(paste0("  ", bad), sep = "\n") }
    invisible(list(worst = worst, mismatches = bad))
}


# Main-grid check: requires one_iter() from 01_run_full_sim.R and
# sim_oracle() in the calling environment (both call the package).
verify_all_fast <- function(param_grid, pilot_frac = NULL, include_rand = TRUE,
                            seeds = NULL, paired_pilot = FALSE, tol = 1e-8, verbose = TRUE) {
    if (is.null(seeds)) seeds <- 2000 + seq_len(nrow(param_grid))
    lab <- all_method_labels(include_rand, pilot_frac)
    worst <- 0; bad <- character(0)
    for (i in seq_len(nrow(param_grid))) {
        set.seed(seeds[i])
        o1 <- suppressWarnings(one_iter(param_grid[i, ], pilot_frac = pilot_frac))
        set.seed(seeds[i])
        o2 <- suppressWarnings(one_iter_all_fast(param_grid[i, ], pilot_frac = pilot_frac,
                                                 include_rand = include_rand,
                                                 paired_pilot = paired_pilot))
        stopifnot(nrow(o1) == nrow(o2), identical(o1$method, lab))
        m1 <- as.matrix(o1[, FAST_COLS[1:7]])
        rel <- abs(m1 - o2[, 1:7]) / pmax(abs(m1), 1e-12)
        w_i <- max(rel, na.rm = TRUE)
        worst <- max(worst, w_i)
        ok_choice <- TRUE
        if (!is.null(pilot_frac)) {
            sel <- !is.na(o2[, 8])
            ok_choice <- identical(CAND_NAMES[o2[sel, 8]], o1$strat_choice[sel])
        }
        if (w_i > tol || !ok_choice)
            bad <- c(bad, sprintf("scenario %d: max rel diff %.2e, choice %s", i, w_i, ok_choice))
        if (verbose && i %% 10 == 0) cat(sprintf("  %d / %d scenarios, worst so far %.2e\n",
                                                 i, nrow(param_grid), worst))
    }
    if (length(bad)) { cat("MISMATCHES:\n"); cat(paste0("  ", bad), sep = "\n") }
    invisible(list(worst = worst, mismatches = bad))
}