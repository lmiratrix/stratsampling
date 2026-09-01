
# Tests for the revised create_stratification_candidates(): within-arm
# quantile cuts, min_cell_size and duplicate screening, and the extended
# constructions (n_groups, interaction_configs, breaks, custom, levels).

make_test_dat <- function(N = 1000, delta = 0, seed = 1) {
    set.seed(seed)
    Z <- rep(0:1, N / 2)
    tibble::tibble(
        Z = Z,
        Yhat = rnorm(N, 3 + delta * Z, 1),
        word_count = rpois(N, 40),
        n_distinct_Yhat = sample(1:3, N, replace = TRUE, prob = c(.5, .35, .15)),
        all_agree = as.integer(n_distinct_Yhat == 1),
        rare = c(rep("a", 3), rep("b", N - 3))
    )
}

# Canonical partition key (same construction as the internal helper)
pkey <- function(s) paste(match(as.character(s), unique(as.character(s))), collapse = ",")


test_that("default call reproduces the legacy candidate names and structure", {
    dat <- make_test_dat()
    cands <- create_stratification_candidates(dat, vars = c("Yhat", "word_count"),
                                              Z = "Z", verbose = FALSE)
    expect_true(all(c("Yhat_tertiles", "Yhat_quartiles", "Yhat_quintiles",
                      "word_count_tertiles", "word_count_quartiles", "word_count_quintiles",
                      "Yhat_word_count_2x2", "Yhat_word_count_2x3",
                      "Yhat_word_count_3x2", "Yhat_word_count_3x3") %in% names(cands)))
    expect_true(all(sapply(cands, length) == nrow(dat)))
    expect_true(all(sapply(cands, is.factor)))
    scr <- attr(cands, "screening")
    expect_s3_class(scr, "tbl_df")
    expect_equal(names(scr), c("candidate", "n_strata", "min_cell", "kept", "reason"))
    expect_true(all(scr$kept[scr$candidate %in% names(cands)]))
})


test_that("quantile cuts are made within arm by default, pooled on request", {
    dat <- make_test_dat(delta = 1.5)   # treated arm shifted up by 1.5 SD
    
    by_arm <- create_stratification_candidates(dat, vars = "Yhat", Z = "Z", verbose = FALSE)
    tab <- table(by_arm$Yhat_tertiles, dat$Z)
    # equal-sized strata within each arm (up to rounding)
    expect_true(all(abs(tab[, "0"] - 500 / 3) <= 1))
    expect_true(all(abs(tab[, "1"] - 500 / 3) <= 1))
    expect_equal(levels(by_arm$Yhat_tertiles), c("Low Yhat", "Med Yhat", "High Yhat"))
    
    pooled <- create_stratification_candidates(dat, vars = "Yhat", Z = "Z",
                                               cut_by_arm = FALSE, verbose = FALSE)
    tab_p <- table(pooled$Yhat_tertiles, dat$Z)
    expect_gt(max(tab_p[, "1"]) - min(tab_p[, "1"]), 100)   # imbalanced under shift
    expect_true(all(abs(rowSums(tab_p) - rowSums(tab)) <= 2))       # both about N/3 overall
    
    # pooled cuts equal the single-sample construction
    expect_equal(as.character(pooled$Yhat_tertiles),
                 as.character(make_quantile_strata(dat$Yhat, 3, "Yhat")))
    
    # without Z, cut_by_arm is irrelevant
    no_z <- create_stratification_candidates(dat, vars = "Yhat", verbose = FALSE)
    expect_equal(as.character(no_z$Yhat_tertiles),
                 as.character(make_quantile_strata(dat$Yhat, 3, "Yhat")))
})


test_that("interactions use within-arm cuts too", {
    dat <- make_test_dat(delta = 1.5)
    cands <- create_stratification_candidates(dat, vars = c("Yhat", "word_count"),
                                              Z = "Z", verbose = FALSE)
    s <- cands$Yhat_word_count_2x2
    tab <- table(s, dat$Z)
    # 4 cells per arm, each about 125 (balanced margins within arm)
    expect_equal(nrow(tab), 4)
    expect_true(all(tab > 80))
})


test_that("min_cell_size screen drops candidates with thin stratum-by-arm cells", {
    dat <- make_test_dat()
    
    c_default <- create_stratification_candidates(dat, vars = c("Yhat", "rare"),
                                                  Z = "Z", verbose = FALSE)
    expect_false("rare_levels" %in% names(c_default))
    scr <- attr(c_default, "screening")
    expect_true(any(grepl("min_cell_size = 2", scr$reason)))
    expect_true(all(sapply(c_default, function(s) min(table(s, dat$Z))) >= 2))
    
    c_off <- create_stratification_candidates(dat, vars = c("Yhat", "rare"),
                                              Z = "Z", min_cell_size = 0, verbose = FALSE)
    expect_true("rare_levels" %in% names(c_off))
    
    c_50 <- create_stratification_candidates(dat, vars = c("Yhat", "rare"),
                                             Z = "Z", min_cell_size = 50, verbose = FALSE)
    expect_true(all(sapply(c_50, function(s) min(table(s, dat$Z))) >= 50))
    
    # Without Z the screen applies to strata rather than cells
    c_noZ <- create_stratification_candidates(dat, vars = "rare",
                                              min_cell_size = 3, verbose = FALSE)
    expect_true("rare_levels" %in% names(c_noZ))
    
    expect_error(create_stratification_candidates(dat, vars = "Yhat", min_cell_size = -1))
})


test_that("a stratum that is empty in one arm fails the screen", {
    dat <- make_test_dat()
    dat$grp <- ifelse(dat$Z == 0 & seq_len(nrow(dat)) <= 10, "only_control", "everyone")
    cands <- create_stratification_candidates(dat, vars = c("Yhat", "grp"), Z = "Z",
                                              verbose = FALSE)
    expect_false(any(grepl("grp", names(cands))))
    scr <- attr(cands, "screening")
    expect_true(all(scr$min_cell[grepl("grp", scr$candidate) & !is.na(scr$min_cell)] == 0))
})


test_that("categorical variables are offered at their native levels", {
    dat <- make_test_dat()
    cands <- create_stratification_candidates(dat, vars = c("Yhat", "n_distinct_Yhat"),
                                              Z = "Z", verbose = FALSE)
    expect_true("n_distinct_Yhat_levels" %in% names(cands))
    expect_equal(nlevels(cands$n_distinct_Yhat_levels), 3)
    expect_true("Yhat_n_distinct_Yhat_2xlevels" %in% names(cands))
    # tertiles of a 3-level variable coincide with its levels: logged as duplicate
    scr <- attr(cands, "screening")
    expect_equal(scr$reason[scr$candidate == "n_distinct_Yhat_tertiles"],
                 "duplicate of 'n_distinct_Yhat_levels'")
    
    # max_levels = 0 switches native levels off for numeric variables
    no_lev <- create_stratification_candidates(dat, vars = c("Yhat", "n_distinct_Yhat"),
                                               Z = "Z", max_levels = 0, verbose = FALSE)
    expect_false("n_distinct_Yhat_levels" %in% names(no_lev))
    expect_true("n_distinct_Yhat_tertiles" %in% names(no_lev))
})


test_that("user-supplied breaks are honoured and participate in crossings", {
    dat <- make_test_dat()
    cands <- create_stratification_candidates(
        dat, vars = c("Yhat", "word_count"), Z = "Z",
        breaks = list(Yhat = c(-Inf, 2, 3, 4, Inf)), verbose = FALSE)
    expect_true("Yhat_cut" %in% names(cands))
    expect_equal(nlevels(cands$Yhat_cut), 4)
    expect_true("Yhat_word_count_cutx2" %in% names(cands))
    expect_true("Yhat_word_count_cutx3" %in% names(cands))
    # break cuts are the same in both arms
    expect_equal(as.character(cands$Yhat_cut),
                 paste0("Yhat=", as.character(cut(dat$Yhat, c(-Inf, 2, 3, 4, Inf),
                                                  include.lowest = TRUE, dig.lab = 4))))
    
    cands2 <- create_stratification_candidates(
        dat, vars = "Yhat",
        breaks = list(Yhat = list(coarse = c(-Inf, 3, Inf), fine = c(-Inf, 2, 3, 4, Inf))),
        verbose = FALSE)
    expect_true(all(c("Yhat_coarse", "Yhat_fine") %in% names(cands2)))
    
    expect_error(create_stratification_candidates(dat, vars = "Yhat",
                                                  breaks = list(word_count = c(0, 40, Inf))))
    
    no_cross <- create_stratification_candidates(
        dat, vars = c("Yhat", "word_count"), Z = "Z",
        breaks = list(Yhat = c(-Inf, 2, 3, 4, Inf)), cross_levels = FALSE, verbose = FALSE)
    expect_false(any(grepl("cutx", names(no_cross))))
})


test_that("custom stratifications are included and crossed", {
    dat <- make_test_dat()
    set.seed(7)
    clus <- sample(paste0("k", 1:4), nrow(dat), replace = TRUE)
    cands <- create_stratification_candidates(
        dat, vars = "Yhat", Z = "Z",
        custom = list(embed = clus, long = function(d) d$word_count > 40),
        verbose = FALSE)
    expect_true(all(c("embed", "long", "embed_x_Yhat_tertiles", "long_x_Yhat_quintiles")
                    %in% names(cands)))
    expect_equal(nlevels(cands$embed), 4)
    
    no_cross <- create_stratification_candidates(
        dat, vars = "Yhat", custom = list(embed = clus), cross_custom = FALSE, verbose = FALSE)
    expect_false(any(grepl("^embed_x_", names(no_cross))))
    
    expect_error(create_stratification_candidates(dat, vars = "Yhat",
                                                  custom = list(bad = 1:5)))
})


test_that("n_groups and interaction_configs are configurable", {
    dat <- make_test_dat()
    cands <- create_stratification_candidates(
        dat, vars = c("Yhat", "word_count"), Z = "Z",
        n_groups = c(2, 6), interaction_configs = list(c(2, 4)), verbose = FALSE)
    expect_true(all(c("Yhat_q2", "Yhat_q6", "Yhat_word_count_2x4") %in% names(cands)))
    expect_false("Yhat_word_count_2x2" %in% names(cands))
    expect_false("Yhat_tertiles" %in% names(cands))
    
    cands_str <- create_stratification_candidates(
        dat, vars = c("Yhat", "word_count"), interaction_configs = "2x2", verbose = FALSE)
    expect_true("Yhat_word_count_2x2" %in% names(cands_str))
    expect_false("Yhat_word_count_3x3" %in% names(cands_str))
    
    cands_none <- create_stratification_candidates(
        dat, vars = c("Yhat", "word_count"), interaction_configs = NULL, verbose = FALSE)
    expect_false(any(grepl("x", sub("^.*_", "", names(cands_none)))))
    
    expect_error(create_stratification_candidates(dat, vars = "Yhat",
                                                  interaction_configs = "2x1"))
})


test_that("duplicate partitions are removed and logged", {
    dat <- make_test_dat()
    # all_agree is exactly n_distinct_Yhat collapsed to two groups
    cands <- create_stratification_candidates(dat, vars = c("n_distinct_Yhat", "all_agree"),
                                              n_groups = c(2, 3), verbose = FALSE)
    keys <- sapply(cands, pkey)
    expect_equal(anyDuplicated(keys), 0)
    scr <- attr(cands, "screening")
    expect_true(any(grepl("^duplicate of", scr$reason)))
    
    keep_dups <- create_stratification_candidates(dat, vars = c("n_distinct_Yhat", "all_agree"),
                                                  n_groups = c(2, 3), drop_duplicates = FALSE,
                                                  verbose = FALSE)
    expect_gt(length(keep_dups), length(cands))
    expect_gt(anyDuplicated(sapply(keep_dups, pkey)), 0)
})


test_that("candidates with missing values are dropped and logged", {
    dat <- make_test_dat()
    dat$Yna <- dat$Yhat
    dat$Yna[1:5] <- NA
    cands <- create_stratification_candidates(dat, vars = c("Yhat", "Yna"), Z = "Z",
                                              verbose = FALSE)
    expect_false(any(grepl("Yna", names(cands))))
    scr <- attr(cands, "screening")
    expect_true(all(scr$reason[grepl("Yna", scr$candidate)] == "contains missing values"))
    expect_true("Yhat_tertiles" %in% names(cands))
})


test_that("verbose message summarises the screen", {
    dat <- make_test_dat()
    expect_message(create_stratification_candidates(dat, vars = c("Yhat", "rare"), Z = "Z"),
                   "constructed")
    expect_silent(create_stratification_candidates(dat, vars = "Yhat", Z = "Z",
                                                   verbose = FALSE))
})