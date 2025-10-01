
# Analyze simulation results ----
library(dplyr)
library(ggplot2)
library(tidyr)


results1=readRDS("~/Dropbox/text and RCT sampling math notes/stratsampling/results/sim_results_balanced_strata.rds")
results2=readRDS("~/Dropbox/text and RCT sampling math notes/stratsampling/results/sim_results_balanced_approx_strata.rds")
results3=readRDS("~/Dropbox/text and RCT sampling math notes/stratsampling/results/sim_results_unbalanced_strata.rds")

results=rbind(results1,results2,results3)

results$method[results$method=="Stratified-Opt"]="Optimal"


# Calculate performance metrics
analysis_results = results %>% filter(!is.na(est)) %>% # need to figure out why some results are NA when using unbalanced stratification with n_code=100
    group_by(N, n_strat, strata_config, tau, bias, residual_var, R2_target, pop.sd, n_code, method) %>%
    mutate(CI_lower = est - 2*SEhat, CI_upper = est+2*SEhat, covers=ifelse(CI_lower<0 & CI_upper>0, 1, 0)) %>%
    summarize(
        # Bias metrics
        mean_est = mean(est),
        bias_est = mean(est - tau),

        # CI Coverage metrics
        coverage = mean(covers),

        # Variance metrics
        emp_var = var(est),
        emp_se = sd(est),
        mean_se = mean(SE),
        mean_sehat = mean(SEhat),
        # MSE metrics
        mse = mean((est - tau)^2)
    ) %>% ungroup()

analysis_results$bias.var.config = sprintf("Bias: %s, Var: %s", analysis_results$bias, analysis_results$residual_var)

res1 = analysis_results %>% filter(strata_config=="balanced-exact")
res2 = analysis_results %>% filter(strata_config=="balanced-approx")
res3 = analysis_results %>% filter(strata_config=="unbalanced")


# Empirical SE across configurations


ggplot(res1[res1$R2_target==0.85,], aes(x = n_code, y = emp_se, color = method)) +
    geom_smooth(se=F) +
    geom_point() +
    facet_wrap(~bias.var.config,ncol=3) +
labs(x = "Proportion of Documents Coded", y = "Empirical SE of impact estimate",
     title = "Efficiency Gains by Bias/Variance Configuration (Balanced Strata, target_R2=0.85)") +
    theme_minimal()
ggsave("figures/empirical_se_balanced_exact_highPredictive.pdf")

ggplot(res1[res1$R2_target==0.4,], aes(x = n_code, y = emp_se, color = method)) +
    geom_smooth(se=F) +
    geom_point() +
    facet_wrap(~bias.var.config,ncol=3) +
    labs(x = "Proportion of Documents Coded", y = "Empirical SE of impact estimate",
         title = "Efficiency Gains by Bias/Variance Configuration (Balanced Strata, target_R2=0.4)") +
    theme_minimal()
ggsave("figures/empirical_se_balanced_exact_lowPredictive.pdf")

# Calculate relative efficiency (based on ratio of empirical SEs across all iterations)
res.oracle=readRDS("~/Dropbox/text and RCT sampling math notes/stratsampling/results/sim_results_oracle.rds")


oracle = res.oracle %>% group_by(N, n_strat, strata_config, tau, bias, residual_var, R2_target, pop.sd) %>%
    summarise(mse.oracle = mean((est-mean(est))^2))



re_analysis <- merge(analysis_results, oracle, by=c("N","n_strat","strata_config","tau","bias","residual_var","R2_target", "pop.sd")) %>%
    mutate(RE_v_oracle=mse/mse.oracle) %>%
    filter(strata_config=="balanced-exact")


# 1. Relative Efficiency across different bias/variance configurations
ggplot(re_analysis[re_analysis$strata_config=="balanced-exact" & re_analysis$R2_target==0.85,], aes(x = n_code, y = RE_v_oracle,col=method)) +
    geom_line() +
    geom_point() +
    facet_grid(bias ~ residual_var,labeller = "label_both") +
    labs(x = "Proportion of Documents Coded", y = "Relative Efficiency vs. Full Coding",
         title = "Efficiency Gains by Bias/Variance Configuration (target_R2=0.85)") +
    theme_minimal()+
    geom_hline(yintercept = 1,lty=1) +
    scale_y_continuous(labels = scales::label_percent(),limits=c(0.95,7),oob=scales::oob_squish)
ggsave("figures/RE_vs_full_coding_balanced_exact_highPredictive.pdf")

ggplot(re_analysis[re_analysis$strata_config=="balanced-exact" & re_analysis$R2_target==0.4,], aes(x = n_code, y = RE_v_oracle,col=method)) +
    geom_line() +
    geom_point() +
    facet_grid(bias ~ residual_var,labeller = "label_both") +
    labs(x = "Proportion of Documents Coded", y = "Relative Efficiency vs. Full Coding",
         title = "Efficiency Gains by Bias/Variance Configuration (target_R2=0.4)") +
    theme_minimal()+
    geom_hline(yintercept = 1,lty=1) +
    scale_y_continuous(labels = scales::label_percent(),limits=c(0.95,7),oob=scales::oob_squish)
ggsave("figures/RE_vs_full_coding_balanced_exact_lowPredictive.pdf")


ggplot(re_analysis[re_analysis$method!="Subset"& re_analysis$R2_target==0.85,], aes(x = n_code, y = RE_v_oracle,col=method)) +
    geom_line() +
    geom_point() +
    facet_grid(bias ~ residual_var,labeller = "label_both") +
    labs(x = "Proportion of Documents Coded", y = "Relative Efficiency vs. Full Coding",
         title = "Efficiency Gains by Bias/Variance Configuration") +
    theme_minimal()+
    geom_hline(yintercept = 1,lty=2) +
    scale_y_continuous(labels = scales::label_percent(),limits=c(0.95,9),oob=scales::oob_squish)


ggplot(re_analysis[!re_analysis$method%in%c("Subset","SRS") & re_analysis$strata_config=="balanced-exact" & re_analysis$R2_target==0.85,],
       aes(x = n_code, y = RE_v_oracle,col=method)) +
    geom_line() +
    geom_point() +
    scale_y_continuous(labels = scales::label_percent(),limits=c(0.95,3),oob=scales::oob_squish)+
    facet_grid(bias ~ residual_var,labeller = "label_both",scales="free_y") +
    labs(x = "Proportion of Documents Coded", y = "Relative Efficiency vs. Full Coding",
         title = "Efficiency Gains by Bias/Variance Configuration") +
    theme_minimal()+
    geom_hline(yintercept = 1,lty=2)


ggplot(re_analysis[!re_analysis$method%in%c("Subset","SRS") & re_analysis$strata_config=="balanced-exact" & re_analysis$R2_target==0.85,],
       aes(x = n_code, y = RE_v_oracle,col=method)) +
    geom_smooth(se=F) +
    geom_point() +
    scale_y_continuous(labels = scales::label_percent(),limits=c(0.95,3),oob=scales::oob_squish)+
    facet_grid(bias ~ residual_var,labeller = "label_both",scales="free_y") +
    labs(x = "Proportion of Documents Coded", y = "Relative Efficiency vs. Full Coding",
         title = "Efficiency Gains by Bias/Variance Configuration (target R2=0.85)") +
    theme_minimal()+
    geom_hline(yintercept = 1,lty=2)

ggsave("figures/RE_optimal_vs_strat_balanced_exact.pdf")





re_analysis2 <- res1 %>% filter(method!="Subset") %>%
    select(N, n_strat, strata_config, tau, bias, residual_var, R2_target, pop.sd, n_code, method, emp_se, mean_se) %>%
    pivot_wider(names_from = method, values_from = c(emp_se, mean_se)) %>%
    select(N:n_code, emp_se_SRS, emp_se_Stratified, emp_se_Optimal) %>%
    pivot_longer(cols=emp_se_Stratified:emp_se_Optimal, names_to=c("method"), names_pattern="emp_se_(.*)", values_to="emp_se") %>%
    mutate(
        RE_vs_SRS = ((emp_se_SRS^2 - emp_se^2)/emp_se_SRS^2),
        var.ratio = emp_se^2/emp_se_SRS^2)



# 1. Relative Efficiency across different bias/variance configurations
ggplot(re_analysis2[re_analysis2$R2_target==0.85,], aes(x = n_code, y = RE_vs_SRS, col=method)) +
    geom_smooth(se=F) +
    geom_point() +
    facet_grid(bias ~ residual_var,labeller = "label_both") +
    labs(x = "Proportion of Documents Coded", y = "% Reduction in Empirical Variance Relative to SRS",
         title = "Efficiency Gains by Bias/Variance Configuration (target R2=0.85)") +
    theme_minimal()+
    geom_hline(yintercept = 0,lty=2) +
    scale_y_continuous(labels = scales::label_percent(),limits=c(-0.2,1),oob=scales::oob_squish)
ggsave("figures/efficiency_vs_SRS_balanced_exact_highPredictive.pdf")


ggplot(re_analysis2[re_analysis2$R2_target==0.4,], aes(x = n_code, y = RE_vs_SRS, col=method)) +
    geom_smooth(se=F) +
    geom_point() +
    facet_grid(bias ~ residual_var,labeller = "label_both") +
    labs(x = "Proportion of Documents Coded", y = "% Reduction in Empirical Variance Relative to SRS",
         title = "Efficiency Gains by Bias/Variance Configuration (target R2=0.4)") +
    theme_minimal()+
    geom_hline(yintercept = 0,lty=2) +
    scale_y_continuous(labels = scales::label_percent(),limits=c(-0.2,1),oob=scales::oob_squish)
ggsave("figures/efficiency_vs_SRS_balanced_exact_lowPredictive.pdf")


# which factors have biggest influence on MSE?

mod.strat = lm(mse ~ strata_config+bias+residual_var+R2_target+n_code, data=analysis_results[analysis_results$method=="Stratified",])

