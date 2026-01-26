library(tidyverse)
library(janitor)
library(broom)
library(scales)
library(gt)


# A consistent color palette for methods
method_pal <- c(
    "Subset"      = "black", #"#9467bd",
    "SRS"         = "#228222",
    "Stratified"  = "#1f77b4",
    "Optimal"     = "#ff6961",
    "Full coding" = "black"
)


# Helper for clean labels
lab_resid <- labeller(residual_var = c(
    "homogeneous" = "Residual variance: homogeneous",
    "heterogeneous" = "Residual variance: heterogeneous",
    "extreme-contrast" = "Residual variance: extreme-contrast"
))
lab_bias <- labeller(bias = c(
    "none" = "Bias: none",
    "small" = "Bias: small",
    "moderate" = "Bias: moderate",
    "large" = "Bias: large",
    "extreme-contrast" = "Bias: extreme-contrast"
))

labs = labeller(residual_var = c(
    "homogeneous" = "Residual variance: homogeneous",
    "heterogeneous" = "Residual variance: heterogeneous",
    "extreme-contrast" = "Residual variance: extreme-contrast"
),
bias = c(
    "none" = "Bias: none",
    "small" = "Bias: small",
    "moderate" = "Bias: moderate",
    "large" = "Bias: large",
    "extreme-contrast" = "Bias: extreme-contrast"
),
R2_target = c("0.4"= "R² target = 0.4", 
              "0.85" = "R² target = 0.85")
)

# -------------------------------
# 1) Figure 1 — Efficiency frontiers (Empirical SE)
#    Two slices: R2 = 0.85 and R2 = 0.40; small-multiples by residual variance.
# -------------------------------
make_frontier_plot <- function(df, r2_val=NULL,facet = "bias~residual_var", bias_vals=NULL, oracle_ref=FALSE,smooth=F) {
    d <- df %>%
        filter(method %in% c("Subset", "SRS", "Stratified", "Optimal"))
    
    if (!is.null(bias_vals)){
        d= d %>% filter(bias %in% bias_vals)
    }
    
    if (!is.null(r2_val)){
        d = d %>% filter(R2_target==r2_val)
    }
    
    if (oracle_ref){
        ref = df %>%
            filter(method=="SRS") %>%
            group_by(any_bias, any_hetero, bias.var.config, R2_target) %>%
            summarise(full_coding_se=mean(emp_se_full))
    }

    if (smooth==T & oracle_ref==T){
        out = ggplot(d, aes(x = coding_frac, y = emp_se, color = method)) + geom_point()+
            geom_smooth(se=F) +  facet_grid(facet, labeller = labs) +
            geom_hline(data=ref, aes(yintercept=full_coding_se), lty=2,col="grey50")

            
    }
    else if (smooth==F & oracle_ref==T){
        out = ggplot(d, aes(x = coding_frac, y = emp_se, color = method)) +geom_point()+
            geom_line(lwd=1.05) + geom_point(size = 2) +
            facet_grid(facet, labeller = labs) +
            geom_hline(data=ref, aes(yintercept=full_coding_se), lty=2,col="grey50")
    }
    else if (smooth==T & oracle_ref==F){
        out = ggplot(d, aes(x = coding_frac, y = emp_se, color = method)) +geom_point()+
            geom_smooth(se=F) +  facet_grid(facet, labeller = labs)
    }
    
    else { 
        out = ggplot(d, aes(x = coding_frac, y = emp_se, color = method)) +geom_point()+
            geom_line() + geom_point(size = 1) +
            facet_grid(facet, labeller = labs) 
        }

    
    out = out + 
        scale_color_manual(name="Method", values = method_pal) +
        scale_x_continuous(breaks=seq(0,1,by=0.2),limits=c(0.05,0.95),labels = percent_format(accuracy = 1)) +
        xlab("Proportion coded")+ylab("Empirical Std. Error")+
        #ggtitle(glue::glue("Efficiency frontiers (R² target = {r2_val})")) +
        theme_minimal() +
        theme(legend.position = "bottom", panel.grid.minor=element_blank() )
    

    
    out
    
}
    




make_delta_plot <- function(df, facet = "bias~residual_var", r2_val=NULL, bias_vals=NULL,smooth=F) {
    base <- df %>%
        filter(method %in% c("Stratified","Optimal")) %>%
        mutate(method=fct_relevel(method, "Stratified","Optimal"))%>%
        select(all_of(cell_keys), method, pct_var_red_vs_SRS)
    
    
    if (!is.null(bias_vals)){
        base = base %>% filter(bias %in% bias_vals)
    }
    
    if (!is.null(r2_val)){
        base = base %>% filter(R2_target==r2_val)
    }
    
    if (smooth){
    g=ggplot(base, aes(x = coding_frac, y = pct_var_red_vs_SRS, color = method)) +
        facet_grid(facet, labeller = labs)+
        geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
        geom_smooth(se=F)
        
    }
    else{
        g=ggplot(base, aes(x = coding_frac, y = pct_var_red_vs_SRS, color = method)) +
            facet_grid(facet, labeller = labs)+
            geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
            geom_line(lwd=1.1) + geom_point(size = 2)
    }
    
    g = g + scale_color_manual(name="Method", 
                               values = method_pal, labels=c("Proportional Allocation", "Optimal Allocation")) +
        scale_x_continuous(breaks=seq(0,1,by=0.2), labels = percent_format(accuracy = 1),limits=c(0.05,0.95)) +
        scale_y_continuous(breaks=seq(-10,80,by=20), labels = percent_format(scale=1),limits=c(-10,85)) +
        xlab("Proportion coded")+ ylab("% Reduction in Empirical Variance vs. SRS")+
        ggtitle(glue::glue("Incremental value of stratified sampling vs SRS (R² target = {r2_val})") ) +
        theme_minimal() +
        theme(legend.position = "bottom")
    return(g)
}



pick_nearest <- function(x, target, tol = 0.02) {
    out <- x[which.min(abs(x - target))]
    if (abs(out - target) > tol) NA_real_ else out
}

make_heatmap = function(dat, r2_val, fixed_budget) {
    
    # Collapse over bias by averaging, but keep bias range for annotation if needed
    heat_df <- dat %>%
        filter(method %in% c("Stratified","Optimal"), R2_target==r2_val) %>%
        group_by(residual_var, bias, coding_frac, method) %>%
        summarise(pct_red = mean(pct_var_red_vs_SRS, na.rm=TRUE), .groups="drop") %>%
        mutate(bias = fct_relevel(bias, "none", "moderate", "large", "extreme-contrast"),
               residual_var=fct_relevel(residual_var, "homogeneous","heterogeneous", "extreme-contrast"))
    
    
    # Choose nearest available grid points to the requested budgets
    avail <- heat_df %>%
        distinct(coding_frac) %>%
        arrange(coding_frac) %>%
        pull(coding_frac)
    
    chosen <- map_dbl(fixed_budget, ~ pick_nearest(avail, .x))
    
    heat_df %>%
        filter(coding_frac == chosen) %>%
        ggplot(aes(x = residual_var, y = bias, fill = pct_red)) +
        geom_tile(color = "white") +
        facet_wrap(~ method) +
        scale_fill_distiller(palette = "Blues", direction = 1, oob = scales::squish, limits=c(0,80),
                             name = "% var. reduction vs SRS", labels=percent_format(scale=1)) +
        xlab("Residual variance pattern")+ ylab("Bias pattern")+
        ggtitle(glue::glue("Heatmap of % variance reduction vs SRS at {percent(fixed_budget)} coded (R² target = {r2_val})")) +
        theme_minimal() + theme(legend.position = "bottom")
    
}



make_full_ratio_plot <- function(df, facet = "bias~residual_var", r2_val=NULL, bias_vals=NULL) {
    d <- df %>%
        filter(method %in% c("Subset","SRS","Stratified","Optimal","Full coding"))
    
    title="Variance inflation vs full coding"
    
    if (!is.null(r2_val)){
        d = d %>% filter(R2_target==r2_val)
        title= paste0(title, glue::glue(" (R² target = {r2_val})"))

    }
    
    if (!is.null(bias_vals)){
        d = d %>% filter(bias %in% bias_vals)
    }
    
    max.y=max(d$re_vs_fullcoding)
    y.ul=ifelse(max.y>10, 10, ceiling(max.y))
    
    ggplot(d, aes(x = coding_frac, y = re_vs_fullcoding, color = method)) +
        geom_hline(yintercept = 1, linetype = 2, color = "grey50") +
        geom_line(lwd=1.1) + geom_point(size = 2) +
        facet_grid(facet,  labeller = labs) +
        scale_color_manual(values = method_pal) +
        scale_x_continuous(breaks=seq(0,1,by=0.2),limits=c(0.05,0.95),labels = percent_format(accuracy = 1)) +
        xlab("Proportion coded")+ylab("Variance Ratio")+
        ggplot2::ggtitle(title)+
        theme_minimal() +
        theme(legend.position = "bottom",panel.grid.minor = element_blank())
}



make_coverage_plot <- function(df, facet = "bias~residual_var", r2_val = NULL, bias_vals = NULL) {
    d <- df %>%
        filter(method %in% c("Subset", "SRS", "Stratified", "Optimal"))
    
    if (!is.null(bias_vals)) {
        d <- d %>% filter(bias %in% bias_vals)
    }
    
    if (!is.null(r2_val)) {
        d <- d %>% filter(R2_target == r2_val)
    }
    
    ggplot(d, aes(x = coding_frac, y = coverage, color = method)) +
        geom_hline(yintercept = 0.95, linetype = 2, color = "grey50") +
        geom_line() + geom_point(size = 1) +
        facet_grid(facet, labeller = labs) +
        scale_color_manual(name = "Method", values = method_pal) +
        scale_x_continuous(breaks = seq(0, 1, by = 0.2), 
                           labels = percent_format(accuracy = 1), 
                           limits = c(0.05, 0.95)) +
        scale_y_continuous(breaks = seq(0.85, 1.0, by = 0.05), 
                           labels = percent_format(accuracy = 1),
                           limits = c(0.85, 1.0)) +
        xlab("Proportion coded") + ylab("95% CI Coverage Rate") +
        ggtitle(glue::glue("Coverage rates across simulation scenarios{ifelse(!is.null(r2_val), paste0(' (R² target = ', r2_val, ')'), '')}")) +
        theme_minimal() +
        theme(legend.position = "bottom", panel.grid.minor = element_blank())
}


make_bias_plot <- function(df, facet = "bias~residual_var", r2_val = NULL, bias_vals = NULL) {
    d <- df %>%
        filter(method %in% c("Subset", "SRS", "Stratified", "Optimal"))
    
    if (!is.null(bias_vals)) {
        d <- d %>% filter(bias %in% bias_vals)
    }
    
    if (!is.null(r2_val)) {
        d <- d %>% filter(R2_target == r2_val)
    }
    
    ggplot(d, aes(x = coding_frac, y = bias_est, color = method)) +
        geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
        geom_line() + geom_point(size = 1) +
        facet_grid(facet, labeller = labs) +
        scale_color_manual(name = "Method", values = method_pal) +
        scale_x_continuous(breaks = seq(0, 1, by = 0.2), 
                           labels = percent_format(accuracy = 1), 
                           limits = c(0.05, 0.95)) +
        scale_y_continuous(limits = c(-0.1, 0.1)) +
        xlab("Proportion coded") + ylab("Monte Carlo Bias") +
        ggtitle(glue::glue("Estimator bias across simulation scenarios{ifelse(!is.null(r2_val), paste0(' (R² target = ', r2_val, ')'), '')}")) +
        theme_minimal() +
        theme(legend.position = "bottom", panel.grid.minor = element_blank())
}
