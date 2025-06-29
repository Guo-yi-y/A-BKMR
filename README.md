
<!-- README.md is generated from README.Rmd. Please edit that file -->

# A-BKMR

<!-- badges: start -->
<!-- badges: end -->

A-BKMR is used to perform advanced BKMR in extremely fast speed.
Additionally, you can use A-BKMR to estimate various effect

## Installation

Attention! The core of A-BKMR is written in **Python**, so Please ensure
that the R package “reticulate” (≥1.41) is installed and that you are
running R version 4.4.3 or later.

We strongly recommend user install python and python packages “numpy”,
“numba”, “scipy”.➡️ [click to view deatiled installation
instruction](docs/python-install.md)

However, if you don’t want to install python, we also provide a virtual
python installation in pairWQS. You can directly
devtools::install_github(“Guo-yi-y/pairWQS”), then a virtual python will
be automatically installed on your computer.

You can install the development version of A-BKMR from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Guo-yi-y/A-BKMR")
```

## Example

This is a basic example which shows you how use A-BKMR:

``` r
library(aBKMR)
library(bkmr)
library(ggplot2)
library(dplyr)
library(glue)
library(ggtext)

dat <- sim_data(n = 500, M = 5, ind=1:3, Zgen="realistic", family = "gaussian"    )

k_id = sam_py_r(R = dat$Z , nd = 20, num_nn = 20, w = F)

km = a_kmbayes(y = dat$y, Z = dat$Z, X = dat$X, knots = dat$Z[k_id,], est.h = T, varsel = T, iter = 10000, family = "gaussian")

PIP_new = ExtractPIPs(km) %>% arrange(desc(PIP))    


ggplot()+
  geom_col(data = PIP_new, aes(x = reorder(variable, PIP, decreasing = T), y = PIP), fill = '#48C5B9')+
  ylab("Estimated PIPs")+
  xlab('')+
  ggtitle('')+
  scale_y_continuous(limits = c(0, 1))+
  theme_bw()
```

1.  You can plot the joint exposure-response curve between various
    exposures and outcome and estimate the joint effect of all
    exposures. For example, all exposures increase simultaneously

``` r
quants = seq(0.25, 0.75, 0.01)
newz = lapply(1:ncol(dat$Z), function(i){
  qf = data.frame(X = quantile(dat$Z[, i], quants))
  names(qf) = names(data.frame(dat$Z))[i]
  return(qf)
}) %>% bind_cols()

overall_res = a_OverallRiskSummaries_vary(km, newz = newz, data.comps = km$data.comps)
overall_res$risk_overall$quants = seq(0.25, 0.75, 0.01)

# By multiplying quants by 100, the estimated JE represents the change in outcome per 1% increase
je_res = com_je(overall_res$preds_e$postmean, overall_res$preds_e$postvar, X = quants*100)

ggplot()+
  geom_hline(yintercept = 0, linetype = 2, color = "red")+
  geom_pointrange(data = overall_res$risk_overall, 
                  aes(quants, est_p, ymin = est_p - 1.96*sd_p, ymax = est_p + 1.96*sd_p), 
                  color = '#48C5B9', 
                  shape = 16
  )+
  annotate("text", x = 0.5, y = 0.4, size = 4,
           label = glue("Joint effect = {sprintf('%.3f', je_res$est)} ({sprintf('%.3f', je_res$lower_quant)}, {sprintf('%.3f', je_res$upper_quant)})") ) + 
  
  theme_bw()
  
```

2.  You can also choose to plot the exposure-response curve with some
    exposure increase but other decrease and estimate the joint effect
    in this scenarios.

``` r

newz = data.frame(z1 = quantile(dat$Z[, 1], seq(0.25, 0.75, 0.01)), z2 = quantile(dat$Z[, 2], seq(0.75, 0.25, -0.01)), z3 = median(dat$Z[, 3]), z4 = median(dat$Z[, 4]), z5 = median(dat$Z[, 5]))

overall_res = a_OverallRiskSummaries_vary(km, newz = newz, data.comps = km$data.comps)
overall_res$risk_overall$quants = seq(0.25, 0.75, 0.01)

je_res = com_je(overall_res$preds_e$postmean, overall_res$preds_e$postvar, X = seq(0.25, 0.75, 0.01)*100)

ggplot()+
  geom_hline(yintercept = 0, linetype = 2, color = "red")+
  geom_pointrange(data = overall_res$risk_overall, 
                  aes(quants, est_p, ymin = est_p - 1.96*sd_p, ymax = est_p + 1.96*sd_p), 
                  color = '#48C5B9', 
                  shape = 16 
  )+
  annotate("text", x = 0.5, y = 0.4, size = 4,
           label = glue("Joint effect = {sprintf('%.3f', je_res$est)} ({sprintf('%.3f', je_res$lower_quant)}, {sprintf('%.3f', je_res$upper_quant)})") ) + 
  
  theme_bw()
  
```

3.  You can also choose to plot the exposure-response curve with all
    exposure increase and fixed the reference point at a specified point
    such as standard level fo exposures.

``` r

newz = lapply(1:ncol(dat$Z), function(i){
  qf = data.frame(X = quantile(dat$Z[, i], quants))
  names(qf) = names(data.frame(dat$Z))[i]
  return(qf)
}) %>% bind_cols()

overall_res = a_OverallRiskSummaries_vary(km, newz = newz, data.comps = km$data.comps, point1 = data.frame(z1 = 1, z2 = 0.5, z3 = 0.1, z4 = 0, z5 = 0))
overall_res$risk_overall$quants = seq(0.25, 0.75, 0.01)

je_res = com_je(overall_res$preds_e$postmean, overall_res$preds_e$postvar, X = seq(0.25, 0.75, 0.01)*100)

ggplot()+
  geom_hline(yintercept = 0, linetype = 2, color = "red")+
  geom_pointrange(data = overall_res$risk_overall, 
                  aes(quants, est_p, ymin = est_p - 1.96*sd_p, ymax = est_p + 1.96*sd_p), 
                  color = '#48C5B9', 
                  shape = 16 
  )+
  annotate("text", x = 0.5, y = 0.4, size = 4,
           label = glue("Joint effect = {sprintf('%.3f', je_res$est)} ({sprintf('%.3f', je_res$lower_quant)}, {sprintf('%.3f', je_res$upper_quant)})") ) + 
  
  theme_bw()
  
```

4.  Below are univariate exposure-response curve and effect

``` r

singlevar_res = a_PredictorResponseUnivar(km, quants = seq(0, 1, 0.1), method = "approx", data.comps = km$data.comps) 

singlevar_data = lapply(1:length(singlevar_res), function(i){
  
  df = singlevar_res[[i]]$risk_overall %>% 
    mutate(var_est = glue("{vars} [{sprintf('%.3f', singlevar_res[[i]]$est$est)} ({sprintf('%.3f', singlevar_res[[i]]$est$lower_quant)}, {sprintf('%.3f', singlevar_res[[i]]$est$upper_quant)}) ]"))
  return(df) 
}) %>% bind_rows() 

ggplot()+
  geom_smooth(data = singlevar_data, 
              aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se), 
              color = '#48C5B9', 
              stat = "identity")+
  facet_wrap(~ var_est)+
  xlab("Exposures")+
  ylab("Outcomes")
```

4.  Interaction effect

``` r

inter <- a_PredictorResponseBivar(fit = km, data.comps = km$data.comps,min.plot.dist = 1)



inter2 <- PredictorResponseBivarLevels(inter, Z = dat$Z, qs = c(0.25, 0.75))

inter2_res = com_2interaction("z1", "z2", inter2)

annotations <- data.frame(
  variable1 = c("z1", "z2"),           
  variable2 = c("z2", "z1"),            
  x = c(0, 0),                          
  y = c(-1, -1),                       
  label = c(glue("<i>P</i><sub>interaction</sub> = {sprintf('%.3f', inter2_res$P[1])}"), glue("<i>P</i><sub>interaction</sub> = {sprintf('%.3f', inter2_res$P[2])}"))                 
)


ggplot(inter2 %>% filter(variable1 %in% c("z1", "z2") & variable2 %in% c("z1", "z2")), 
       aes(z1, est)) + 
  geom_smooth(aes(col = quantile), stat = "identity") + 
  facet_grid(variable2 ~ variable1) +
  ggtitle("h(expos1 | quantiles of expos2)") +
  xlab("expos1") +geom_richtext(
    data = annotations, 
    aes(x = x, y = y, label = label), 
    inherit.aes = FALSE,         
    fill = NA,                   
    label.color = NA,            
    color = "red",               
    size = 5,                    
    hjust = 0,                   
    vjust = 0                    
  ) 


multi_inter_res <- a_SingVarRiskSummaries(km, y = dat$y, Z = dat$Z, data.comps = km$data.comps,
                                          qs.diff = c(0.25, 0.75), 
                                          q.fixed = c(0.25, 0.50, 0.75),
                                          method = "approx")

multi_2inter = multi_inter_res %>% filter(variable == "z1" | variable == "z2")

multi_inter1 = sprintf("%.3f", trend(multi_2inter$est_p[1:3], multi_2inter$sd_p[1:3])$p_trend)

multi_inter2 = sprintf("%.3f", trend(multi_2inter$est_p[4:6], multi_2inter$sd_p[4:6])$p_trend)

ggplot(multi_2inter, aes(variable, est_p, ymin = est_p - 1.96*sd_p, 
                         ymax = est_p + 1.96*sd_p, col = q.fixed)) + 
  geom_pointrange(position = position_dodge(width = 0.75)) + 
  coord_flip() + annotate("text", x = 2.4, y = 0.8, size = 6,
                          label = glue("~italic(P)[trend] == '{multi_inter1}' "),
                          parse = TRUE) + annotate("text", x = 1.4, y = 0.7, size = 6,
                                                   label = glue("~italic(P)[trend] == '{multi_inter2}'"),
                                                   parse = TRUE) 
```

\##Reference

Yi Guo, Huixun Jia, Ziwei Peng, Xinming Xu, Zhicheng Zhang, Keyu Pan,
Yuqin Zhou, Haidong Kan, Zhenyu Wu, Cong Liu. Advanced Bayesian Kernel
Machine Regression for Large-Scale Exposome Studies: Making the
Impossible Possible.
