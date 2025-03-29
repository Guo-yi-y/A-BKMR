
<!-- README.md is generated from README.Rmd. Please edit that file -->

# A-BKMR

<!-- badges: start -->
<!-- badges: end -->

A-BKMR is used to perform advanced BKMR in extremely fast speed.
Additionally, you can use A-BKMR to estimate various effect that other
machine learning method can not output

## Installation

You can install the development version of pairWQS from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Guo-yi-y/A-BKMR")
```

## Example

This is a basic example which shows you how use A-BKMR:

``` r
library(aBKMR)
library(ggplot2)
k_id = sam_py_r(R = dat$Z , nd = 50, num_nn = 50, w = F)

km = a_kmbayes(y = dat$y, Z = dat$Z, X = dat$X, knots = dat$Z[k_id,], est.h = T, varsel = T, iter = 10000, family = "gaussian")

PIP_new = ExtractPIPs(km) %>% arrange(desc(PIP))    


ggplot()+
  geom_col(data = PIP_new, aes(x = reorder(variable, PIP, decreasing = T), y = PIP), fill = '#48C5B9')+
  ylab("Estimated PIPs")+
  xlab('')+
  ggtitle('')+
  scale_y_continuous(limits = c(0, 1))+
  theme_bw()

overall_res = a_OverallRiskSummaries(km, data.comps = km$data.comps)
je_res = com_je(overall_res, s=500)


ggplot()+
  geom_hline(yintercept = 0, linetype = 2, color = "red")+
  geom_pointrange(data = risks.overall, 
                  aes(quant, est_p, ymin = est_p - 1.96*sd_p, ymax = est_p + 1.96*sd_p), 
                  color = '#48C5B9', 
                  shape = 16, 
                  position = position_nudge(x = -0.01, y = 0), 
  )+
  annotate("text", x = 0.35, y = 0.4, size = 6,
           label = glue("Joint effect = {je_res$est} (je_res$upper, je_res$lower)") ) + 
  
  theme_bw()

singlevar_res = a_SingVarRiskSummaries(fit = km, y = dat$y, Z = dat$Z, X = dat$X, 
                                 qs.diff = c(0.25, 0.75), 
                                 q.fixed = c(0.25, 0.50, 0.75),
                                 method = "approx", data.comps = km$data.comps)

ue_res = com_ue(singlevar_res, s = 500)


ggplot()+
  geom_smooth(data = singlevar_res, 
              aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se), 
              color = '#48C5B9', 
              stat = "identity")+
  facet_wrap(~ variable)+
  xlab("Exposures")+
  ylab("Outcomes")



inter <- a_PredictorResponseBivar(fit = km, data.comps = fitkm$data.comps,min.plot.dist = 1)



inter2 <- PredictorResponseBivarLevels(inter2, Z = dat$Z, qs = c(0.25, 0.75))

inter2_res = com_2interaction(var1, var2, inter2)

annotations <- data.frame(
  variable1 = c("var1", "var2"),            # 对应的 variable1 值
  variable2 = c("var2", "var1"),            # 对应的 variable2 值
  x = c(-2, 1),                          # 注释在 z1 轴上的位置
  y = c(0.7, 0.7),                       # 注释在 est 轴上的位置
  label = c(glue("<i>P</i><sub>interaction</sub> = {inter2_res$P[1]}"), glue("<i>P</i><sub>interaction</sub> = {inter2_res$P[2]}"))                  # 注释文本
)


ggplot(inter2 %>% filter(variable1 %in% c(var1, var2) & variable2 %in% c(var1, var2)), 
  aes(z1, est)) + 
  geom_smooth(aes(col = quantile), stat = "identity") + 
  facet_grid(variable2 ~ variable1) +
  ggtitle("h(expos1 | quantiles of expos2)") +
  xlab("expos1") +geom_richtext(
    data = annotations, 
    aes(x = x, y = y, label = label), 
    inherit.aes = FALSE,         # 不继承主图层的美学映射
    fill = NA,                   # 去掉背景填充
    label.color = NA,            # 去掉标签边框
    color = "red",               # 注释文本颜色
    size = 5,                    # 注释文本大小
    hjust = 0,                   # 水平对齐方式
    vjust = 0                    # 垂直对齐方式
  ) 


multi_inter_res <- a_SingVarRiskSummaries(km, y = dat$y, Z = dat$Z, data.comps = km$data.comps,
                                      qs.diff = c(0.25, 0.75), 
                                      q.fixed = c(0.25, 0.50, 0.75),
                                      method = "approx")

multi_2inter = risks.singvar %>% filter(variable == var1 | variable == var2)

multi_inter1 = trend(multi_2inter$est[1:3], multi_2inter$sd[1:3])$P

multi_inter2 = trend(multi_2inter$est[4:6], multi_2inter$sd[4:6])$P

ggplot(multi_2inter, aes(variable, est, ymin = est - 1.96*sd, 
                         ymax = est + 1.96*sd, col = q.fixed)) + 
  geom_pointrange(position = position_dodge(width = 0.75)) + 
  coord_flip() + annotate("text", x = 2.4, y = 0.8, size = 6,
                          label = glue("~italic(P)[trend] == '{multi_inter1}' "),
                          parse = TRUE) + annotate("text", x = 1.4, y = -0.2, size = 6,
                                                   label = glue("~italic(P)[trend] == '{multi_inter2}'"),
                                                   parse = TRUE) 


                                      




```

Reference
