#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# plot of bayesian LCA sens/ spec
#######################################
rm(list=ls())
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)

d = read.csv("~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/bayesian_LCA.csv", header=TRUE, sep=",")

sens_kk = d[,c("label","sens_kk", "sens_kk_lb", "sens_kk_ub")] %>%
  mutate(diagnostic = "Kato-Katz", measure = "Sensitivity") %>%
  rename(est = sens_kk, lb = sens_kk_lb, ub = sens_kk_ub)

sens_q = d[,c("label","sens_q", "sens_q_lb", "sens_q_ub")] %>%
  mutate(diagnostic = "qPCR", measure = "Sensitivity") %>%
  rename(est = sens_q, lb = sens_q_lb, ub = sens_q_ub)

spec_kk = d[,c("label","spec_kk", "spec_kk_lb", "spec_kk_ub")] %>%
  mutate(diagnostic = "Kato-Katz", measure = "Specificity") %>%
  rename(est = spec_kk, lb = spec_kk_lb, ub = spec_kk_ub)

spec_q = d[,c("label","spec_q", "spec_q_lb", "spec_q_ub")] %>%
  mutate(diagnostic = "qPCR", measure = "Specificity") %>%
  rename(est = spec_q, lb = spec_q_lb, ub = spec_q_ub)


dl = rbind(sens_kk, sens_q, spec_kk, spec_q)

ggplot(dl, aes(x = label ,y = est, group=diagnostic)) +
  geom_point(aes(col = diagnostic), position = position_dodge(width=0.4)) + 
  geom_errorbar(aes(ymin = lb, ymax = ub, col = diagnostic), 
                position = position_dodge(width=0.4), width=0.2) + 
  facet_wrap(~measure) + 
  theme_bw()  +
  xlab("Organism") + ylab("Estimate (95% CI)") +
  theme(legend.position="bottom")

