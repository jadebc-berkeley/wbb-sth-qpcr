#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# Check whether there is a correlation
# between the KK counter and false positive for Ascaris

# by Jade Benjamin-Chung
# jadebc@berkeley.edu
#######################################
rm(list=ls())

# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

# load data
load(paste0(data_dir,"qdata.RData"))

# create indicator for false positive Ascaris in KK
qdata = qdata %>% mutate(
  false_al = ifelse(alkk==1 & positive.Al==0,1,0),
  counter = as.factor(counter)
)

# plot P(False positive | counter)
false_pos_counter = lapply(as.list(levels(qdata$counter)), function(x) mean_se(qdata$false_al[qdata$counter==x]))
false_pos_counter = bind_rows(false_pos_counter) %>%
  mutate(Counter = levels(qdata$counter))

false_plot = ggplot(false_pos_counter, aes(x = Counter, y = y)) + 
  geom_bar(stat="identity",width=0.5, aes(fill=Counter)) +
  geom_errorbar(aes(ymin=ymin, ymax=ymax), width=0.2) +
  xlab("Kato-Katz Technician Number") + 
  ylab("Probability of false positive") +
  theme_bw() +
  theme(legend.position="none")

ggsave(false_plot, file = paste0(fig_dir, "wbb-qpcr-kk-falsepos-counter.pdf"),
       width = 6, height = 3)

# test association
glm_fit = glm(false_al ~ counter, data = qdata, family = "binomial")
summary(glm_fit)

x = qdata %>% group_by(labdate, counter) %>%
  summarise(falsepos = mean(false_al))

ggplot(x, aes(x = labdate, y = falsepos, group=counter))+
  # geom_point(aes(col=counter)) +
  geom_line(aes(col=counter))




