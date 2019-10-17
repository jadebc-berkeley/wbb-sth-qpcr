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
  false_al = ifelse(alkk==1 & positive.Al==0,1,0)
)

# since each slide was read by two counters, 
# need to transform data so that it is in long 
# format
counter_data = qdata %>% dplyr::select(counter1, counter2, false_al) %>%
  gather("counter","name",-false_al) %>%
  mutate(name = as.factor(name))

# plot P(False positive | counter)
false_pos_counter = lapply(as.list(levels(counter_data$name)), function(x)
  mean_se(counter_data$false_al[counter_data$name==x]))
false_pos_counter = bind_rows(false_pos_counter) %>%
  mutate(Counter = levels(counter_data$name))


false_plot = ggplot(false_pos_counter, aes(x = Counter, y = y)) + 
  geom_bar(stat="identity",width=0.5, aes(fill=Counter)) +
  geom_errorbar(aes(ymin=ymin, ymax=ymax), width=0.2) +
  xlab("Kato-Katz Technician Number") + 
  ylab("Probability of false positive") +
  theme_bw() +
  theme(legend.position="none")

ggsave(false_plot, file = paste0(fig_dir, "wbb-qpcr-kk-falsepos-counter.pdf"),
       width = 4, height = 3)

# test association
glm_fit = glm(false_al ~ counter, data = counter_data, family = "binomial")
summary(glm_fit)

# x = qdata %>% group_by(labdate, counter) %>%
#   summarise(falsepos = mean(false_al))
# 
# ggplot(x, aes(x = labdate, y = falsepos, group=counter))+
#   # geom_point(aes(col=counter)) +
#   geom_line(aes(col=counter))
# 



