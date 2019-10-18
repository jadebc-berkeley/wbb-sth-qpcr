#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# false positive for Ascaris by epg distribution

# by Jade Benjamin-Chung
# jadebc@berkeley.edu
#######################################
rm(list=ls())

# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

# load data
load(paste0(data_dir,"qdata.RData"))

qdata = qdata %>% mutate(kkpos_concordance = case_when(
  alkk == 1 & positive.Al == 0 ~ "KK positive,\nqPCR negative",
  alkk == 1 & positive.Al == 1 ~ "KK positive,\nqPCR positive",
  alkk == 0 ~ "KK -"
)) %>%
  mutate(kkpos_concordance = as.factor(kkpos_concordance))


epg_fp = ggplot(qdata %>% filter(alkk==1 & !is.na(kkpos_concordance)), 
       aes(x = kkpos_concordance, y = alepg)) + 
  geom_boxplot(aes(fill = kkpos_concordance), width = 0.5) +
  scale_y_continuous(trans = "log10", 
                     labels = c(10, 10, 100, 1000, 10000, 50000),
                     breaks = c(10, 10, 100, 1000, 10000, 50000)) + 
  ylab(expression(paste(italic("A. lumbricoides"), " eggs per gram"))) +
  xlab("") +
  scale_fill_manual(values = c("#FC7171", "#80D96D")) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 11, color = "black"))


ggsave(epg_fp, file = paste0(fig_dir, "wbb-qpcr-fp-epg.pdf"),
       width = 4, height = 4)


# how many above moderate intensity of 5000/12 
qdata %>% filter(alkk==1 & kkpos_concordance == "KK positive, qPCR negative") %>%
  filter(alepg > 5000) %>% nrow()


