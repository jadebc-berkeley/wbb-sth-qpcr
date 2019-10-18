#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# cumulative probability of false pos/neg over time

# by Jade Benjamin-Chung
# jadebc@berkeley.edu
#######################################
rm(list=ls())

# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

# load data
load(paste0(data_dir,"qdata.RData"))

qdata = qdata %>% mutate(date = as.Date(labdate, format = "%d%b%Y"))

#------------------------------------------------------
# ascaris false positive
#------------------------------------------------------

N_fp = nrow(qdata %>% filter(positive.Al == 0) )

fp = qdata %>% 
  filter(positive.Al == 0) %>%
  mutate(date = as.Date(labdate, format = "%d%b%Y")) %>%
  group_by(date) %>%
  summarise(sumalkk = sum(alkk)) %>%
  mutate(cumsumalkk = cumsum(sumalkk)) %>%
  mutate(cuminc_alkk_pos = cumsumalkk / N_fp) 
  
date_plot_fp = ggplot(fp, aes(x = date, y = cuminc_alkk_pos)) + 
  geom_line() +
  scale_x_date(breaks = "1 month", date_labels = "%b %Y") +
  theme_bw() + 
  ylab("Cumulative probability of false positive") +
  xlab("Kato-Katz date")


ggsave(date_plot_fp, file = paste0(fig_dir, "wbb-qpcr-kk-falsepos-kk-date.pdf"),
       width = 6, height = 3)

##########################################
# False negatives by KK 
##########################################

#-----------------------------------------
# Ascaris
#-----------------------------------------
N_fn_al = nrow(qdata %>% filter(positive.Al == 1) )

fn_al = qdata %>% 
  filter(positive.Al == 1) %>%
  mutate(rev_alkk = -(alkk - 1)) %>%
  group_by(date) %>%
  summarise(sumkk = sum(rev_alkk)) %>%
  mutate(cumsumkk = cumsum(sumkk)) %>%
  mutate(cuminc_kk_fn = cumsumkk / N_fn_al) %>%
  mutate(outcome = "Ascaris")

#-----------------------------------------
# Hookworm
#-----------------------------------------
N_fn_hw = nrow(qdata %>% filter(positive.Hw == 1) )

fn_hw = qdata %>% 
  filter(positive.Hw == 1) %>%
  mutate(rev_hwkk = -(hwkk - 1)) %>%
  group_by(date) %>%
  summarise(sumkk = sum(rev_hwkk)) %>%
  mutate(cumsumkk = cumsum(sumkk)) %>%
  mutate(cuminc_kk_fn = cumsumkk / N_fn_hw) %>%
  mutate(outcome = "Hookworm")

#-----------------------------------------
# Trichuris
#-----------------------------------------
N_fn_hw = nrow(qdata %>% filter(positive.Tt == 1) )

fn_tt = qdata %>% 
  filter(positive.Tt == 1) %>%
  mutate(rev_ttkk = -(ttkk - 1)) %>%
  group_by(date) %>%
  summarise(sumkk = sum(rev_ttkk)) %>%
  mutate(cumsumkk = cumsum(sumkk)) %>%
  mutate(cuminc_kk_fn = cumsumkk / N_fn_hw) %>%
  mutate(outcome = "Trichuris")

fn_df = bind_rows(fn_al, fn_hw, fn_tt)

date_plot_fn = ggplot(fn_df, aes(x = date, y = cuminc_kk_fn)) + 
  geom_line(aes(col = outcome)) +
  scale_x_date(breaks = "1 month", date_labels = "%b %Y") + 
  theme_bw() + 
  ylab("Cumulative probability of false negative") +
  xlab("Kato-Katz date") +
  theme(legend.title = element_blank(),
        legend.position = "bottom")  + 
  expand_limits(x = as.Date("5/1/2016", format = "%m/%d/%Y"))

ggsave(date_plot_fn, file = paste0(fig_dir, "wbb-qpcr-kk-falseneg-kk-date.pdf"),
       width = 7, height = 4)






