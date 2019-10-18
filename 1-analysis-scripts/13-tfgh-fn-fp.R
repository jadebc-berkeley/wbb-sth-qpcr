#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# false positive and false negative probabilities
# by date of DNA extraction

# by Jade Benjamin-Chung
# jadebc@berkeley.edu
#######################################
rm(list=ls())

# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

# load qpcr data
load(paste0(data_dir,"qdata.RData"))

# date of assay 
d = read.csv("~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/Extraction date tracking_Bangladesh.csv")

d = d %>% rename(sampleid = Barcode) %>%
  mutate(Process.Date = as.character(Process.Date)) %>%
  # impute dates for last batch that are missing dates
  mutate(Process.Date = ifelse(is.na(Process.Date), "2018-03-02", Process.Date)) %>%
  mutate(date = as.Date(Process.Date, "%m/%d/%y")) %>%
  dplyr::select(sampleid, date) 

# merge datasets
qdate = left_join(qdata, d, by = "sampleid") %>%
  filter(!is.na(positive.Al))

assert_that(nrow(qdate) == 2799)


#-----------------------------------------
# Ascaris false positives by KK 
#-----------------------------------------
N_fp = nrow(qdate %>% filter(positive.Al == 0) )

fp = qdate %>% 
  filter(positive.Al == 0) %>%
  group_by(date) %>%
  summarise(sumalkk = sum(alkk)) %>%
  mutate(cumsumalkk = cumsum(sumalkk)) %>%
  mutate(cuminc_alkk_pos = cumsumalkk / N_fp) 

date_plot_al_fp = ggplot(fp, aes(x = date, y = cuminc_alkk_pos)) + 
  geom_line() +
  scale_x_date(breaks = "1 month", date_labels = "%b %Y") + 
  theme_bw() + 
  ylab("Cumulative probability of false positive") +
  xlab("DNA extraction date")

ggsave(date_plot_al_fp, file = paste0(fig_dir, "wbb-qpcr-kk-falsepos-dna-ext-date.pdf"),
       width = 6, height = 3)

##########################################
# False negatives by KK 
##########################################

#-----------------------------------------
# Ascaris
#-----------------------------------------
N_fn_al = nrow(qdate %>% filter(positive.Al == 1) )

fn_al = qdate %>% 
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
N_fn_hw = nrow(qdate %>% filter(positive.Hw == 1) )

fn_hw = qdate %>% 
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
N_fn_hw = nrow(qdate %>% filter(positive.Tt == 1) )

fn_tt = qdate %>% 
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
  xlab("DNA extraction date") +
  theme(legend.title = element_blank(),
        legend.position = "bottom")

ggsave(date_plot_fn, file = paste0(fig_dir, "wbb-qpcr-kk-falseneg-dna-ext-date.pdf"),
       width = 6, height = 4)

#-----------------------------------------
# Ascaris false positives by KK 
#-----------------------------------------
# 
# # indicator for Al false pos
# al_fp = qdate %>% 
#   filter(positive.Al == 0) %>%
#   mutate(false_pos_al = ifelse(alkk==1 & positive.Al ==0, 1, 0)) %>%
#   mutate(date = as.Date(date))
# 
# # calculate mean and CI by date
# mean_al_fp= matrix(NA, length(unique(al_fp$date)), 3)
# for(i in 1:length(unique(al_fp$date))){
#   values = al_fp$false_pos_al[al_fp$date == unique(al_fp$date)[i]]
#   values_complete = values[!is.na(values)]
#   mean_al_fp[i,1] = mean_se(values_complete)$y
#   mean_al_fp[i,2] = mean_se(values_complete)$ymin
#   mean_al_fp[i,3] = mean_se(values_complete)$ymax
# }
# mean_al_fp= data.frame(mean_al_fp)
# colnames(mean_al_fp) = c("y", "ymin","ymax")
# mean_al_fp= mean_al_fp%>% mutate(date = unique(al_fp$date))
# 
# # all dates
# date_al_fp = data.frame(date = seq.Date(from = as.Date(min(mean_al_fp$date, na.rm = TRUE)),
#                     to = as.Date(max(mean_al_fp$date, na.rm = TRUE)), 
#                     by = "days"))
# 
# plot_data_al_fp = full_join(date_al_fp, mean_al_fp, by = "date")
# 
# date_plot_al_fp = ggplot(plot_data_al_fp, aes(x = date, y = y)) + 
#   geom_bar(stat = "identity") + 
#   # geom_linerange(aes(ymin = ymin, ymax = ymax)) + 
#   scale_x_date(breaks = "1 month", date_labels = "%b %Y") +
#   xlab("Date of DNA extraction") + 
#   ylab("Probability of false positive") +
#   theme_bw()

# ggsave(date_plot_al_fp, file = paste0(fig_dir, "wbb-qpcr-kk-falsepos-dna-ext-date.pdf"),
#   width = 10, height = 3)

# ##########################################
# # False negatives by KK 
# ##########################################
# 
# #-----------------------------------------
# # Ascaris
# #-----------------------------------------
# # indicator for Al false negative
# al_fn = qdate %>% 
#   filter(positive.Al == 1) %>%
#   mutate(false_neg_al = ifelse(alkk==0 & positive.Al ==1, 1, 0))
# 
# # calculate mean and CI by date
# mean_al_fn= matrix(NA, length(unique(al_fn$date)), 3)
# for(i in 1:length(unique(al_fn$date))){
#   values = al_fn$false_neg_al[al_fn$date == unique(al_fn$date)[i]]
#   values_complete = values[!is.na(values)]
#   mean_al_fn[i,1] = mean_se(values_complete)$y
#   mean_al_fn[i,2] = mean_se(values_complete)$ymin
#   mean_al_fn[i,3] = mean_se(values_complete)$ymax
# }
# 
# mean_al_fn= data.frame(mean_al_fn)
# colnames(mean_al_fn) = c("y", "ymin","ymax")
# mean_al_fn= mean_al_fn%>% mutate(date = unique(al_fn$date))
# 
# # all dates
# date_al_fn = data.frame(date = seq.Date(from = min(mean_al_fn$date, na.rm = TRUE),
#                                       to = max(mean_al_fn$date, na.rm = TRUE), 
#                                       by = "days"))
# 
# plot_data_al_fn = full_join(date_al_fn, mean_al_fn, by = "date") %>%
#   mutate(outcome = "Ascaris")
# 
# 
# #-----------------------------------------
# # Ascaris
# #-----------------------------------------
# # indicator for Hw false negative
# hw_fn = qdate %>% 
#   filter(positive.Al == 1) %>%
#   mutate(false_neg_hw = ifelse(hwkk==0 & positive.Al ==1, 1, 0))
# 
# # calculate mean and CI by date
# mean_hw_fn= matrix(NA, length(unique(hw_fn$date)), 3)
# for(i in 1:length(unique(hw_fn$date))){
#   values = hw_fn$false_neg_hw[hw_fn$date == unique(hw_fn$date)[i]]
#   values_complete = values[!is.na(values)]
#   mean_hw_fn[i,1] = mean_se(values_complete)$y
#   mean_hw_fn[i,2] = mean_se(values_complete)$ymin
#   mean_hw_fn[i,3] = mean_se(values_complete)$ymax
# }
# 
# mean_hw_fn= data.frame(mean_hw_fn)
# colnames(mean_hw_fn) = c("y", "ymin","ymax")
# mean_hw_fn= mean_hw_fn%>% mutate(date = unique(hw_fn$date))
# 
# # hwl dates
# date_hw_fn = data.frame(date = seq.Date(from = min(mean_hw_fn$date, na.rm = TRUE),
#                                         to = max(mean_hw_fn$date, na.rm = TRUE), 
#                                         by = "days"))
# 
# plot_data_hw_fn = full_join(date_hw_fn, mean_hw_fn, by = "date") %>%
#   mutate(outcome = "Hookworm")
# 
# 
# #-----------------------------------------
# # Trichuris
# #-----------------------------------------
# # indicator for Hw false negative
# tt_fn = qdate %>% 
#   filter(positive.Al == 1) %>%
#   mutate(false_neg_tt = ifelse(ttkk==0 & positive.Al ==1, 1, 0))
# 
# # calculate mean and CI by date
# mean_tt_fn= matrix(NA, length(unique(tt_fn$date)), 3)
# for(i in 1:length(unique(tt_fn$date))){
#   values = tt_fn$false_neg_tt[tt_fn$date == unique(tt_fn$date)[i]]
#   values_complete = values[!is.na(values)]
#   mean_tt_fn[i,1] = mean_se(values_complete)$y
#   mean_tt_fn[i,2] = mean_se(values_complete)$ymin
#   mean_tt_fn[i,3] = mean_se(values_complete)$ymax
# }
# 
# mean_tt_fn= data.frame(mean_tt_fn)
# colnames(mean_tt_fn) = c("y", "ymin","ymax")
# mean_tt_fn= mean_tt_fn%>% mutate(date = unique(tt_fn$date))
# 
# # ttl dates
# date_tt_fn = data.frame(date = seq.Date(from = min(mean_tt_fn$date, na.rm = TRUE),
#                                         to = max(mean_tt_fn$date, na.rm = TRUE), 
#                                         by = "days"))
# 
# plot_data_tt_fn = full_join(date_tt_fn, mean_tt_fn, by = "date") %>%
#   mutate(outcome = "Trichuris")
# 
# 
# #-----------------------------------------
# # Plot
# #-----------------------------------------
# all_fn = bind_rows(plot_data_al_fn, plot_data_hw_fn, plot_data_tt_fn)
# 
# date_plot_fn = ggplot(all_fn, aes(x = date, y = y)) + 
#   geom_bar(stat = "identity") + 
#   # geom_linerange(aes(ymin = ymin, ymax = ymax)) + 
#   scale_x_date(breaks = "1 month", date_labels = "%b %Y") +
#   xlab("Date of DNA extraction") + 
#   ylab("Probability of false negative") +
#   theme_bw() +
#   facet_wrap(~outcome, ncol = 1)
# 
# ggsave(date_plot_fn, file = paste0(fig_dir, "wbb-qpcr-kk-falseneg-dna-ext-date.pdf"),
#        width = 10, height = 8)
