#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# compare individual LL counter accuracy against expert

# by Jade Benjamin-Chung
# jadebc@berkeley.edu
#######################################

# load raw kk data
d = read.csv("~/Dropbox/WASHB Parasites/Data/Year 2/Temp/2-WASHB-P-kk-temp.csv")
d = d %>% mutate(dataid = as.character(dataid))

# load qpcr kk data
load(paste0(data_dir,"qdata.RData"))

# filter to ids in this analysis
qdata_ids = qdata %>% dplyr::select(dataid, personid) %>% 
  distinct() %>%
  mutate(inqpcr = 1)

d = full_join(d, qdata_ids, by = c("dataid", "personid")) %>%
  filter(inqpcr == 1)

# fix counter typo
d = d %>% 
  mutate(counter = as.character(counter)) %>%
  mutate(counter = ifelse(counter=="HKC ", "HKC", counter))

d = d %>%
  # create pos/neg for each 
  mutate(
    alpos_orig = ifelse(originalAL > 0, 1, 0),
    alpos_recount = ifelse(recountAL > 0, 1, 0),
    alpos_expert = ifelse(expertAL > 0, 1, 0),
    
    hwpos_orig = ifelse(originalHW > 0, 1, 0),
    hwpos_recount = ifelse(recountHW > 0, 1, 0),
    hwpos_expert = ifelse(expertHW > 0, 1, 0),
    
    ttpos_orig = ifelse(originalTT > 0, 1, 0),
    ttpos_recount = ifelse(recountTT > 0, 1, 0),
    ttpos_expert = ifelse(expertTT > 0, 1, 0)
  
  ) %>%
  # indicator for concordance with expert
  mutate(
    al_orig_match =    ifelse(alpos_orig == alpos_expert, 1, 0),
    al_recount_match = ifelse(alpos_recount == alpos_expert, 1, 0),
    
    hw_orig_match =    ifelse(hwpos_orig == hwpos_expert, 1, 0),
    hw_recount_match = ifelse(hwpos_recount == hwpos_expert, 1, 0),
    
    tt_orig_match =    ifelse(ttpos_orig == ttpos_expert, 1, 0),
    tt_recount_match = ifelse(ttpos_recount == ttpos_expert, 1, 0),
    
    counter = as.factor(counter)
  ) 

########################################################
# Kappa statistics
########################################################

get_kappa = function(counter_name, origvar, expertvar){
  kappa_df = d %>% 
    filter(counter == counter_name) %>%
    filter(!is.na(!!sym(expertvar))) %>%
    dplyr::select(!!sym(origvar), !!sym(expertvar))
  
  kappa_test = kappa2(kappa_df, "unweighted")
  kappa = kappa_test$value
  pval = kappa_test$p.value
  
  out = list(kappa = kappa, pval = pval)
  
  return(out)
}

al_kappa = list()
for(i in 1:length(unique(d$counter))){
  al_kappa[[i]] = get_kappa(
    counter_name = unique(d$counter)[i],
    origvar = "alpos_orig",
    expertvar = "alpos_expert"
  )
}

al_kappa_df = bind_rows(al_kappa) %>% 
  mutate(outcome = "Ascaris",
         counter = unique(d$counter))

hw_kappa = list()
for(i in 1:length(unique(d$counter))){
  hw_kappa[[i]] = get_kappa(
    counter = unique(d$counter)[i],
    origvar = "hwpos_orig",
    expertvar = "hwpos_expert"
  )
}
hw_kappa_df = bind_rows(hw_kappa) %>% 
  mutate(outcome = "Hookworm",
         counter = unique(d$counter))

tt_kappa = list()
for(i in 1:length(unique(d$counter))){
  tt_kappa[[i]] = get_kappa(
    counter = unique(d$counter)[i],
    origvar = "ttpos_orig",
    expertvar = "ttpos_expert"
  )
}
tt_kappa_df = bind_rows(tt_kappa) %>% 
  mutate(outcome = "Trichuris",
         counter = unique(d$counter))

kappa_all = bind_rows(al_kappa_df, hw_kappa_df, tt_kappa_df) %>%
  mutate(countern = case_when(
    counter == "HKC" ~ 1,
    counter == "MHR" ~ 2, 
    counter == "RK" ~ 3,
    counter == "SNJ" ~ 4
  )) %>% 
  mutate(pval_f = ifelse(pval < 0.001, "<0.001", sprintf("%0.03f", pval))) %>%
  dplyr::select(-c(pval, counter)) 

kappa_wide_kappa = dcast(kappa_all %>% dplyr::select(-pval_f), 
                         countern ~ outcome, id.var = "countern", value.var = "kappa")
kappa_wide_pval= dcast(kappa_all %>% dplyr::select(-kappa), 
                       countern ~ outcome, id.var = "countern", value.var = "pval_f") %>%
  rename(
    al_pval = Ascaris,
    hw_pval = Hookworm,
    tt_pval = Trichuris
  )

kappa_result = full_join(kappa_wide_kappa, kappa_wide_pval, by = "countern") %>%
  dplyr::select(countern, Ascaris, al_pval, 
         Hookworm, hw_pval, 
         Trichuris, tt_pval)

write.csv(kappa_result, paste0(tab_dir, "kk_tech_kappa.csv"))

########################################################
# Calculate probability of concordance and SE 
########################################################
# Ascaris --------------------------------------- 
diff_alpos_orig = lapply(as.list(levels(d$counter)), function(x)
  mean_se(d$al_orig_match[d$counter==x]))
diff_alpos_orig = bind_rows(diff_alpos_orig) %>%
  mutate(Counter = levels(d$counter),
         outcome = "Ascaris",
         count = "Original count")

alpos_recount = d %>% filter(!is.na(recountAL))
diff_alpos_recount = lapply(as.list(levels(alpos_recount$recounter)), function(x)
  mean_se(alpos_recount$al_orig_match[alpos_recount$recounter==x]))
diff_alpos_recount = bind_rows(diff_alpos_recount) %>%
  mutate(Counter = levels(alpos_recount$recounter),
         outcome = "Ascaris",
         count = "Recount") %>%
  filter(Counter != "")

# Hookworm --------------------------------------- 
diff_hwpos_orig = lapply(as.list(levels(d$counter)), function(x)
  mean_se(d$hw_orig_match[d$counter==x]))
diff_hwpos_orig = bind_rows(diff_hwpos_orig) %>%
  mutate(Counter = levels(d$counter),
         outcome = "Hookworm",
         count = "Original count")

hwpos_recount = d %>% filter(!is.na(recountHW))
diff_hwpos_recount = lapply(as.list(levels(hwpos_recount$recounter)), function(x)
  mean_se(hwpos_recount$al_orig_match[hwpos_recount$recounter==x]))
diff_hwpos_recount = bind_rows(diff_hwpos_recount) %>%
  mutate(Counter = levels(hwpos_recount$recounter),
         outcome = "Hookworm",
         count = "Recount") %>%
  filter(Counter != "")

# Trichuris --------------------------------------- 
diff_ttpos_orig = lapply(as.list(levels(d$counter)), function(x)
  mean_se(d$tt_orig_match[d$counter==x]))
diff_ttpos_orig = bind_rows(diff_ttpos_orig) %>%
  mutate(Counter = levels(d$counter),
         outcome = "Trichuris",
         count = "Original count")

ttpos_recount = d %>% filter(!is.na(recountTT))
diff_ttpos_recount = lapply(as.list(levels(ttpos_recount$recounter)), function(x)
  mean_se(ttpos_recount$al_orig_match[ttpos_recount$recounter==x]))
diff_ttpos_recount = bind_rows(diff_ttpos_recount) %>%
  mutate(Counter = levels(ttpos_recount$recounter),
         outcome = "Trichuris",
         count = "Recount") %>%
  filter(Counter != "")



########################################################
# plot
########################################################

palette = brewer.pal(n = 4, name = "Set1")

plotpos_data = bind_rows(
  diff_alpos_orig, diff_alpos_recount,
  diff_hwpos_orig, diff_hwpos_recount,
  diff_ttpos_orig, diff_ttpos_recount
) %>%
  mutate(countern = case_when(
    Counter == "HKC" ~ 1,
    Counter == "MHR" ~ 2, 
    Counter == "RK" ~ 3,
    Counter == "SNJ" ~ 4
  )) %>%
  mutate(countern = as.factor(countern))

diffpos_plot = ggplot(plotpos_data, aes(x = countern, y = y)) + 
  geom_bar(stat="identity",width=0.5, aes(fill=countern), 
           position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = ymin, ymax = ymax, fill=countern), 
                 position = position_dodge(width = 0.5)) + 
  facet_grid( outcome ~ count, scales = "free") +
  geom_hline(yintercept = 0) + 
  scale_fill_manual(values = palette) + 
  theme_bw() + 
  ylab("Percent agreement with expert technician") +
  xlab("Kato-Katz Technician Number") +
  theme(legend.title = element_blank(),
        legend.position = "bottom")

ggsave(diffpos_plot, file = paste0(fig_dir, "wbb-qpcr-kk-counter-posneg-qa.pdf"),
       width = 6, height = 6)



########################################################
# calculate difference between original counts
# and recounts and expert EPG among positive
# by original, recount, or expert
########################################################
# Ascaris --------------------------------------- 
d_al = d %>% 
  filter(originalAL > 0 | recountAL > 0 | expertAL >0) %>%
  dplyr::select(counter, recounter, expert, 
                originalAL, recountAL, expertAL) %>%
  
  # drop if no expert measurement
  filter(!is.na(expertAL)) %>%
  
  # calculate difference
  mutate(
    diff_al_orig = originalAL - expertAL,
    diff_al_recount = recountAL - expertAL,
    
    counter = as.factor(counter),
    recounter = as.factor(recounter)
  ) 

# Hookworm --------------------------------------- 
d_hw = d %>% 
  filter(originalHW > 0 | recountHW > 0 | expertHW >0) %>%
  dplyr::select(counter, recounter, expert, 
                originalHW, recountHW, expertHW) %>%
  
  # drop if no expert measurement
  filter(!is.na(expertHW)) %>%
  
  # calculate difference
  mutate(
    diff_hw_orig = originalHW - expertHW,
    diff_hw_recount = recountHW - expertHW,
    counter = as.factor(counter),
    recounter = as.factor(recounter)
  ) 

# Trichuris --------------------------------------- 
d_tt = d %>% 
  filter(originalTT > 0 | recountTT > 0 | expertTT >0) %>%
  dplyr::select(counter, recounter, expert, 
                originalTT, recountTT, expertTT) %>%
  
  # drop if no expert measurement
  filter(!is.na(expertTT)) %>%
  
  # calculate difference
  mutate(
    diff_tt_orig = originalTT - expertTT,
    diff_tt_recount = recountTT - expertTT,
    counter = as.factor(counter),
    recounter = as.factor(recounter)
  ) 


########################################################
# estimate the mean diff and SE for each counter
########################################################

# Ascaris --------------------------------------- 
diff_al_orig = lapply(as.list(levels(d_al$counter)), function(x)
  mean_se(d_al$diff_al_orig[d_al$counter==x]))
diff_al_orig = bind_rows(diff_al_orig) %>%
  mutate(Counter = levels(d_al$counter),
         outcome = "Ascaris",
         count = "Original count")

al_recount = d_al %>% filter(!is.na(recountAL))
diff_al_recount = lapply(as.list(levels(al_recount$recounter)), function(x)
  mean_se(al_recount$diff_al_recount[al_recount$recounter==x]))
diff_al_recount = bind_rows(diff_al_recount) %>%
  mutate(Counter = levels(al_recount$recounter),
         outcome = "Ascaris",
         count = "Recount") %>%
  filter(Counter != "")

# Hookworm --------------------------------------- 
diff_hw_orig = lapply(as.list(levels(d_hw$counter)), function(x)
  mean_se(d_hw$diff_hw_orig[d_hw$counter==x]))
diff_hw_orig = bind_rows(diff_hw_orig) %>%
  mutate(Counter = levels(d_hw$counter),
         outcome = "Hookworm",
         count = "Original count")

# manually add back HKC for original HW since they
# only found negatives
diff_hw_orig = bind_rows(diff_hw_orig,
                         data.frame(y = 0.1, ymin = 0, ymax = 0,
                                    Counter = "HKC", outcome = "Hookworm",
                                    count = "Original count"))

hw_recount = d_hw %>% filter(!is.na(recountHW))
diff_hw_recount = lapply(as.list(levels(hw_recount$recounter)), function(x)
  mean_se(hw_recount$diff_hw_recount[hw_recount$recounter==x]))
diff_hw_recount = bind_rows(diff_hw_recount) %>%
  mutate(Counter = levels(hw_recount$recounter),
         outcome = "Hookworm",
         count = "Recount") %>%
  filter(Counter != "")


# Trichuris --------------------------------------- 
diff_tt_orig = lapply(as.list(levels(d_tt$counter)), function(x)
  mean_se(d_tt$diff_tt_orig[d_tt$counter==x]))
diff_tt_orig = bind_rows(diff_tt_orig) %>%
  mutate(Counter = levels(d_tt$counter),
         outcome = "Trichuris",
         count = "Original count")

tt_recount = d_tt %>% filter(!is.na(recountTT))
diff_tt_recount = lapply(as.list(levels(tt_recount$recounter)), function(x)
  mean_se(tt_recount$diff_tt_recount[tt_recount$recounter==x]))
diff_tt_recount = bind_rows(diff_tt_recount) %>%
  mutate(Counter = levels(tt_recount$recounter),
         outcome = "Trichuris",
         count = "Recount") %>%
  filter(Counter != "")


########################################################
# plot
########################################################

palette = brewer.pal(n = 4, name = "Set1")

plot_data = bind_rows(
  diff_al_orig, diff_al_recount,
  diff_hw_orig, diff_hw_recount,
  diff_tt_orig, diff_tt_recount
) %>%
  mutate(countern = case_when(
    Counter == "HKC" ~ 1,
    Counter == "MHR" ~ 2, 
    Counter == "RK" ~ 3,
    Counter == "SNJ" ~ 4
  )) %>%
  mutate(countern = as.factor(countern))

diff_plot = ggplot(plot_data, aes(x = countern, y = y)) + 
  geom_bar(stat="identity",width=0.5, aes(fill=countern), 
           position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = ymin, ymax = ymax, fill=countern), 
                 position = position_dodge(width = 0.5)) + 
  facet_grid( outcome ~ count, scales = "free") +
  geom_hline(yintercept = 0) + 
  scale_fill_manual(values = palette) + 
  theme_bw() + 
  ylab("Mean difference in single slide egg count") +
  xlab("Kato-Katz Technician Number") +
  theme(legend.title = element_blank(),
        legend.position = "bottom")

ggsave(diff_plot, file = paste0(fig_dir, "wbb-qpcr-kk-counter-epg-qa.pdf"),
       width = 6, height = 6)


