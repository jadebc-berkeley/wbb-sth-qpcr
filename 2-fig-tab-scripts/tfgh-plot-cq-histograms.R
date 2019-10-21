#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# histograms of Cq values

# by Jade Benjamin-Chung
# jadebc@berkeley.edu
#######################################

x = qdata %>% dplyr::select(CTmean.Na, CTmean.Al, CTmean.Ac, 
                     CTmean.Tt) %>%
  melt() %>%
  filter(!is.na(value)) %>%
  mutate(sth = case_when(
    variable == "CTmean.Na" ~ "N. americanus",
    variable == "CTmean.Al" ~ "A. lumbricoides",
    variable == "CTmean.Ac" ~ "A. ceylanicum",
    variable == "CTmean.Tt" ~ "T. trichiura"
  ))

cq_bimodal = ggplot(x, aes(x = value)) + 
  geom_histogram(bins = 50) +
  facet_wrap(~ sth, ncol = 1, scale = "free") + 
  theme_bw() +
  xlab("Mean Cq value") + 
  ylab("Number of samples")


ggsave(cq_bimodal, file = paste0(fig_dir, "wbb-qpcr-kk-cq-bimodal.pdf"),
       width = 5, height = 8)

