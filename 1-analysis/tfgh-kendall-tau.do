capture log close
log using "~/Documents/CRG/wash-benefits/bangladesh/src/wbb-sth-qpcr/1-analysis/tfgh-kendall-tau.log", replace
**************************************
* WASH Benefits Bangladesh STH KK qPCR validation
* calculate correlation between epg and CT value
* using Kendall's tau beta
**************************************

insheet using "~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/qdata.csv", clear

foreach var of varlist ctmeanac ctmeanad ctmeanal ctmeanna ctmeantt{
	replace `var' = "0" if `var'=="NA"
	destring `var', replace
}


ktau alepg ctmeanal, stats(taub)
ktau hwepg ctmeanna, stats(taub)
ktau hwepg ctmeanac, stats(taub)
ktau ttepg ctmeantt, stats(taub)


bootstrap r(tau_b), reps(1000) seed(1234) cluster(block): ktau alepg ctmeanal, stats(taub)
bootstrap r(tau_b), reps(1000) seed(1234) cluster(block): ktau hwepg ctmeanna, stats(taub)
bootstrap r(tau_b), reps(1000) seed(1234) cluster(block): ktau hwepg ctmeanac, stats(taub)
bootstrap r(tau_b), reps(1000) seed(1234) cluster(block): ktau ttepg ctmeantt, stats(taub)

log close
