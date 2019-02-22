capture log close
log using "~/Documents/CRG/wash-benefits/bangladesh/src/wbb-sth-qpcr/1-analysis/7-tfgh-kendall-tau.log", replace
**************************************
* WASH Benefits Bangladesh STH KK qPCR validation
* calculate correlation between epg and CT value
* using Kendall's tau beta
**************************************

insheet using "~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/qdata.csv", clear

foreach var of varlist ctmeanac ctmeanad ctmeanal2 ctmeanna ctmeantt positivena positiveac positivead{
	replace `var' = "" if `var'=="NA"
	destring `var', replace
}

gen ctmeanhw = .
replace ctmeanhw = ctmeanna if positivena==1 & positiveac==0 & positivead==0
replace ctmeanhw = ctmeanac if positivena==0 & positiveac==1 & positivead==0
replace ctmeanhw = ctmeanad if positivena==0 & positiveac==0 & positivead==1
replace ctmeanhw = (ctmeanna + ctmeanac)/2 if positivena==1 & positiveac==1 & positivead==0
replace ctmeanhw = (ctmeanna + ctmeanad)/2 if positivena==1 & positiveac==0 & positivead==1
replace ctmeanhw = (ctmeanna + ctmeanac)/2 if positivena==1 & positiveac==1 & positivead==0
replace ctmeanhw = (ctmeanac + ctmeanad)/2 if positivena==0 & positiveac==1 & positivead==1


ktau alepg ctmeanal2, stats(taub)
ktau hwepg ctmeanna, stats(taub)
ktau hwepg ctmeanac, stats(taub)
ktau hwepg ctmeanhw, stats(taub)
ktau ttepg ctmeantt, stats(taub)


bootstrap r(tau_b), reps(1000) seed(1234) cluster(block): ktau alepg ctmeanal2, stats(taub)
bootstrap r(tau_b), reps(1000) seed(1234) cluster(block): ktau hwepg ctmeanna, stats(taub)
bootstrap r(tau_b), reps(1000) seed(1234) cluster(block): ktau hwepg ctmeanac, stats(taub)
bootstrap r(tau_b), reps(1000) seed(1234) cluster(block): ktau ttepg ctmeantt, stats(taub)

log close
