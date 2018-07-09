capture log close
log using "~/Documents/CRG/wash-benefits/bangladesh/src/wbb-sth-qpcr/1-analysis/tfgh-kendall-tau.log", replace
**************************************
* WASH Benefits Bangladesh STH KK qPCR validation
* calculate correlation between epg and DNA copies
* using Kendall's tau beta
**************************************

insheet using "~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/concentration.csv", clear

replace alepg="." if alepg=="NA"
replace hwepg="." if hwepg=="NA"
replace ttepg="." if ttepg=="NA"
destring alepg hwepg ttepg, replace

foreach var of varlist copiesal copiesna copiesad copiesac copiestt{
	replace `var' = "0" if `var'=="NA"
	destring `var', replace
}


ktau alepg copiesal, stats(taub)
ktau hwepg copiesna, stats(taub)
ktau hwepg copiesac, stats(taub)
ktau ttepg copiestt, stats(taub)


bootstrap r(tau_b), reps(1000) seed(1234) cluster(block): ktau alepg copiesal, stats(taub)
bootstrap r(tau_b), reps(1000) seed(1234) cluster(block): ktau hwepg copiesna, stats(taub)
bootstrap r(tau_b), reps(1000) seed(1234) cluster(block): ktau hwepg copiesac, stats(taub)
bootstrap r(tau_b), reps(1000) seed(1234) cluster(block): ktau ttepg copiestt, stats(taub)

log close
