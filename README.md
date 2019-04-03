# wbb-sth-qpcr

### Comparison of multi-parallel qPCR and Kato-Katz for detection of soil-transmitted helminth infection among children in rural Bangladesh

This is a repository for the above-mentioned manuscript comparing the performance of Kato-Katz and qPCR in a sample of 2,800 children who enrolled in the WASH Benefits Bangladesh trial [NCT01590095](https://clinicaltrials.gov/ct2/show/NCT01590095). This analysis was conducted by Jade Benjamin-Chung (jadebc@berkeley.edu).

### Associated protocols and datasets

All data required to run the analyses will be made available in concert with the publication of the article through the Open Science Framework: [https://osf.io/agk6w/](https://osf.io/agk6w/). Once data is downloaded from OSF, the data directory for the user must be changed in `0-config.R`. This will allow for replication of study findings using scripts in this repository. 

### Directory structure

**`1-analysis`** : analysis scripts

* `0-run-analysis.sh`: bash script that reproduces entire R workflow in this subdirectory
* `1-tfgh-preliminary.R`: R script that reads in raw data and prepares the analysis data (which is shared on OSF)
* `2-tfgh-descriptive.R`: R script that summarizes characteristics of study participants
* `3-tfgh-prev-mean.R`: R script that calculates prevalence and geometric mean of STH using each diagnostic
* `4-tfgh-perm-test.R`: R script that performs a clustered permutation test to assess whether prevalence is different between kk and qPCR
* `5-tfgh-kendall-tau.do`: Stata script that performs kendall's Tau test 
* `6-tfgh-sens-qpcr-gold.R`: R script that estimates sensitivity, specificity of kk using qPCR as the gold standard
* `7-tfgh-sens-pooled-gold.R`: R script that estimates sensitivity, specificity of kk using either diagnostic as the gold standard
* `8-tfgh-lca-inputs.R`: R script that produces the initial values for Bayesian latent class analysis models
* `9-tfgh-lca-al.odc`: WinBugs script that runs Bayesian LCA model for *A. lumbricoides*
* `9-tfgh-lca-hw.odc`: WinBugs script that runs Bayesian LCA model for hookworm
* `9-tfgh-lca-tt.odc`: WinBugs script that runs Bayesian LCA model for *T. trichiura*
* `10-tfgh-kappa-test`: R script that performs kappa test of agreement between KK and qPCR classification


**`2-fig-tab`** :  tables and figures

* `0-base-plot-functions.R`: R script that includes functions used for plotting
* `0-base-table-functions.R`: R script that includes functions used for table generation
* `tfgh-boxplot-al.R`: R script that creates a box plot comparing the qPCR Cq values for *A. lumbricoides* within levels of Kato-Katz infection intensity
* `tfgh-boxplot-al.R`: R script that creates a box plot comparing the qPCR Cq values for *A. lumbricoides* within levels of Kato-Katz infection intensity (manuscript Figure 4)
* `tfgh-plot-bar.R`: R script that creates a bar plot comparing the prevalence of each STH using Kato-Katz vs. qPCR (manuscript Figure 1)
* `tfgh-plot-scatter.R`: R script that creates a scatter plot comparing qPCR Cq values vs. Kato-Katz eggs per gram values (manuscript Figure 3)
* `tfgh-plot-venn-diagram.R`: R script that creates a venn diagram of co-infection (manuscript Figure 2)
* `tfgh-table-2x2.R`: R script that creates 2x2 tables of concordance and discordance between Kato-Katz and qPCR (manuscript Table 3)
* `tfgh-table-prevint.R`: R script that creates a table comparing the prevalence from each diagnostic, geometric mean of Kato-Katz eggs per gram, and median Cq values from qPCR (manuscript Table 2)
* `tfgh-table-sens-spec-pooled-gold.R`: R script that creates a table comparing estimates of sensitivity, specificity using either method as gold standard (manuscript Table 4, Panel 2)
* `tfgh-table-sens-spec-qgold.R`: R script that creates a table comparing estimates of sensitivity, specificity using pooled kk qPCR as gold standard (not included in manuscript)









