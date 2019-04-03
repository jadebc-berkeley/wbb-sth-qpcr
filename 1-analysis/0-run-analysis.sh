
#!/bin/bash

rm -f 1-tfgh-preliminary.Rout
rm -f 2-tfgh-descriptive.Rout
rm -f 3-tfgh-prev-mean.Rout
rm -f 4-tfgh-perm-test.Rout
rm -f 6-tfgh-sens-qpcr-gold.Rout
rm -f 7-tfgh-sens-pooled-gold.Rout
rm -f 8-tfgh-lca-inputs.Rout

R CMD BATCH 1-tfgh-preliminary.R
R CMD BATCH 2-tfgh-descriptive.R
R CMD BATCH 3-tfgh-prev-mean.R
R CMD BATCH 4-tfgh-perm-test.R
R CMD BATCH 6-tfgh-sens-qpcr-gold.R
R CMD BATCH 7-tfgh-sens-pooled-gold.R
R CMD BATCH 8-tfgh-lca-inputs.R
