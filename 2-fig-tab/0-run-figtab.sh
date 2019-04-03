
#!/bin/bash

rm -f tfgh-plot-boxplot-al.Rout
rm -f tfgh-plot-bar.Rout
rm -f tfgh-plot-scatter.Rout
rm -f tfgh-plot-venn-diagram.Rout
rm -f tfgh-table-2x2.Rout
rm -f tfgh-table-prevint.Rout
rm -f tfgh-table-sens-spec-pooled-gold.Rout
rm -f tfgh-table-sens-spec-qgold.Rout

R CMD BATCH tfgh-plot-boxplot-al.R
R CMD BATCH tfgh-plot-bar.R
R CMD BATCH tfgh-plot-scatter.R
R CMD BATCH tfgh-plot-venn-diagram.R
R CMD BATCH tfgh-table-2x2.R
R CMD BATCH tfgh-table-prevint.R
R CMD BATCH tfgh-table-sens-spec-pooled-gold.R
R CMD BATCH tfgh-table-sens-spec-qgold.R
