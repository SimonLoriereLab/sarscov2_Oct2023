devtools::install_github("outbreak-info/R-outbreak-info")
library(outbreakinfo)
library(dplyr)
library(tidyr)
library(reshape2)
outbreakinfo::authenticateUser()

lineages_of_interest <- c("BA.1","BA.5","BQ.1.1","XBB.1","XBB.1.16.1","XBB.1.5","XBB.1.9.1","XBF","EG.5.1.3","BA.2.86")
mutations = getMutationsByLineage(pangolin_lineage=lineages_of_interest, frequency=0.75, logInfo = FALSE)
mutS=mutations[mutations$gene=="S",]
mutS=mutS[mutS$prevalence>0.1,]
widemutS=dcast(mutS, lineage ~ mutation, value.var="prevalence")

svg(file="mutations.svg",width = 40,height=4)
plotMutationHeatmap(mutS, title = "S-gene mutations in lineages")
dev.off()
