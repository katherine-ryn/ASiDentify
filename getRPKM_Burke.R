# ---
# Title: getRPKM_Burke.R
# Purpose: This script generates Burke_timecourse.csv and Burke_metadata.csv
#           for use in `make_ASiD_dataset_PC.R`, `make_ASiD_dataset_RBP.R` and
#           `make_ASiD_dataset_CR.R`.
# This uses Github code from Burke et al. 2020: https://github.com/LieberInstitute/libd_stem_timecourse/blob/master/tc_analysis/time_course/run_voom_gene_counts.R
#          
# ---

library(SummarizedExperiment)

setwd("") # Change this to the path that contains the ASiD folder

load("./data/libd_stemcell_timecourse_rseGene_n157.rda")


dropInd = which(rse_gene$Class != "Naked genomes" | rse_gene$CONDITION %in% c("RENEW","NEURONS_ALONE"))

rse_gene = rse_gene[,-dropInd] 

gene <- getRPKM(rse_gene)
gene <- as.data.frame(gene)
gene <- tibble::rownames_to_column(gene, var = "rowname")
colnames(gene)[1] <- "rn"

write.csv(gene, file = "./data/Burke_timecourse.csv", row.names = F)

#Get metadata
pd = colData(rse_gene)
pd$REP = as.numeric(pd$BIO_REP)
pd$REP[is.na(pd$REP)] = 1
pd$DAY = ordered(pd$DAY, levels=sort(as.numeric(levels(as.factor(pd$DAY)))) )
pd$DIV = as.numeric(as.character(pd$DAY))
pd$EXPERIMENT[pd$EXPERIMENT=="TIMECOURSE1"] = "TIMECOURSE"
pd$LINE = gsub('-','_',pd$LINE)

pd_df <- as.data.frame(pd)
write.csv(pd_df, file = "./data/Burke_metadata.csv", row.names = F)



