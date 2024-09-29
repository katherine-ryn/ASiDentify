# ---
# Title: make_ASiD_dataset_PC.R
# Purpose: This script generates the columns in ASiD (for all protein-coding genes) from various databases.
#          
# ---

library(data.table)
library(dplyr)
library(biomaRt)
library(readxl)
library(SummarizedExperiment)
library(stringr)

setwd("") # Change this to the path that contains the ASiD folder
source("./nested_cv_funcs.R")

genes <- fread("./data/hgnc_complete_set_2024-04-01.tsv", data.table = F)
genes <- genes %>% filter(locus_group == 'protein-coding gene')

PCDHA_cluster <- fread("./data/HGNC_PCDHA_cluster.txt", data.table = F)
colnames(PCDHA_cluster) <- c("hgnc_id", "status", "symbol", "name", "ensembl_gene_id")

genes <- dplyr::bind_rows(genes, PCDHA_cluster) 
genes <- genes %>% filter(!ensembl_gene_id=='') 
genes <- genes[,c(1,2,3,6,20)]
genes <- genes[,c(1,4, 2, 3, 5)]
genes <- genes[order(genes$symbol),]



###1###
#Add snRNA-seq data (multiple cortical areas):
snHuman.path <- "./data/multiple_cortical_trimmed_means.csv"
snHuman <- fread(snHuman.path, data.table = F)
exc_indices <- grep("Exc", colnames(snHuman))
inh_indices <- grep("Inh", colnames(snHuman))
astro_indices <- grep("Astro", colnames(snHuman))
micro_indices <- grep("Micro", colnames(snHuman))
oligo_indices <- grep("Oligo", colnames(snHuman))

i = 1
for (i in 1:nrow(snHuman)) {
  exc_mean <- mean(as.numeric(snHuman[i, exc_indices]))
  snHuman$ExcExpression[i] <- exc_mean
  inh_mean <- mean(as.numeric(snHuman[i, inh_indices]))
  snHuman$InhExpression[i] <- inh_mean
  astro_mean <- mean(as.numeric(snHuman[i, astro_indices]))
  snHuman$AstroExpression[i] <- astro_mean
  micro_mean <- mean(as.numeric(snHuman[i, micro_indices]))
  snHuman$MicroExpression[i] <- micro_mean
  oligo_mean <- mean(as.numeric(snHuman[i, oligo_indices]))
  snHuman$OligoExpression[i] <- oligo_mean
}

#Use Ensembl Version 111 (Jan 2024)
mart <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", version = 111, dataset = "hsapiens_gene_ensembl")

gene_info <- read.csv("./data/human_MTG_2018-06-14_genes-rows.csv", header = T)
snHuman <- merge(snHuman, gene_info[, c("gene", "entrez_id")], by.x = "feature", by.y = "gene", all.x = TRUE) 

gene_df_snHuman_NCBI <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "ensembl_gene_id_version", "gene_biotype", "entrezgene_id"),
  filters = c("chromosome_name", "entrezgene_id", "biotype"),
  values = list(c(1:22, "X", "Y", "MT"), snHuman$entrez_id, c("protein_coding", "lncRNA")),
  mart = mart
)

#Only keep unique Ensembl rows. 
gene_df_distinct <- distinct_at(gene_df_snHuman_NCBI, vars(ensembl_gene_id,
                                                           external_gene_name,
                                                           ensembl_gene_id_version,
                                                           gene_biotype), .keep_all = TRUE)

snHuman <- merge(snHuman, gene_df_distinct[, c("ensembl_gene_id", "entrezgene_id")], by.x = "entrez_id", by.y = "entrezgene_id", all.x = TRUE)

#Merge.
genes <- merge(genes, snHuman[, c("ensembl_gene_id",
                                  "ExcExpression",
                                  "InhExpression",
                                  "AstroExpression",
                                  "MicroExpression",
                                  "OligoExpression")],
               by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all.x = TRUE)


###2###
#Add brain region gene expression
GTEx_data.path <- "./data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct"
GTEx_data <- fread(GTEx_data.path, data.table = FALSE)
Amygdala_index <- c(10)
Basal_ganglia_indices <- c(12, 19, 20, 22)
Cerebellum_indices <- c(13, 14)
Cortex_indices <- c(11, 15, 16)
Hippocampus_index <- c(17)
Hypothalamus_index <- c(18)
brain_indices <- grep("Brain", colnames(GTEx_data))
Non_brain_indices <- setdiff(c(3:56), brain_indices)

i = 1
for (i in 1:nrow(GTEx_data)) {
  Amygdala_expr <- as.numeric(GTEx_data[i,Amygdala_index])
  Basal_ganglia_median <- median(as.numeric(GTEx_data[i,Basal_ganglia_indices]))
  Cerebellum_median <- median(as.numeric(GTEx_data[i, Cerebellum_indices]))
  Cortex_median <- median(as.numeric(GTEx_data[i,Cortex_indices]))
  Hippocampus_expr <- as.numeric(GTEx_data[i,Hippocampus_index])
  Hypothalamus_expr <- as.numeric(GTEx_data[i,Hypothalamus_index])
  Non_brain_median <- median(as.numeric(GTEx_data[i,Non_brain_indices]))
  
  GTEx_data$Amygdala_RNA[i] <- log2(Amygdala_expr + 1)
  GTEx_data$Basal_ganglia_RNA[i] <- log2(Basal_ganglia_median + 1)
  GTEx_data$Cerebellum_RNA[i] <- log2(Cerebellum_median + 1)
  GTEx_data$Cortex_RNA[i] <- log2(Cortex_median + 1)
  GTEx_data$Hippocampus_RNA[i] <- log2(Hippocampus_expr + 1)
  GTEx_data$Hypothalamus_RNA[i] <- log2(Hypothalamus_expr + 1)
  GTEx_data$Non_brain_RNA[i] <- log2(Non_brain_median + 1)
}

#Used Ensembl ID converter web tool to update necessary ENSG ID to v111.
GTEx_data$Ensembl_v111 <- GTEx_data$Name
GTEx_data$Ensembl_v111 <- gsub("\\.\\d+", "", GTEx_data$Ensembl_v111)

GTEx_data <- GTEx_data %>% mutate(Ensembl_v111 = recode(Ensembl_v111, 'ENSG00000274267' = 'ENSG00000286522',
                                                        'ENSG00000278272' = 'ENSG00000287080',
                                                        'ENSG00000244693' = 'ENSG00000289604',
                                                        'ENSG00000274791' = 'ENSG00000288709',
                                                        'ENSG00000277203' = 'ENSG00000288722'))

genes <- merge(genes, GTEx_data[, c("Ensembl_v111",
                                    "Amygdala_RNA", 
                                    "Basal_ganglia_RNA",
                                    "Cerebellum_RNA",
                                    "Cortex_RNA",
                                    "Hippocampus_RNA",
                                    "Hypothalamus_RNA",
                                    "Non_brain_RNA")],
               by.x = "ensembl_gene_id", by.y = "Ensembl_v111", all.x = TRUE)


###3###
#Add iPSC developmental time course data
geneRpkm <- read.csv("./data/Burke_timecourse.csv", header = T)
cell_metadata <- read.csv( "./data/Burke_metadata.csv", header = T)
cell_metadata$SampleID <- gsub('R16-', 'R16.', cell_metadata$SampleID)
cell_metadata <- cell_metadata %>% relocate(c(LINE, DIV, REP), .after = SampleID)
cell_metadata <- arrange(cell_metadata, DIV, LINE)

#Get the average of (biological) reps from the same line, on the same day treated with the same conditions
#Beginning of accelerated dorsal period: Day 2, NPC stage: Day 15, Neural rosette: day 21, End of neuron maturation: Day 77
i = 1
for (i in 1:nrow(geneRpkm)) {
  geneRpkm$Day2_165_B_3[i] <- mean(as.numeric(geneRpkm[i,cell_metadata[c(2,3),1]]))
  geneRpkm$Day2_165_B_6X[i] <- mean(as.numeric(geneRpkm[i,cell_metadata[c(5,6),1]]))
  geneRpkm$Day2_21_B_8[i] <- mean(as.numeric(geneRpkm[i,cell_metadata[c(10,11),1]]))
  geneRpkm$Day2_66_A_3[i] <- mean(as.numeric(geneRpkm[i,cell_metadata[c(16,17),1]]))
  geneRpkm$Day2_90_A_5[i] <- mean(as.numeric(geneRpkm[i,cell_metadata[c(21,22),1]]))
  geneRpkm$Day15_165_B_6X[i] <- mean(as.numeric(geneRpkm[i,cell_metadata[c(61,62),1]]))
  geneRpkm$Day15_21_B_8[i] <- mean(as.numeric(geneRpkm[i,cell_metadata[c(65,66),1]]))
  geneRpkm$Day15_66_A_3[i] <- mean(as.numeric(geneRpkm[i,cell_metadata[c(68,69),1]]))
  geneRpkm$Day15_90_A_5[i] <- mean(as.numeric(geneRpkm[i,cell_metadata[c(72,73),1]]))
  geneRpkm$Day21_165_B_8X[i] <- mean(as.numeric(geneRpkm[i,cell_metadata[c(76,77),1]]))
  geneRpkm$Day21_66_A_9[i] <- mean(as.numeric(geneRpkm[i,cell_metadata[c(82,83),1]]))
  geneRpkm$Day21_90_A_10[i] <- mean(as.numeric(geneRpkm[i,cell_metadata[c(84,85),1]]))
  geneRpkm$Day77_165_B_6X[i] <- mean(as.numeric(geneRpkm[i,cell_metadata[c(106,107),1]]))
  geneRpkm$Day77_21_B_8[i] <- mean(as.numeric(geneRpkm[i,cell_metadata[c(109,110),1]]))
  geneRpkm$Day77_66_A_3[i] <- mean(as.numeric(geneRpkm[i,cell_metadata[c(112,113),1]]))
  geneRpkm$Day77_90_A_5[i] <- mean(as.numeric(geneRpkm[i,cell_metadata[c(116,117),1]]))
}

#Take the median and log2 of the mean (of bio reps) and others from the same time point:
i = 1
for (i in 1:nrow(geneRpkm)) {
  Day2_median <-
    median(as.numeric(geneRpkm[i, c(
      "R16.191_H5TKNBBXX",
      "R16.293_H5TKNBBXX",
      "R16.697_H7V3JBBXX",
      "R16.036_H5H27BBXX",
      "R16.855_H7VHFBBXX",
      "R16.693_H7V3JBBXX",
      "R16.033_H5H27BBXX",
      "R16.212_H5TKNBBXX",
      "R16.218_H5TKNBBXX",
      "R16.685_H7V3JBBXX",
      "R16.039_H5H27BBXX",
      "R16.042_H5H27BBXX",
      "R16.689_H7V3JBBXX",
      "Day2_165_B_3",
      "Day2_165_B_6X",
      "Day2_21_B_8",
      "Day2_66_A_3",
      "Day2_90_A_5"
    )]))
  Day15_median <-
    median(as.numeric(geneRpkm[i, c(
      "R16.853_H7V3TBBXX",
      "R16.077_H5H27BBXX",
      "R16.864_H7VHFBBXX",
      "R16.079_H5H27BBXX" ,
      "R16.083_H5H27BBXX",
      "R16.081_H5H27BBXX",
      "Day15_165_B_6X",
      "Day15_21_B_8",
      "Day15_66_A_3",
      "Day15_90_A_5"
    )]))
  Day21_median <-
    median(as.numeric(geneRpkm[i, c(
      "R16.854_H7V3TBBXX",
      "R16.100_H5TJTBBXX" ,
      "R16.865_H7VHFBBXX" ,
      "R16.098_H5TJTBBXX" ,
      "R16.096_H5TJTBBXX",
      "R16.099_H5TJTBBXX",
      "R16.097_H5TJTBBXX",
      "Day21_165_B_8X",
      "Day21_66_A_9",
      "Day21_90_A_10"
    )]))
  Day77_median <-
    median(as.numeric(geneRpkm[i, c(
      "R16.922_H7V3TBBXX",
      "R16.917_H7V3TBBXX" ,
      "R16.916_H7V3TBBXX" ,
      "R16.923_H7V3TBBXX",
      "R16.911_H7V3TBBXX",
      "Day77_165_B_6X",
      "Day77_21_B_8",
      "Day77_66_A_3",
      "Day77_90_A_5"
    )]))
  
  geneRpkm$Day2_median_RNA[i] <- log2(Day2_median + 1)
  geneRpkm$Day15_median_RNA[i] <- log2(Day15_median + 1)
  geneRpkm$Day21_median_RNA[i] <- log2(Day21_median + 1)
  geneRpkm$Day77_median_RNA[i] <- log2(Day77_median + 1)
}


#Used Ensembl ID converter web tool to update necessary ENSG ID to v111.
geneRpkm$Ensembl_v111 <- geneRpkm$rn
geneRpkm$Ensembl_v111 <- gsub("\\.\\d+", "", geneRpkm$Ensembl_v111)

geneRpkm <- geneRpkm %>% mutate(Ensembl_v111 = recode(Ensembl_v111, 'ENSG00000273547' = 'ENSG00000284662',
                                                      'ENSG00000278566' = 'ENSG00000284733',
                                                      'ENSG00000274267' = 'ENSG00000286522',
                                                      'ENSG00000278272' = 'ENSG00000287080',
                                                      'ENSG00000244693' = 'ENSG00000289604',
                                                      'ENSG00000274791' = 'ENSG00000288709',
                                                      'ENSG00000277203' = 'ENSG00000288722'))


genes <- merge(genes, geneRpkm[, c("Ensembl_v111",
                                   "Day2_median_RNA",
                                   "Day15_median_RNA",
                                   "Day21_median_RNA", 
                                   "Day77_median_RNA")],
               by.x = "ensembl_gene_id", by.y = "Ensembl_v111", all.x = TRUE)

###4###
#Add Protein expression --> the roughly 5000 reliably quantified proteins (minus the ambiguous ones)
protein_expression.path <- "./data/41593_2017_11_MOESM7_ESM_Carlyle_Table_S5.xlsx"
protein_expression <- read_excel(protein_expression.path, sheet = 1)

BrainSpanID <- read_excel("./data/BrainSpan_RNASeq_Specimen_IDs.xlsx")
BrainSpanID <- data.frame(BrainSpanID)
Carlyle_sample_metadata <- read.table("./data/Carlyle_TableS1.txt", header = T)

BrainSpan_RNA_expr <- read.csv("./data/BrainSpan_Gencode_v10_genes_matrix/expression_matrix.csv", header = F, row.names = 1)

BrainSpan_RNA_cols <- read.csv("./data/BrainSpan_Gencode_v10_genes_matrix/columns_metadata.csv", header = T)
BrainSpan_RNA_cols <- cbind(BrainSpan_RNA_cols, concat_name = NA, external_name = NA, new_col_name = NA)
i = 1
for (i in 1:nrow(BrainSpan_RNA_cols)) {
  BrainSpan_RNA_cols$concat_name[i] <- paste(BrainSpan_RNA_cols[i,3], BrainSpan_RNA_cols[i,7], sep = ".")
}

i = 1
for (i in 1:nrow(BrainSpan_RNA_cols)) {
  row_match = which(grepl(BrainSpan_RNA_cols[i,9], BrainSpanID$Specimen.name))
  if (isTRUE(length(row_match) == 1)) {
    BrainSpan_RNA_cols[i,10] <- BrainSpanID[row_match, 2]
  }
} 

i = 1
for (i in 1:nrow(BrainSpan_RNA_cols)) {
  BrainSpan_RNA_cols$new_col_name[i] <- paste(BrainSpan_RNA_cols[i,10], BrainSpan_RNA_cols[i,7], sep = "_")
}

#Set col names
colnames(BrainSpan_RNA_expr) <- BrainSpan_RNA_cols$new_col_name

#Set gene names
BrainSpan_RNA_genes <- read.csv("./data/BrainSpan_Gencode_v10_genes_matrix/rows_metadata.csv", header = T)
BrainSpan_RNA_expr <- cbind(BrainSpan_RNA_expr, BrainSpan_RNA_genes$ensembl_gene_id, BrainSpan_RNA_genes$gene_symbol, BrainSpan_RNA_genes$entrez_id)
colnames(BrainSpan_RNA_expr)[525:527] <- c("ensembl_gene_id", "gene_symbol", "entrez_id")
BrainSpan_RNA_expr <- BrainSpan_RNA_expr %>% relocate(c(ensembl_gene_id, gene_symbol, entrez_id), .before = 1)

BrainSpan_RNA_expr <- BrainSpan_RNA_expr[,-c(grep("^NA_", colnames(BrainSpan_RNA_expr)))]

BrainSpan_RNA_expr$CBC_sum <- rowSums(BrainSpan_RNA_expr[,c(grep("_CBC$", colnames(BrainSpan_RNA_expr)))])
BrainSpan_RNA_expr$MD_sum <- rowSums(BrainSpan_RNA_expr[,c(grep("_MD$", colnames(BrainSpan_RNA_expr)))])
BrainSpan_RNA_expr$STR_sum <- rowSums(BrainSpan_RNA_expr[,c(grep("_STR$", colnames(BrainSpan_RNA_expr)))])
BrainSpan_RNA_expr$AMY_sum <- rowSums(BrainSpan_RNA_expr[,c(grep("_AMY$", colnames(BrainSpan_RNA_expr)))])
BrainSpan_RNA_expr$HIP_sum <- rowSums(BrainSpan_RNA_expr[,c(grep("_HIP$", colnames(BrainSpan_RNA_expr)))])
BrainSpan_RNA_expr$V1C_sum <- rowSums(BrainSpan_RNA_expr[,c(grep("_V1C$", colnames(BrainSpan_RNA_expr)))])
BrainSpan_RNA_expr$DFC_sum <- rowSums(BrainSpan_RNA_expr[,c(grep("DFC$", colnames(BrainSpan_RNA_expr)))])

#Update protein ENSG.
protein_expression$Ensembl_v111 <- protein_expression$EnsemblID
protein_expression$Ensembl_v111 <- gsub("\\.\\d+", "", protein_expression$Ensembl_v111)

protein_expression <- protein_expression %>% mutate(Ensembl_v111 = recode(Ensembl_v111, 'ENSG00000274267' = 'ENSG00000286522',
                                                                          'ENSG00000274791' = 'ENSG00000288709',
                                                                          'ENSG00000277203' = 'ENSG00000288722',
                                                                          'ENSG00000278272' = 'ENSG00000287080'))

###Update necessary RNA ENSG ID. 
BrainSpan_RNA_expr <- BrainSpan_RNA_expr %>% mutate(ensembl_gene_id = recode(ensembl_gene_id, 'ENSG00000124529' = 'ENSG00000278705', 
                                                                             'ENSG00000170616' = 'ENSG00000261678',
                                                                             'ENSG00000170727' = 'ENSG00000261236',
                                                                             'ENSG00000172660' = 'ENSG00000270647',
                                                                             'ENSG00000174595' = 'ENSG00000266265',
                                                                             'ENSG00000188987' = 'ENSG00000277157',
                                                                             'ENSG00000196176' = 'ENSG00000278637',
                                                                             'ENSG00000197914' = 'ENSG00000273542',
                                                                             'ENSG00000198327' = 'ENSG00000274618',
                                                                             'ENSG00000198558' = 'ENSG00000275126',
                                                                             'ENSG00000221995' = 'ENSG00000196535',
                                                                             'ENSG00000227184' = 'ENSG00000261150',
                                                                             'ENSG00000229583' = 'ENSG00000268994'))

#Remove any rows with ambiguous protein mapping. 
protein_expression <- protein_expression[!grepl(";",protein_expression$GeneSymbol),]


#Grep in the proteins by ENSG. 
ens_gene = c()
i = 1
for (i in 1:nrow(protein_expression)) {
  ens_id = protein_expression[i,21]
  row_match = which(grepl(ens_id, genes$ensembl_gene_id))
  if (length(row_match) > 1) {
    print(ens_id)
  }
  if (!identical(row_match, integer(0))) {
    genes[row_match, c(22:28)] <- protein_expression[i,7:13]
  }
} 


###If = NA, check if RNA sum = 0. If yes, then add a 0. Do for each brain region.
#0 RNA in CBC:
i = 1
for (i in 1:nrow(genes)){
  if (is.na(genes[i,22])) {
    CBC_RNA = which(grepl(genes[i,1], BrainSpan_RNA_expr$ensembl_gene_id))
    if (length(CBC_RNA) > 1) {
      print(CBC_RNA)
    }
    if (isTRUE(BrainSpan_RNA_expr[CBC_RNA, 217] == 0)) {
      genes[i, 22] <- 0
    }
  }
}
sum(is.na(genes$`log10(CBC)`)) #13573

#0 RNA in MD:
i = 1
for (i in 1:nrow(genes)){
  if (is.na(genes[i,23])) {
    MD_RNA = which(grepl(genes[i,1], BrainSpan_RNA_expr$ensembl_gene_id))
    if (length(MD_RNA) > 1) {
      print(MD_RNA)
    }
    if (isTRUE(BrainSpan_RNA_expr[MD_RNA, 218] == 0)) {
      genes[i, 23] <- 0
    }
  }
}
sum(is.na(genes$`log10(MD)`)) #13603


#0 RNA in STR:
i = 1
for (i in 1:nrow(genes)){
  if (is.na(genes[i,24])) {
    STR_RNA = which(grepl(genes[i,1], BrainSpan_RNA_expr$ensembl_gene_id))
    if (length(STR_RNA) > 1) {
      print(STR_RNA)
    }
    if (isTRUE(BrainSpan_RNA_expr[STR_RNA, 219] == 0)) {
      genes[i, 24] <- 0
    }
  }
}
sum(is.na(genes$`log10(STR)`)) #13445

#0 RNA in AMY:
i = 1
for (i in 1:nrow(genes)){
  if (is.na(genes[i,25])) {
    AMY_RNA = which(grepl(genes[i,1], BrainSpan_RNA_expr$ensembl_gene_id))
    if (length(AMY_RNA) > 1) {
      print(AMY_RNA)
    }
    if (isTRUE(BrainSpan_RNA_expr[AMY_RNA, 220] == 0)) {
      genes[i, 25] <- 0
    }
  }
}
sum(is.na(genes$`log10(AMY)`)) #13773

#0 RNA in HIP:
i = 1
for (i in 1:nrow(genes)){
  if (is.na(genes[i,26])) {
    HIP_RNA = which(grepl(genes[i,1], BrainSpan_RNA_expr$ensembl_gene_id))
    if (length(HIP_RNA) > 1) {
      print(HIP_RNA)
    }
    if (isTRUE(BrainSpan_RNA_expr[HIP_RNA, 221] == 0)) {
      genes[i, 26] <- 0
    }
  }
}
sum(is.na(genes$`log10(HIP)`)) #13668

#0 RNA in V1C:
i = 1
for (i in 1:nrow(genes)){
  if (is.na(genes[i,27])) {
    V1C_RNA = which(grepl(genes[i,1], BrainSpan_RNA_expr$ensembl_gene_id))
    if (length(V1C_RNA) > 1) {
      print(V1C_RNA)
    }
    if (isTRUE(BrainSpan_RNA_expr[V1C_RNA, 222] == 0)) {
      genes[i, 27] <- 0
    }
  }
}
sum(is.na(genes$`log10(V1C)`)) #13553

#0 RNA in DFC:
i = 1
for (i in 1:nrow(genes)){
  if (is.na(genes[i,28])) {
    DFC_RNA = which(grepl(genes[i,1], BrainSpan_RNA_expr$ensembl_gene_id))
    if (length(DFC_RNA) > 1) {
      print(DFC_RNA)
    }
    if (isTRUE(BrainSpan_RNA_expr[DFC_RNA, 223] == 0)) {
      genes[i, 28] <- 0
    }
  }
}
sum(is.na(genes$`log10(DFC)`)) #13535


###5###
#Mutational constraint (LOEUF)
LOEUF.path <- "./data/gnomad.v2.1.1.lof_metrics.by_gene.txt"
LOEUF <- fread(LOEUF.path, data.table = F)
LOEUF_no_dup <- LOEUF[-c(11645, 7027, 4975, 19208, 5970,
                         2067, 7784, 10148, 7528, 2640,
                         18271, 15770, 19490, 15735, 2990,
                         14210, 16427, 11384, 5939, 13660,
                         5310, 8785, 18617, 19570, 15149,
                         11652, 11325, 2672, 12947, 8351,
                         6237, 10921, 16410, 11992, 12743,
                         17709, 15573, 18069, 15244, 7007), ] # remove duplicated gene rows

LOEUF_no_dup <- LOEUF_no_dup %>% mutate(gene_id = recode(gene_id,  'ENSG00000170616' = 'ENSG00000261678',
                                                         'ENSG00000170727' = 'ENSG00000261236',
                                                         'ENSG00000227184' = 'ENSG00000261150',
                                                         'ENSG00000124529' = 'ENSG00000278705',
                                                         'ENSG00000174595' = 'ENSG00000266265',
                                                         'ENSG00000188987' = 'ENSG00000277157',
                                                         'ENSG00000196176' = 'ENSG00000278637',
                                                         'ENSG00000197914' = 'ENSG00000273542',
                                                         'ENSG00000198327' = 'ENSG00000274618',
                                                         'ENSG00000198558' = 'ENSG00000275126',
                                                         'ENSG00000255800' = 'ENSG00000277556', 
                                                         'ENSG00000257019' = 'ENSG00000276119'))

genes <- merge(genes, LOEUF_no_dup[, c("gene_id", "oe_lof_upper")], by.x = "ensembl_gene_id", by.y = "gene_id", all.x = TRUE) 

###6###
#Add synapse localization data
synapse.path <-"./data/syngo_genes_release_2023.xlsx"
synapse <- read_excel(synapse.path)
synapse <- data.frame(synapse)

genes$Synapse <- as.numeric(genes$ensembl_gene_id %in% synapse$ensembl_id)

###7###
#ASD susceptibility genes
sfari.path <- "./data/SFARI-Gene_genes_01-16-2024release_02-25-2024export.csv"
sfari <- read.csv(sfari.path, header = T, na.strings=c("","NA"))

#Manually add missing ENSG IDs. 
sfari[13,4] <- "ENSG00000196839"
sfari[72,4] <- "ENSG00000169083"
sfari[603,4] <- "ENSG00000105976"

genes$Autism_susceptibility <- as.numeric(genes$ensembl_gene_id %in% sfari$ensembl.id)

#Add additional MSSNG (From Trost et al.) ASD genes
MSSNG <- read_excel("./data/1-s2.0-S0092867422013241-mmc2_MSSNG.xlsx", sheet = 5, skip = 1)
MSSNG[7,1] <- "KMT5B"
MSSNG[49,1] <- "SRPRA"

gene_df_MSSNG <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype", "entrezgene_id"),
  filters = c("chromosome_name", "external_gene_name", "biotype"),
  values = list(c(1:22, "X", "Y", "MT"), MSSNG$Gene, c("protein_coding", "lncRNA")),
  mart = mart
)

MSSNG <- merge(MSSNG, gene_df_MSSNG[, c("ensembl_gene_id",
                                        "external_gene_name")],
               by.x = "Gene", by.y = "external_gene_name", all.x = TRUE)

#remove rows if they're in SFARI
MSSNG <- MSSNG[!(MSSNG$ensembl_gene_id %in% sfari$ensembl.id),] 

ens_gene = c()
i = 1
for (i in 1:nrow(MSSNG)) {
  ens_id = MSSNG[i,14]
  row_match = which(grepl(ens_id, genes$ensembl_gene_id))
  if (length(row_match) > 1) { print(ens_id) }
  if (!identical(row_match, integer(0))) {
    genes[row_match, 31] <- 1
  } 
}


#Prep final df
names(genes) <- c("Ensembl ID",
                  "HGNC ID",
                  "HGNC status",
                  "Gene name",
                  "Gene description", 
                  "Excitatory Expression",
                  "Inhibitory Expression",
                  "Astrocyte Expression", 
                  "Microglia Expression", 
                  "Oligodendrocyte Expression",
                  "Amygdala expression",
                  "Basal ganglia expression",
                  "Cerebellum expression",
                  "Cortex expression",
                  "Hippocampus expression",
                  "Hypothalamus expression",
                  "Non-brain expression",
                  "Accelerated dorsal expression (D2)",
                  "NPC expression (D15)",
                  "Neural rosette expression (D21)",
                  "Neuron expression (D77)",
                  "CBC protein expression", 
                  "MD protein expression",
                  "STR protein expression",
                  "AMY protein expression",
                  "HIP protein expression",
                  "V1C protein expression",
                  "DFC protein expression",
                  "LOEUF",
                  "Synapse localization",
                  "Autism susceptibility")

genes <- impute(genes, "mean")


#Export.
outfile.path <- "./model_input/PC/PC_genes_input.tsv"
write.table(genes, outfile.path, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
