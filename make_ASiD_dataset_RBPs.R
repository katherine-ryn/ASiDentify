# ---
# Title: make_ASiD_dataset_RBPs.R
# Purpose: This script generates the columns in ASiD from various databases for all 4 lists of RBPs.
#          
# ---

library(data.table)
library(dplyr)
library(biomaRt)
library(readxl)
library(stringr)


setwd("") # Change this to the path that contains the ASiD folder

rbps.path <- "./RBP_lists/Comp_RBPs.csv"
rbps <- fread(rbps.path, data.table = FALSE)
names(rbps) <- gsub("\\.| ", "_", names(rbps))


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
  mart = mart)

#Only keep unique Ensembl rows. 
gene_df_distinct <- distinct_at(gene_df_snHuman_NCBI, vars(ensembl_gene_id,
                                                           external_gene_name,
                                                           ensembl_gene_id_version,
                                                           gene_biotype),
                                .keep_all = TRUE)

snHuman <- merge(snHuman, gene_df_distinct[, c("ensembl_gene_id", "entrezgene_id")],
                 by.x = "entrez_id", by.y = "entrezgene_id", all.x = TRUE)

#Using grep to match (not merge) because 'rbps' df has >1 Ensembl ID per row
rbps <- cbind(rbps, ExcExpression = NA,
              InhExpression = NA,
              AstroExpression = NA,
              MicroExpression = NA,
              OligoExpression = NA)
ens_gene = c()
i = 1
for (i in 1:nrow(snHuman)) {
  ens_id = snHuman[i,129]
  row_match = which(grepl(ens_id, rbps$Ensembl_ID))
  if (length(row_match) > 1) {
    print(ens_id)
  }
  if (!identical(row_match, integer(0))) {
    rbps[row_match, c(8:12)] <- snHuman[i,124:128]
  }
} 

###2###
#Add brain region gene expression.
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

rbps <- cbind(rbps, Amygdala_RNA = NA,
              Basal_ganglia_RNA = NA,
              Cerebellum_RNA = NA, 
              Cortex_RNA = NA, 
              Hippocampus_RNA = NA,
              Hypothalamus_RNA = NA,
              Non_brain_RNA = NA)

ens_gene = c()
i = 1
for (i in 1:nrow(GTEx_data)) {
  ens_id = GTEx_data[i,64]
  row_match = which(grepl(ens_id, rbps$Ensembl_ID))
  if (length(row_match) > 1) {
    print(ens_id)
  }
  if (!identical(row_match, integer(0))) {
    rbps[row_match, c(13:19)] <- GTEx_data[i,57:63]
  }
} 
sum(is.na(rbps$Non_brain_RNA)) #7

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
      "Day2_90_A_5")]))
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
      "Day15_90_A_5")]))
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
      "Day21_90_A_10")]))
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
      "Day77_90_A_5")]))
  
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
                                                     'ENSG00000239665' = 'ENSG00000165630', 
                                                     'ENSG00000244693' = 'ENSG00000289604', 
                                                     'ENSG00000274791' = 'ENSG00000288709', 
                                                     'ENSG00000277203' = 'ENSG00000288722'))


rbps <- cbind(rbps, Day2_median_RNA = NA, Day15_median_RNA = NA, Day21_median_RNA = NA, Day77_median_RNA = NA)
    
ens_gene = c()
i = 1
for (i in 1:nrow(geneRpkm)) {
   ens_id = geneRpkm[i,139]
   row_match = which(grepl(ens_id, rbps$Ensembl_ID))
   if (length(row_match) > 1) {
     print(ens_id)
   }
   if (!identical(row_match, integer(0))) {
     rbps[row_match, c(20:23)] <- geneRpkm[i,135:138]
   }
} 
sum(is.na(rbps$Day2_median_RNA)) #0



###4###
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

LOEUF_no_dup = LOEUF_no_dup %>% mutate(gene_id = recode(gene_id,  'ENSG00000170616' = 'ENSG00000261678',
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

ens_gene = c()
i = 1
for (i in 1:nrow(LOEUF_no_dup)) {
  ens_id = LOEUF_no_dup[i,64]
  row_match = which(grepl(ens_id, rbps$Ensembl_ID))
  if (length(row_match) > 1) {
    print(ens_id)
  }
  if (!identical(row_match, integer(0))) {
    rbps[row_match, c(24)] <- LOEUF_no_dup[i,30]
  }
} 
colnames(rbps)[24] <- c("LOEUF")
sum(is.na(rbps$LOEUF)) #60

###5###
#Add synapse localization data
synapse.path <-"./data/syngo_genes_release_2023.xlsx"
synapse <- read_excel(synapse.path)
synapse <- data.frame(synapse)

ens_gene = c()
i = 1
for (i in 1:nrow(synapse)) {
  ens_id = synapse[i,5]
  row_match = which(grepl(ens_id, rbps$Ensembl_ID))
  if (length(row_match) > 1) { print(ens_id) }
  if (!identical(row_match, integer(0))) {
    rbps[row_match, c(25)] <- 1
  } 
}
colnames(rbps)[25] <- c("Synapse")
rbps <- rbps %>% mutate(Synapse = ifelse(is.na(Synapse), 0, Synapse))


###6###
#ASD susceptibility genes
sfari.path <- "./data/SFARI-Gene_genes_01-16-2024release_02-25-2024export.csv"
sfari <- read.csv(sfari.path, header = T, na.strings=c("","NA"))

#Add SFARI genes. Including any gene which has Score 1,2 or 3 or Syndromic.
ens_gene = c()
i = 1
for (i in 1:nrow(sfari)) {
  ens_id = sfari[i,4]
  row_match = which(grepl(ens_id, rbps$Ensembl_ID))
    if (length(row_match) > 1) { print(ens_id) }
    if (!identical(row_match, integer(0))) {
      rbps[row_match, c(26)] <- 1
    } 
}
colnames(rbps)[26] <- c("Autism_susceptibility") 

#Add additional MSSNG (from Trost et al.) ASD genes
MSSNG <- read_excel("./data/1-s2.0-S0092867422013241-mmc2_MSSNG.xlsx", sheet = 5, skip = 1)
MSSNG[7,1] <- "KMT5B"
MSSNG[49,1] <- "SRPRA"

gene_df_MSSNG <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype", "entrezgene_id"),
  filters = c("chromosome_name", "external_gene_name", "biotype"),
  values = list(c(1:22, "X", "Y", "MT"), MSSNG$Gene, c("protein_coding", "lncRNA")),
  mart = mart)

MSSNG <- merge(MSSNG, gene_df_MSSNG[, c("ensembl_gene_id", "external_gene_name")],
               by.x = "Gene", by.y = "external_gene_name", all.x = TRUE)

#remove rows if they're in SFARI
MSSNG <- MSSNG[!(MSSNG$ensembl_gene_id %in% sfari$ensembl.id),] 

ens_gene = c()
i = 1
for (i in 1:nrow(MSSNG)) {
  ens_id = MSSNG[i,14]
  row_match = which(grepl(ens_id, rbps$Ensembl_ID))
  if (length(row_match) > 1) { print(ens_id) }
  if (!identical(row_match, integer(0))) {
    rbps[row_match, c(26)] <- 1
  } 
}

rbps <- rbps %>% mutate(Autism_susceptibility = ifelse(is.na(Autism_susceptibility), 0, Autism_susceptibility))

###Prep final df
names(rbps) <- c("Gene name",
                 "Ensembl ID",
                 "Gene description",
                 "Canonical RBD", 
                 "Pfam ID", 
                 "Protein Domains",
                 "Pfam Description",
                 "Excitatory Expression",
                 "Inhibitory Expression",
                 "Astrocyte Expression", 
                 "Microglia Expression",
                 "Oligodendrocyte expression",
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
                 "LOEUF", 
                 "Synapse localization", 
                 "Autism susceptibility")


###Make ASiD table for all 4 RBP lists. Export. 
filenames <- list.files("./RBP_lists", pattern="*.csv", full.names = T)
ldf <- lapply(filenames, read.csv)
df_names <- gsub("./RBP_lists/|.csv", "", filenames)
names(ldf) <- df_names

i = 1
for(i in names(ldf)){
  ldf[[i]] <- merge(ldf[[i]], rbps[, c("Gene name",
                                       "Excitatory Expression",
                                       "Inhibitory Expression",
                                       "Astrocyte Expression",
                                       "Microglia Expression",
                                       "Oligodendrocyte expression",
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
                                       "LOEUF",
                                       "Synapse localization",
                                       "Autism susceptibility")],
                    by.x = "Gene.name", by.y = "Gene name", all.x = TRUE)
  
  names(ldf[[i]]) <- c("Gene name",
                       "Ensembl ID",
                       "Gene description",
                       "Canonical RBD",
                       "Pfam ID", 
                       "Protein Domains",
                       "Pfam Description",
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
                       "LOEUF", 
                       "Synapse localization",
                       "Autism susceptibility")
  ldf[[i]] <-  ldf[[i]] %>% filter(!(if_any(`Excitatory Expression`:`LOEUF`, is.na))) 
  
}

#Export
i = 1
for(i in names(ldf)){
   write.table(ldf[[i]], str_glue("./model_input/RBP/{i}.tsv"), 
               quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}

