# ---
# Title: nested_cv.R
# Purpose: This script generates the columns in ASiD from various databases.
# ---

library(data.table)
library(dplyr)
library(glmnet)
library(ROCR)

setwd("") # Change this to the path that contains the ASiD folder
source("./nested_cv_funcs.R")

#Each model run is below:

###1### -----------------------------------------------------------------
#All protein-coding genes:
genes <- fread("./model_input/PC/PC_genes_input.tsv", data.table = FALSE) 


# Train machine learning model
autism_res <- trainCV(genes, K = 5, indep_var = "autism", standardized_col = c(6:22), feature_col = c(6:23))
autism_pred <- do.call("rbind", autism_res$predictions)
autism_pred <- autism_pred %>% arrange(-predictions)



#Export table with original (unstandardized) values. 
genes <- merge(genes, autism_pred[, c("Gene name", "predictions")],
               by.x = "Gene name", by.y = "Gene name", all.x = TRUE)
genes <- genes %>% arrange(-predictions)
genes <- genes %>% relocate(c(`Gene name`), .after = `HGNC status`)

#Export in analyzing_PC_output.R


###2### -----------------------------------------------------------------
#For the 4 RBP models:
#For autism_res_CompRBP
CompRBP <- fread("./model_input/RBP/Comp_RBPs.tsv", data.table = FALSE)  
autism_res_CompRBP <- trainCV(CompRBP, K = 5, indep_var = "autism", standardized_col = c(8:24), feature_col = c(8:25))

autism_pred_CompRBP <- do.call("rbind", autism_res_CompRBP$predictions)
autism_pred_CompRBP <- autism_pred_CompRBP %>% arrange(-predictions)

#Export table with original (unstandardized) values. 
CompRBP <- merge(CompRBP, autism_pred_CompRBP[, c("Gene name", "predictions")],
               by.x = "Gene name", by.y = "Gene name", all.x = TRUE)
CompRBP <- CompRBP %>% arrange(-predictions)



#For autism_res_HCRBP
HCRBP <- fread("./model_input/RBP/Highconfidence_RBPs.tsv", data.table = FALSE)  
autism_res_HCRBP <- trainCV(HCRBP, K = 5, indep_var = "autism", standardized_col = c(8:24), feature_col = c(8:25))

autism_pred_HCRBP <- do.call("rbind", autism_res_HCRBP$predictions)
autism_pred_HCRBP <- autism_pred_HCRBP %>% arrange(-predictions)

#Export table with original (unstandardized) values. 
HCRBP <- merge(HCRBP, autism_pred_HCRBP[, c("Gene name", "predictions")],
                by.x = "Gene name", by.y = "Gene name", all.x = TRUE)
HCRBP <- HCRBP %>% arrange(-predictions)


#For autism_res_CRBP
CRBP <- fread("./model_input/RBP/Canonical_RBPs.tsv", data.table = FALSE)  
autism_res_CRBP <- trainCV(CRBP, K = 5, indep_var = "autism", standardized_col = c(8:24), feature_col = c(8:25))

autism_pred_CRBP <- do.call("rbind", autism_res_CRBP$predictions)
autism_pred_CRBP <- autism_pred_CRBP %>% arrange(-predictions)

#Export table with original (unstandardized) values. 
CRBP <- merge(CRBP, autism_pred_CRBP[, c("Gene name", "predictions")],
                by.x = "Gene name", by.y = "Gene name", all.x = TRUE)
CRBP <- CRBP %>% arrange(-predictions)



#For autism_res_NCRBP
NCRBP <- fread("./model_input/RBP/Noncanonical_RBPs.tsv", data.table = FALSE)  
autism_res_NCRBP <- trainCV(NCRBP, K = 5, indep_var = "autism", standardized_col = c(8:24), feature_col = c(8:25))

autism_pred_NCRBP <- do.call("rbind", autism_res_NCRBP$predictions)
autism_pred_NCRBP <- autism_pred_NCRBP %>% arrange(-predictions)

#Export table with original (unstandardized) values. 
NCRBP <- merge(NCRBP, autism_pred_NCRBP[, c("Gene name", "predictions")],
                by.x = "Gene name", by.y = "Gene name", all.x = TRUE)
NCRBP <- NCRBP %>% arrange(-predictions)


###3### -----------------------------------------------------------------
#Chromatin regulators:
genes <- fread("./model_input/CR/CR_input.tsv", data.table = FALSE)  

# Train machine learning model
autism_res <- trainCV(genes, K = 5, indep_var = "autism", standardized_col = c(20:36), feature_col = c(20:37))
autism_pred <- do.call("rbind", autism_res$predictions)
autism_pred <- autism_pred %>% arrange(-predictions)

#Export table with original (unstandardized) values. 
genes <- merge(genes, autism_pred[, c("Gene name", "predictions")],
               by.x = "Gene name", by.y = "Gene name", all.x = TRUE)
genes <- genes %>% arrange(-predictions)
genes <- genes %>% relocate(c(`Gene name`), .after = `HGNC ID`)

#Export in analyzing_CR_output.R


