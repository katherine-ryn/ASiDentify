# ---
# Title: nested_cv_random.R
#   For running 10 random sets of 819 genes.  
# ---

library(data.table)
library(dplyr)
library(glmnet)
library(ROCR)
library(stringr)

setwd("") # Change this to the path that contains the ASiD folder
source("./nested_cv_funcs.R")

genes <- fread("./model_input/PC/PC_genes_input.tsv", data.table = FALSE) 

#Randomly pick 10 sets of 82 ASD genes and 719 non-ASD genes.
set.seed(1567)
genes_set1_ASD <- genes %>% filter(`Autism susceptibility` == 1) %>% slice_sample(n = 82, replace = FALSE)
genes_set1_nonASD <- genes %>% filter(`Autism susceptibility` == 0) %>% slice_sample(n = 719, replace = FALSE)
genes_set1 <- genes_set1_ASD %>% bind_rows(genes_set1_nonASD)

set.seed(15)
genes_set2_ASD <- genes %>% filter(`Autism susceptibility` == 1) %>% slice_sample(n = 82, replace = FALSE)
genes_set2_nonASD <- genes %>% filter(`Autism susceptibility` == 0) %>% slice_sample(n = 719, replace = FALSE)
genes_set2 <- genes_set2_ASD %>% bind_rows(genes_set2_nonASD)

set.seed(159)
genes_set3_ASD <- genes %>% filter(`Autism susceptibility` == 1) %>% slice_sample(n = 82, replace = FALSE)
genes_set3_nonASD <- genes %>% filter(`Autism susceptibility` == 0) %>% slice_sample(n = 719, replace = FALSE)
genes_set3 <- genes_set3_ASD %>% bind_rows(genes_set3_nonASD)

set.seed(478)
genes_set4_ASD <- genes %>% filter(`Autism susceptibility` == 1) %>% slice_sample(n = 82, replace = FALSE)
genes_set4_nonASD <- genes %>% filter(`Autism susceptibility` == 0) %>% slice_sample(n = 719, replace = FALSE)
genes_set4 <- genes_set4_ASD %>% bind_rows(genes_set4_nonASD)

set.seed(2677)
genes_set5_ASD <- genes %>% filter(`Autism susceptibility` == 1) %>% slice_sample(n = 82, replace = FALSE)
genes_set5_nonASD <- genes %>% filter(`Autism susceptibility` == 0) %>% slice_sample(n = 719, replace = FALSE)
genes_set5 <- genes_set5_ASD %>% bind_rows(genes_set5_nonASD)

set.seed(984)
genes_set6_ASD <- genes %>% filter(`Autism susceptibility` == 1) %>% slice_sample(n = 82, replace = FALSE)
genes_set6_nonASD <- genes %>% filter(`Autism susceptibility` == 0) %>% slice_sample(n = 719, replace = FALSE)
genes_set6 <- genes_set6_ASD %>% bind_rows(genes_set6_nonASD)

set.seed(502)
genes_set7_ASD <- genes %>% filter(`Autism susceptibility` == 1) %>% slice_sample(n = 82, replace = FALSE)
genes_set7_nonASD <- genes %>% filter(`Autism susceptibility` == 0) %>% slice_sample(n = 719, replace = FALSE)
genes_set7 <- genes_set7_ASD %>% bind_rows(genes_set7_nonASD)

set.seed(111)
genes_set8_ASD <- genes %>% filter(`Autism susceptibility` == 1) %>% slice_sample(n = 82, replace = FALSE)
genes_set8_nonASD <- genes %>% filter(`Autism susceptibility` == 0) %>% slice_sample(n = 719, replace = FALSE)
genes_set8 <- genes_set8_ASD %>% bind_rows(genes_set8_nonASD)

set.seed(345)
genes_set9_ASD <- genes %>% filter(`Autism susceptibility` == 1) %>% slice_sample(n = 82, replace = FALSE)
genes_set9_nonASD <- genes %>% filter(`Autism susceptibility` == 0) %>% slice_sample(n = 719, replace = FALSE)
genes_set9 <- genes_set9_ASD %>% bind_rows(genes_set9_nonASD)

set.seed(6101)
genes_set10_ASD <- genes %>% filter(`Autism susceptibility` == 1) %>% slice_sample(n = 82, replace = FALSE)
genes_set10_nonASD <- genes %>% filter(`Autism susceptibility` == 0) %>% slice_sample(n = 719, replace = FALSE)
genes_set10 <- genes_set10_ASD %>% bind_rows(genes_set10_nonASD)


set_list <- list(genes_set1, genes_set2, genes_set3, genes_set4, genes_set5, 
                 genes_set6, genes_set7, genes_set8, genes_set9, genes_set10)
i = 1
for (i in 1:10) {
  autism_res <- trainCV(set_list[[i]], K = 5, indep_var = "autism", standardized_col = c(6:22), feature_col = c(6:23))
  autism_pred <- do.call("rbind", autism_res$predictions)
  autism_pred <- autism_pred %>% arrange(-predictions)
  
  
  assign(paste0(str_glue("genes_set{i}"), "_autism_res"), autism_res)
  assign(paste0(str_glue("genes_set{i}"), "_autism_pred"), autism_pred)
}

#Export in analyzing_random_output.R


