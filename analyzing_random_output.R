# ---
# Title: analyzing_random_output.R
# Purpose: This script uses the nested_cv_random.R output from 10 set of random 801 genes
#           and generates the odds ratios. 
# ---

library(ggplot2)
library(scales)
library(cowplot)
library(ggbeeswarm)
library(openxlsx)


setwd("") # Change this to the path that contains the ASiD folder


# SETUP -------------------------------------------------------------------
#Extract beta coefficients from the 5 folds. 
get_beeswarm_tb <- function(res, K) {
  
  df <- as.data.frame(res[["betaCoefs"]][1])
  names(df) <- "fold1"
  df$variable <- rownames(df)
  df <- df[, c(2, 1)]
  
  i = 2
  for (i in 2:K) {
    df2 <- as.data.frame(res[["betaCoefs"]][i])
    names(df2) <- paste0("fold", i)
    df2$variable <- rownames(df2)
    df2 <- df2[, c(2:1)]
    df <- merge(df, df2, 
                by = "variable",
                all = T)
  }
  
  df[is.na(df)] <- 0
  
  coefs <- numeric((nrow(df)-1) * 5)
  variables <- character((nrow(df)-1) * 5)
  
  i = 2
  for (i in (2:nrow(df))) {
    ifelse(i == 2, 
           coefs <- as.numeric(df[c(i), c(2:ncol(df))]),
           coefs <- c(coefs, as.numeric(df[c(i), c(2:ncol(df))])))
    
    ifelse(i == 2,
           variables <- rep(df[i, 1], 5),
           variables <- c(variables, rep(df[i, 1], 5)))
  }
  
  tb <- data.frame(coef = coefs,
                   variable = variables)
  
  return(tb)
}

random_names <- c("set1", "set2", "set3", "set4", "set5", "set6", "set7", "set8", "set9", "set10")
random_runs <- list(genes_set1_autism_res = genes_set1_autism_res,
                    genes_set2_autism_res = genes_set2_autism_res,
                    genes_set3_autism_res = genes_set3_autism_res,
                    genes_set4_autism_res = genes_set4_autism_res,
                    genes_set5_autism_res = genes_set5_autism_res,
                    genes_set6_autism_res = genes_set6_autism_res,
                    genes_set7_autism_res = genes_set7_autism_res,
                    genes_set8_autism_res = genes_set8_autism_res,
                    genes_set9_autism_res = genes_set9_autism_res,
                    genes_set10_autism_res = genes_set10_autism_res)

# Run nested_cv.R to get autism_res1-10.
n = 1
for (n in 1:length(random_names)) {
  autism_tb <- get_beeswarm_tb(res = random_runs[[n]], K = 5)
  
  autism_tb$OR <- exp(autism_tb$coef)
  
  autism_CIs <- autism_tb %>% group_by(variable) %>% summarize(Mean = mean(OR)) %>% ungroup()
  
  groups <- seq(1, nrow(autism_CIs)*5, by = 5)
  CI_low <- numeric(0)
  CI_high <- numeric(0)
  
  for (i in groups) {
    members <- autism_tb$OR[i:(i+5-1)]
    if (sum(members) == 5) { #This is if they are all = 1 OR
      CI_low <- c(CI_low, 0)
      CI_high <- c(CI_high, 0)
      next
    }
    ci <- as.numeric(t.test(members, conf.level = 0.95)$conf.int)
    CI_low <- c(CI_low, ci[1])
    CI_high <- c(CI_high, ci[2])
  }
  
  
  autism_CIs$CI_low <- CI_low
  autism_CIs$CI_high <- CI_high
  
  
  autism_CIs <- as.data.frame(autism_CIs)
  autism_CIs$group <- random_names[n]
  autism_tb$group <- random_names[n]
  autism_CIs$group <- factor(autism_CIs$group, levels = random_names[n])
  autism_tb$group <- factor(autism_tb$group, levels = random_names[n])
  
  assign(paste0(random_names[n], "_autism_tb"), autism_tb)
  assign(paste0(random_names[n], "_autism_CIs"), autism_CIs)
}


#Export data.
wb <- createWorkbook()
addWorksheet(wb = wb, sheetName = "CI Random Set1", gridLines = T)
writeDataTable(wb = wb, sheet = 1, x = set1_autism_CIs, withFilter = F, tableStyle = "None")
addWorksheet(wb = wb, sheetName = "CI Random Set2", gridLines = T)
writeData(wb = wb, sheet = 2, x = set2_autism_CIs)
addWorksheet(wb = wb, sheetName = "CI Random Set3", gridLines = T)
writeData(wb = wb, sheet = 3, x = set3_autism_CIs)
addWorksheet(wb = wb, sheetName = "CI Random Set4", gridLines = T)
writeData(wb = wb, sheet = 4, x = set4_autism_CIs)
addWorksheet(wb = wb, sheetName = "CI Random Set5", gridLines = T)
writeData(wb = wb, sheet = 5, x = set5_autism_CIs)
addWorksheet(wb = wb, sheetName = "CI Random Set6", gridLines = T)
writeData(wb = wb, sheet = 6, x = set6_autism_CIs)
addWorksheet(wb = wb, sheetName = "CI Random Set7", gridLines = T)
writeData(wb = wb, sheet = 7, x = set7_autism_CIs)
addWorksheet(wb = wb, sheetName = "CI Random Set8", gridLines = T)
writeData(wb = wb, sheet = 8, x = set8_autism_CIs)
addWorksheet(wb = wb, sheetName = "CI Random Set9", gridLines = T)
writeData(wb = wb, sheet = 9, x = set9_autism_CIs)
addWorksheet(wb = wb, sheetName = "CI Random Set10", gridLines = T)
writeData(wb = wb, sheet = 10, x = set10_autism_CIs)

saveWorkbook(wb, "./tables/Random_output.xlsx", overwrite = TRUE)

#Also save as csv.
write.csv(set1_autism_CIs, "./tables/csv/random_set1_CIs.csv", 
          row.names = F)
write.csv(set2_autism_CIs, "./tables/csv/random_set2_CIs.csv", 
          row.names = F)
write.csv(set3_autism_CIs, "./tables/csv/random_set3_CIs.csv", 
          row.names = F)
write.csv(set4_autism_CIs, "./tables/csv/random_set4_CIs.csv", 
          row.names = F)
write.csv(set5_autism_CIs, "./tables/csv/random_set5_CIs.csv", 
          row.names = F)
write.csv(set6_autism_CIs, "./tables/csv/random_set6_CIs.csv", 
          row.names = F)
write.csv(set7_autism_CIs, "./tables/csv/random_set7_CIs.csv", 
          row.names = F)
write.csv(set8_autism_CIs, "./tables/csv/random_set8_CIs.csv", 
          row.names = F)
write.csv(set9_autism_CIs, "./tables/csv/random_set9_CIs.csv", 
          row.names = F)
write.csv(set10_autism_CIs, "./tables/csv/random_set10_CIs.csv", 
          row.names = F)



