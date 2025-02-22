# ---
# Title: analyzing_PC_output.R
# Purpose: This script uses the nested_cv.R section 1 output from all PC genes and generates the AUROC/AUPRC curves and output supp. tables. 
#          
# ---

library(ggplot2)
library(scales)
library(cowplot)
library(ggbeeswarm)
library(openxlsx)
library(readxl)
library(pROC)
library(ggpubr)

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

# Run nested_cv.R to get autism_res
autism_tb <- get_beeswarm_tb(autism_res, K = 5)

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
autism_CIs$group <- "Protein-coding genes"
autism_tb$group <- "Protein-coding genes"
autism_CIs$group <- factor(autism_CIs$group, levels = c("Protein-coding genes"))
autism_tb$group <- factor(autism_tb$group, levels = c("Protein-coding genes"))

#To make OR figure: go to ./making_figs.R

AUROC_AUPRC <- as.data.frame(unlist(autism_res[["aucs"]]))
AUROC_AUPRC <- cbind(AUROC_AUPRC ,unlist(autism_res[["aucprs"]]))
colnames(AUROC_AUPRC) <- c("AUROC", "AUPRC")

#Label 303 new candidate ASD risk genes, not found by 4 other papers or in bona fide ASD risk gene list:
Krishnan <- read_excel("./data/Krishnan.xlsx", col_names = T, sheet = 1, range = "A1:O25826")
Krishnan <- Krishnan[-12506,]#remove duplicate
Krishnan <- Krishnan %>% mutate(K_quantile = ntile(desc(probability), 10)) 

Lin <- read_excel("./data/Lin.xlsx", skip = 1)
Lin <- Lin %>% mutate(L_quantile = ntile(desc(Score), 10))

Duda <- read_excel("./data/Duda.xlsx")
Duda <- Duda %>% mutate(D_quantile = ntile(desc(PredictionScore), 10))

Brueggeman <- read.csv("./data/Brueggeman.csv")
Brueggeman <- Brueggeman %>% mutate(B_quantile = ntile(desc(forecASD), 10))

genes <- genes %>% mutate(decile = ntile(desc(predictions), 10))
PC_top <- genes %>% filter(decile == 1) 
PC_top <- PC_top$`Gene name`

Krishnan_top <- Krishnan %>% filter(K_quantile == 1)
Krishnan_top <- Krishnan_top$symbol
Krishnan_top <- Krishnan_top[!is.na(Krishnan_top)]

Lin_top <- Lin %>% filter(L_quantile == 1)
Lin_top <- Lin_top$Gene

Duda_top <- Duda %>% filter(D_quantile == 1)
Duda_top <- Duda_top$Gene

Brueggeman_top <- Brueggeman %>% filter(B_quantile == 1)
Brueggeman_top <- Brueggeman_top$symbol

#ASiD top decile vs. union of all others' top deciles:
union_top <- c(Krishnan_top, Lin_top, Duda_top, Brueggeman_top)
union_top <- unique(union_top)
PC_only <- setdiff(PC_top, union_top) 

genes$newASDriskgene <- NA
r = 1
for (r in 1:nrow(genes)) {
  if (genes[r,26] == 1 & genes[r,24] == 0 & !(genes[r,4] %in% union_top)) {
    genes[r,27] <- 1
  } else if (genes[r,26] != 1) {
    genes[r,27] <- NA
  } else (genes[r,27] <- 0)
}

#Extract lambda.1se from each fold. 
lambda <- data.frame(Fold = integer(), lambda.1se = numeric())
i = 1
for (i in 1:5) {
  lambda[i,1] <- i
  lambda[i,2] <- autism_res[["models"]][[i]][["lambda.1se"]]
}


#Export data.
wb <- createWorkbook()
addWorksheet(wb = wb, sheetName = "Predictions", gridLines = T)
writeDataTable(wb = wb, sheet = 1, x = genes, withFilter = F, tableStyle = "None")
addWorksheet(wb = wb, sheetName = "Odds Ratios", gridLines = T)
writeData(wb = wb, sheet = 2, x = autism_tb)
addWorksheet(wb = wb, sheetName = "CIs", gridLines = T)
writeData(wb = wb, sheet = 3, x = autism_CIs)
addWorksheet(wb = wb, sheetName = "AUROC_AUPRC", gridLines = T)
writeData(wb = wb, sheet = 4, x = AUROC_AUPRC)
addWorksheet(wb = wb, sheetName = "lambda", gridLines = T)
writeData(wb = wb, sheet = 5, x = lambda)

saveWorkbook(wb, "./tables/PC_genes_output.xlsx", overwrite = TRUE)

#Also save files as csv's:
write.csv(genes, "./tables/csv/PC_Predictions.csv", 
          row.names = F)
write.csv(autism_tb, "./tables/csv/PC_Odds_Ratios.csv", 
          row.names = F)
write.csv(autism_CIs, "./tables/csv/PC_CIs.csv", 
          row.names = F)
write.csv(AUROC_AUPRC, "./tables/csv/PC_AUROC_AUPRC.csv", 
          row.names = F)
write.csv(genes, "./tables/csv/PC_lambda.csv", 
          row.names = F)


#Make AUROC and AUPRC figure.
# ROC
plotROC <- function(results, dependentVarIndex, predictionIndex) {
  
  # Performing checks of user input
  if (class(results) != "trainCV") {
    stop("results should be an S3 object of class trainCV.")
  }
  
  if (is.numeric(dependentVarIndex) == FALSE) {
    stop("dependentVarIndex should be a positive interger.")
  }
  
  if (is.numeric(dependentVarIndex) == TRUE && dependentVarIndex < 1) {
    stop("dependentVarIndex should be a positive interger.")
  }
  
  # Initializing variables
  K <- length(results$models)
  allpf <- data.frame(fpr = numeric(),
                      tpr = numeric(),
                      fold = numeric())
  aucs <- numeric()
  
  # Adding predictions and auc values from each model to allpf
  for (i in 1:K) {
    pred <- ROCR::prediction(results$predictions[[i]][[predictionIndex]],
                             results$predictions[[i]][[dependentVarIndex]])
    perf <- ROCR::performance(pred, "tpr", "fpr")
    
    area <- pROC::auc(results$predictions[[i]][[dependentVarIndex]],
                      results$predictions[[i]][[predictionIndex]])
    aucs <- append(aucs, area)
    
    pf <- data.frame(fpr = perf@x.values,
                     tpr = perf@y.values,
                     fold = i)
    names(pf) <- c("fpr", "tpr", "fold")
    
    allpf <- rbind(allpf, pf)
  }
  
  # Formatting the legend labels
  legendLabels <- character()
  for (i in 1:K) {
    content <- paste("Fold ", as.character(i), " (AUC = ",
                     as.character(format(round(aucs[i], 2), nsmall = 2)), ")", sep = "")
    legendLabels <- append(legendLabels, content)
  }
  
  # Generating the ROC curves
  rocCurve <-
    ggplot2::ggplot(data = allpf, aes(x = fpr, y = tpr, color = as.factor(fold))) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::theme_bw() + 
    ggplot2::theme(axis.text.x = element_text(size = 11, colour = "black"),
                  axis.text.y = element_text(size = 11, colour = "black"),
                  axis.title = element_text(size = 14, colour = "black"),
                  axis.line = element_line(size = 1),
                  legend.text = element_text(size = 12),
                  legend.title = element_blank(),
                  legend.position = c(0.65, 0.35),
                  legend.key.size = unit(0.70, "cm"),
                  legend.box.background = element_blank(), 
                  panel.border = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank()) +
    ggsci::scale_color_lancet(labels = legendLabels) +
    ggplot2::labs(x = "False positive rate", y = "True positive rate")

  return(rocCurve)
}

# PR Curves
plotPR <- function(results, dependentVarIndex, predictionIndex) {
  # Performing checks of user input
  if (class(results) != "trainCV") {
    stop("results should be an S3 object of class trainCV.")
  }
  
  if (is.numeric(dependentVarIndex) == FALSE) {
    stop("dependentVarIndex should be a positive interger.")
  }
  
  if (is.numeric(dependentVarIndex) == TRUE && dependentVarIndex < 1) {
    stop("dependentVarIndex should be a positive interger.")
  }
  
  
  # Initializing variables
  K <- length(results$models)
  allpf <- data.frame(fpr = numeric(),
                      tpr = numeric(),
                      fold = numeric())
  aucprs <- numeric()
  
  # Adding predictions and aucpr values from each model to allpf
  for (i in 1:K) {
    pred <- ROCR::prediction(results$predictions[[i]][[predictionIndex]],
                             results$predictions[[i]][[dependentVarIndex]])
    perf <- ROCR::performance(pred, "prec", "rec")
    
    aucprObj <- ROCR::performance(pred, "aucpr")
    aucpr <- aucprObj@y.values[[1]]
    aucprs <- append(aucprs, aucpr)
    
    pf <- data.frame(rec = perf@x.values,
                     prec = perf@y.values,
                     fold = i)
    data.frame(rec = perf@x.values, prec = perf@y.values)
    names(pf) <- c("rec", "prec", "fold")
    
    allpf <- rbind(allpf, pf)
  }
  
  # Formatting the legend labels
  legendLabels <- character()
  for (i in 1:K) {
    content <- paste("Fold ", as.character(i), " (AUPR = ",
                     as.character(round(aucprs[i], 2)), ")", sep = "")
    legendLabels <- append(legendLabels, content)
  }
  
  # Generating the precision-recall curves
  prCurve <-
    ggplot2::ggplot(data = allpf, aes(x = rec, y = prec, colour = as.factor(fold))) +
    ggplot2::geom_line(linewidth = 0.6, na.rm = T) +
    ggplot2::theme_bw() + 
    ggplot2::theme(axis.text.x = element_text(size = 11, colour = "black"),
                   axis.text.y = element_text(size = 11, colour = "black"),
                   axis.title = element_text(size = 14, colour = "black"),
                   axis.line = element_line(size = 1),
                   legend.text = element_text(size = 12),
                   legend.title = element_blank(),
                   legend.position = c(0.65, 0.72),
                   legend.key.size = unit(0.70, "cm"),
                   legend.box.background = element_blank(), 
                   panel.border = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank()) +
    ggplot2::geom_hline(yintercept = ((sum(genes$`Autism susceptibility` == 1)) / (sum(genes$`Autism susceptibility` == 1 | genes$`Autism susceptibility` == 0))),
                                      linetype='dashed', col = "black") + #Adding PR baseline = 0.064
    ggsci::scale_color_lancet(labels = legendLabels) +
    ggplot2::labs(x = "Recall", y = "Precision")

  return(prCurve)
}

# Run nested_cv.R first to get autism_res
asd_roc <- plotROC(autism_res, 24, 25)
asd_roc
asd_pr <- plotPR(autism_res, 24, 25)
asd_pr

#Organize panels:
Figure_3 <- ggdraw() +
  draw_plot(asd_roc, x = 0, y = 0, width = .5, height = 1) +
  draw_plot(asd_pr, x = .5, y = 0, width = .5, height = 1) +
  draw_plot_label(label = c("a", "b"), size = 14, x = c(0, 0.5), y = c(1, 1))

svg(paste0("./figures/Figure_3.svg"), height = 3.5, width = 7.2)
print(Figure_3)
dev.off()




