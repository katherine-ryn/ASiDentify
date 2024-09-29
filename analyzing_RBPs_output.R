# ---
# Title: analyzing_RBPs_output.R
# Purpose: This script uses the nested_cv.R output from all RBP genes (including all 4 lists of RBPs)
#           and generates the AUROC/AUPRC curves and output supp. tables. 
#          
# ---

library(ggplot2)
library(scales)
library(cowplot)
library(ggbeeswarm)
library(openxlsx)
library(pROC)
library(ggpubr)
library(stringr)

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

RBP_names <- c("CompRBP", "HCRBP", "CRBP", "NCRBP")
RBP_runs <- list(autism_res_CompRBP = autism_res_CompRBP,
                 autism_res_HCRBP = autism_res_HCRBP, 
                 autism_res_CRBP = autism_res_CRBP, 
                 autism_res_NCRBP = autism_res_NCRBP)

# Run nested_cv.R to get autism_res_CompRBP, autism_res_HCRBP, autism_res_CRBP, autism_res_NCRBP

n = 1
for (n in 1:length(RBP_names)) {
autism_tb <- get_beeswarm_tb(res = RBP_runs[[n]], K = 5)

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
autism_CIs$group <- RBP_names[n]
autism_tb$group <- RBP_names[n]
autism_CIs$group <- factor(autism_CIs$group, levels = RBP_names[n])
autism_tb$group <- factor(autism_tb$group, levels = RBP_names[n])

#To make OR figure: go to ./making_figs.R

AUROC_AUPRC <- as.data.frame(unlist(RBP_runs[[n]][["aucs"]]))
AUROC_AUPRC <- cbind(AUROC_AUPRC ,unlist(RBP_runs[[n]][["aucprs"]]))
colnames(AUROC_AUPRC) <- c("AUROC", "AUPRC")

assign(paste0(RBP_names[n], "_autism_tb"), autism_tb)
assign(paste0(RBP_names[n], "_autism_CIs"), autism_CIs)
assign(paste0(RBP_names[n], "_AUROC_AUPRC"), AUROC_AUPRC)

}


#Export data.
wb <- createWorkbook()
addWorksheet(wb = wb, sheetName = "Comprehensive RBP Predictions", gridLines = T)
writeDataTable(wb = wb, sheet = 1, x = CompRBP, withFilter = F, tableStyle = "None")
addWorksheet(wb = wb, sheetName = "Comprehensive RBP Odds Ratios", gridLines = T)
writeData(wb = wb, sheet = 2, x = CompRBP_autism_tb)
addWorksheet(wb = wb, sheetName = "Comprehensive RBP CIs", gridLines = T)
writeData(wb = wb, sheet = 3, x = CompRBP_autism_CIs)
addWorksheet(wb = wb, sheetName = "Comprehensive RBP AUROC_AUPRC", gridLines = T)
writeData(wb = wb, sheet = 4, x = CompRBP_AUROC_AUPRC)

addWorksheet(wb = wb, sheetName = "HC RBP Predictions", gridLines = T)
writeData(wb = wb, sheet = 5, x = HCRBP)
addWorksheet(wb = wb, sheetName = "HC Odds Ratios", gridLines = T)
writeData(wb = wb, sheet = 6, x = HCRBP_autism_tb)
addWorksheet(wb = wb, sheetName = "HC CIs", gridLines = T)
writeData(wb = wb, sheet = 7, x = HCRBP_autism_CIs)
addWorksheet(wb = wb, sheetName = "HC AUROC_AUPRC", gridLines = T)
writeData(wb = wb, sheet = 8, x = HCRBP_AUROC_AUPRC)

addWorksheet(wb = wb, sheetName = "C RBP Predictions", gridLines = T)
writeData(wb = wb, sheet = 9, x = CRBP)
addWorksheet(wb = wb, sheetName = "C RBP Odds Ratios", gridLines = T)
writeData(wb = wb, sheet = 10, x = CRBP_autism_tb)
addWorksheet(wb = wb, sheetName = "C RBP CIs", gridLines = T)
writeData(wb = wb, sheet = 11, x = CRBP_autism_CIs)
addWorksheet(wb = wb, sheetName = "C RBP AUROC_AUPRC", gridLines = T)
writeData(wb = wb, sheet = 12, x = CRBP_AUROC_AUPRC)

addWorksheet(wb = wb, sheetName = "NC RBP Predictions", gridLines = T)
writeData(wb = wb, sheet = 13, x = NCRBP)
addWorksheet(wb = wb, sheetName = "NC RBP Odds Ratios", gridLines = T)
writeData(wb = wb, sheet = 14, x = NCRBP_autism_tb)
addWorksheet(wb = wb, sheetName = "NC RBP CIs", gridLines = T)
writeData(wb = wb, sheet = 15, x = NCRBP_autism_CIs)
addWorksheet(wb = wb, sheetName = "NC RBP AUROC_AUPRC", gridLines = T)
writeData(wb = wb, sheet = 16, x = NCRBP_AUROC_AUPRC)

saveWorkbook(wb, "./tables/RBP_output.xlsx", overwrite = TRUE)


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
                  legend.position = c(0.65, 0.30),
                  legend.key.size = unit(0.6, "cm"),
                  legend.box.background = element_blank(), 
                  panel.border = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank()) +
    ggplot2::coord_cartesian(ylim = c(0, 1.1)) +
    ggplot2::scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    ggsci::scale_color_lancet(labels = legendLabels) +
    ggplot2::labs(x = "False positive rate", y = "True positive rate")

  return(rocCurve)
}

# PR Curves
plotPR <- function(results, dependentVarIndex, predictionIndex, ASDRBP, CompRBP) {
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
                     as.character(format(round(aucprs[i], 2), nsmall = 2)), ")", sep = "")
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
                   legend.position = c(0.675, 0.8),
                   legend.key.size = unit(0.6, "cm"),
                   legend.box.background = element_blank(), 
                   panel.border = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank()) +
    ggplot2::coord_cartesian(ylim = c(0, 1.1)) +
    ggplot2::scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    ggplot2::geom_hline(yintercept = (ASDRBP/CompRBP), linetype='dashed', col = "black") + #Adding PR baseline
    ggsci::scale_color_lancet(labels = legendLabels) +
    ggplot2::labs(x = "Recall", y = "Precision")

  return(prCurve)
}

# Run nested_cv.R first to get autism_res
asd_roc_CompRBP <- plotROC(autism_res_CompRBP, 33, 34)
asd_roc_CompRBP
asd_pr_CompRBP <- plotPR(autism_res_CompRBP, 33, 34, 293, 3300)
asd_pr_CompRBP

asd_roc_HCRBP <- plotROC(autism_res_HCRBP, 33, 34)
asd_roc_HCRBP
asd_pr_HCRBP <- plotPR(autism_res_HCRBP, 33, 34, 81, 1016)
asd_pr_HCRBP

asd_roc_CRBP <- plotROC(autism_res_CRBP, 33, 34)
asd_roc_CRBP
asd_pr_CRBP <- plotPR(autism_res_CRBP, 33, 34, 82, 819)
asd_pr_CRBP

asd_roc_NCRBP <- plotROC(autism_res_NCRBP, 33, 34)
asd_roc_NCRBP
asd_pr_NCRBP <- plotPR(autism_res_NCRBP, 33, 34, 211, 2481)
asd_pr_NCRBP

#Organize panels:
Figure_S3 <- ggdraw() +
  draw_plot(asd_roc_CompRBP, x = 0, y = 0.5, width = .5, height = .5) +
  draw_plot(asd_pr_CompRBP, x = .5, y = .5, width = .5, height = .5) +
  draw_plot(asd_roc_HCRBP, x = 0, y = 0, width = .5, height = .5) +
  draw_plot(asd_pr_HCRBP, x = .5, y = 0, width = .5, height = .5) +
  draw_plot_label(label = c("a", "b", "c", "d"), size = 14, x = c(0, 0.5, 0, 0.5), y = c(1, 1, 0.5, 0.5))

svg(paste0("./figures/Figure_S3.svg"), height = 7, width = 7.2)
print(Figure_S3)
dev.off()


Figure_S4 <- ggdraw() +
  draw_plot(asd_roc_CRBP, x = 0, y = 0.5, width = .5, height = .5) +
  draw_plot(asd_pr_CRBP, x = .5, y = .5, width = .5, height = .5) +
  draw_plot(asd_roc_NCRBP, x = 0, y = 0, width = .5, height = .5) +
  draw_plot(asd_pr_NCRBP, x = .5, y = 0, width = .5, height = .5) +
draw_plot_label(label = c("a", "b", "c", "d"), size = 14, x = c(0, 0.5, 0, 0.5), y = c(1, 1, 0.5, 0.5))

svg(paste0("./figures/Figure_S4.svg"), height = 7, width = 7.2)
print(Figure_S4)
dev.off()




