# ---
# Title: analyzing_CR_output.R
# Purpose: This script uses the nested_cv.R output from all CR genes and generates the AUROC/AUPRC curves and output supp. tables. 
#          
# ---

library(ggplot2)
library(scales)
library(cowplot)
library(ggbeeswarm)
library(openxlsx)
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
autism_CIs$group <- "Chromatin regulators"
autism_tb$group <- "Chromatin regulators"
autism_CIs$group <- factor(autism_CIs$group, levels = c("Chromatin regulators"))
autism_tb$group <- factor(autism_tb$group, levels = c("Chromatin regulators"))

#To make OR figure: go to ./making_figs.R

AUROC_AUPRC <- as.data.frame(unlist(autism_res[["aucs"]]))
AUROC_AUPRC <- cbind(AUROC_AUPRC ,unlist(autism_res[["aucprs"]]))
colnames(AUROC_AUPRC) <- c("AUROC", "AUPRC")

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

saveWorkbook(wb, "./tables/CR_output.xlsx", overwrite = TRUE)


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
    ggplot2::coord_cartesian(ylim = c(0, 1.15)) +
    ggplot2::scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
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
                   legend.position = c(0.66, 0.815),
                   legend.key.size = unit(0.6, "cm"),
                   legend.box.background = element_blank(), 
                   panel.border = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank()) +
    ggplot2::coord_cartesian(ylim = c(0, 1.15)) +
    ggplot2::scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    ggplot2::geom_hline(yintercept = ((sum(genes$`Autism susceptibility` == 1)) / (sum(genes$`Autism susceptibility` == 1 | genes$`Autism susceptibility` == 0))),
                        linetype='dashed', col = "black") + #Adding PR baseline
    ggsci::scale_color_lancet(labels = legendLabels) +
    ggplot2::labs(x = "Recall", y = "Precision")

  return(prCurve)
}

# Run nested_cv.R first to get autism_res
asd_roc <- plotROC(autism_res, 45, 46)
asd_roc
asd_pr <- plotPR(autism_res, 45, 46)
asd_pr

#Organize panels:
Figure_S5 <- ggdraw() +
  draw_plot(asd_roc, x = 0, y = 0, width = .5, height = 1) +
  draw_plot(asd_pr, x = .5, y = 0, width = .5, height = 1) +
  draw_plot_label(label = c("a", "b"), size = 14, x = c(0, 0.5), y = c(1, 1))

svg(paste0("./figures/Figure_S5.svg"), height = 3.5, width = 7.2)
print(Figure_S5)
dev.off()




