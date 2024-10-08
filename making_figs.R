# ---
# Title: making_figs.R
# Purpose: This script generates the figures in the paper. 
#          
# ---

library(biomaRt)
library(cowplot)
library(data.table)
library(dplyr)
library(ggbeeswarm)
library(ggplot2)
library(ggpubr)
library(VennDiagram)
library(gprofiler2)
library(readxl)
library(rstatix)
library(scales)
library(stats)
library(stringr)
library(wesanderson)


setwd("") # Change this to the path that contains the ASiD folder


palette1 <- wes_palette("Darjeeling1")
palette2 <- wes_palette("Darjeeling2")
palette <- append(palette1, palette2)
palette3 <- wes_palette("Royal2")


# FIGURE 2 ----------------------------------------------------------------

#figure generated in 'analyzing_PC_output.R'



# FIGURE 3 ----------------------------------------------------------------
###FIGURE 3A:
PC_pred <- read_excel("./tables/PC_genes_output.xlsx", sheet = 1)
PC_pred <- as.data.frame(PC_pred)
PC_pred <- PC_pred[,-c(33,34)] #Remove these columns so the code runs with the correct column indices

PC_pred <- arrange(PC_pred, desc(predictions))
PC_pred$Rank <- rank(desc(PC_pred[,'predictions']))
PC_pred <- PC_pred %>% mutate(quantile = ntile(desc(predictions), 10))

PC_pred$`Autism susceptibility` <- as.factor(PC_pred$`Autism susceptibility`)
PC_pred$`Autism susceptibility` <- factor(PC_pred$`Autism susceptibility`,levels=c("1","0"))

names(PC_pred) <- gsub("\\.| ", "_", names(PC_pred))
PC_pred$Autism_susceptibility <- as.factor(PC_pred$Autism_susceptibility)

stat.test.3a <- PC_pred %>% wilcox_test(predictions ~ Autism_susceptibility)
stat.test.3a <- stat.test.3a %>% add_xy_position(x = "Autism_susceptibility") 
stat.test.3a <- stat.test.3a %>% add_significance("p")

Fig_3a <- ggplot(PC_pred, aes(x = Autism_susceptibility, y=predictions)) +
    geom_boxplot(aes(fill = as.factor(Autism_susceptibility)), outlier.size = 1) +
    labs(x ="ASD Gene", y = "ASiD Score") +
    coord_cartesian(ylim = c(0, 1.2)) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_x_discrete(labels=c("1" = "Yes", "0" = "No")) +
    theme_bw() +
    theme(legend.position="none",
          axis.text.x.bottom = element_text(size = 11, color = "black"),
          axis.text.y.left = element_text(size = 11, color = "black"),
          axis.title = element_text(size = 11, color = "black")) +
    scale_fill_manual(values = palette[c(1,4)]) +
    stat_pvalue_manual(stat.test.3a, bracket.nudge.y = 0.06, label = "p", y.position = 1.04, label.size = 3, tip.length =  0.04)
Fig_3a


###FIGURE 3B:
SFARI <- read.csv("./data/SFARI-Gene_genes_01-16-2024release_02-25-2024export.csv", header = T)

PC_pred <- merge(PC_pred, SFARI[, c("ensembl.id", "gene.score")], 
                 by.x = "Ensembl_ID", by.y = "ensembl.id", all.x = TRUE)

PC_pred_scores <- PC_pred %>% filter(!is.na(gene.score)) #remove non-SFARI genes, or SFARI Syndromic genes
PC_pred_scores$gene.score <- as.factor(PC_pred_scores$gene.score)

stat.test.3b <- PC_pred_scores %>% wilcox_test(predictions ~ gene.score, p.adjust.method = "bonferroni")
stat.test.3b <- stat.test.3b %>% add_xy_position(x = "gene.score") 

Fig_3b <- ggplot(PC_pred_scores, aes(x = gene.score, y=predictions)) +
    geom_boxplot(aes(fill = as.factor(gene.score)), outlier.size = 1) +
    labs(x="SFARI Score", y = "ASiD Score") +
    coord_cartesian(ylim = c(0, 1.2)) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) + 
    theme_bw() +
    theme(legend.position="none",
          axis.text.x.bottom = element_text(size = 11, color = "black"),
          axis.text.y.left = element_text(size = 11, color = "black"),
          axis.title = element_text(size = 11, color = "black")) +
    scale_fill_manual(values = palette[c(7, 5, 9)]) +
    stat_pvalue_manual(stat.test.3b[c(1,2),], label = "p.adj", y.position = c(1.04, 1.15), tip.length = 0.02, label.size = 3) +
    stat_pvalue_manual(stat.test.3b[3,], label = "p.adj.signif", y.position = 0.97, tip.length = 0.02, label.size = 3)
Fig_3b


###FIGURE 3C:
SFARI_new <- read.csv("./data/SFARI-Gene_genes_08-19-2024release_08-22-2024export.csv", header = T)

new_genes <- setdiff(SFARI_new$gene.symbol, SFARI$gene.symbol)
SFARI_new_only <- filter(SFARI_new, gene.symbol  %in% new_genes)

i = 1
for (i in 1:nrow(SFARI_new_only)) {
  ENSG <- SFARI_new_only[i,4]
  row_match <- which(grepl(ENSG, PC_pred$Ensembl_ID))
  SFARI_new_only[i,c(11:13)] <- PC_pred[row_match,c(32:34)]
}


#Assume 10% chance base probability
new_SFARI_binom <- binom.test(11, 29, p = 0.1, alternative = c("greater")) #11 out of 29 in top decile
p.val3c <- new_SFARI_binom[["p.value"]]


Fig_3c <- ggplot(SFARI_new_only, aes(quantile, after_stat(density))) +
    geom_histogram(fill = palette[2], color = "black", binwidth = 1) +
    labs(y = paste0("Fraction of", "\n", "genes"), x = "Decile", title = "29 New SFARI Genes") +
    coord_cartesian(ylim = c(0, 0.72)) +
    scale_x_continuous(breaks = c(1:10)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_bw() +
    theme(legend.position="none",
          axis.text.x.bottom = element_text(size = 11, color = "black"), 
          axis.text.y.left = element_text(size = 11, color = "black"),
          axis.title = element_text(size = 11, color = "black"), 
          panel.grid.major.x = element_blank(), 
          plot.title = element_text(hjust = 0.95, vjust = -7, size = 12), 
          plot.margin =  unit(c(-13,5.5,5.5,5.5), "points")) +
    annotate('segment', x = 1, xend = 1, y = 0.4, yend = 0.5, linewidth = 0.7, colour = "black") +
    annotate('segment', x = 1, xend = 2, y = 0.5, yend = 0.5, linewidth = 0.7, colour = "black") +
    annotate("text", label = "italic(P) == 6.2e-05", parse = TRUE, x = 4, y = 0.5)
Fig_3c


#Calculate enrichment in top 20%:
new_SFARI_binom20 <- binom.test(23, 29, p = 0.2, alternative = c("greater")) #23 out of 29 in top 20%
new_SFARI_binom20[["p.value"]]


###FIGURE 3D:
AD <- read_excel("./data/Bellenguez_2022.xlsx", sheet = 21, skip = 1)
AD <- filter(AD, `Gene Prioritization Tier` %in% c("Tier 1" , "Tier 2"))
AD_genes <- AD$`ENSG ID`
AD_genes <- append(AD_genes, "ENSG00000130203") #add APOE. Total of 56 genes

AD <- rbind(AD, c("APOE", "", "ENSG00000130203", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", ""))


i = 1
for (i in 1:nrow(AD)) {
  ENSG <- AD[i, 3]
  row_match <- which(grepl(ENSG, PC_pred$Ensembl_ID))
  AD[i,c(23:25)] <- PC_pred[row_match,c(32:34)]
}

AD_binom <- binom.test(5, 56, p = 0.1, alternative = c("greater"))
AD_binom[["p.value"]]


Fig_3d <- ggplot(AD, aes(quantile, after_stat(density))) +
  geom_histogram(fill = palette[2], color = "black", binwidth = 1) +
  labs(y = paste0("Fraction of", "\n", "genes"), x = "Decile", title = "56 Alzheimer's Disease Genes") +
  coord_cartesian(ylim = c(0, 0.4)) +
  scale_x_continuous(breaks = c(1:10)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(legend.position="none", 
        axis.text.x.bottom = element_text(size = 11, color = "black"),
        axis.text.y.left = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 11, color = "black"),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(hjust = 0.95, vjust = -7, size = 12), 
        plot.margin =  unit(c(-13,5.5,5.5,5.5), "points")) +
  annotate('segment', x = 1, xend = 1, y = 0.1, yend = 0.25, linewidth = 0.7, colour = "black") +
  annotate('segment', x = 1, xend = 2, y = 0.25, yend = 0.25, linewidth = 0.7, colour = "black") +
  annotate("text", label = "italic(P) == 0.67", parse = TRUE, x = 3.5, y = 0.25)
Fig_3d


###FIGURE 3E:
PD <- read_excel("./data/Kim_2024.xlsx", sheet = 12)

#Only keep unique rows. 
PD <- PD[!duplicated(PD$`Gene Name`), ]
PD <- PD[-c(4, 13, 15, 28), -c(1:3, 6:26)]
PD[23,1] <- "ENSG00000100147"
PD[13,1] <- "ENSG00000166582"
PD[4,1] <- "ENSG00000087269"
PD <- PD[-6, ]  #remove ENSG not in PC_pred

i = 1
for (i in 1:nrow(PD)) {
  ENSG <- PD[i, 1]
  row_match <- which(grepl(ENSG, PC_pred$Ensembl_ID))
  PD[i,c(3:5)] <- PC_pred[row_match,c(32:34)]
}

PD_binom <- binom.test(3, 23, p = 0.1, alternative = c("greater"))
PD_binom[["p.value"]]


Fig_3e <- ggplot(PD, aes(quantile, after_stat(density))) +
  geom_histogram(fill = palette[2], color = "black", binwidth = 1) + 
  labs(y = paste0("Fraction of", "\n", "genes"), x = "Decile", title = "23 Parkinson's Disease Genes") +
  coord_cartesian(ylim = c(0, 0.4)) +
  scale_x_continuous(breaks = c(1:10)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(legend.position="none", 
        axis.text.x.bottom = element_text(size = 11, color = "black"),
        axis.text.y.left = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 11, color = "black"),
        panel.grid.major.x = element_blank(), 
        plot.title = element_text(hjust = 0.95, vjust = -7, size = 12), 
        plot.margin =  unit(c(-13,5.5,5.5,5.5), "points")) +
  annotate('segment', x = 1, xend = 1, y = 0.15, yend = 0.225, linewidth = 0.7, colour = "black") +
  annotate('segment', x = 1, xend = 2, y = 0.225, yend = 0.225, linewidth = 0.7, colour = "black") +
  annotate("text", label = "italic(P) == 0.41", parse = TRUE, x = 3.5, y = 0.225)
Fig_3e


#Organize panels:
Figure3 <- ggdraw() +
  draw_plot(Fig_3a, x = 0, y = .645, width = .5, height = .355) +
  draw_plot(Fig_3b, x = .5, y = .645, width = .5, height = .355) +
  draw_plot(Fig_3c, x = 0, y = 0.43, width = 1, height = 0.215) +
  draw_plot(Fig_3d, x = 0, y = .215, width = 1, height = 0.215) +
  draw_plot(Fig_3e, x = 0, y = 0, width = 1, height = 0.215) +
  draw_plot_label(label = c("a", "b", "c", "d", "e"), 
                  size = 14, 
                  x = c(0, 0.5, 0, 0, 0),
                  y = c(1, 1, 0.67, 0.45, 0.235))


svg(paste0("./figures/Figure_3.svg"), height = 7, width = 3.5)
print(Figure3)
dev.off()



# FIGURE 4 ----------------------------------------------------------------
Krishnan <- read_excel("./data/Krishnan.xlsx", col_names = T, sheet = 1, range = "A1:O25826")
Krishnan <- Krishnan[-12506,]#remove duplicate

Lin <- read_excel("./data/Lin.xlsx", skip = 1)

Duda <- read_excel("./data/Duda.xlsx")

Brueggeman <- read.csv("./data/Brueggeman.csv")

PC_pred <- merge(PC_pred, Krishnan[, c("symbol", "probability")],
                 by.x = "Gene_name", by.y = "symbol", all.x = TRUE)
PC_pred <- merge(PC_pred, Lin[, c("Gene", "Score")], 
                 by.x = "Gene_name", by.y = "Gene", all.x = TRUE)
PC_pred <- merge(PC_pred, Duda[, c("Gene", "PredictionScore")],
                 by.x = "Gene_name", by.y = "Gene", all.x = TRUE)
PC_pred <- merge(PC_pred, Brueggeman[, c("symbol", "forecASD")],
                 by.x = "Gene_name", by.y = "symbol", all.x = TRUE)

names(PC_pred)[names(PC_pred) == 'probability'] <- 'Krishnan'
names(PC_pred)[names(PC_pred) == 'Score'] <- 'Lin'
names(PC_pred)[names(PC_pred) == 'PredictionScore'] <- 'Duda'
names(PC_pred)[names(PC_pred) == 'forecASD'] <- 'Brueggeman'


#Make 4 new dataframes, with genes present in both analyses:
PC_pred_Krishnan <- PC_pred[!is.na(PC_pred$Krishnan), ] 
PC_pred_Lin <- PC_pred[!is.na(PC_pred$Lin), ]
PC_pred_Duda <- PC_pred[!is.na(PC_pred$Duda), ]
PC_pred_Brueggeman <- PC_pred[!is.na(PC_pred$Brueggeman), ]

#Re-calculate the top decile according to these genes:
PC_pred_Krishnan <- PC_pred_Krishnan %>%
  mutate(quantile_K = ntile(desc(predictions), 10)) 
PC_pred_Krishnan <- PC_pred_Krishnan %>%
  mutate(quantile_Krishnan = ntile(desc(Krishnan), 10)) 

PC_pred_Lin <- PC_pred_Lin %>%
  mutate(quantile_K = ntile(desc(predictions), 10)) 
PC_pred_Lin <- PC_pred_Lin %>%
  mutate(quantile_Lin = ntile(desc(Lin), 10)) 

PC_pred_Duda <- PC_pred_Duda %>%
  mutate(quantile_K = ntile(desc(predictions), 10)) 
PC_pred_Duda <- PC_pred_Duda %>%
  mutate(quantile_Duda = ntile(desc(Duda), 10)) 

PC_pred_Brueggeman <- PC_pred_Brueggeman %>% 
  mutate(quantile_K = ntile(desc(predictions), 10)) 
PC_pred_Brueggeman <- PC_pred_Brueggeman %>%
  mutate(quantile_Brueggeman = ntile(desc(Brueggeman), 10)) 

#Krishnan
#matrix(c(intersection, (in B not in A), (in A not in B), (total - union)) 
matrix_Krishnan <- matrix(c(length(intersect(PC_pred_Krishnan[PC_pred_Krishnan$quantile_K == 1, 1], PC_pred_Krishnan[PC_pred_Krishnan$quantile_Krishnan == 1, 1])),
                            length(setdiff(PC_pred_Krishnan[PC_pred_Krishnan$quantile_Krishnan == 1, 1], PC_pred_Krishnan[PC_pred_Krishnan$quantile_K == 1, 1])), 
                            length(setdiff(PC_pred_Krishnan[PC_pred_Krishnan$quantile_K == 1, 1], PC_pred_Krishnan[PC_pred_Krishnan$quantile_Krishnan == 1, 1])),
                            nrow(PC_pred_Krishnan) - length(union(PC_pred_Krishnan[PC_pred_Krishnan$quantile_K == 1, 1], PC_pred_Krishnan[PC_pred_Krishnan$quantile_Krishnan == 1, 1]))), nrow = 2)
Fisher_Krishnan <- fisher.test(matrix_Krishnan, alternative = "greater") #7.964607e-298

Krishnan_venn_list <- list(`Top 10%` = PC_pred_Krishnan[PC_pred_Krishnan$quantile_K == 1, 1],
                           `Top 10% Krishnan` = PC_pred_Krishnan[PC_pred_Krishnan$quantile_Krishnan == 1, 1])

Kvenn <- venn.diagram(Krishnan_venn_list, category.names = c("ASiD", paste0("Krishnan", "\n", "et al. 2016")), filename = NULL,
                           scaled = T,
                           sub = expression(paste(italic("P"), " = 7.96e-298")),
                           sub.pos = c(0.5,0.11),
                           fontfamily  = "sans",
                           cat.fontfamily = "sans", 
                           sub.fontfamily = "sans",
                           fill = c("#00A08A", "#9986A5"), 
                           alpha = 0.2, 
                           col = c("#00A08A", "#9986A5"), 
                           cex = 1.3, 
                           cat.cex = 1.1, 
                           cat.pos = c(-20, 20), 
                           cat.dist = c(0.03, 0.07)) 



#Lin
matrix_Lin <- matrix(c(length(intersect(PC_pred_Lin[PC_pred_Lin$quantile_K == 1, 1], PC_pred_Lin[PC_pred_Lin$quantile_Lin == 1, 1])),
                       length(setdiff(PC_pred_Lin[PC_pred_Lin$quantile_Lin == 1, 1], PC_pred_Lin[PC_pred_Lin$quantile_K == 1, 1])), 
                       length(setdiff(PC_pred_Lin[PC_pred_Lin$quantile_K == 1, 1], PC_pred_Lin[PC_pred_Lin$quantile_Lin == 1, 1])),
                       nrow(PC_pred_Lin) - length(union(PC_pred_Lin[PC_pred_Lin$quantile_K == 1, 1], PC_pred_Lin[PC_pred_Lin$quantile_Lin == 1, 1]))), nrow = 2)
Fisher_Lin <- fisher.test(matrix_Lin, alternative = "greater") #0

Lin_venn_list <- list(`Top 10%` = PC_pred_Lin[PC_pred_Lin$quantile_K == 1, 1],
                      `Top 10% Lin` = PC_pred_Lin[PC_pred_Lin$quantile_Lin == 1, 1])

Lvenn <- venn.diagram(Lin_venn_list, category.names = c("ASiD", paste0("Lin", "\n", "et al. 2020")), filename = NULL,
                      scaled = T,
                      sub = expression(paste(italic("P"), " = 0.00")),
                      sub.pos = c(0.5,0.10),
                      fontfamily  = "sans",
                      cat.fontfamily = "sans", 
                      sub.fontfamily = "sans",
                      fill = c("#00A08A", "#F98400"), 
                      alpha = 0.2, 
                      col = c("#00A08A", "#F98400"), 
                      cex = 1.3, 
                      cat.cex = 1.1, 
                      cat.pos = c(-20, 20), 
                      cat.dist = c(0.04, 0.07)) 



#Duda
matrix_Duda <- matrix(c(length(intersect(PC_pred_Duda[PC_pred_Duda$quantile_K == 1, 1], PC_pred_Duda[PC_pred_Duda$quantile_Duda == 1, 1])),
                        length(setdiff(PC_pred_Duda[PC_pred_Duda$quantile_Duda == 1, 1], PC_pred_Duda[PC_pred_Duda$quantile_K == 1, 1])), 
                        length(setdiff(PC_pred_Duda[PC_pred_Duda$quantile_K == 1, 1], PC_pred_Duda[PC_pred_Duda$quantile_Duda == 1, 1])),
                        nrow(PC_pred_Duda) - length(union(PC_pred_Duda[PC_pred_Duda$quantile_K == 1, 1], PC_pred_Duda[PC_pred_Duda$quantile_Duda == 1, 1]))), nrow = 2)
Fisher_Duda <- fisher.test(matrix_Duda, alternative = "greater") #0


Duda_venn_list <- list(`Top 10%` = PC_pred_Duda[PC_pred_Duda$quantile_K == 1, 1],
                       `Top 10% Duda` = PC_pred_Duda[PC_pred_Duda$quantile_Duda == 1, 1])

Dvenn <- venn.diagram(Duda_venn_list, category.names = c("ASiD", paste0("Duda", "\n", "et al. 2018")), filename = NULL,
                      scaled = T,
                      sub = expression(paste(italic("P"), " = 0.00")),
                      sub.pos = c(0.5,0.10),
                      fontfamily  = "sans",
                      cat.fontfamily = "sans", 
                      sub.fontfamily = "sans",
                      fill = c("#00A08A", "#F8AFA8"), 
                      alpha = 0.2, 
                      col = c("#00A08A", "#F8AFA8"), 
                      cex = 1.3, 
                      cat.cex = 1.1, 
                      cat.pos = c(-20, 20), 
                      cat.dist = c(0.04, 0.07)) 



#Brueggeman
matrix_Brueggeman <- matrix(c(length(intersect(PC_pred_Brueggeman[PC_pred_Brueggeman$quantile_K == 1, 1], PC_pred_Brueggeman[PC_pred_Brueggeman$quantile_Brueggeman == 1, 1])),
                              length(setdiff(PC_pred_Brueggeman[PC_pred_Brueggeman$quantile_Brueggeman == 1, 1], PC_pred_Brueggeman[PC_pred_Brueggeman$quantile_K == 1, 1])), 
                              length(setdiff(PC_pred_Brueggeman[PC_pred_Brueggeman$quantile_K == 1, 1], PC_pred_Brueggeman[PC_pred_Brueggeman$quantile_Brueggeman == 1, 1])),
                              nrow(PC_pred_Brueggeman) - length(union(PC_pred_Brueggeman[PC_pred_Brueggeman$quantile_K == 1, 1], PC_pred_Brueggeman[PC_pred_Brueggeman$quantile_Brueggeman == 1, 1]))), nrow = 2)
Fisher_Brueggeman <- fisher.test(matrix_Brueggeman, alternative = "greater") #4.144624e-305


Brueggeman_venn_list <- list(`Top 10%` = PC_pred_Brueggeman[PC_pred_Brueggeman$quantile_K == 1, 1],
                             `Top 10% Brueggeman` = PC_pred_Brueggeman[PC_pred_Brueggeman$quantile_Brueggeman == 1, 1])


Bvenn <- venn.diagram(Brueggeman_venn_list, category.names = c("ASiD", paste0("Brueggeman", "\n", "et al. 2020")), filename = NULL,
                      scaled = T,
                      sub = expression(paste(italic("P"), " = 4.14e-305")),
                      sub.pos = c(0.5,0.11),
                      fontfamily  = "sans",
                      cat.fontfamily = "sans", 
                      sub.fontfamily = "sans",
                      fill = c("#00A08A", "#FF0000"), 
                      alpha = 0.2, 
                      col = c("#00A08A", "#FF0000"), 
                      cex = 1.3, 
                      cat.cex = 1.1, 
                      cat.pos = c(-20, 20), 
                      cat.dist = c(0.04, 0.07)) 


#Organize panels:
Figure_4 <- ggdraw() +
  draw_plot(Bvenn, x = 0, y = 0.5, width = .5, height = .5) +
  draw_plot(Dvenn, x = .5, y = .5, width = .5, height = .5) +
  draw_plot(Kvenn, x = 0, y = 0, width = 0.5, height = 0.5) +
  draw_plot(Lvenn, x = 0.5, y = 0, width = 0.5, height = 0.5) +
  draw_plot_label(label = c("a", "b", "c", "d"), 
                  size = 14, 
                  x = c(0, 0.5, 0, 0.5),
                  y = c(1, 1, 0.5, 0.5))

svg(paste0("./figures/Figure_4.svg"), height = 6, width = 6.5)
print(Figure_4)
dev.off()





# FIGURE 5 ----------------------------------------------------------------
PC_ORs <- read_excel("./tables/PC_genes_output.xlsx", sheet = 2)
PC_CIs <- read_excel("./tables/PC_genes_output.xlsx", sheet = 3)


PC_CIs <- PC_CIs %>% filter(!(CI_low < 1 & CI_high > 1))
PC_feat <- PC_CIs$variable

PC_ORs <- PC_ORs %>% filter(variable %in% PC_feat)

PC_CIs$group <- "Protein-coding genes"
PC_ORs$group <- "Protein-coding genes"
PC_CIs$group <- factor(PC_CIs$group, levels = c("Protein-coding genes"))
PC_ORs$group <- factor(PC_ORs$group, levels = c("Protein-coding genes"))


Figure_5 <- ggplot(PC_CIs, aes(x = variable, color = as.factor(group))) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), position = position_dodge(0.6), width = 0, linewidth = 0.75) +
  geom_point(aes(y = Mean), position = position_dodge(0.6), shape = 3, size = 6, stroke = 0.75) +
  geom_beeswarm(data = PC_ORs, aes(x = variable, y = OR, fill = as.factor(group)), size = 2.5, dodge.width = 0.6, shape = 21, color = "black", alpha = 0.85) +
  ylab("Odds Ratio") +
  scale_color_manual(values = "darkgrey") + 
  scale_fill_manual(values = palette[7]) + 
  scale_x_discrete(labels = label_wrap(10)) +
  scale_y_continuous(breaks = seq(0.0, 2.0, by = 0.2),
                     labels =  label_number(accuracy = 0.1)) +
  geom_hline(yintercept=1, linetype='dotted', col = "black", linewidth = 1) +
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14, colour = "black"),
        axis.text.x.bottom = element_text(size = 11, colour = "black"),
        axis.text.y.left = element_text(size = 11, colour = "black"),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        panel.grid.major.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "none")
Figure_5


svg(paste0("./figures/Figure_5.svg"), height = 3, width = 5)
print(Figure_5)
dev.off()


# FIGURE 6 ----------------------------------------------------------------
#STEP 1: Get all paralogs for all ENSG IDs. 
mart <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", version = 111, dataset = "hsapiens_gene_ensembl")

ensembl_id <- PC_pred$Ensembl_ID

paralogs <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","hsapiens_paralog_ensembl_gene","hsapiens_paralog_associated_gene_name"),
                 filters = "ensembl_gene_id",
                 values = ensembl_id,
                 mart = mart)


#STEP 2: Aggregate all paralogs together 
merged_copy <- paralogs
dat.merged <- merged_copy %>%
  dplyr::group_by(ensembl_gene_id) %>%
  dplyr::summarise(hsapiens_paralog_associated_gene_name = paste(hsapiens_paralog_associated_gene_name, collapse = ";"))

dat.merged$hsapiens_paralog_associated_gene_name = sub("^$", "N/A", dat.merged$hsapiens_paralog_associated_gene_name)

PC_pred <- merge(PC_pred, dat.merged[, c("ensembl_gene_id", "hsapiens_paralog_associated_gene_name")],
                 by.x = "Ensembl_ID", by.y = "ensembl_gene_id", all.x = TRUE)
names(PC_pred)[names(PC_pred) == "hsapiens_paralog_associated_gene_name"] <- "All_Paralogs"

#STEP 3: Search ASD genes to check if present or not
merged_copy <- cbind(merged_copy, ASD=NA)
merged_copy$hsapiens_paralog_ensembl_gene = sub("^$", "N/A", merged_copy$hsapiens_paralog_ensembl_gene)
merged_copy$hsapiens_paralog_associated_gene_name = sub("^$", "N/A", merged_copy$hsapiens_paralog_associated_gene_name)

PC_pred_ASD <- PC_pred %>% filter(Autism_susceptibility == 1)

i = 1
for (i in 1:nrow(merged_copy)) {
  if (any(grepl(merged_copy[i,3], PC_pred_ASD$Ensembl_ID))) {
    merged_copy[i,5] = "Yes"
  } else {
    merged_copy[i,5] = "No"
  }
}

#STEP 4: Keep only genes that are present, then aggregate them 
merged_copy_only_present <- merged_copy
merged_copy_only_present <- merged_copy_only_present[!grepl("No", merged_copy_only_present$ASD),]

ASD.merged <- merged_copy_only_present %>%
  dplyr::group_by(ensembl_gene_id) %>%
  dplyr::summarise(hsapiens_paralog_associated_gene_name = paste(hsapiens_paralog_associated_gene_name, collapse = ";"))


PC_pred <- merge(PC_pred, ASD.merged[, c("ensembl_gene_id", "hsapiens_paralog_associated_gene_name")],
                 by.x = "Ensembl_ID", by.y = "ensembl_gene_id", all.x = TRUE)
names(PC_pred)[names(PC_pred) == "hsapiens_paralog_associated_gene_name"] <- "ASD_Paralogs"

###STEP 5: PARALOGY vs PREDICTION SCORE 
#New col w/ each of the 4 classifications
paralogy_classes <- PC_pred
paralogy_classes <- cbind(paralogy_classes, Paralogy_Class = NA)

i = 1
for (i in 1:nrow(paralogy_classes)) {
  if (paralogy_classes[i,31] == 1 & (!is.na(paralogy_classes[i,"ASD_Paralogs"]))) {
    paralogy_classes[i,"Paralogy_Class"] = "ASD Gene with ASD Paralog"
  } else if (paralogy_classes[i,31] == 1 & (is.na(paralogy_classes[i,"ASD_Paralogs"]))){
    paralogy_classes[i,"Paralogy_Class"] = "ASD Gene without ASD Paralog"
  } else if (paralogy_classes[i,31] == 0 & (!is.na(paralogy_classes[i,"ASD_Paralogs"]))){
    paralogy_classes[i,"Paralogy_Class"] = "Non-ASD Gene with ASD Paralog"
  } else if (paralogy_classes[i,31] == 0 & (is.na(paralogy_classes[i,"ASD_Paralogs"]))) {
    paralogy_classes[i,"Paralogy_Class"] = "Non-ASD Gene without ASD Paralog"
  }
}

#Remove genes w/out any paralogs:
paralogy_classes <- paralogy_classes[!grepl("N/A", paralogy_classes$All_Paralogs),]

#Non-ASD gene w/ ASD Paralog vs Non-ASD gene w/o ASD paralog:
match <- c("Non-ASD Gene with ASD Paralog", "Non-ASD Gene without ASD Paralog")
match_df <- filter(paralogy_classes, Paralogy_Class %in% match)
match_df$Paralogy_Class <- as.factor(match_df$Paralogy_Class) 

stat.test.6 <- match_df %>% wilcox_test(predictions ~ Paralogy_Class) %>% add_significance("p")
stat.test.6 <- stat.test.6 %>% add_xy_position(x = "Paralogy_Class")
stat.test.6
Fig_6 <- ggplot(match_df, aes(x=Paralogy_Class, y=predictions)) +
    geom_boxplot(aes(fill = as.factor(Paralogy_Class)), outlier.size = 1) +  
    labs(x = element_blank(), y = "ASiD Score") + 
    coord_cartesian(ylim = c(0, 1.1)) + 
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) + 
    theme_bw() +
    theme(legend.position="none", 
          axis.text.x.bottom = element_text(size = 11, color = "black"), 
          axis.text.y.left = element_text(size = 11, color = "black"),
          axis.title.x = element_text(size = 11, color = "black"), 
          axis.title.y = element_text(size = 14, color = "black")) +
    scale_fill_manual(values = palette[c(2,3)]) +
    stat_pvalue_manual(stat.test.6, label = "p", y.position = 1.04, label.size = 4) + 
    scale_x_discrete(labels = label_wrap(13)) 
Fig_6


svg(paste0("./figures/Figure_6.svg"), height = 3.5, width = 3.5)
print(Fig_6)
dev.off()






# FIGURE 7 ----------------------------------------------------------------
#All RBPs
AllRBP_ORs <- read_excel("./tables/RBP_output.xlsx", sheet = 2)
AllRBP_CIs <- read_excel("./tables/RBP_output.xlsx", sheet = 3)

AllRBP_CIs <- AllRBP_CIs %>% filter(!(CI_low < 1 & CI_high > 1))
AllRBP_feat <- AllRBP_CIs$variable

AllRBP_ORs <- AllRBP_ORs %>% filter(variable %in% AllRBP_feat)

AllRBP_CIs$group <- "RBP"
AllRBP_ORs$group <- "RBP"
AllRBP_CIs$group <- factor(AllRBP_CIs$group, levels = c("RBP"))
AllRBP_ORs$group <- factor(AllRBP_ORs$group, levels = c("RBP"))


Figure_7a <- ggplot(AllRBP_CIs, aes(x = variable, color = as.factor(group))) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), position = position_dodge(0.6), width = 0, linewidth = 0.75) +
  geom_point(aes(y = Mean), position = position_dodge(0.6), shape = 3, size = 6, stroke = 0.75) +
  geom_beeswarm(data = AllRBP_ORs, aes(x = variable, y = OR, fill = as.factor(group)), size = 2.5, dodge.width = 0.6, shape = 21, color = "black", alpha = 0.85) +
  ylab("Odds Ratio") +
  scale_color_manual(values = "darkgrey") + 
  scale_fill_manual(values = palette[4]) + 
  scale_x_discrete(labels = label_wrap(10)) +
  scale_y_continuous(breaks = seq(0.0, 2.0, by = 0.2),
                     labels =  label_number(accuracy = 0.1)) +
  coord_cartesian(ylim = c(0.2, 1.5)) +
  geom_hline(yintercept=1, linetype='dotted', col = "black", linewidth = 1) +
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14, colour = "black"),
        axis.text.x.bottom = element_text(size = 11, colour = "black"),
        axis.text.y.left = element_text(size = 11, colour = "black"),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        panel.grid.major.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "none")
Figure_7a

#High-confidence RBPs
HCRBP_ORs <- read_excel("./tables/RBP_output.xlsx", sheet = 6)
HCRBP_CIs <- read_excel("./tables/RBP_output.xlsx", sheet = 7)

HCRBP_CIs <- HCRBP_CIs %>% filter(!(CI_low < 1 & CI_high > 1))
HCRBP_feat <- HCRBP_CIs$variable

HCRBP_ORs <- HCRBP_ORs %>% filter(variable %in% HCRBP_feat)

HCRBP_CIs$group <- "RBP"
HCRBP_ORs$group <- "RBP"
HCRBP_CIs$group <- factor(HCRBP_CIs$group, levels = c("RBP"))
HCRBP_ORs$group <- factor(HCRBP_ORs$group, levels = c("RBP"))


Figure_7b <- ggplot(HCRBP_CIs, aes(x = variable, color = as.factor(group))) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), position = position_dodge(0.6), width = 0, linewidth = 0.75) +
  geom_point(aes(y = Mean), position = position_dodge(0.6), shape = 3, size = 6, stroke = 0.75) +
  geom_beeswarm(data = HCRBP_ORs, aes(x = variable, y = OR, fill = as.factor(group)), size = 2.5, dodge.width = 0.6, shape = 21, color = "black", alpha = 0.85) +
  ylab("Odds Ratio") +
  scale_color_manual(values = "darkgrey") + 
  scale_fill_manual(values = palette[4]) + 
  scale_x_discrete(labels = label_wrap(10)) +
  scale_y_continuous(breaks = seq(0.0, 2.0, by = 0.2),
                     labels =  label_number(accuracy = 0.1)) +
  coord_cartesian(ylim = c(0.2, 1.5)) +
  geom_hline(yintercept=1, linetype='dotted', col = "black", linewidth = 1) +
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14, colour = "black"),
        axis.text.x.bottom = element_text(size = 11, colour = "black"),
        axis.text.y.left = element_text(size = 11, colour = "black"),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        panel.grid.major.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "none")
Figure_7b

#Noncanonical RBPs
NCRBP_ORs <- read_excel("./tables/RBP_output.xlsx", sheet = 14)
NCRBP_CIs <- read_excel("./tables/RBP_output.xlsx", sheet = 15)

NCRBP_CIs <- NCRBP_CIs %>% filter(!(CI_low < 1 & CI_high > 1))
NCRBP_feat <- NCRBP_CIs$variable

NCRBP_ORs <- NCRBP_ORs %>% filter(variable %in% NCRBP_feat)

NCRBP_CIs$group <- "RBP"
NCRBP_ORs$group <- "RBP"
NCRBP_CIs$group <- factor(NCRBP_CIs$group, levels = c("RBP"))
NCRBP_ORs$group <- factor(NCRBP_ORs$group, levels = c("RBP"))


Figure_7c <- ggplot(NCRBP_CIs, aes(x = variable, color = as.factor(group))) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), position = position_dodge(0.6), width = 0, linewidth = 0.75) +
  geom_point(aes(y = Mean), position = position_dodge(0.6), shape = 3, size = 6, stroke = 0.75) +
  geom_beeswarm(data = NCRBP_ORs, aes(x = variable, y = OR, fill = as.factor(group)), size = 2.5, dodge.width = 0.6, shape = 21, color = "black", alpha = 0.85) +
  ylab("Odds Ratio") +
  scale_color_manual(values = "darkgrey") + 
  scale_fill_manual(values = palette[4]) + 
  scale_x_discrete(labels = label_wrap(10)) +
  scale_y_continuous(breaks = seq(0.0, 2.0, by = 0.2),
                     labels =  label_number(accuracy = 0.1)) +
  coord_cartesian(ylim = c(0.2, 1.5)) +
  geom_hline(yintercept=1, linetype='dotted', col = "black", linewidth = 1) +
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14, colour = "black"),
        axis.text.x.bottom = element_text(size = 11, colour = "black"),
        axis.text.y.left = element_text(size = 11, colour = "black"),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        panel.grid.major.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "none")
Figure_7c

#Canonical RBPs
CRBP_ORs <- read_excel("./tables/RBP_output.xlsx", sheet = 10)
CRBP_CIs <- read_excel("./tables/RBP_output.xlsx", sheet = 11)

CRBP_CIs <- CRBP_CIs %>% filter(!(CI_low < 1 & CI_high > 1))
CRBP_feat <- CRBP_CIs$variable

CRBP_ORs <- CRBP_ORs %>% filter(variable %in% CRBP_feat)

CRBP_CIs$group <- "RBP"
CRBP_ORs$group <- "RBP"
CRBP_CIs$group <- factor(CRBP_CIs$group, levels = c("RBP"))
CRBP_ORs$group <- factor(CRBP_ORs$group, levels = c("RBP"))


Figure_7d <- ggplot(CRBP_CIs, aes(x = variable, color = as.factor(group))) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), position = position_dodge(0.6), width = 0, linewidth = 0.75) +
  geom_point(aes(y = Mean), position = position_dodge(0.6), shape = 3, size = 6, stroke = 0.75) +
  geom_beeswarm(data = CRBP_ORs, aes(x = variable, y = OR, fill = as.factor(group)), size = 2.5, dodge.width = 0.6, shape = 21, color = "black", alpha = 0.85) +
  ylab("Odds Ratio") +
  scale_color_manual(values = "darkgrey") + 
  scale_fill_manual(values = palette[4]) + 
  scale_x_discrete(labels = label_wrap(10)) +
  scale_y_continuous(breaks = seq(0.0, 2.0, by = 0.2),
                     labels =  label_number(accuracy = 0.1)) +
  coord_cartesian(ylim = c(0.2, 1.5)) +
  geom_hline(yintercept=1, linetype='dotted', col = "black", linewidth = 1) +
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14, colour = "black"),
        axis.text.x.bottom = element_text(size = 11, colour = "black"),
        axis.text.y.left = element_text(size = 11, colour = "black"),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        panel.grid.major.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "none")
Figure_7d




#Organize panels:
Figure_7 <- ggdraw() +
  draw_plot(Figure_7a, x = 0, y = 0.5, width = .55, height = .5) +
  draw_plot(Figure_7b, x = .55, y = .527, width = .45, height = .473) +
  draw_plot(Figure_7c, x = 0, y = 0, width = 0.55, height = 0.5) +
  draw_plot(Figure_7d, x = 0.55, y = 0.027, width = 0.45, height = 0.473) +
  draw_plot_label(label = c("a", "b", "c", "d"), 
                  size = 14, 
                  x = c(0, 0.55, 0, 0.55),
                  y = c(1, 1, 0.5, 0.5))

svg(paste0("./figures/Figure_7.svg"), height = 6, width = 7.2)
print(Figure_7)
dev.off()



# FIGURE 8 ----------------------------------------------------------------
CR_ORs <- read_excel("./tables/CR_output.xlsx", sheet = 2)
CR_CIs <- read_excel("./tables/CR_output.xlsx", sheet = 3)

CR_CIs <- CR_CIs %>% filter(!(CI_low < 1 & CI_high > 1))
CR_feat <- CR_CIs$variable

CR_ORs <- CR_ORs %>% filter(variable %in% CR_feat)

CR_CIs$group <- "CR"
CR_ORs$group <- "CR"
CR_CIs$group <- factor(CR_CIs$group, levels = c("CR"))
CR_ORs$group <- factor(CR_ORs$group, levels = c("CR"))


Figure_8 <- ggplot(CR_CIs, aes(x = variable, color = as.factor(group))) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), position = position_dodge(0.6), width = 0, linewidth = 0.75) +
  geom_point(aes(y = Mean), position = position_dodge(0.6), shape = 3, size = 6, stroke = 0.75) +
  geom_beeswarm(data = CR_ORs, aes(x = variable, y = OR, fill = as.factor(group)), size = 2.5, dodge.width = 0.6, shape = 21, color = "black", alpha = 0.85) +
  ylab("Odds Ratio") +
  scale_color_manual(values = "darkgrey") + 
  scale_fill_manual(values = palette[2]) + 
  scale_x_discrete(labels = label_wrap(10)) +
  scale_y_continuous(breaks = seq(0.0, 2.0, by = 0.2),
                     labels =  label_number(accuracy = 0.1)) +
  geom_hline(yintercept=1, linetype='dotted', col = "black", linewidth = 1) +
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14, colour = "black"),
        axis.text.x.bottom = element_text(size = 11, colour = "black"),
        axis.text.y.left = element_text(size = 11, colour = "black"),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        panel.grid.major.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "none")
Figure_8


svg(paste0("./figures/Figure_8.svg"), height = 3, width = 4)
print(Figure_8)
dev.off()


### SUPPLEMENTARY ###
# FIGURE S1 -----------------------------------------------------------------
#PC_pred_ASD <- PC_pred %>% filter(Autism_susceptibility == 1)

table(PC_pred_ASD$quantile)
#Assume 10% chance base probability
PC_pred_ASD_binom <- binom.test(491, 1163, p = 0.1, alternative = c("greater")) #491 out of 1163 bona fide ASD risk genes
p.val1s <- PC_pred_ASD_binom[["p.value"]]


S1_title <- expression(paste("1,163 ", italic("bona fide"), " ASD risk genes"))

Fig_S1 <- ggplot(PC_pred_ASD, aes(quantile, after_stat(density))) +
  geom_histogram(fill = palette[2], color = "black", binwidth = 1) +
  labs(y = paste0("Fraction of", "\n", "genes"), x = "Decile", title = S1_title) +
  coord_cartesian(ylim = c(0, 0.7)) +
  scale_x_continuous(breaks = c(1:10)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(legend.position="none",
        axis.text.x.bottom = element_text(size = 11, color = "black"), 
        axis.text.y.left = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 11, color = "black"), 
        panel.grid.major.x = element_blank(), 
        plot.title = element_text(hjust = 0.95, vjust = -8, size = 12), 
        plot.margin =  unit(c(-13,5.5,5.5,5.5), "points")) +
  annotate('segment', x = 1, xend = 1, y = 0.45, yend = 0.5, linewidth = 0.7, colour = "black") +
  annotate('segment', x = 1, xend = 2, y = 0.5, yend = 0.5, linewidth = 0.7, colour = "black") +
  annotate("text", label = "italic(P) == 4.5e-180", parse = TRUE, x = 3.75, y = 0.5)
Fig_S1


svg(paste0("./figures/Figure_S1.svg"),  height = 2, width = 4)
print(Fig_S1)
dev.off()


# FIGURE S2 & TABLE S4 ---------------------------------------------------------------
GO_decile <- gost((filter(PC_pred, quantile == 1) %>% pull(Gene_name)), 
                  organism = "hsapiens", 
                  significant = T,
                  exclude_iea = F, 
                  user_threshold = 0.05,
                  correction_method = "fdr", 
                  custom_bg = PC_pred$Gene_name, 
                  sources = c("GO:BP", "GO:MF", "GO:CC"))

#Filter for term_size <2000:
GO_decile_2000 <- GO_decile[["result"]] %>% filter(term_size <2000)

#Keep top 5 of each GO term category:
GO_decile_2000 <- GO_decile_2000[with(GO_decile_2000, order(GO_decile_2000$source, GO_decile_2000$p_value)),]
topBP <- GO_decile_2000 %>% filter(source == "GO:BP") %>% arrange(p_value) %>% dplyr::slice(1:5)
topCC <- GO_decile_2000 %>% filter(source == "GO:CC") %>% arrange(p_value) %>% dplyr::slice(1:5)
topMF <- GO_decile_2000 %>% filter(source == "GO:MF") %>% arrange(p_value) %>% dplyr::slice(1:5)

top5 <- rbind(topBP, topCC, topMF)
top5$term_name <- factor(top5$term_name, levels = rev(top5$term_name))

Figure_S2 <- ggplot(top5, aes(x = term_name, y = -log10(p_value), fill = source)) +
  geom_col(position = "dodge") +
  labs(x = element_blank(), y = expression(paste(-log[10], "(adjusted ", italic("P"), ")"))) + 
  coord_flip() +
  theme_bw() +
  theme(legend.position="none",
        axis.text.x.bottom = element_text(size = 11, color = "black"),
        axis.text.y.left = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        strip.text = element_text(size = 11)) +
  scale_fill_manual(values = palette[c(1,2,3)])
Figure_S2

Figure_S2 <- Figure_S2 +
  theme(strip.placement = "outside",
        panel.spacing = unit(0, "lines")) +
  facet_wrap(~factor(source, c("GO:BP", "GO:CC", "GO:MF")), nrow=3, strip.position = "right", scales = "free_y")
Figure_S2


svg(paste0("./figures/Figure_S2.svg"),  height = 4, width = 7.2)
print(Figure_S2)
dev.off()


#Export Table S4:
GO_decile_2000$parents <- sapply(GO_decile_2000$parents, paste, collapse = ", ")

write.csv(GO_decile_2000, file = "./tables/Table_S4.csv", row.names = F)


# FIGURE S3 ---------------------------------------------------------------
#figure generated in 'analyzing_RBPs_output.R'

# FIGURE S4 ---------------------------------------------------------------
#figure generated in 'analyzing_RBPs_output.R'


# FIGURE S5 ---------------------------------------------------------------
#figure generated in 'analyzing_CR_output.R'



# TABLE S3 ----------------------------------------------------------------
fisher_df <- data.frame(Paper = NA,
                        Decile = NA,
                        Genes_in_ASiD_decile = NA,
                        Genes_in_Paper_decile = NA,
                        Intersection = NA,
                        pvalue = NA)

#Krishnan
i = 1
for (i in 1:10) {
  matrix_Krishnan <- matrix(c(length(intersect(PC_pred_Krishnan[PC_pred_Krishnan$quantile_K == i, 1], PC_pred_Krishnan[PC_pred_Krishnan$quantile_Krishnan == i, 1])),
                              length(setdiff(PC_pred_Krishnan[PC_pred_Krishnan$quantile_Krishnan == i, 1], PC_pred_Krishnan[PC_pred_Krishnan$quantile_K == i, 1])), 
                              length(setdiff(PC_pred_Krishnan[PC_pred_Krishnan$quantile_K == i, 1], PC_pred_Krishnan[PC_pred_Krishnan$quantile_Krishnan == i, 1])),
                              nrow(PC_pred_Krishnan) - length(union(PC_pred_Krishnan[PC_pred_Krishnan$quantile_K == i, 1], PC_pred_Krishnan[PC_pred_Krishnan$quantile_Krishnan == i, 1]))), nrow = 2)
  Fisher_Krishnan <- fisher.test(matrix_Krishnan, alternative = "greater") 
  
  assign(paste0("Fisher_K_", i), Fisher_Krishnan)
  fisher_df[i,] <- c("Krishnan", i, (matrix_Krishnan[1,1]+matrix_Krishnan[1,2]), (matrix_Krishnan[1,1]+matrix_Krishnan[1,2]),  matrix_Krishnan[1,1], Fisher_Krishnan$p.value)
}


#Lin
i = 1
for (i in i:10) {
  matrix_Lin <- matrix(c(length(intersect(PC_pred_Lin[PC_pred_Lin$quantile_K == i, 1], PC_pred_Lin[PC_pred_Lin$quantile_Lin == i, 1])),
                         length(setdiff(PC_pred_Lin[PC_pred_Lin$quantile_Lin == i, 1], PC_pred_Lin[PC_pred_Lin$quantile_K == i, 1])), 
                         length(setdiff(PC_pred_Lin[PC_pred_Lin$quantile_K == i, 1], PC_pred_Lin[PC_pred_Lin$quantile_Lin == i, 1])),
                         nrow(PC_pred_Lin) - length(union(PC_pred_Lin[PC_pred_Lin$quantile_K == i, 1], PC_pred_Lin[PC_pred_Lin$quantile_Lin == i, 1]))), nrow = 2)
  Fisher_Lin <- fisher.test(matrix_Lin, alternative = "greater") 
  
  assign(paste0("Fisher_L_", i), Fisher_Lin)
  fisher_df[(i+10),] <- c("Lin", i, (matrix_Lin[1,1]+matrix_Lin[1,2]), (matrix_Lin[1,1]+matrix_Lin[1,2]),  matrix_Lin[1,1], Fisher_Lin$p.value)
}



#Duda
i = 1
for (i in 1:10) {
  matrix_Duda <- matrix(c(length(intersect(PC_pred_Duda[PC_pred_Duda$quantile_K == i, 1], PC_pred_Duda[PC_pred_Duda$quantile_Duda == i, 1])),
                          length(setdiff(PC_pred_Duda[PC_pred_Duda$quantile_Duda == i, 1], PC_pred_Duda[PC_pred_Duda$quantile_K == i, 1])), 
                          length(setdiff(PC_pred_Duda[PC_pred_Duda$quantile_K == i, 1], PC_pred_Duda[PC_pred_Duda$quantile_Duda == i, 1])),
                          nrow(PC_pred_Duda) - length(union(PC_pred_Duda[PC_pred_Duda$quantile_K == i, 1], PC_pred_Duda[PC_pred_Duda$quantile_Duda == i, 1]))), nrow = 2)
  Fisher_Duda <- fisher.test(matrix_Duda, alternative = "greater") 
  
  assign(paste0("Fisher_D_", i), Fisher_Duda)
  fisher_df[(i+20),] <- c("Duda", i, (matrix_Duda[1,1]+matrix_Duda[1,2]), (matrix_Duda[1,1]+matrix_Duda[1,2]),  matrix_Duda[1,1], Fisher_Duda$p.value)}

#Brueggeman
i = 1
for (i in 1:10) {
  matrix_Brueggeman <- matrix(c(length(intersect(PC_pred_Brueggeman[PC_pred_Brueggeman$quantile_K == i, 1], PC_pred_Brueggeman[PC_pred_Brueggeman$quantile_Brueggeman == i, 1])),
                                length(setdiff(PC_pred_Brueggeman[PC_pred_Brueggeman$quantile_Brueggeman == i, 1], PC_pred_Brueggeman[PC_pred_Brueggeman$quantile_K == i, 1])), 
                                length(setdiff(PC_pred_Brueggeman[PC_pred_Brueggeman$quantile_K == i, 1], PC_pred_Brueggeman[PC_pred_Brueggeman$quantile_Brueggeman == i, 1])),
                                nrow(PC_pred_Brueggeman) - length(union(PC_pred_Brueggeman[PC_pred_Brueggeman$quantile_K == i, 1], PC_pred_Brueggeman[PC_pred_Brueggeman$quantile_Brueggeman == i, 1]))), nrow = 2)
  Fisher_Brueggeman <- fisher.test(matrix_Brueggeman, alternative = "greater") 
  
  assign(paste0("Fisher_B_", i), Fisher_Brueggeman)
  fisher_df[(i+30),] <- c("Brueggeman", i, (matrix_Brueggeman[1,1]+matrix_Brueggeman[1,2]), (matrix_Brueggeman[1,1]+matrix_Brueggeman[1,2]),  matrix_Brueggeman[1,1], Fisher_Brueggeman$p.value)
  }

fisher_df$pvalue <- as.numeric(fisher_df$pvalue)
fisher_df$p_adjusted <- p.adjust(fisher_df$pvalue, method = "bonferroni")

#Export table
write.csv(fisher_df, "./tables/Table_S3.csv", row.names = F)

