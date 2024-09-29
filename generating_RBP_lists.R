# ---
# Title: generating_RBP_lists.R
# Purpose: This script uses the COMPILED_TABLE from RBPbase (re-run with adjustments described in Materials & Methods)
#         to generate the 4 RBP lists (comprehensive, high-confidence, canonical and non-canonical)
#          
# ---

#If you ran make_ASiD_dataset_PC.R in the same R session, detach and reload the dplyr package. 
# detach("package:dplyr", unload=TRUE)

library(stringr)
library(dplyr)

setwd("") # Change this to the path that contains the ASiD folder



load("./data/COMPILED_TABLE.Rda")
collapse_df <- COMPILED_TABLE[["Hs"]]
collapse_df <- as.data.frame(collapse_df)


# UNION OF EXPERIMENTS ----------------------------------------------------
#Combine genes into one column if they were from similar experiments. This counts as 1 RIC. 
#Perez-Perri 34 and 35 are 1 experiment
collapse_df$sum_PP <- collapse_df %>% select(RBPBASE000000034.1, RBPBASE000000035.1) %>% rowSums(na.rm = TRUE)
collapse_df <- mutate(collapse_df, sum_PP = if_else(condition = sum_PP > 0, true = 1, false = 0))
  nrow(filter(collapse_df, sum_PP == 1)) #690
collapse_df$sum_PP <- as.logical(collapse_df$sum_PP)

#Backlund 38, 39, 40 and 41 are 1 experiment
collapse_df$sum_Backlund <- collapse_df %>% select(RBPBASE000000038.1, RBPBASE000000039.1, RBPBASE000000040.1, RBPBASE000000041.1) %>% rowSums(na.rm = TRUE)
collapse_df <- mutate(collapse_df, sum_Backlund = if_else(condition = sum_Backlund > 0, true = 1, false = 0))
  nrow(filter(collapse_df, sum_Backlund == 1)) #845
collapse_df$sum_Backlund <- as.logical(collapse_df$sum_Backlund)

#Urdaneta 49, 50, and 51 are 1 experiment
collapse_df$sum_Urdaneta <- collapse_df %>% select(RBPBASE000000049.1, RBPBASE000000050.1, RBPBASE000000051.1) %>% rowSums(na.rm = TRUE)
collapse_df <- mutate(collapse_df, sum_Urdaneta = if_else(condition = sum_Urdaneta > 0, true = 1, false = 0))
  nrow(filter(collapse_df, sum_Urdaneta == 1)) #3030
collapse_df$sum_Urdaneta <- as.logical(collapse_df$sum_Urdaneta)

#Trendel 59 and 62 are 1 experiment (60, and 61 will each count as their own)
collapse_df$sum_Trendel <- collapse_df %>% select(RBPBASE000000059.1, RBPBASE000000062.1) %>% rowSums(na.rm = TRUE)
collapse_df <- mutate(collapse_df, sum_Trendel = if_else(condition = sum_Trendel > 0, true = 1, false = 0))
  nrow(filter(collapse_df, sum_Trendel == 1)) #1208
collapse_df$sum_Trendel <- as.logical(collapse_df$sum_Trendel)

#Garcia-Moreno 66 and 67 will count as 1 experiment
collapse_df$sum_GarM <- collapse_df %>% select(RBPBASE000000066.1, RBPBASE000000067.1) %>% rowSums(na.rm = TRUE)
collapse_df <- mutate(collapse_df, sum_GarM = if_else(condition = sum_GarM > 0, true = 1, false = 0))
  nrow(filter(collapse_df, sum_GarM == 1)) #752
collapse_df$sum_GarM <- as.logical(collapse_df$sum_GarM)


# CALCULATE NUMBER OF RIC EXPERIMENTS USING ABOVE COLLAPSED EXPERIMENTS----------------------------------------
row_list = c()
i = 1
for (i in 1:nrow(collapse_df)) {
  if (isTRUE(collapse_df[i, "sum_PP"] |
             collapse_df[i, "sum_Backlund"] |
             collapse_df[i, "sum_Urdaneta"] |
             collapse_df[i, "sum_Trendel"] |
             collapse_df[i, "sum_GarM"] |
             collapse_df[i,"RBPBASE000000007.1"] |
             collapse_df[i,"RBPBASE000000008.1"] |
             collapse_df[i,"RBPBASE000000009.1"] | 
             collapse_df[i, "RBPBASE000000010.1"] |
             collapse_df[i, "RBPBASE000000011.1"] |
             collapse_df[i,"RBPBASE000000012.1"] |
             collapse_df[i,"RBPBASE000000013.1"] |
             collapse_df[i, "RBPBASE000000032.1"] |
             collapse_df[i, "RBPBASE000000033.1"] |
             collapse_df[i, "RBPBASE000000036.1"] |
             collapse_df[i,"RBPBASE000000037.1"] | 
             collapse_df[i,"RBPBASE000000046.1"] | 
             collapse_df[i,"RBPBASE000000047.1"] | 
             collapse_df[i,"RBPBASE000000048.1"] | 
             collapse_df[i, "RBPBASE000000052.1"] |
             collapse_df[i,"RBPBASE000000060.1"] |
             collapse_df[i,"RBPBASE000000061.1"] |
             collapse_df[i, "RBPANNO000000001.1"] |
             collapse_df[i, "RBPANNO000000002.1"] |
             collapse_df[i, "RBPANNO000000017.1"] |
             collapse_df[i, "RBPANNO000000023.1"] |
             collapse_df[i, "RBPANNO000000028.1"] |
             collapse_df[i, "RBPANNO000000038.1"] == "TRUE")) {
    row_list = append(row_list, i)
  }
}

new_df <- collapse_df %>% filter(row_number() %in% row_list) 

new_df$sum_hits <- new_df %>% select(sum_PP, sum_Backlund,
                                     sum_Urdaneta, sum_Trendel, 
                                     sum_GarM, RBPBASE000000007.1,
                                     RBPBASE000000008.1, RBPBASE000000009.1,
                                     RBPBASE000000010.1, RBPBASE000000012.1,
                                     RBPBASE000000013.1, RBPBASE000000032.1,
                                     RBPBASE000000033.1,  RBPBASE000000036.1,
                                     RBPBASE000000037.1, RBPBASE000000046.1,
                                     RBPBASE000000047.1, RBPBASE000000048.1,
                                     RBPBASE000000060.1, RBPBASE000000061.1) %>% 
  rowSums(na.rm = TRUE) #max number is 20


#Simplify columns in df. 
new_df <- select(new_df, UNIQUE,
                 ID, Description,
                 RBPBASE000000007.1, RBPBASE000000008.1,
                 RBPBASE000000009.1, RBPBASE000000010.1,
                 RBPBASE000000011.1, RBPBASE000000012.1,
                 RBPBASE000000013.1, RBPBASE000000032.1,
                 RBPBASE000000033.1, RBPBASE000000034.1,
                 RBPBASE000000035.1, RBPBASE000000036.1,
                 RBPBASE000000037.1, RBPBASE000000038.1,
                 RBPBASE000000039.1, RBPBASE000000040.1,
                 RBPBASE000000041.1, RBPBASE000000046.1,
                 RBPBASE000000047.1, RBPBASE000000048.1,
                 RBPBASE000000049.1, RBPBASE000000050.1,
                 RBPBASE000000051.1, RBPBASE000000052.1,
                 RBPBASE000000059.1, RBPBASE000000060.1,
                 RBPBASE000000061.1, RBPBASE000000062.1,
                 RBPBASE000000066.1, RBPBASE000000067.1,
                 RBPANNO000000001.1, RBPANNO000000002.1,
                 RBPANNO000000005.1, RBPANNO000000008.1,
                 RBPANNO000000011.1, RBPANNO000000014.1,
                 RBPANNO000000017.1, RBPANNO000000023.1,
                 RBPANNO000000028.1, RBPANNO000000034.1,
                 RBPANNO000000038.1, RBPANNO000000043.1,
                 RBPANNO000000046.1, RBPANNO000000047.1,
                 RBPANNO000000048.1, RBPANNO000000057.1,
                 RBPANNO000000058.1, RBPANNO000000059.1,
                 RBPANNO000000060.1, RBPANNO000000061.1,
                 RBPANNO000000062.1, RBPANNO000000063.1,
                 RBPANNO000000064.1, RBPANNO000000065.1,
                 RBPANNO000000070.1, RBPANNO000000071.1,
                 RBPANNO000000072.1, RBPANNO000000074.1,
                 RBPANNO000000077.1, RBPANNO000000078.1,
                 any_Hs, hits_Hs, 
                 sum_PP, sum_Backlund, 
                 sum_Urdaneta, sum_Trendel,
                 sum_GarM, sum_hits)


# PRODUCING 4 RBP LISTS -------------------------------------------------
#>= 2 RIC
RIC2 = c()
i = 1
for (i in 1:nrow(new_df)) {
    if (new_df[i,71] >= 2) {
    RIC2 <- append(RIC2, new_df[i,1])
  }
} #2776


#Annotated. Includes any RBP from RBPANNO000000001.1, RBPANNO000000023.1 or RBPANNO000000028.1
Anno = c()
i = 1
for (i in 1:nrow(new_df)) {
  if (isTRUE(new_df[i, "RBPANNO000000001.1"] |
             new_df[i, "RBPANNO000000023.1"] |
             new_df[i, "RBPANNO000000028.1"] == "TRUE")) {
    Anno <- append(Anno, new_df[i,1])
  }
} #1540


Comp_RBPs_vec <- union(RIC2, Anno)
Comp_RBPs <- new_df[new_df$UNIQUE %in% Comp_RBPs_vec, ]

Canonical_RBPs <- Comp_RBPs[Comp_RBPs$RBPANNO000000017.1, ]

Noncanonical_RBPs_vec <- setdiff(Comp_RBPs_vec, Canonical_RBPs$UNIQUE)
Noncanonical_RBPs <- new_df[new_df$UNIQUE %in% Noncanonical_RBPs_vec, ]

Highconfidence_RBPs_vec <- intersect(RIC2, Anno)
Highconfidence_RBPs <- new_df[new_df$UNIQUE %in% Highconfidence_RBPs_vec, ]

df_list <- list(Comp_RBPs = Comp_RBPs,
                Canonical_RBPs = Canonical_RBPs,
                Noncanonical_RBPs = Noncanonical_RBPs,
                Highconfidence_RBPs = Highconfidence_RBPs)


#Export 4 dfs. 
i = 1
for(i in 1:length(df_list)){
  temp_df <- as.data.frame(df_list[[i]])
  temp_df <- temp_df[,c(1:3, 40, 46:48)] 
  temp_df[,"RBPANNO000000017.1"] <- as.numeric(temp_df[,"RBPANNO000000017.1"])
  names(temp_df) <- c("Gene name",
                           "Ensembl ID",
                           "Gene description",
                           "Contains canonical RBD",
                           "Pfam ID",
                           "Protein domains",
                           "Pfam Description")

  write.csv(temp_df, str_glue("./RBP_lists/{names(df_list)[i]}.csv"), row.names = F)
}

#Export details of 3,300 RBPs from RBPbase:
#write.csv(Comp_RBPs[,c(1:35, 40:42, 44)], "./tables/RBP_RIC_Annotation.csv", row.names = F)


