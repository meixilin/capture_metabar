# Title: Plot metabarcoding results of sites --------
# Author: Meixi Lin
# Date: Tue Oct 30 21:28:48 2018
# Author:
# Date:
# Modification:

# preparation --------
rm(list = ls())
cat("\014")

setwd("/Users/linmeixi/UCLA/Lab/Capture Array/metabar/seq_analysis/result_table/")
library(dplyr)
date() # the execution date

cutoff <- 0.95

# read in the data 12SV5 ---- 
asv.12SV5 <- read.delim("nochim_merged12SV5.txt"); 
colnames(asv.12SV5)[1] <- "id"
taxo.12SV5 <- read.delim("nochim_merged12SV5.tag.tsv")
joined.12SV5 <- left_join(asv.12SV5, taxo.12SV5, by = "id")

forplot.12SV5 <- joined.12SV5 %>% 
    filter(best_identity.12SV5_db > cutoff) %>% 
    select(id, starts_with("X12SV5"), starts_with("best_identity"), scientific_name, species_list.12SV5_db)

write.table(forplot.12SV5, file = "forplot.12SV5.tsv", quote = F, sep = "\t", row.names = T, col.names = NA)

# read in the data Mamm2 ----
asv.MAMM2 <- read.delim("nochim_mergedMAMM2.txt"); 
colnames(asv.MAMM2)[1] <- "id"
taxo.MAMM2 <- read.delim("nochim_mergedMAMM2.tag.tsv")
joined.MAMM2 <- left_join(asv.MAMM2, taxo.MAMM2, by = "id")

forplot.MAMM2 <- joined.MAMM2 %>% 
    filter(best_identity.MAMM2_db > cutoff) %>% 
    select(id, starts_with("MAMM2"), starts_with("best_identity"), scientific_name, species_list.MAMM2_db)

write.table(forplot.MAMM2, file = "forplot.MAMM2.tsv", quote = F, sep = "\t", row.names = T, col.names = NA)

