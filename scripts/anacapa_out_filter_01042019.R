# Title: read in and subset anacapa output from capture metabar --------
# Author: Meixi Lin
# Date: Wed Nov 28 10:41:30 2018
# Author:
# Date:
# Modification:

# preparation --------
rm(list = ls())
cat("\014")

setwd("~/UCLA/Lab/Capture Array/metabar/seq_analysis")
library(ranacapa)
library(dplyr)
library(reshape2)
library(phyloseq)
library(vegan)
date() # the execution date

primers <- c("12SV5", "MAMM2", "PITS", "18S", "CO1")
samples <- c("FOD.DNA.1", "MVP.DNA.2", "MVP.DNA.3", "pcr.blank", "positive.mammal", "PPB.DNA.1", "SMM.DNA.1", "SMV.DNA.1", "SSA.DNA.1", "extract.blank")

# read in the data --------
# read raw asv
asv <- read.delim(file = "anacapa_out/12SV5_taxonomy_tables/12SV5_ASV_taxonomy_detailed.txt", stringsAsFactors = F)
samplename <- colnames(asv)[-(1:4)] %>%
    strsplit(split = "_") %>%
    lapply(function(xx) {xx <- xx[2]}) %>%
    unlist

# subset the asv by 
# 1. length 
# 2. taxonomy confidence 
asvnew <- asv %>%
    filter(input_sequence_length > 30)
    

# subset the data fromwe 

names(sumlist) <- primers
# read metadata 
biom <- read.csv(file = "anacapa_out/metadata.csv", stringsAsFactors = F)

# build phyloseq from sum.taxonomy --------
phylist <- lapply(sumlist, function(xx) {
    physeq <- convert_anacapa_to_phyloseq(xx, biom)
})
phylist

# perform decontamination --------
# filter out singletons and doubletons 
phylist3 <- lapply(phylist, function(xx) {
    physeq <- filter_taxa(xx, function(x) sum(x) > 2, T)
    return(physeq)
})
phylist3
lapply(phylist3, function(xx) {sample_sums(xx)})

# filter out 
# 1. any taxonomical entries if they appeared in the blanks more than 3 times (total)
phylist3d <- lapply(phylist3, function(xx) {
    # check if contain negative control
    xxmeta <- data.frame(sample_data(xx))[,"sample_type"]
    if ("neg_control" %in% xxmeta) {
        blanksum <- xx %>%
            subset_samples(sample_type == "neg_control") %>%
            taxa_sums
        goodtaxa <- names(blanksum[blanksum <= 2])
        xxd <- prune_taxa(goodtaxa, xx)
    } else {
        xxd <- xx
    }
    samplesum <- sample_sums(xxd)
    goodsample <- names(samplesum[samplesum > 0])
    xxds <- prune_samples(goodsample, xxd)
    return(xxds)
})
phylist3d
lapply(phylist3d, function(xx) {sample_sums(xx)})
save(phylist3d, file = "phylist3d.RData")


# phylist3d --------
# plot_barplot 
plot_bar(phylist3d$`12SV5`, fill = "Species")
plot_bar(phylist3d$`MAMM2`, fill = "Species")
plot_bar(phylist3d$`PITS`, fill = "Species")
physeq <- phylist3d$PITS
TopNOTUs <- names(sort(taxa_sums(physeq), TRUE)[1:10])
phy10   <- prune_taxa(TopNOTUs, physeq)
# you can also convert to relative abundance 
phy10re <- transform_sample_counts(phy10, function(xx){xx / sum(xx)})
pp <- plot_bar(phy10re, fill = "Species")
plot_bar(phy10, fill = "Family") 
plot_bar(phylist3d$`CO1`, fill = "Genus")

# plot ordination for 18S samples 
physeq1 <- phylist3d$`18S`
sample_sums(physeq1)
set.seed(17)
ord.res <- ordinate(physeq1, method = "NMDS", distance = "jaccard", binary = T)
pp <- plot_ordination(physeq1, ord.res, type = "samples", color = "loc") +
    scale_color_brewer(type = "qual", palette = 1) +
    ggtitle("NMDS on 18S, jaccard dist")
ord.res.2 <- ordinate(physeq1, method = "NMDS", distance = "bray")
pp <- plot_ordination(physeq1, ord.res.2, type = "samples", color = "loc") +
    geom_point(shape = 1) + 
    scale_color_brewer(type = "qual", palette = 1) +
    ggtitle("NMDS on 18S, bray dist")

# plot qq plot for 18S samples (mvp)
forplot <- as.data.frame(otu_table(physeq1)) %>%
    tibble::rownames_to_column(var = "sum.taxonomy") %>%
    glom_tax_df(taxlevel = "Family")
forplot <- forplot[-which(rownames(forplot) == "unknown"),]
xylim <- range(forplot[,2], forplot[,3])
plot(forplot[,2], forplot[,3], xlim = xylim, ylim = xylim)
abline(a = 0, b = 1, col = "red")


# perform rarefaction --------
load(file = "phylist3d.RData")
# settings #####
rarefaction_depth <- 10000
rarefaction_reps  <- 10
seed <- 7

# perform ####
phylist3dr <- lapply(phylist3d, function(xx) {
    xxr <- custom_rarefaction_2(xx, sample_size = rarefaction_depth, replicates = rarefaction_reps, myseed = seed)
    return(xxr)
})
xxr.18S <- custom_rarefaction_2(phylist3d, sample_size = 10000, replicates = rarefaction_reps, myseed = seed)
xxr.12SV5 <- custom_rarefaction_2(phylist3d$`12SV5`, sample_size = 30000, replicates = rarefaction_reps, myseed = seed)


#Check rarefaction, CHANGE THIS LINE BY OBJECT ####
physeq_obj <- phylist3d$`12SV5`
physeq_obj_rare <- xxr.12SV5
sample_sums(physeq_obj_rare)
hist(sample_sums(physeq_obj_rare))
dim(otu_table(physeq_obj_rare))
dim(otu_table(physeq_obj))

# OPTIONAL: plot rarefaction
p <- ggrare(physeq_obj_rare,step = 100, color = "index", se=FALSE) + theme_bw() + theme_ranacapa()
p <- p + ggtitle('Rarified Samples') +
    theme(plot.title = element_text(hjust = 0.5)) + theme_ranacapa() +
    theme(panel.grid.minor.y=element_blank(),panel.grid.minor.x=element_blank(),panel.grid.major.y=element_blank(),panel.grid.major.x=element_blank(), legend.position = "none")
p

# continue analysis on phylist3d --------
load("phylist3d.RData")
# glom to family level and convert to boolean 
dffam <- lapply(phylist3d, function(xx) {
    xx <- as.data.frame(otu_table(xx)) %>%
        tibble::rownames_to_column(var = "sum.taxonomy") %>%
        glom_tax_df(., taxlevel = "Family")
    xx <- (xx > 0) * 1 
})
# plot the family level result 
phyfam <- lapply(dffam, function(xx) {
    OTU <- otu_table(xx, taxa_are_rows = TRUE)
    taxmat <- matrix(rownames(xx), nrow = nrow(xx), ncol = 1)
    colnames(taxmat) <- "Family"
    rownames(taxmat) <- rownames(xx)
    TAX <- tax_table(taxmat)
    phy <- phyloseq(OTU, TAX)
})
phyloseq::plot_heatmap(phyfam$`18S`, trans = NULL, method = NULL)
# only take the mvp 
dfmvp <- lapply(dffam, function(xx) {
    xx <- cbind(xx[,c("MVP.DNA.2", "MVP.DNA.3")])
})
gplots::venn(dfmvp$`18S`)
