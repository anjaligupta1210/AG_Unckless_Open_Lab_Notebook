##All plots and data analysis for Evol2024 talk

#Set working directory
setwd("/Users/anjaligupta/Library/CloudStorage/GoogleDrive-anjaligupta1210@gmail.com/.shortcut-targets-by-id/0B-qq8zQGecn_Y0Y5ZENITnlUMVU/UncklessLab/Personnel/Anjali/Talks/Evol2024/DataAndCode")

#Synteny plots
source("asynt.R")

#import alignments
ST_SR1 <- import.paf("map_DaffSTfemale5.1.main_to_DaffSR.X.remapped.qry.paf")
ST_SR2 <- import.paf("map_DaffSTfemale5.1.main_to_DaffSR2.STmapped.qry.paf")
SR1_SR2 <- import.paf("map_DaffSR2.STmapped.qry_to_DaffSR.X.remapped.qry.paf")

#Next we import load scaffold length data which is necessary to plot the scaffolds
ST <- import.genome(fai_file="DaffSTfemale5.1.main.fasta.fai")
SR1 <- import.genome(fai_file="DaffSR.X.remapped.qry.fasta.fai")
SR2 <- import.genome(fai_file = "DaffSR2.STmapped.qry.fasta.fai")

#now define the scaffolds we're interested in
reference_scafs <- "ChrX_MullerAD"
query_scafs <- "ChrX_MullerAD"

#keep only alignments involving these scaffolds
ST_SR1 <- subset(ST_SR1, 
                     query %in% query_scafs & 
                       reference %in% reference_scafs)
ST_SR2 <- subset(ST_SR2, 
                 query %in% query_scafs & 
                   reference %in% reference_scafs)
SR1_SR2 <- subset(SR1_SR2, 
                 query %in% query_scafs & 
                   reference %in% reference_scafs)


#multiple scaffold plot of synetny blocks 
synblocks_ST_SR1 <- get.synteny.blocks.multi(ST_SR1, min_subblock_size=5000)
synblocks_ST_SR2 <- get.synteny.blocks.multi(ST_SR2, min_subblock_size=5000)
synblocks_SR1_SR2 <- get.synteny.blocks.multi(SR1_SR2, min_subblock_size=5000)

plot.alignments.multi(synblocks_ST_SR1, 
                      reference_lens=SR1$seq_len, 
                      query_lens=ST$seq_len, 
                      sigmoid=T,
                      reference_above = TRUE, 
                      show_labels = FALSE, 
                      show_outline = FALSE, 
                      show_alignment_tracts=TRUE)

plot.alignments.multi(synblocks_ST_SR2, 
                      reference_lens=SR2$seq_len, 
                      query_lens=ST$seq_len, 
                      sigmoid=T,
                      reference_above = TRUE, 
                      show_labels = FALSE, 
                      show_outline = FALSE, 
                      show_alignment_tracts=TRUE)

plot.alignments.multi(synblocks_SR1_SR2, 
                      reference_lens=SR1$seq_len, 
                      query_lens=SR2$seq_len, 
                      sigmoid=T,
                      reference_above = TRUE, 
                      show_labels = FALSE, 
                      show_outline = FALSE, 
                      show_alignment_tracts=TRUE)


### Trees

#Load required packages
library(ape)
library(phangorn)
library(spatstat)
library(spatstat.geom)
library(phytools)


EntireX <- phytools::read.newick("GATKfilteredSNPchrX_NoAztNoBlanks.min4.phy.treefile")

#Root using pseudoobscura
EntireX<-root(EntireX,"Dpse")

##make ultrametric

EntireX<-chronos(EntireX)

plot(EntireX, type="phylogram",nodes="intermediate", alpha = .05)
nodelabels(EntireX$node.label, adj = c(1.5, -0.5), frame = "none", cex = 3, pos=1, col = "red")

##Autosomes
Autosomes <- phytools::read.newick("GATKfilteredSNPAutosome_NoAztNoBlanks.min4.phy.treefile")

#Root using pseudoobscura
Autosomes<-root(Autosomes,"Dpse")

##make ultrametric

Autosomes<-chronos(Autosomes)

plotTree(Autosomes, type="cladogram",nodes="intermediate", alpha = .05)
nodelabels(Autosomes$node.label, adj = c(1.5, -0.5), frame = "none", cex = 3, pos=2, col = "red")


##XL vs XR

XL <- phytools::read.newick("chrX_XL.min4.phy.treefile")

#Root using pseudoobscura
XL<-root(XL,"Dpse")

##make ultrametric

XL<-chronos(XL)

plotTree(XL, type="cladogram",nodes="intermediate", alpha = .05)
nodelabels(XL$node.label, adj = c(1.5, -0.5), frame = "none", cex = 3, pos=1, col = "red")


XR <- phytools::read.newick("chrX_XR.min4.phy.treefile")

#Root using pseudoobscura
XR<-root(XR,"Dpse")

##make ultrametric

XR<-chronos(XR)

plotTree(XR, type="cladogram",nodes="intermediate", alpha = .05)
nodelabels(XR$node.label, adj = c(1.5, -0.5), frame = "none", cex = 3, pos=1, col = "red")


#Inversions

##XL 

XL_PreInv <- phytools::read.newick("chrX_XL_PreInv.min4.phy.treefile")

#Root using pseudoobscura
XL_PreInv<-root(XL_PreInv,"Dpse")

##make ultrametric

XL_PreInv<-chronos(XL_PreInv)

plotTree(XL_PreInv, type="cladogram",nodes="intermediate", alpha = .05)
nodelabels(XL_PreInv$node.label, adj = c(1.5, -0.5), frame = "none", cex = 3, pos=1, col = "red")


##
XL_Inv <- phytools::read.newick("chrX_XL_Inv.min4.phy.treefile")

#Root using pseudoobscura
XL_Inv<-root(XL_Inv,"Dpse")

##make ultrametric

XL_Inv<-chronos(XL_Inv)

plotTree(XL_Inv, type="cladogram",nodes="intermediate", alpha = .05)
nodelabels(XL_Inv$node.label, adj = c(1.5, -0.5), frame = "none", cex = 3, pos=1, col = "red")

##
XL_PostInv <- phytools::read.newick("chrX_XL_PostInv.min4.phy.treefile")

#Root using pseudoobscura
XL_PostInv<-root(XL_PostInv,"Dpse")

##make ultrametric

XL_PostInv<-chronos(XL_PostInv)

plotTree(XL_PostInv, type="cladogram",nodes="intermediate", alpha = .05)
nodelabels(XL_PostInv$node.label, adj = c(1.5, -0.5), frame = "none", cex = 3, pos=1, col = "red")



##XR 

XR_PreInv <- phytools::read.newick("chrX_XR_PreInv.min4.phy.treefile")

#Root using pseudoobscura
XR_PreInv<-root(XR_PreInv,"Dpse")

##make ultrametric

XR_PreInv<-chronos(XR_PreInv)

plotTree(XR_PreInv, type="cladogram",nodes="intermediate", alpha = .05)
nodelabels(XR_PreInv$node.label, adj = c(1.5, -0.5), frame = "none", cex = 3, pos=1, col = "red")


##
XR_Inv <- phytools::read.newick("chrX_XR_Inv.min4.phy.treefile")

#Root using pseudoobscura
XR_Inv<-root(XR_Inv,"Dpse")

##make ultrametric

XR_Inv<-chronos(XR_Inv)

plotTree(XR_Inv, type="cladogram",nodes="intermediate", alpha = .05)
nodelabels(XR_Inv$node.label, adj = c(1.5, -0.5), frame = "none", cex = 3, pos=1, col = "red")

##
XR_PostInv <- phytools::read.newick("chrX_XR_PostInv.min4.phy.treefile")

#Root using pseudoobscura
XR_PostInv<-root(XR_PostInv,"Dpse")

##make ultrametric

XR_PostInv<-chronos(XR_PostInv)

plotTree(XR_PostInv, type="cladogram",nodes="intermediate", alpha = .05)
nodelabels(XR_PostInv$node.label, adj = c(1.5, -0.5), frame = "none", cex = 3, pos=1, col = "red")



###3MB WINDOW Trees
library(readr)

# List of treefiles
treefiles <- list.files(path = "3mbTreesX/",
                        pattern = "ChrX_MullerAD_[0-9]+_[0-9]+\\.min4\\.phy\\.treefile$",
                        full.names = TRUE)

# Initialize an empty list to store data frames for each treefile
WindowTrees_list <- list()

# Root and label trees
for (treefile in treefiles) {
  # Root tree
  tree <- phytools::read.newick(treefile)
  rooted_tree <- root(tree,"Dpse")
  rooted_tree <- chronos(rooted_tree)
  
  rooted_tree$node.label <- NULL
  rooted_tree$edge.length <- NULL
  
  # Extract file name
  file_name <- basename(treefile)
  # Extract topology
  topology <- write.tree(rooted_tree)
  # Extract Position
  positions <- as.numeric(unlist(strsplit(gsub("ChrX_MullerAD_|\\.min4\\.phy\\.treefile", "", file_name), "_")))
  Start <- positions[1]
  End <- positions[2]
  # Create a data frame for the current treefile
  tree_data <- data.frame(File = file_name, topology = topology, Start = Start, End = End)
  # Append the data frame to the list
  WindowTrees_list[[length(WindowTrees_list) + 1]] <- tree_data
}

# Combine all data frames into one
WindowTrees <- data.frame()
WindowTrees <- do.call(rbind, WindowTrees_list)

WindowTrees$topology <- as.factor(WindowTrees$topology)

Log <- read_csv("3mbTreeslogLikelihood.csv")

WindowTrees <- merge(WindowTrees,Log, by = "File")

# Define monophyletic and not monophyletic topologies
monophyletic_topologies <- c(
  "((((Daff_SR1,Daff_SR2),Daff_ST),Dalg),(Datha_ea,Datha_eb),Dpse);",
  "((((Daff_SR1,Daff_ST),Daff_SR2),(Datha_ea,Datha_eb)),Dalg,Dpse);",
  "((((Daff_SR1,Daff_ST),Daff_SR2),Dalg),(Datha_ea,Datha_eb),Dpse);",
  "(((Daff_SR1,Daff_SR2),Daff_ST),(Dalg,(Datha_ea,Datha_eb)),Dpse);",
  "(((Daff_SR1,Daff_ST),Daff_SR2),(Dalg,(Datha_ea,Datha_eb)),Dpse);",
  "((Daff_SR1,(Daff_SR2,Daff_ST)),(Dalg,(Datha_ea,Datha_eb)),Dpse);"
)

not_monophyletic_topologies <- c(
)

Mono_SR1SR2_close <- c(
  "((((Daff_SR1,Daff_SR2),Daff_ST),Dalg),(Datha_ea,Datha_eb),Dpse);",
  "(((Daff_SR1,Daff_SR2),Daff_ST),(Dalg,(Datha_ea,Datha_eb)),Dpse);"
  
)

Mono_SR1ST_close <- c(
  "((((Daff_SR1,Daff_ST),Daff_SR2),(Datha_ea,Datha_eb)),Dalg,Dpse);",
  "((((Daff_SR1,Daff_ST),Daff_SR2),Dalg),(Datha_ea,Datha_eb),Dpse);",
  "(((Daff_SR1,Daff_ST),Daff_SR2),(Dalg,(Datha_ea,Datha_eb)),Dpse);"
)

Mono_SR2ST_close <- c(
  "((Daff_SR1,(Daff_SR2,Daff_ST)),(Dalg,(Datha_ea,Datha_eb)),Dpse);"
)


# Add monophyly column based on topology
WindowTrees$monophyly <- NA
WindowTrees$monophyly[WindowTrees$topology %in% monophyletic_topologies] <- "D. affinis is monophyletic"
WindowTrees$monophyly[WindowTrees$topology %in% not_monophyletic_topologies] <- "D. affinis is polyphyletic"


WindowTrees$Mono <- NA
WindowTrees$Mono[WindowTrees$topology %in% not_monophyletic_topologies] <- "D. affinis is polyphyletic"
WindowTrees$Mono[WindowTrees$topology %in% Mono_SR1SR2_close] <- "((Daff_SR1,Daff_SR2),Daff_ST)"
WindowTrees$Mono[WindowTrees$topology %in% Mono_SR1ST_close] <- "((Daff_SR1,Daff_ST),Daff_SR2)"
WindowTrees$Mono[WindowTrees$topology %in% Mono_SR2ST_close] <- "((Daff_SR1,(Daff_SR2,Daff_ST))"





WindowTrees$Region <- ifelse(WindowTrees$Start<5490855,"XL_PreInv",
                             ifelse(WindowTrees$Start<10233152 &
                                      WindowTrees$Start>5490855, "XL_Inv",
                                    ifelse(WindowTrees$Start>10233152 &
                                             WindowTrees$Start<30000000, "XL_PostInv",
                                           ifelse(WindowTrees$Start>30000000 &
                                                    WindowTrees$Start<44023827, "XR_PreInv",
                                                  ifelse(WindowTrees$Start<62734159 &
                                                           WindowTrees$Start>44023827, "XR_Inv",
                                                         ifelse(WindowTrees$Start>62734159, "XR_PostInv", NA))))))

WindowTrees$Region <- as.factor(WindowTrees$Region)


levels(WindowTrees$Region) <- c("XL_PreInv",
                                "XL_Inv",
                                "XL_PostInv",
                                "XR_PreInv",
                                "XR_Inv",
                                "XR_PostInv")

library(ggplot2)
ggplot(WindowTrees, aes(x = as.numeric(Start),
                        y = as.numeric(`Log-likelihood of the tree`),
                        fill = Mono)) +
  geom_bar(stat = "identity") +
  labs(x = "X chromosome position",
       y = "Log-likelihood of the tree",
       fill = "") +
  theme_minimal() +
  theme(legend.position = "top",
        title = element_text(size = 14, face = "bold"),
        text = element_text(size = 14, face = "bold")) +
  geom_vline(xintercept = c(5490855,
                            10233152,
                            30000000,
                            44023827,
                            62734159),
             linetype = "dashed",
             size = 1.6) +
  ggtitle("3MB windows")


## PCA

#Get Sample ID's from vcf
snps <- vcfR::read.vcfR("chrX_XR_PostInv.vcf", 
                        convertNA  = TRUE)
snps_num <- vcfR::extract.gt(snps, 
                             element = "GT",
                             IDtoRowNames  = F,
                             convertNA = T,
                             return.alleles = F)
snps_num_t <- t(snps_num)
snps_num_df <- data.frame(snps_num_t) 
sample_id <- row.names(snps_num_df)

file_plot_pairs <- list(
  list("pcangsd_filtered_SNPs_GATK_PopGen_ChrX_haploid.cov", "EntireX"),
  list("pcangsdchrX_XL_Inv.cov", "XL_Inv"),
  list("pcangsdchrX_XL_PostInv.cov", "XL_PostInv"),
  list("pcangsdchrX_XL_PreInv.cov", "XL_PreInv"),
  list("pcangsdchrX_XL.cov", "XL"),
  list("pcangsdchrX_XR_Inv.cov", "XR_Inv"),
  list("pcangsdchrX_XR_PostInv.cov", "XR_PostInv"),
  list("pcangsdchrX_XR_PreInv.cov", "XR_PreInv"),
  list("pcangsdchrX_XR.cov", "XR"),
  list("pcangsd_filtered_SNPs_GATK_PopGen_Autosomes_diploid_Chr2_MullerE.cov", "Chr2_MullerE"),
  list("pcangsd_filtered_SNPs_GATK_PopGen_Autosomes_diploid_Chr2.group4_MullerE.cov", "Chr2.group4_MullerE"),
  list("pcangsd_filtered_SNPs_GATK_PopGen_Autosomes_diploid_Chr3_MullerC.cov", "Chr3_MullerC"),
  list("pcangsd_filtered_SNPs_GATK_PopGen_Autosomes_diploid_Chr4_MullerB.cov", "Chr4_MullerB"),
  list("pcangsd_filtered_SNPs_GATK_PopGen_Autosomes_diploid_Chr5_MullerF.cov", "Chr5_MullerF"),
  list("pcangsd_filtered_SNPs_GATK_PopGen_Autosomes_diploid_mtDNA.cov", "mtDNA"),
  list("pcangsd_filtered_SNPs_GATK_PopGen_Autosomes_diploid_Unknown_69.cov", "Unkown69")
)

for(pair in file_plot_pairs){

CovMat <- pair[[1]]
  
X <- as.matrix(read.table(CovMat)) # Reads estimated covariance matrix

# Plot PCA plot
e <- eigen(X)
PCA <- as.data.frame(e$vectors)
PCA$ind <- c(sample_id)

#Get sex ratio information
library(readr)
AffinisPopGenData <- read_csv("AffinisPopGenData.csv")
popdata <- read_csv("affinis_pop_samples - Sheet1 (1).csv")

popdata <- popdata[,c(1,4)]

AffinisPopGenData <- merge(AffinisPopGenData,
                           popdata,
                           by="Sample")

AffinisPopGenData$ind <- AffinisPopGenData$IlluminaPoolName
New <- merge(PCA, AffinisPopGenData, by="ind", all.x = TRUE)

New$pop <- gsub("[01-9]*",
                "",
                New$ind)

library(ggplot2)
library(ggrepel)

New[,81:83][is.na(New[,81:83])] <- 0

# Initialize a vector to store the p-values
p_values <- numeric(nrow(New))

# Perform the binomial test for each row
for (i in 1:nrow(New)) {
  female_counts <- New[i, 81]
  male_counts <- New[i, 82]
  total_count <- New[i, 83]
  if(male_counts+female_counts>0){
    total_count <- female_counts + male_counts
    # We can use binom.test to check if the proportion of females (successes) is significantly different from 0.5
    binom_test <- binom.test(female_counts, total_count, p = 0.5)
    p_values[i] <- binom_test$p.value
  }
  else
    p_values[i] <- NA
}

# Add the p-values to the dataframe
New$p_value <- p_values

# Define significance levels
New$significance <- ifelse(New$p_value < 0.001, "***", 
                           ifelse(New$p_value < 0.01, "**", 
                                  ifelse(New$p_value < 0.05, "*", "ns")))


New$SR <- ifelse(New$significance=="ns",
                 "50/50",
                 ifelse(!New$significance=="ns" & 
                          New$PropFemale<0.5,
                        "50/50",
                        ifelse(New$significance%in%c("**",
                                                     "***") &
                                 New$PropFemale>0.5,
                               "SR",
                               "NA")))

IDaccX <- read_csv("IDaccX.csv")

New <- merge(New, IDaccX, by="ind", all.x = TRUE)

pca.eigenval.sum = sum(abs(e$values)) #sum of eigenvalues

varPC1 <- round((e$values[1]/pca.eigenval.sum)*100,
                digits = 2) 
varPC2 <- round((e$values[2]/pca.eigenval.sum)*100,
                digits = 2) 


plot <- ggplot(New, aes(x=V1,
                y=-V2,
                col=as.factor(Collection))) +
  geom_label_repel(data = subset(New,
                                 ind%in%c("Daff_SR1",
                                          "Daff_SR2",
                                          "Daff_ST")),
                   aes(label=ind),
                   col="blue",
                   hjust=-0.5,
                   vjust=0.8) +
  geom_point(aes(pch=X),
             size=3,
             alpha=0.7) +
  theme_bw() +
  theme(legend.position = "top") +
  ggtitle(paste(pair[[2]])) +
  labs(x=paste("PC1=",varPC1,"%"),
       y=paste("PC2=",varPC2,"%"))

assign(paste(pair[[2]]), plot)

}

library(ggpubr)
ggarrange(XL_PreInv,XL_Inv,XL_PostInv,
          XR_PreInv,XR_Inv,XR_PostInv, 
          nrow = 2, ncol = 3, 
          common.legend = TRUE)

ggarrange(XL, XR, 
          ncol = 2, nrow = 1, 
          common.legend = TRUE)

EntireX

ggarrange(Chr2.group4_MullerE,Chr2_MullerE,
          Chr3_MullerC,Chr4_MullerB,
          Chr5_MullerF,Unkown69,mtDNA,
          nrow=2, ncol = 4, 
          common.legend = TRUE)


## HeatMap from vcf

library(vcfR)

vcf_file <- "filtered_SNPs_GATK_PopGen_ChrX_haploid.vcf.gz"

vcf <- read.vcfR(vcf_file, verbose = FALSE)

# Extract the fixed information (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO)
fixed <- as.data.frame(getFIX(vcf))

# Extract the genotype information
genotype <- extract.gt(vcf)

# Combine fixed and genotype information into a single data frame
vcf_table <- cbind(fixed, 
                   genotype)

library(dplyr)

vcf_table <- vcf_table %>%
  mutate(across(8:85, ~ as.numeric(as.character(.))))

vcf_table[, 8:85] <- apply(vcf_table[, 8:85], 2, 
                           function(x) ifelse(!x %in% c(0, NA), 1, x))

vcf_table <- vcf_table[,c(1:2,4:85)]

num_samples <- ncol(vcf_table) - 6
threshold <- num_samples * 0.8

vcf_table <- vcf_table[rowSums(!is.na(vcf_table[, 7:ncol(vcf_table)])) >= threshold, ]

Daff_ST <- read.table("Daff_ST.txt")
ST_samples <- levels(as.factor(Daff_ST$V1))

Daff_SR1 <- read.table("Daff_SR1.txt")
SR1_samples <- levels(as.factor(Daff_SR1$V1))

Daff_SR2 <- read.table("Daff_SR2.txt")
SR2_samples <- levels(as.factor(Daff_SR2$V1))

ST <- vcf_table[, c(names(vcf_table)[1:6], ST_samples)]
SR1 <- vcf_table[, c(names(vcf_table)[1:6], SR1_samples)]
SR2 <- vcf_table[, c(names(vcf_table)[1:6], SR2_samples)]


library(pheatmap)

data_list <- list(ST = ST, SR1 = SR1, SR2 = SR2)

# Loop through the list of dataframes
for (df_name in names(data_list)) {
  df <- data_list[[df_name]]

# Extract position information for X-axis
positions <- df$POS

# Create a matrix for the heatmap
genotype_matrix <- as.matrix(df[, 7:ncol(df)])

# Transpose the matrix to switch X and Y axes
genotype_matrix <- t(genotype_matrix)

# Add column names (positions) and row names (individual IDs)
colnames(genotype_matrix) <- positions
rownames(genotype_matrix) <- colnames(df)[7:ncol(df)]

# Sort the columns of the genotype matrix based on positions
sorted_indices <- order(positions)
genotype_matrix <- genotype_matrix[, sorted_indices]

# Create a discrete color scale
discrete_colors <- c("lightblue", "yellow")

# Plot the heatmap
Plot <- pheatmap(genotype_matrix, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE, 
         color = discrete_colors,
         show_rownames = FALSE,
         show_colnames = FALSE,  # Hide column names (positions)
         na_col = "white",
         legend_breaks = c(0, 1),
         legend_labels = c("Consensus", "Variant"),
         na.rm=TRUE,
         legend = FALSE)

assign(paste(df_name,"Plot"), Plot)

}


grob_ST <- `ST Plot`$gtable
grob_SR1 <- `SR1 Plot`$gtable
grob_SR2 <- `SR2 Plot`$gtable

n_samples <- c(ST = 52, SR1 = 20, SR2 = 6)
total_samples <- sum(n_samples)
height_proportions <- n_samples / total_samples

library(gridExtra)

grid.arrange(grob_ST,
             grob_SR1,
             grob_SR2,
             nrow = 3,
             heights = height_proportions)

## PopGen stats - 100kb

library(readr)
PopGenStats <- read_csv("SNPs_GATK_PopGen_ChrX_diploid.Fst.Dxy.pi.csv")

# Set up the plotting layout
par(mfrow = c(1, 3))

# Plot each graph
plot(PopGenStats$start, 
     PopGenStats$pi_ST, 
     xlab="position", 
     ylab="diversity", 
     main="pi - ST - 100kb")

plot(PopGenStats$start, 
     PopGenStats$pi_SR1, 
     xlab="position", 
     ylab="diversity", 
     main="pi - SR1 - 100kb")

plot(PopGenStats$start, 
     PopGenStats$pi_SR2, 
     xlab="position", 
     ylab="diversity", 
     main="pi - SR2 - 100kb")

# Reset the layout
par(mfrow = c(1, 1))

## Fst and dxy

# Set up the plotting layout
par(mfrow = c(2, 3))

# Plot each graph
plot(PopGenStats$start, 
     PopGenStats$Fst_ST_SR1, 
     xlab="position", 
     ylab="Fst", 
     main="Fst - ST & SR1 - 100kb")

plot(PopGenStats$start, 
     PopGenStats$Fst_ST_SR2, 
     xlab="position", 
     ylab="Fst", 
     main="Fst - ST & SR2 - 100kb")

plot(PopGenStats$start, 
     PopGenStats$Fst_SR1_SR2, 
     xlab="position", 
     ylab="Fst", 
     main="Fst - SR1 & SR2 - 100kb")

plot(PopGenStats$start, 
     PopGenStats$dxy_ST_SR1, 
     xlab="position", 
     ylab="dxy", 
     main="dxy - ST & SR1 - 100kb")

plot(PopGenStats$start, 
     PopGenStats$dxy_ST_SR2, 
     xlab="position", 
     ylab="dxy", 
     main="dxy - ST & SR2 - 100kb")

plot(PopGenStats$start, 
     PopGenStats$dxy_SR1_SR2, 
     xlab="position", 
     ylab="dxy", 
     main="dxy - SR1 & SR2 - 100kb")


## PopGen stats 100kb window

# Read the data
pi.all <- read.table("pi_PopGen_100kb.windowed.pi", header=T)
Tajima <- read.table("TajimaDPopGen_100kb.Tajima.D", header=T)
Fst <- read.table("FstPopGenAll_100kb.windowed.weir.fst", header=T)
FstSTSR1 <- read.table("FstPopGen_STSR1_100kb.windowed.weir.fst", header=T)
FstSTSR2 <- read.table("FstPopGen_STSR2_100kb.windowed.weir.fst", header=T)
FstSR1SR2 <- read.table("FstPopGen_SR1SR2_100kb.windowed.weir.fst", header=T)

# Set up the plotting layout
par(mfrow = c(2, 3))

# Plot each graph
plot(pi.all$BIN_START, pi.all$PI, xlab="position", ylab="diversity", main="pi - all samples - 100kb")
plot(Tajima$BIN_START, Tajima$TajimaD, xlab="position", ylab="TajimaD", main="TajimaD - all samples - 100kb")
plot(Fst$BIN_START, Fst$MEAN_FST, xlab="position", ylab="Fst", main="Fst - all samples - 100kb")
plot(FstSTSR1$BIN_START, FstSTSR1$MEAN_FST, xlab="position", ylab="Fst", main="Fst - ST and SR1 - 100kb")
plot(FstSTSR2$BIN_START, FstSTSR2$MEAN_FST, xlab="position", ylab="Fst", main="Fst - ST and SR2 - 100kb")
plot(FstSR1SR2$BIN_START, FstSR1SR2$MEAN_FST, xlab="position", ylab="Fst", main="Fst - SR1 and SR2 - 100kb")

# Reset the layout
par(mfrow = c(1, 1))


## PopGen stats - 1Mb

library(readr)
PopGenStats <- read_csv("SNPs_GATK_PopGen_ChrX_diploid.Fst.Dxy.pi_1mb.csv")

# Set up the plotting layout
par(mfrow = c(1, 3))

# Plot each graph
plot(PopGenStats$start, 
     PopGenStats$pi_ST, 
     xlab="position", 
     ylab="diversity", 
     main="pi - ST - 1Mb")

plot(PopGenStats$start, 
     PopGenStats$pi_SR1, 
     xlab="position", 
     ylab="diversity", 
     main="pi - SR1 - 1Mb")

plot(PopGenStats$start, 
     PopGenStats$pi_SR2, 
     xlab="position", 
     ylab="diversity", 
     main="pi - SR2 - 1Mb")

# Reset the layout
par(mfrow = c(1, 1))

## Fst and dxy

# Set up the plotting layout
par(mfrow = c(2, 3))

# Plot each graph
plot(PopGenStats$start, 
     PopGenStats$Fst_ST_SR1, 
     xlab="position", 
     ylab="Fst", 
     main="Fst - ST & SR1 - 1Mb")

plot(PopGenStats$start, 
     PopGenStats$Fst_ST_SR2, 
     xlab="position", 
     ylab="Fst", 
     main="Fst - ST & SR2 - 1Mb")

plot(PopGenStats$start, 
     PopGenStats$Fst_SR1_SR2, 
     xlab="position", 
     ylab="Fst", 
     main="Fst - SR1 & SR2 - 1Mb")

plot(PopGenStats$start, 
     PopGenStats$dxy_ST_SR1, 
     xlab="position", 
     ylab="dxy", 
     main="dxy - ST & SR1 - 1Mb")

plot(PopGenStats$start, 
     PopGenStats$dxy_ST_SR2, 
     xlab="position", 
     ylab="dxy", 
     main="dxy - ST & SR2 - 1Mb")

plot(PopGenStats$start, 
     PopGenStats$dxy_SR1_SR2, 
     xlab="position", 
     ylab="dxy", 
     main="dxy - SR1 & SR2 - 1Mb")



## PopGen stats 1Mb window

# Read the data
pi.all <- read.table("pi_PopGen_1Mb.windowed.pi", header=T)
Tajima <- read.table("TajimaDPopGen_1Mb.Tajima.D", header=T)
Fst <- read.table("FstPopGenAll_1Mb.windowed.weir.fst", header=T)
FstSTSR1 <- read.table("FstPopGen_STSR1_1Mb.windowed.weir.fst", header=T)
FstSTSR2 <- read.table("FstPopGen_STSR2_1Mb.windowed.weir.fst", header=T)
FstSR1SR2 <- read.table("FstPopGen_SR1SR2_1Mb.windowed.weir.fst", header=T)

# Set up the plotting layout
par(mfrow = c(2, 3))

# Plot each graph
plot(pi.all$BIN_START, pi.all$PI, xlab="position", ylab="diversity", main="pi - all samples - 1Mb")
plot(Tajima$BIN_START, Tajima$TajimaD, xlab="position", ylab="TajimaD", main="TajimaD - all samples - 1Mb")
plot(Fst$BIN_START, Fst$MEAN_FST, xlab="position", ylab="Fst", main="Fst - all samples - 1Mb")
plot(FstSTSR1$BIN_START, FstSTSR1$MEAN_FST, xlab="position", ylab="Fst", main="Fst - ST and SR1 - 1Mb")
plot(FstSTSR2$BIN_START, FstSTSR2$MEAN_FST, xlab="position", ylab="Fst", main="Fst - ST and SR2 - 1Mb")
plot(FstSR1SR2$BIN_START, FstSR1SR2$MEAN_FST, xlab="position", ylab="Fst", main="Fst - SR1 and SR2 - 1Mb")

# Reset the layout
par(mfrow = c(1, 1))



## SNPeff output

library(vcfR)

vcf_file <- "Snpeffoutput.vcf"

vcf <- read.vcfR(vcf_file, verbose = FALSE)

# Extract the fixed information (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO)
fixed <- as.data.frame(getFIX(vcf))

# Extract the genotype information
genotype <- extract.gt(vcf)

library(tidyr)
# Extract the INFO field
info <- as.data.frame(vcf@fix[, "INFO"], stringsAsFactors = FALSE)
colnames(info) <- "INFO"

# Extract the ANN field from the INFO column
info$ANN <- sapply(strsplit(info$INFO, ";"), function(x) {
  ann_field <- grep("^ANN=", x, value = TRUE)
  if (length(ann_field) > 0) {
    return(sub("^ANN=", "", ann_field))
  } else {
    return(NA)
  }
})

# Remove the original INFO column to avoid redundancy
info$INFO <- NULL

# Separate the ANN field into its components
info <- info %>%
  separate(ANN, into = c("Effect", "Impact", "Gene", "Gene_ID", "Feature_Type", "Feature_ID", "Transcript_BioType", "Rank", "HGVS.c", "HGVS.p", "cDNA.pos", "CDS.pos", "Protein.pos", "Distance", "ERRORS"), sep = "\\|", fill = "right")


# Combine fixed and genotype information into a single data frame
vcf_table <- cbind(fixed, info, genotype)

Data <- vcf_table[,c(1:2,4:7,9,13,15,23:100)]

##Now lets analyze Data
levels(as.factor(Data$Impact))

Coding <- subset(Data,
                 Rank=="protein_coding")

levels(as.factor(Coding$Impact))

library(dplyr)

Coding <- Coding %>%
  mutate(across(10:87, ~ as.numeric(as.character(.))))

Coding[, 10:87] <- apply(Coding[, 10:87], 2, 
                         function(x) ifelse(!x %in% c(0, NA), 1, x))

DeleterCoding <- subset(Coding, Impact %in% c("initiator_codon_variant",                          
                                              "missense_variant",
                                              "start_lost",            
                                              "stop_gained",
                                              "stop_lost",
                                              "synonymous_variant"))


library(dplyr)

# Select columns 10 to 87
columns_of_interest <- DeleterCoding %>% select(10:87)

# Calculate the proportion of "1" in each column
proportion_of_ones <- columns_of_interest %>%
  summarise(across(everything(), ~ mean(. == 1, na.rm = TRUE)))


library(reshape2)

DelCod <- melt(proportion_of_ones)

library(readr)

Id <- read_csv("IDaccX.csv")

colnames(DelCod)[1] <- "ind"

Data2 <- merge(DelCod, Id, by="ind", all.x = TRUE)

library(ggplot2)

ggplot(Data2, aes(x=X,
                  y=value,
                  col=X)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(),
             alpha=0.5) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x="",
       y="Proportion of deleterious variants in translating region")

model <- aov(value ~ X, data = Data2)
summary(model)
