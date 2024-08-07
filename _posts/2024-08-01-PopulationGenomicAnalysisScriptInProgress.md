---
layout: post
title: Population Genomics Analysis - in progress
date: 01 August 2024    
category: [ Computational Pipelines ]
tags: [ PCangsd, GATK, Population genomics ]
---

# Population Genomics Analysis for *D. affinis* SR1, SR2, and ST

I have Illumina short reads from ~80 male *D. affinis* samples.

First I will concatenate and trim the reads and align them to my reference genome


```python
## Load required modules
ml bwa

## Index reference genome fasta file
bwa index Daffinis_STfemale_v5.1.masked.fasta
```

Concatenate reads --> Run fastp --> Align them to reference genome & index it


```python
#!/bin/bash

# Specify the directory containing your local files
File_dir="/work/unckless/a948g501/PopGen/"

for i in "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 74 75 76 77 78 79 80 87 88 89 90 91"
do

    # Copy files from remote server to local directory
    scp a948g501@dtn.ku.edu:/resfs/GROUPS/MB/NIH0072614/Daffinis/RawReads/Popgen/pool_${i}_*.fq.gz ${File_dir} || { echo "Error: Failed to copy files for pool $i"; exit 1; }

    chmod 777 "${File_dir}pool_${i}_"*.fq.gz

    # Concatenate *1.fq.gz files
    cat "${File_dir}pool_${i}_"*1.fq.gz > "${File_dir}pool_${i}_1.fq.gz"

    # Concatenate *2.fq.gz files
    cat "${File_dir}pool_${i}_"*2.fq.gz > "${File_dir}pool_${i}_2.fq.gz"

    # Run fastp
    module load conda
    conda activate fastp
    fastp -f 5 -l 5 -i "${File_dir}pool_${i}_1.fq.gz" -I "${File_dir}pool_${i}_2.fq.gz" -o "${File_dir}fil.pool${i}_1.fq.gz" -O "${File_dir}fil.pool${i}_2.fq.gz"
    conda deactivate

    # Remove concatenated reads - only keep filtered reads
    rm "${File_dir}pool_${i}_"*.fq.gz 

    # Map to genome
    module load bwa
    module load samtools
    bwa mem Daffinis_STfemale_v5.1.masked.fasta "${File_dir}fil.pool${i}_1.fq.gz" "${File_dir}fil.pool${i}_2.fq.gz" | samtools view -hb -F 4 - | samtools sort - > "${File_dir}pool${i}.bam"
    samtools index "${File_dir}pool${i}.bam" 

    # Remove reads
    rm "${File_dir}fil.pool${i}_"*.fq.gz

done
```

### Check coverage of samples and get rid of low quality samples

This is a quality control step and you can skip this step if you have high coverage for all your samples


```python
samtools coverage <bam-file-name>
```

If you have a few very high coverage samples and all other samples are low coverage, downsize coverage for the few high coverage samples or 


if you have mostly high coverage samples and a few low coverage samples, get rid of your low coverage samples of possible

### Variant Calling Using GATK

This is a two step process -

1. First I'll AddOrReplaceReadGroups using Picard in the alignments
2. Next I'll call Variants separately for Chromosome X (haploid) and Autosomes (diploid)

### (i) Picard AddOrReplaceReadGroups


```python
## Load required modules
ml picard

## List my raw reads files
SamplesArray="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 74 75 76 77 78 79 80 87 88 89 90 91"

## Run Picard
for file in "pool${SamplesArray}"
do
picard AddOrReplaceReadGroups I=${file}.bam \
    O=${file}.picard.bam \
        RGLB=${file} \
            RGID=${file} \
                RGPL=illumina \
                    RGPU=unit1 \
                        RGSM=${file}
done

#I=String                      Input file (BAM or SAM or a GA4GH url).  Required.
#O=File                        Output file (BAM or SAM).  Required. 
#RGID=String                   Read-Group ID  Default value: 1. This option can be set to 'null' to clear the default value. 
#RGLB=String                   Read-Group library  Required. 
#RGPL=String                   Read-Group platform (e.g. ILLUMINA, SOLID)  Required.
#RGPU=String                   Read-Group platform unit (eg. run barcode)  Required. 
#RGSM=String                   Read-Group sample name  Required. 


## Index Picard alignments
for file in "pool${SamplesArray}"
do
samtools index ${file}.picard.bam
done
```

### (ii) Calling Variants using GATK

GATK requires indexing reference fasta file using a .dict extension using gatk CreateSequenceDictionary


```python
gatk CreateSequenceDictionary -R Daffinis_STfemale_v5.1.masked.fasta
```

I’m going to use gatk HaplotypeCaller to call SNPs separately for ChrX (haploid) and autosomes (diploid)

quick things:

- -L is to specify only region to include
- -XL is to specify only region to exclude

GATK is slow and will take a long time to run


```python
## Make a list of picard bam files

ls *.picard.bam > DaffPopGenBams.list
```

The list file should have a ".list" extension for GATK to recognise it


```python
## Load required modules
ml gatk

## Run HaplotypeCaller

## ChrX
gatk HaplotypeCaller \
     -I DaffPopGenBams.list \
     -O GATK_PopGen_ChrX_haploid.vcf \
     -R Daffinis_STfemale_v5.1.masked.fasta \
     --add-output-vcf-command-line true \
     -ploidy 1 \
     -L ChrX_MullerAD \
     --native-pair-hmm-threads 16 

## Autosomes
module load gatk
gatk HaplotypeCaller \
     -I DaffPopGenBams.list \
     -O GATK_PopGen_Autosomes_diploid.vcf \
     -R Daffinis_STfemale_v5.1.masked.fasta \
     --add-output-vcf-command-line true \
     -ploidy 2 \
     -XL ChrX_MullerAD \
     --native-pair-hmm-threads 16 
```

### Filtering SNPs from GATK output --> Convert to Beagle --> PCA

#### Chromosome X


```python
pre=GATK_PopGen_ChrX_haploid
VCF=GATK_PopGen_ChrX_haploid.vcf

module load bcftools

bcftools filter -i 'TYPE="snp" && FORMAT/DP > 3 && FORMAT/DP < 500' $VCF > "filtered_SNPs_"$VCF


module load vcftools

vcftools --vcf "filtered_SNPs_"$VCF --chr ChrX_MullerAD --BEAGLE-PL --out "filtered_SNPs_"$pre


gzip "filtered_SNPs_"$pre".BEAGLE.PL"


module load conda

conda activate pcangsd

pcangsd -b "filtered_SNPs_"$pre".BEAGLE.PL.gz" -e 2 -t 64 -o "pcangsd_filtered_SNPs_"$pre

conda deactivate
```

#### Autosomes


```python
pre=GATK_PopGen_Autosomes_diploid
VCF=GATK_PopGen_Autosomes_diploid.vcf

module load bcftools

bcftools filter -i 'TYPE="snp" && FORMAT/DP > 3 && FORMAT/DP < 500' $VCF > "filtered_SNPs_"$VCF


module load vcftools

module load conda

conda activate pcangsd

for chrom in Chr4_MullerB Chr2_MullerE Chr2.group4_MullerE Chr3_MullerC Chr5_MullerF Unknown_69 mtDNA
do
vcftools --vcf "filtered_SNPs_"$VCF --chr $chrom --BEAGLE-PL --out "filtered_SNPs_"$pre"_"$chrom
gzip "filtered_SNPs_"$pre"_"$chrom".BEAGLE.PL"
pcangsd -b "filtered_SNPs_"$pre"_"$chrom".BEAGLE.PL.gz" -e 2 -t 64 -o "pcangsd_filtered_SNPs_"$pre"_"$chrom
done


conda deactivate
```

### Analysis using R

PCangsd will output .cov files --> Analyze it in R using [MyRscript](https://github.com/anjaligupta1210/AG_Unckless_Open_Lab_Notebook/blob/master/rscripts/Evol2024DataAndCode/ALLplots.R)



Basic R code for PC analysis using .cov file -

```R
C <- as.matrix(read.table("pcangsd.cov")) # Reads estimated covariance matrix

# Plot PCA plot
e <- eigen(C)
plot(e$vectors[,1:2], xlab="PC1", ylab="PC2", main="PCAngsd")

```

### To verify PCA you could make a heatmap from your vcf file in R using the code below -

```R
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

## I'm using the cluster IDs from ST, SR1 and SR2 as I got it from PCanmgsd output PCA clusters

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
```

### Do a PCA in R on your vcf

This was really slow since it was a large vcf file

I used the code from [Part1](https://rpubs.com/madisondougherty/980777) & [Part2](https://rpubs.com/torris459/981399)

### Some further analysis that you might want to do are - 

#### 1. Selection scans

Use scripts from - https://speciationgenomics.github.io/haplotypes/

#### 2. Coalescent analysis

Use scripts from - https://speciationgenomics.github.io/easysfs/ & https://speciationgenomics.github.io/fastsimcoal2/

#### 3. Retracing ancestral recombination graphs

Try RELEARN & SINGER -- yet to try!!
