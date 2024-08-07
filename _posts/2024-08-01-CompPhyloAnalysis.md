---
layout: post
title: Comparative Phylogenetic Anlaysis Full Script
date: 01 August 2024   
category: [ Computational Pipelines ]
tags: [ GATK, Sliding window trees, SNP calling ]
---

# Comparative phylogenetic analysis for *D. affinis* SR1, SR2, and ST

First I will navigate to my working directory and then copy all my raw Illumina short reads to this folder


```python
## Specify your working directory

dir=/work/unckless/a948g501/CompPhylo
```

## 1. Run fastqc and fastp on your raw reads

This is the quality control step for your raw reads,

First I'll run fastqc - if you want to combine multiple fastqc file reports into a single file for ease of viewing - you can use multiQC - https://multiqc.info/docs/

Based on your fastqc reports, you can decide whether you need to trim the reads!


```python
# Navigate to my working directory
cd $dir

# List my raw reads files
SamplesRaw="Daff_ST_1.fastq Daff_ST_2.fastq \
    Daff_SR1_1.fastq Daff_SR1_2.fastq \
        Daff_SR2_1.fastq Daff_SR2_2.fastq \
            Dalg_1.fastq Dalg_2.fastq \
                Datha_ea_1.fastq Datha_ea_2.fastq \
                    Datha_eb_1.fastq Datha_eb_2.fastq \
                        Dpse_1.fastq Dpse_2.fastq \
                            Dazt_1.fastq Dazt_2.fastq \
                                Dnrg_1.fastq Dnrg_2.fastq"

# Install fastqc and fastp if you do not have it installed on the cluster

## You can use multiqc to compile fastqc reports into a single file if you have a large number of samples

## load required modules
ml fastqc
ml conda
conda activate fastp

## Running fastqc on all my raw reads
for file in SamplesRaw
do
fastqc $file
done

### Look at your fastqc reports and see if you need to trim any reads using fastp

### I had to trim the first and last 5 bases

## Trimming reads using fastp
for file in SamplesRaw
do
fastp -i $file \
    -f 5 -t 5 \
        -o $file".gz"
done

## deactivate the conda environment
conda deactivate
```

## 2. Align to reference genome using bwa mem

First copy the reference genome to this folder and index it.

I'm using *D. affinis* masked female ST genome 

I will align the reads to the reference genome and index the alignment file


```python
## Load required modules
ml bwa
ml samtools

# List my raw reads files
Samples="Daff_ST \
    Daff_SR1 \
        Daff_SR2 \
            Dalg \
                Datha_ea \
                    Datha_eb \
                        Dpse \
                            Dazt \
                                Dnrg"


## Align using bwa
for file in Samples
do
bwa mem Daffinis_STfemale_v5.1.masked.fasta $file"_1.fastq.gz" $file"_2.fastq.gz" \
    | samtools view -hb -F 4 - | \
        samtools sort - > $file".bam"
done


# Indexing the aligned file
for file in Samples
do
samtools index $file".bam"
done

```

## 3. Calling variants using GATK

This is a two step process -


    1. First I'll AddOrReplaceReadGroups using Picard in the alignments
    2. Next I'll call Variants separately for Chromosome X (haploid) and Autosomes (diploid) 

### (i) Picard AddOrReplaceReadGroups


```python
## Load required modules
ml picard

## List my raw reads files
Samples="Daff_ST \
    Daff_SR1 \
        Daff_SR2 \
            Dalg \
                Datha_ea \
                    Datha_eb \
                        Dpse \
                            Dazt \
                                Dnrg"

## Run Picard
for file in Samples
do
picard AddOrReplaceReadGroups I=${file}.bam \
    O=${file}.picard.bam \
        RGLB=${file} \
            RGID=${file} \
                RGPL=illumina \
                    RGPU=unit1 \
                        RGSM=${file}
done

## Index Picard alignments
for file in Samples
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

## load required modules
ml gatk

## Run GATK HaplotypeCaller
### ChrX
gatk HaplotypeCaller \
     -I Daff_SR1.picard.bam \
     -I Daff_SR2.picard.bam \
     -I Daff_ST.picard.bam \
     -I Dalg.picard.bam \
     -I Datha_ea.picard.bam \
     -I Datha_eb.picard.bam \
     -I Dazt.picard.bam \
     -I Dpse.picard.bam \
     -I Dnrg.picard.bam \
     -O GATK_ChrX_haploid.vcf \
     -R Daffinis_STfemale_v5.1.masked.fasta \
     --add-output-vcf-command-line true \
     -ploidy 1 \
     -L ChrX_MullerAD

### Autosomes
gatk HaplotypeCaller \
     -I Daff_SR1.picard.bam \
     -I Daff_SR2.picard.bam \
     -I Daff_ST.picard.bam \
     -I Dalg.picard.bam \
     -I Datha_ea.picard.bam \
     -I Datha_eb.picard.bam \
     -I Dazt.picard.bam \
     -I Dpse.picard.bam \
     -I Dnrg.picard.bam \
     -O GATK_Autosomes_diploid.vcf \
     -R Daffinis_STfemale_v5.1.masked.fasta \
     --add-output-vcf-command-line true \
     -ploidy 2 \
     -XL ChrX_MullerAD

## 4. Filtering SNPs from the GATK output using bcftools

I'm going to use only SNPs with quality > 40 for downstream analysis


```python
## Load required modules
ml bcftools
ml samtools

## Filter using bcftools

### ChrX
bcftools filter \
    -i 'QUAL > 40 && TYPE="snp"' \
    GATK_ChrX_haploid.vcf \
    -o GATK_filtered_ChrX_snp_haploid.vcf

### Autosomes
bcftools filter \
    -i 'QUAL > 40 && TYPE="snp"' \
    GATK_Autosomes_diploid.vcf \
    -o GATK_filtered_Autosome_snp_diploid.vcf

### bgzip and index vcf file
bgzip GATK_filtered_ChrX_snp_haploid.vcf
tabix GATK_filtered_ChrX_snp_haploid.vcf.gz

bgzip GATK_filtered_Autosome_snp_diploid.vcf
tabix GATK_filtered_Autosome_snp_diploid.vcf.gz
```

## 5. Making trees for phylogenetic analysis

Making trees will be done in a 3 step process -
1. Convert vcf to phylip using [vcf2phylip.py](https://raw.githubusercontent.com/edgardomortiz/vcf2phylip/master/vcf2phylip.py)
2. Make trees using iqtree
3. Make trees using RaxML

I'll make trees for -
1. All autosomes
2. Entire X
3. XL and XR arm 
4. XL and XR arm: pre-inversion, inversion, and post-inversion
5. Sliding windows of 3Mb

### (i) Partitioning vcf files using bcftools

#### A. Partitioning XL and XR arm


```python
## Load required modules
ml bcftools

## Partition using bcftools
bcftools view \
    -r ChrX_MullerAD:1-30000000 \
    GATK_filtered_ChrX_snp_haploid.vcf.gz \
    -o chrX_XL.vcf

bcftools view \
    -r ChrX_MullerAD:30000000- \
    GATK_filtered_ChrX_snp_haploid.vcf.gz \
    -o chrX_XR.vcf
```

#### B. Partitioning XL and XR arm: pre-inversion, inversion, and post-inversion


```python
## Load required modules
ml bcftools

## Partition using bcftools
bcftools view \
    -r ChrX_MullerAD:1-5490855 \
    GATK_filtered_ChrX_snp_haploid.vcf.gz \
    -o chrX_XL_PreInv.vcf

bcftools view \
    -r ChrX_MullerAD:5490855-10233152 \
    GATK_filtered_ChrX_snp_haploid.vcf.gz \
    -o chrX_XL_Inv.vcf

bcftools view \
    -r ChrX_MullerAD:10233152-30000000 \
    GATK_filtered_ChrX_snp_haploid.vcf.gz \
    -o chrX_XL_PostInv.vcf

bcftools view \
    -r ChrX_MullerAD:30000000-44023827 \
    GATK_filtered_ChrX_snp_haploid.vcf.gz \
    -o chrX_XR_PreInv.vcf

bcftools view \
    -r ChrX_MullerAD:44023827-62734159 \
    GATK_filtered_ChrX_snp_haploid.vcf.gz \
    -o chrX_XR_Inv.vcf

bcftools view \
    -r ChrX_MullerAD:62734159- \
    GATK_filtered_ChrX_snp_haploid.vcf.gz \
    -o chrX_XR_PostInv.vcf
```

#### C. Partitioning for sliding windows of 3Mb

First I’m going to make a bed file with information about each partition


```python
j=1
for i in {1..23}
do
    echo -e "ChrX_MullerAD\t${j}\t$((j+2999999))" >> partitions_X_3mb.bed
    j=$((j+3000000))
done
```

Next I’ll run bcftools view looping over each line in my partitions file


```python
# Define the input BED file
bed_file="partitions_X_3mb.bed"

# Define the input VCF file
vcf_file="GATK_filtered_ChrX_snp_haploid.vcf.gz"

## Load required modules
ml bcftools

# Loop through each line of the BED file
while IFS=$'\t' read -r chromosome start end; do
    # Define the output file prefix
    output_prefix="${chromosome}_${start}_${end}"

    # Define the output VCF file name
    output_vcf="${output_prefix}.vcf"

    # Subset the VCF file using bcftools view
    bcftools view -r "${chromosome}:${start}-${end}" "$vcf_file" -o "$output_vcf"

    echo "Subset VCF file generated: $output_vcf"
done < "$bed_file"
```

### (ii) Convert VCFs to phylip

#### Autosomes


```python
## Load required modules
ml python

## Convert vcf to phylip
python vcf2phylip.py -i GATK_filtered_Autosome_snp_diploid.vcf.gz
```

#### Entire X


```python
## Load required modules
ml python

## Convert vcf to phylip
python vcf2phylip.py -i GATK_filtered_ChrX_snp_haploid.vcf.gz
```

#### Partitioned XL and XR


```python
## Load required modules
ml python

## Convert vcf to phylip
for file in chrX_*.vcf
do
python vcf2phylip.py -i $file
done
```

#### Sliding windows of 3Mb


```python
bed_file="partitions_X_3mb.bed"

## Load required modules
ml python

# Loop through each VCF file matching the pattern
while IFS=$'\t' read -r chromosome start end; do
    # Define the output file prefix
    output_prefix="${chromosome}_${start}_${end}"

    for vcf_file in "${output_prefix}.vcf"; do
    python vcf2phylip.py -i "$vcf_file"
    done
done < "$bed_file"
```

### (iii) Make trees using iqtree2 and RaxML

#### Trees using iqtree2


```python
## Load required modules
ml iqtree2

## Make trees using ML and bootstrap replicates
for file in *.min4.phy
do
iqtree2 -s $file \
    -bb 1000 \
        -wbt \
            -nt AUTO \
                -o Dpse
done
```

#### Trees using RaxML

RaxML will take a lot longer to run compared to iqtree.



iqtree will take a few seconds to a minute to run while RaxML will end up taking a few to several hours


```python
## Load required modules
ml raxml-ng

## Make trees using ML and bootstrap replicates
for file in *.min4.phy
do
raxml-ng --all \
    --msa $file \
    --msa-format PHYLIP \
    --outgroup Dpse \
    --model TVMef \
    --threads 8 \
    --bs-metric fbp,tbe \
    --force
done
```





You can visulaise newick trees on the terminal using the following command from newick_utils


```python
ml conda
conda activate newick_utils
nw_display <Newick-tree-file>
```
