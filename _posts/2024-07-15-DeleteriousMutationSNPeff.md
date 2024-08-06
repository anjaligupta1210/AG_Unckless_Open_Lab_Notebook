---
layout: post
title: Deleterious mutations on the X
date: 15 July 2024  
category: [ Computational Pipelines ]
tags: [ EnsemblVEP, SNPeff, GenomeDelta ]
---

## Estimating whether SR1 and SR2 are accumulating deleterious mutations

### 1. Using Ensembl-VEP

### Getting my files ready


```python
mkdir -p /work/unckless/a948g501/Evol2024/EnsemblVEP

cd /work/unckless/a948g501/Evol2024/EnsemblVEP

## scp Helixer ST genome gff to this folder

cp /work/unckless/a948g501/PAML/Daffinis_STfemale_v5.1.masked.fasta .

module load samtools
samtools faidx Daffinis_STfemale_v5.1.masked.fasta


module load bcftools

bcftools filter -i 'TYPE="snp" && FORMAT/DP > 3 && FORMAT/DP < 500' /work/unckless/a948g501/Evol2024/PopGenPCA/GATK_PopGen_ChrX_haploid.vcf > "filtered_SNPs_GATK_PopGen_ChrX_haploid.vcf"

module load samtools
bgzip filtered_SNPs_GATK_PopGen_ChrX_haploid.vcf
tabix filtered_SNPs_GATK_PopGen_ChrX_haploid.vcf.gz

bgzip Daffinis_STfemale_v5.1.masked.fasta

samtools faidx Daffinis_STfemale_v5.1.masked.fasta.gz

module load conda
conda activate gffread_env

# Convert GFF to GTF
gffread STgenome_helixer.gff -T -o STgenome_helixer.gtf

# Sort the GTF file
sort -k1,1 -k4,4n STgenome_helixer.gtf > STgenome_helixer.sorted.gtf

# Compress the GTF file
bgzip STgenome_helixer.sorted.gtf

# Index the GTF file
tabix -p gff STgenome_helixer.sorted.gtf.gz

conda deactivate

```


```python
#!/bin/bash
#SBATCH --job-name=EnsemblVEP  # Job name
#SBATCH --partition=sixhour      # Partition Name (Required)
#SBATCH --mail-type=END,FAIL     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=anjaligupta@ku.edu   # Where to send mail
#SBATCH --chdir=/work/unckless/a948g501/Evol2024/EnsemblVEP/
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1          # Run on a single CPU
#SBATCH --mem=64gb           # Job memory request
#SBATCH --time=06:00:00        # Time limit days-hrs:min:sec
#SBATCH --output=EnsemblVEP_%j.log  # Standard output and error log


module load conda
conda activate vep


vep -i filtered_SNPs_GATK_PopGen_ChrX_haploid.vcf.gz \
-o ChrX_EnsemblVEP.txt \
--format vcf \
--fasta Daffinis_STfemale_v5.1.masked.fasta.gz \
--gtf STgenome_helixer.sorted.gtf.gz

conda deactivate
```

### 2. Using GenomeDelta to annotate TEs


```python
#!/bin/bash
#SBATCH --job-name=GenomeDelta  # Job name
#SBATCH --partition=unckless      # Partition Name (Required)
#SBATCH --mail-type=END,FAIL     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=anjaligupta@ku.edu   # Where to send mail
#SBATCH --chdir=/work/unckless/a948g501/Evol2024/GenomeDelta
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1          # Run on a single CPU
#SBATCH --mem=64gb           # Job memory request
#SBATCH --time=60-00:00:00        # Time limit days-hrs:min:sec
#SBATCH --output=GenomeDelta_%j.log  # Standard output and error log

module load bedtools
module load blast+


for bam_file in pool10.picard.bam
do
    base_name=$(basename "$bam_file" .picard.bam)
    bash /work/unckless/a948g501/Evol2024/GenomeDelta/GenomeDelta/linux/main.sh \
        --bam "$bam_file" \
        --fa Daffinis_STfemale_v5.1.masked.fasta \
        --of /work/unckless/a948g501/Evol2024/GenomeDelta/"$base_name" \
        --prefix "$base_name" \
        --t 20
done

```

### 3. Premature stop codons in protein coding genes


```python
bedtools intersect -a $bam_file".bam" -b $gtf_file -wa | \
    samtools fasta -N - > $bam_file".fa"

transeq -sequence $bam_file".fa" -outseq "CDS"$bam_file".fa"

##look for stop at stop codon in transeq and then use bioawk
```

### 4. Using SNPeff


```python
module load java
module load openjdk
source ~/.bashrc

java -jar /work/unckless/a948g501/Evol2024/SNPeff/snpEff/snpEff.jar databases
```


```python
module load conda
conda activate gffread_env
gffread --gtf /work/unckless/a948g501/PAML/Daffinis_STfemale_v5.1.ago2.gtf -L > Daff_genes.gff
gffread -w Daff_CDS.fa -g /work/unckless/a948g501/PAML/Daffinis_STfemale_v5.1.masked.fasta /work/unckless/a948g501/PAML/Daffinis_STfemale_v5.1.ago2.gtf
conda deactivate
```


```python
DBNAME="Drosophila_affinis"
GFF3="/work/unckless/a948g501/Evol2024/SNPeff/AnnotationInputFiles/STgenome_helixer.gff"
FASTA="/work/unckless/a948g501/Evol2024/SNPeff/AnnotationInputFiles/Daffinis_STfemale_v5.1.masked.fasta"


#Go into the snpEff directory and create a directory for your files
cd /work/unckless/a948g501/Evol2024/SNPeff/snpEff
mkdir -p data/$DBNAME

#Copy the files into snpEff's directory structure
cp $GFF3 data/$DBNAME/genes.gff
cp $FASTA data/$DBNAME/sequences.fa

#Edit snpEff.config and insert your specific database information:
echo "# $DBNAME
$DBNAME.genome : $DBNAME
$DBNAME.reference : data/$DBNAME/sequences.fa" >> snpEff.config

#Build the database
module load java
module load openjdk
source ~/.bashrc

java -jar /work/unckless/a948g501/Evol2024/SNPeff/snpEff/snpEff.jar build -gff3 -v $DBNAME
```


```python
module load java
module load openjdk
source ~/.bashrc

java -Xmx4g \
-jar /work/unckless/a948g501/Evol2024/SNPeff/snpEff/snpEff.jar \
ann \
-configOption 'Drosophila_affinis'.genome='Drosophila_affinis' -configOption 'Drosophila_affinis'.codonTable='Standard' 'Drosophila_affinis' \
/work/unckless/a948g501/Evol2024/PopGenPCA/filtered_SNPs_GATK_PopGen_ChrX_haploid.vcf.gz \
> SNPeff_ChrX_SNP_annotations.vcf
```
