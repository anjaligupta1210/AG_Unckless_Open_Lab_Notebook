---
layout: post
title: PCA using pcangsd and VC using GATK and Popgen stats
date: 31 May 2024  
category: [ Computational Pipelines ]
tags: [ PCangsd, GATK, Popgen stats on vcf, Selection scan, Population genomics ]
---

## PCA following Variant calling using GATK 

Run an interactive session


```python
srun --time=6:00:00 --ntasks=1 --nodes=1 --partition=sixhour --pty /bin/bash -l
```


```python
mkdir /work/unckless/a948g501/Evol2024/PopGenPCA
```

Copy *D. affinis* ST female genome reference assembly to this folder


```python
module load GATK
gatk CreateSequenceDictionary -R Daffinis_STfemale_v5.1.masked.fasta
```

Make affinis sample ID list - ``affinisPoolList.txt``


```python
1
2
3
4
5
6
7
8
9
10
11
12
13
14
15
16
17
18
19
20
21
22
23
24
26
27
28
29
30
31
32
33
34
35
36
37
38
39
40
41
42
43
44
45
46
47
48
57
58
59
60
61
62
63
64
65
66
67
68
69
70
71
72
74
75
76
77
78
79
80
87
88
89
90
91
```

Make an array list


```python
array=$(paste -sd, affinisPoolList.txt)
```

### Picard AddOrReplaceReadGroups


```python
#!/bin/bash
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=06:00:00
#SBATCH --mem=10GB
#SBATCH --partition=sixhour
#SBATCH --chdir=/work/unckless/a948g501/Evol2024/PopGenPCA
#SBATCH --array=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,74,75,76,77,78,79,80,87,88,89,90,91


module load picard

picard AddOrReplaceReadGroups \
 I=/work/unckless/a948g501/PopGen/pool$SLURM_ARRAY_TASK_ID.bam \
 O=pool$SLURM_ARRAY_TASK_ID.picard.bam \
 RGLB=pool$SLURM_ARRAY_TASK_ID \
 RGID=pool$SLURM_ARRAY_TASK_ID \
 RGPL=illumina \
 RGPU=unit1 \
 RGSM=pool$SLURM_ARRAY_TASK_ID

#I=String                      Input file (BAM or SAM or a GA4GH url).  Required.
#O=File                        Output file (BAM or SAM).  Required. 
#RGID=String                   Read-Group ID  Default value: 1. This option can be set to 'null' to clear the default value. 
#RGLB=String                   Read-Group library  Required. 
#RGPL=String                   Read-Group platform (e.g. ILLUMINA, SOLID)  Required.
#RGPU=String                   Read-Group platform unit (eg. run barcode)  Required. 
#RGSM=String                   Read-Group sample name  Required. 
```

#### Index Picard alignment files


```python
#!/bin/bash
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=06:00:00
#SBATCH --mem=10GB
#SBATCH --partition=sixhour
#SBATCH --chdir=/work/unckless/a948g501/Evol2024/PopGenPCA
#SBATCH --array=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,74,75,76,77,78,79,80,87,88,89,90,91


module load samtools 

samtools index pool$SLURM_ARRAY_TASK_ID.picard.bam
```

### Copy Daff ST, SR1, and SR2 files to this folder


```python
cp /work/unckless/a948g501/SlidingTrees/Daff_*.picard.* /work/unckless/a948g501/PopGenGATK/ &
```

### List of all alignment files to use for Variant calling


```python
Daff_ST.picard.bam
Daff_SR1.picard.bam
Daff_SR2.picard.bam
pool1.picard.bam
pool2.picard.bam
pool3.picard.bam
pool4.picard.bam
pool5.picard.bam
pool6.picard.bam
pool7.picard.bam
pool8.picard.bam
pool9.picard.bam
pool10.picard.bam
pool11.picard.bam
pool12.picard.bam
pool13.picard.bam
pool14.picard.bam
pool15.picard.bam
pool16.picard.bam
pool17.picard.bam
pool18.picard.bam
pool19.picard.bam
pool20.picard.bam
pool21.picard.bam
pool22.picard.bam
pool23.picard.bam
pool24.picard.bam
pool26.picard.bam
pool27.picard.bam
pool28.picard.bam
pool29.picard.bam
pool30.picard.bam
pool31.picard.bam
pool32.picard.bam
pool33.picard.bam
pool34.picard.bam
pool35.picard.bam
pool36.picard.bam
pool37.picard.bam
pool38.picard.bam
pool39.picard.bam
pool40.picard.bam
pool41.picard.bam
pool42.picard.bam
pool43.picard.bam
pool44.picard.bam
pool45.picard.bam
pool46.picard.bam
pool47.picard.bam
pool48.picard.bam
pool57.picard.bam
pool58.picard.bam
pool59.picard.bam
pool60.picard.bam
pool61.picard.bam
pool62.picard.bam
pool63.picard.bam
pool64.picard.bam
pool65.picard.bam
pool66.picard.bam
pool67.picard.bam
pool68.picard.bam
pool69.picard.bam
pool70.picard.bam
pool71.picard.bam
pool72.picard.bam
pool74.picard.bam
pool75.picard.bam
pool76.picard.bam
pool77.picard.bam
pool78.picard.bam
pool79.picard.bam
pool80.picard.bam
pool87.picard.bam
pool88.picard.bam
pool89.picard.bam
pool90.picard.bam
pool91.picard.bam
```

### Index reference D. affinis ST female genome file


```python
module load samtools
samtools faidx Daffinis_STfemale_v5.1.masked.fasta
```

### GATK Variant calls over ChrX


```python
#!/bin/bash
#SBATCH --nodes=1 --ntasks=1
#SBATCH --cpus-per-task=16 
#SBATCH --time=60-00:00:00
#SBATCH --mail-type=END,FAIL     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=anjaligupta@ku.edu   # Where to send mail	
#SBATCH --mem=100GB
#SBATCH --partition=unckless
#SBATCH --chdir=/work/unckless/a948g501/Evol2024/PopGenPCA


module load gatk
gatk HaplotypeCaller \
     -I DaffPopGenBams.list \
     -O GATK_PopGen_ChrX_haploid.vcf \
     -R Daffinis_STfemale_v5.1.masked.fasta \
     --add-output-vcf-command-line true \
     -ploidy 1 \
     -L ChrX_MullerAD \
     --native-pair-hmm-threads 16 
```

### GATK Variant calls over Autosomes


```python
#!/bin/bash
#SBATCH --nodes=1 --ntasks=1
#SBATCH --cpus-per-task=16 
#SBATCH --time=60-00:00:00
#SBATCH --mail-type=END,FAIL     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=anjaligupta@ku.edu   # Where to send mail	
#SBATCH --mem=100GB
#SBATCH --partition=unckless
#SBATCH --chdir=/work/unckless/a948g501/Evol2024/PopGenPCA


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

### Filter VCF --> convert to Beagle --> Do PCA 

#### ChrX


```python
#!/bin/bash
#SBATCH --job-name=VCFToBeagleToPCA_X  # Job name
#SBATCH --partition=sixhour      # Partition Name (Required)
#SBATCH --mail-type=END,FAIL     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=anjaligupta@ku.edu   # Where to send mail	
#SBATCH --ntasks=8
#SBATCH --chdir=/work/unckless/a948g501/Evol2024/PopGenPCA
#SBATCH --cpus-per-task=1          # Run on a single CPU
#SBATCH --mem=64gb           # Job memory request
#SBATCH --time=06:00:00        # Time limit days-hrs:min:sec
#SBATCH --output=VCFToBeagleToPCA_X_%j.log  # Standard output and error log


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
#!/bin/bash
#SBATCH --job-name=VCFToBeagleToPCA_auto  # Job name
#SBATCH --partition=unckless      # Partition Name (Required)
#SBATCH --mail-type=END,FAIL     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=anjaligupta@ku.edu   # Where to send mail	
#SBATCH --ntasks=8
#SBATCH --chdir=/work/unckless/a948g501/Evol2024/PopGenPCA
#SBATCH --cpus-per-task=1          # Run on a single CPU
#SBATCH --mem=64gb           # Job memory request
#SBATCH --time=06-00:00:00        # Time limit days-hrs:min:sec
#SBATCH --output=VCFToBeagleToPCA_auto_%j.log  # Standard output and error log


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

### Partition ChrX filtered VCF aligning to inversion breakpoints

First I need to zip and index my vcf using samtools bgzip and tabix


```python
VCF=GATK_PopGen_ChrX_haploid.vcf

module load samtools

bgzip "filtered_SNPs_"$VCF

tabix "filtered_SNPs_"$VCF".gz"
```

### PCA for left arm and right arm of ChrX & for pre and post inversion breakpoints

Now I will partition my vcf file into left and right arm sections assuming the centromere is at 30 Mb

And partitioning vcf for the two single inversions using inversion breakpoint positions that Rob gave me


```python
VCF=filtered_SNPs_GATK_PopGen_ChrX_haploid.vcf.gz

module load bcftools

bcftools view -r ChrX_MullerAD:1-30000000 $VCF -o chrX_XL.vcf

bcftools view -r ChrX_MullerAD:30000000- $VCF -o chrX_XR.vcf

bcftools view -r ChrX_MullerAD:1-5490855 $VCF -o chrX_XL_PreInv.vcf

bcftools view -r ChrX_MullerAD:5490855-10233152 $VCF -o chrX_XL_Inv.vcf

bcftools view -r ChrX_MullerAD:10233152-30000000 $VCF -o chrX_XL_PostInv.vcf

bcftools view -r ChrX_MullerAD:30000000-44023827 $VCF -o chrX_XR_PreInv.vcf

bcftools view -r ChrX_MullerAD:44023827-62734159 $VCF -o chrX_XR_Inv.vcf

bcftools view -r ChrX_MullerAD:62734159- $VCF -o chrX_XR_PostInv.vcf
```

Convert VCFs to Beagle format


```python
module load vcftools

for VCF in "chrX_XL" "chrX_XR" "chrX_XL_PreInv" "chrX_XL_Inv" "chrX_XL_PostInv" "chrX_XR_PreInv" "chrX_XR_Inv" "chrX_XR_PostInv"
do
vcftools --vcf $VCF".vcf" --chr ChrX_MullerAD --BEAGLE-PL --out $VCF
gzip $VCF".BEAGLE.PL"
done
```

PCA


```python
#!/bin/bash
#SBATCH --job-name=PCA_X_part  # Job name
#SBATCH --partition=sixhour      # Partition Name (Required)
#SBATCH --mail-type=END,FAIL     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=anjaligupta@ku.edu   # Where to send mail	
#SBATCH --ntasks=8
#SBATCH --chdir=/work/unckless/a948g501/Evol2024/PopGenPCA
#SBATCH --cpus-per-task=1          # Run on a single CPU
#SBATCH --mem=64gb           # Job memory request
#SBATCH --time=06:00:00        # Time limit days-hrs:min:sec
#SBATCH --output=PCA_X_part_%j.log  # Standard output and error log


module load conda

conda activate pcangsd

for VCF in "chrX_XL" "chrX_XR" "chrX_XL_PreInv" "chrX_XL_Inv" "chrX_XL_PostInv" "chrX_XR_PreInv" "chrX_XR_Inv" "chrX_XR_PostInv"
do
pcangsd -b $VCF".BEAGLE.PL.gz" -e 2 -t 64 -o "pcangsd"$VCF
done


conda deactivate
```

## Nucleotide diversity from vcf


```python
mkdir -p /work/unckless/a948g501/Evol2024/NucleotideDiversityPopGen

cd /work/unckless/a948g501/Evol2024/NucleotideDiversityPopGen
```


```python
cp /work/unckless/a948g501/Evol2024/PopGenPCA/filtered_SNPs_GATK_PopGen_ChrX_haploid.vcf.gz .
```


```python
bcftools view filtered_SNPs_GATK_PopGen_ChrX_haploid.vcf.gz | sed -r -e 's@\t([a-Z0-9.]):@\t\1/\1:@g' | bcftools view -Oz -o filtered_SNPs_GATK_PopGen_ChrX_hap2dip.vcf.gz

bcftools +fixploidy filtered_SNPs_GATK_PopGen_ChrX_haploid.vcf.gz -- -f 2 > filtered_SNPs_GATK_PopGen_ChrX_diploid.vcf

```


```python
ml vcftools

vcftools --vcf filtered_SNPs_GATK_PopGen_ChrX_diploid.vcf --TajimaD 100000 --out TajimaDPopGen_100kb

vcftools --vcf filtered_SNPs_GATK_PopGen_ChrX_diploid.vcf --window-pi 100000 --window-pi-step 100000 --out pi_PopGen_100kb

vcftools --vcf filtered_SNPs_GATK_PopGen_ChrX_diploid.vcf --weir-fst-pop Daff_ST.txt --weir-fst-pop Daff_SR1.txt --weir-fst-pop Daff_SR2.txt --fst-window-size 100000 --fst-window-step 100000 --out FstPopGenAll_100kb 

vcftools --vcf filtered_SNPs_GATK_PopGen_ChrX_diploid.vcf --weir-fst-pop Daff_ST.txt --weir-fst-pop Daff_SR1.txt --fst-window-size 100000 --fst-window-step 100000 --out FstPopGen_STSR1_100kb

vcftools --vcf filtered_SNPs_GATK_PopGen_ChrX_diploid.vcf --weir-fst-pop Daff_ST.txt --weir-fst-pop Daff_SR2.txt --fst-window-size 100000 --fst-window-step 100000 --out FstPopGen_STSR2_100kb

vcftools --vcf filtered_SNPs_GATK_PopGen_ChrX_diploid.vcf --weir-fst-pop Daff_SR2.txt --weir-fst-pop Daff_SR1.txt --fst-window-size 100000 --fst-window-step 100000 --out FstPopGen_SR1SR2_100kb

```


```python
ml vcftools

vcftools --vcf filtered_SNPs_GATK_PopGen_ChrX_diploid.vcf --TajimaD 1000000 --out TajimaDPopGen_1Mb

vcftools --vcf filtered_SNPs_GATK_PopGen_ChrX_diploid.vcf --window-pi 1000000 --window-pi-step 1000000 --out pi_PopGen_1Mb

vcftools --vcf filtered_SNPs_GATK_PopGen_ChrX_diploid.vcf --weir-fst-pop Daff_ST.txt --weir-fst-pop Daff_SR1.txt --weir-fst-pop Daff_SR2.txt --fst-window-size 1000000 --fst-window-step 1000000 --out FstPopGenAll_1Mb 

vcftools --vcf filtered_SNPs_GATK_PopGen_ChrX_diploid.vcf --weir-fst-pop Daff_ST.txt --weir-fst-pop Daff_SR1.txt --fst-window-size 1000000 --fst-window-step 1000000 --out FstPopGen_STSR1_1Mb

vcftools --vcf filtered_SNPs_GATK_PopGen_ChrX_diploid.vcf --weir-fst-pop Daff_ST.txt --weir-fst-pop Daff_SR2.txt --fst-window-size 1000000 --fst-window-step 1000000 --out FstPopGen_STSR2_1Mb

vcftools --vcf filtered_SNPs_GATK_PopGen_ChrX_diploid.vcf --weir-fst-pop Daff_SR2.txt --weir-fst-pop Daff_SR1.txt --fst-window-size 1000000 --fst-window-step 1000000 --out FstPopGen_SR1SR2_1Mb

```

### Following speciation genomics

https://speciationgenomics.github.io/sliding_windows/

#### Convert vcf file to geno.gz


```python
cd /work/unckless/a948g501/Evol2024/PopulationGenomicsStats
```


```python
ml python
ml samtools

python /work/unckless/a948g501/Evol2024/PopulationGenomicsStats/genomics_general/VCF_processing/parseVCF.py \
    -i /work/unckless/a948g501/Evol2024/NucleotideDiversityPopGen/filtered_SNPs_GATK_PopGen_ChrX_diploid.vcf | bgzip > filtered_SNPs_GATK_PopGen_ChrX_diploid.geno.gz
```


```python
python /work/unckless/a948g501/Evol2024/PopulationGenomicsStats/genomics_general/popgenWindows.py \
    -g filtered_SNPs_GATK_PopGen_ChrX_diploid.geno.gz -o SNPs_GATK_PopGen_ChrX_diploid.Fst.Dxy.pi.csv.gz \
   -f phased -w 100000 -m 100 -s 100000 \
   -p ST -p SR1 -p SR2 \
   --popsFile ind.info
```


```python
python /work/unckless/a948g501/Evol2024/PopulationGenomicsStats/genomics_general/popgenWindows.py \
    -g filtered_SNPs_GATK_PopGen_ChrX_diploid.geno.gz -o SNPs_GATK_PopGen_ChrX_diploid.Fst.Dxy.pi_1mb.csv.gz \
   -f phased -w 1000000 -m 100 -s 1000000 \
   -p ST -p SR1 -p SR2 \
   --popsFile ind.info
```

#### Selection scans

https://speciationgenomics.github.io/phasing/


```python
cd /work/unckless/a948g501/Evol2024/PopulationGenomicsStats/SelectiveSweep
```


```python
wget https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r904.glibcv2.17.linux.tar.gz

tar -zxvf shapeit.v2.r904.glibcv2.17.linux.tar.gz 
```


```python
SHAPEIT_EXEC=/work/unckless/a948g501/Evol2024/PopulationGenomicsStats/SelectiveSweep/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit

$SHAPEIT_EXEC --help
```


```python
cp /work/unckless/a948g501/Evol2024/NucleotideDiversityPopGen/filtered_SNPs_GATK_PopGen_ChrX_diploid.vcf .

ml samtools
bgzip filtered_SNPs_GATK_PopGen_ChrX_diploid.vcf

bcftools view -m2 -M2 -v snps filtered_SNPs_GATK_PopGen_ChrX_diploid.vcf.gz -Oz -o filtered_SNPs_GATK_PopGen_ChrX_diploid_biallelic.vcf.gz

VCF="filtered_SNPs_GATK_PopGen_ChrX_diploid_biallelic.vcf.gz"

tabix $VCF

ml bcftools

bcftools query -l $VCF | wc -l
bcftools index -n $VCF


OUTPUT="filtered_SNPs_GATK_PopGen_ChrX_diploid_biallelic_phased"


```


```python
#!/bin/bash
#SBATCH --job-name=Shapeit  # Job name
#SBATCH --partition=unckless      # Partition Name (Required)
#SBATCH --mail-type=END,FAIL     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=anjaligupta@ku.edu   # Where to send mail	
#SBATCH --ntasks=8
#SBATCH --chdir=/work/unckless/a948g501/Evol2024/PopulationGenomicsStats/SelectiveSweep
#SBATCH --cpus-per-task=1          # Run on a single CPU
#SBATCH --mem=64gb           # Job memory request
#SBATCH --time=10-00:00:00        # Time limit days-hrs:min:sec
#SBATCH --output=Shapeit_%j.log  # Standard output and error log


SHAPEIT_EXEC=/work/unckless/a948g501/Evol2024/PopulationGenomicsStats/SelectiveSweep/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit
VCF="filtered_SNPs_GATK_PopGen_ChrX_diploid_biallelic.vcf.gz"
OUTPUT="filtered_SNPs_GATK_PopGen_ChrX_diploid_biallelic_phased"

$SHAPEIT_EXEC --input-vcf $VCF \
 -O $OUTPUT \
--window 1 -T 8 \
--force
```


```python
SHAPEIT_EXEC=/work/unckless/a948g501/Evol2024/PopulationGenomicsStats/SelectiveSweep/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit
OUTPUT="filtered_SNPs_GATK_PopGen_ChrX_diploid_biallelic_phased"

$SHAPEIT_EXEC -convert \
--input-haps ${OUTPUT} \
--output-vcf ${OUTPUT}.vcf

ml samtools
ml bcftools

bgzip ${OUTPUT}.vcf
bcftools index ${OUTPUT}.vcf.gz
```

Looking at how this is different from our standard vcf


```python
bcftools view -H ${VCF} | head | cut -f 1-12
bcftools view -H ${OUTPUT}.vcf.gz | head | cut -f 1-12
```


```python
# set a new variable
VCF=/work/unckless/a948g501/Evol2024/PopulationGenomicsStats/SelectiveSweep/filtered_SNPs_GATK_PopGen_ChrX_diploid_biallelic_phased.vcf.gz

# look at the sample names
bcftools query -l $VCF


# extract sample names
cp *.txt /work/unckless/a948g501/Evol2024/PopulationGenomicsStats/SelectiveSweep/
```


```python
VCF=/work/unckless/a948g501/Evol2024/PopulationGenomicsStats/SelectiveSweep/filtered_SNPs_GATK_PopGen_ChrX_diploid_biallelic_phased.vcf.gz

bcftools view -S Daff_SR1.txt -O z -o Daff_SR1.vcf.gz $VCF
bcftools view -S Daff_SR2.txt -O z -o Daff_SR2.vcf.gz $VCF
bcftools view -S Daff_ST.txt -O z -o Daff_ST.vcf.gz $VCF
```

#### Coalescent fastsimcoal2

https://speciationgenomics.github.io/easysfs/


```python
git clone https://github.com/isaacovercast/easySFS.git
```


```python
#!/bin/bash
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=06:00:00
#SBATCH --mem=20GB
#SBATCH --partition=sixhour
#SBATCH --chdir=/work/unckless/a948g501/Evol2024/PopulationGenomicsStats

ml python

python /work/unckless/a948g501/Evol2024/PopulationGenomicsStats/easySFS/easySFS.py \
-i /work/unckless/a948g501/Evol2024/PopGenPCA/filtered_SNPs_GATK_PopGen_ChrX_haploid.vcf \
-p /work/unckless/a948g501/Evol2024/PopulationGenomicsStats/ind.info \
--preview -a \
--ploidy 1 

```


```python
SR1
(2, 155834.0)	(3, 232452.0)	(4, 283855.0)	(5, 322519.0)	(6, 353219.0)	(7, 378315.0)	(8, 399100.0)	(9, 416230.0)	(10, 429821.0)	(11, 439499.0)	(12, 444775.0)	(13, 444571.0)	(14, 435751.0)	(15, 414920.0)	(16, 376460.0)	(17, 314415.0)	(18, 227191.0)	(19, 125185.0)	(20, 39676.0)	

SR2
(2, 102591.0)	(3, 151589.0)	(4, 183594.0)	(5, 178655.0)	(6, 106011.0)	

ST
(2, 242029.0)	(3, 362262.0)	(4, 448978.0)	(5, 519242.0)	(6, 579217.0)	(7, 631967.0)	(8, 679429.0)	(9, 722659.0)	(10, 762377.0)	(11, 799065.0)	(12, 833112.0)	(13, 864724.0)	(14, 894134.0)	(15, 921400.0)	(16, 946809.0)	(17, 970445.0)	(18, 992102.0)	(19, 1012002.0)	(20, 1029689.0)	(21, 1045504.0)	(22, 1059350.0)	(23, 1070761.0)	(24, 1079637.0)	(25, 1086065.0)	(26, 1089433.0)	(27, 1090039.0)	(28, 1087252.0)	(29, 1080243.0)	(30, 1069337.0)	(31, 1053526.0)	(32, 1032175.0)	(33, 1004375.0)	(34, 969609.0)	(35, 926451.0)	(36, 874425.0)	(37, 812325.0)	(38, 739724.0)	(39, 656665.0)	(40, 564918.0)	(41, 468340.0)	(42, 369614.0)	(43, 275959.0)	(44, 192496.0)	(45, 124192.0)	(46, 72630.0)	(47, 37899.0)	(48, 17598.0)	(49, 7166.0)	(50, 2543.0)	(51, 909.0)	(52, 313.0)	
```


```python
#!/bin/bash
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=06:00:00
#SBATCH --mem=20GB
#SBATCH --partition=sixhour
#SBATCH --chdir=/work/unckless/a948g501/Evol2024/PopulationGenomicsStats

ml python

python ceasySFS/easySFS.py \
-i /work/unckless/a948g501/Evol2024/PopGenPCA/filtered_SNPs_GATK_PopGen_ChrX_haploid.vcf \
-p /work/unckless/a948g501/Evol2024/PopulationGenomicsStats/ind.info -a -f \
--proj 12,4,27 \
--ploidy 1 \
-o outputSFS

```


```python
ml conda
conda activate fastsimcoal2

fsc27093 --help

conda deactivate
```

Follow https://speciationgenomics.github.io/fastsimcoal2/ for coalescent analysis
