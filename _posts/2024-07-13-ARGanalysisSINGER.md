---
layout: post
title: ARG analysis using SINGER - doesn't work yet
date: 13 July 2024  
category: [ Computational Pipelines ]
tags: [ ARG, SINGER ]
---

### Ancestral Recombination Graphs

Doesn't work yet !!!

### Using SINGER

#### Installation


```python
mkdir /work/unckless/a948g501/Evol2024/Singer

cd /work/unckless/a948g501/Evol2024/Singer
```


```python
git clone https://github.com/popgenmethods/SINGER.git
```


```python
cp /work/unckless/a948g501/Evol2024/PopGenPCA/filtered_SNPs_GATK_PopGen_ChrX_haploid.vcf.gz .

module load samtools

bgzip -d filtered_SNPs_GATK_PopGen_ChrX_haploid.vcf.gz 
```


```python


ml python
ml conda
conda activate tskit


python /work/unckless/a948g501/Evol2024/Singer/SINGER/SINGER/SINGER/parallel_singer \
    -vcf filtered_SNPs_GATK_PopGen_ChrX_haploid \
    -Ne 2e5 \
        -m 1e-9 \
              -n 10000 \
                  -thin 50 \
                    -output ARGsingerOut

conda deactivate
```
