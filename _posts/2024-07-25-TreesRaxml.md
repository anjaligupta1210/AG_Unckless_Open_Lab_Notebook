---
layout: post
title: Trees using RaxML
date: 25 July 2024  
category: [ Computational Pipelines ]
tags: [ Sliding window trees, Comparative genomics, RaxML ]
---



## Sliding Trees using iqtree and raxML


```python
/work/unckless/a948g501/Evol2024/SlidingTrees
```

File list


```python
/work/unckless/a948g501/SlidingTrees/KUCGposter/chrX_XL_Inv.min4.phy
/work/unckless/a948g501/SlidingTrees/KUCGposter/chrX_XL.min4.phy
/work/unckless/a948g501/SlidingTrees/KUCGposter/chrX_XL_PostInv.min4.phy
/work/unckless/a948g501/SlidingTrees/KUCGposter/chrX_XL_PreInv.min4.phy
/work/unckless/a948g501/SlidingTrees/KUCGposter/chrX_XR_Inv.min4.phy
/work/unckless/a948g501/SlidingTrees/KUCGposter/chrX_XR.min4.phy
/work/unckless/a948g501/SlidingTrees/KUCGposter/chrX_XR_PostInv.min4.phy
/work/unckless/a948g501/SlidingTrees/KUCGposter/chrX_XR_PreInv.min4.phy
/work/unckless/a948g501/SlidingTrees/KUCGposter/GATKfilteredSNPAutosome_NoAztNoBlanks.min4.phy
/work/unckless/a948g501/SlidingTrees/KUCGposter/GATKfilteredSNPchrX_NoAztNoBlanks.min4.phy
```


```python
iqtree2 -s /work/unckless/a948g501/SlidingTrees/KUCGposter/chrX_XL_Inv.min4.phy -bb 1000 -wbt -nt AUTO -o Dpse -redo
iqtree2 -s /work/unckless/a948g501/SlidingTrees/KUCGposter/chrX_XL.min4.phy -bb 1000 -wbt -nt AUTO -o Dpse -redo
iqtree2 -s /work/unckless/a948g501/SlidingTrees/KUCGposter/chrX_XL_PostInv.min4.phy -bb 1000 -wbt -nt AUTO -o Dpse -redo
iqtree2 -s /work/unckless/a948g501/SlidingTrees/KUCGposter/chrX_XL_PreInv.min4.phy -bb 1000 -wbt -nt AUTO -o Dpse -redo
iqtree2 -s /work/unckless/a948g501/SlidingTrees/KUCGposter/chrX_XR_Inv.min4.phy -bb 1000 -wbt -nt AUTO -o Dpse -redo
iqtree2 -s /work/unckless/a948g501/SlidingTrees/KUCGposter/chrX_XR.min4.phy -bb 1000 -wbt -nt AUTO -o Dpse -redo
iqtree2 -s /work/unckless/a948g501/SlidingTrees/KUCGposter/chrX_XR_PostInv.min4.phy -bb 1000 -wbt -nt AUTO -o Dpse -redo
iqtree2 -s /work/unckless/a948g501/SlidingTrees/KUCGposter/chrX_XR_PreInv.min4.phy -bb 1000 -wbt -nt AUTO -o Dpse -redo
iqtree2 -s /work/unckless/a948g501/SlidingTrees/KUCGposter/GATKfilteredSNPAutosome_NoAztNoBlanks.min4.phy -bb 1000 -wbt -nt AUTO -o Dpse -redo
iqtree2 -s /work/unckless/a948g501/SlidingTrees/KUCGposter/GATKfilteredSNPchrX_NoAztNoBlanks.min4.phy -bb 1000 -wbt -nt AUTO -o Dpse -redo

```

## Trees using RaxML


```python
for phylip in "chrX_XL_Inv.min4.phy" "chrX_XL.min4.phy" "chrX_XL_PostInv.min4.phy" "chrX_XL_PreInv.min4.phy" "chrX_XR_Inv.min4.phy" "chrX_XR.min4.phy" "chrX_XR_PostInv.min4.phy" "chrX_XR_PreInv.min4.phy" "GATKfilteredSNPAutosome_NoAztNoBlanks.min4.phy" "GATKfilteredSNPchrX_NoAztNoBlanks.min4.phy"
do
sbatch RaxmlSupport.sh $phylip
done
```


```python
#!/bin/bash
#SBATCH --job-name=Raxml  # Job name
#SBATCH --partition=sjmac      # Partition Name (Required)
#SBATCH --mail-type=END,FAIL     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=anjaligupta@ku.edu   # Where to send mail
#SBATCH --chdir=/work/unckless/a948g501/Evol2024/SlidingTrees/Raxml_BootstrapSupport
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1          # Run on a single CPU
#SBATCH --mem=4gb           # Job memory request
#SBATCH --time=10-00:00:00        # Time limit days-hrs:min:sec
#SBATCH --output=SlidingTrees_%j.log  # Standard output and error log

phylip=$1

ml raxml-ng

raxml-ng --all \
--msa "/work/unckless/a948g501/Evol2024/SlidingTrees/Raxml_ML/"$phylip \
--msa-format PHYLIP \
--outgroup Dpse \
--model TVMef \
--threads 8 \
--bs-metric fbp,tbe \
--force
```


```python
ml conda
conda activate newick_utils
nw_display <Newick-tree-file>
```

## Trees for Runt by blasting it to Dpse 


```python
for genome in "Daffinis" "Dalgonquin" "DathabascaEA" "DathabascaEB" "Dpseudoobscura"
do
sbatch blast.sh $genome
done
```


```python
#!/bin/bash
#SBATCH --job-name=Blast  # Job name
#SBATCH --partition=sixhour     # Partition Name (Required)
#SBATCH --mail-type=END,FAIL     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=anjaligupta@ku.edu   # Where to send mail
#SBATCH --chdir=/work/unckless/a948g501/Evol2024/SlidingTrees/BlastRunt
#SBATCH --mem=4gb           # Job memory request
#SBATCH --time=06:00:00        # Time limit days-hrs:min:sec
#SBATCH --output=Blast_%j.log  # Standard output and error log

genome=$1

ml blast+s

makeblastdb -in $genome".fna" \
-dbtype nucl \
-out $genome"_db"

blastn -query Runt_xa.fna \
-db $genome"_db" \
-out $genome"_results.txt" \
-outfmt 6
```


```python
ml samtools

for genome in "Daffinis" "Dalgonquin" "DathabascaEA" "DathabascaEB" "Dpseudoobscura"
do
samtools faidx $genome".fna"
done
```


```python
for genome in "Daffinis" "Dalgonquin" "DathabascaEA" "DathabascaEB" "Dpseudoobscura"
do
sort -k 3,3gr -k 11,11g $genome"_results.txt" > sorted_tmp.txt
mv sorted_tmp.txt $genome"_results.txt"
done
```


```python
ml samtools

for genome in "Daffinis" "Dalgonquin" "DathabascaEA" "DathabascaEB" "Dpseudoobscura"
do
best_hit=$(head -n 1 $genome"_results.txt")
contig=$(echo "$best_hit" | awk '{print $2}')
output_file="Runt_contig_"$genome".fasta"
samtools faidx $genome".fna" "${contig}" > "$output_file"
done
```


```python
cat Runt_*fasta > Runt_Multi.fasta
cat Runt_xa.fna >> Runt_Multi.fasta
```


```python
#!/bin/bash
#SBATCH --job-name=RuntAlign  # Job name
#SBATCH --partition=sixhour     # Partition Name (Required)
#SBATCH --mail-type=END,FAIL     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=anjaligupta@ku.edu   # Where to send mail
#SBATCH --chdir=/work/unckless/a948g501/Evol2024/SlidingTrees/BlastRunt
#SBATCH --mem=4gb           # Job memory request
#SBATCH --time=06:00:00        # Time limit days-hrs:min:sec
#SBATCH --output=RuntAlign_%j.log  # Standard output and error log


ml muscle
muscle -in Runt_Multi.fasta -out "alignment_Runt".phy -maxiters 1 -diags1 
```


```python
ml iqtree
iqtree2 -s "alignment_Runt".phy
```
