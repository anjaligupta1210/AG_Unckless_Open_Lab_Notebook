---
layout: post
title: Repeat Masker and modeler on long reads assembly
date: 05 July 2024  
category: [ Computational Pipelines ]
tags: [ RepeatMasker, RepeatModeler, Mummer ]
---

# Running RepeatModeler and RepeatMasker on *D. affinis* SR2 genome 

### RepeatModeler


```python
#!/bin/bash
#SBATCH --job-name=RMSR2  # Job name
#SBATCH --partition=unckless      # Partition Name (Required)
#SBATCH --mail-type=END,FAIL     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=anjaligupta@ku.edu   # Where to send mail	
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1          # Run on a single CPU
#SBATCH --mem=64gb           # Job memory request
#SBATCH --time=60-00:00:00        # Time limit days-hrs:min:sec
#SBATCH --output=RMSR2_%j.log  # Standard output and error log
#SBATCH --chdir=/work/unckless/a948g501/SR2_Affinis

module load repeatmodeler
module load repeatmasker/4.0.9

echo "STARTING"

BuildDatabase -name Daff_SR2 -engine ncbi SR2ONT.pilon3.fasta

RepeatModeler -engine ncbi -pa 8 -database Daff_SR2

echo "DONE"

```

### RepeatMasker


```python
#!/bin/bash
#SBATCH --job-name=RMSR2  # Job name
#SBATCH --partition=unckless      # Partition Name (Required)
#SBATCH --mail-type=END,FAIL     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=anjaligupta@ku.edu   # Where to send mail	
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1          # Run on a single CPU
#SBATCH --mem=64gb           # Job memory request
#SBATCH --time=60-00:00:00        # Time limit days-hrs:min:sec
#SBATCH --output=RMSR2_%j.log  # Standard output and error log
#SBATCH --chdir=/work/unckless/a948g501/SR2_Affinis

module load repeatmodeler
module load repeatmasker/4.0.9

echo "STARTING"

RepeatMasker -pa 8 -gff -lib Daff_SR2-families.fa -dir MaskerOutput.DaffSR2 SR2ONT.pilon3.fasta

echo "DONE"

```


```python
cp SR2ONT.pilon3.fasta.masked SR2ONT.pilon3.masked.fasta
```

## Mummer


```python
#!/bin/bash
#SBATCH --job-name=MummerSR2  # Job name
#SBATCH --partition=unckless      # Partition Name (Required)
#SBATCH --mail-type=END,FAIL     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=anjaligupta@ku.edu   # Where to send mail	
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1          # Run on a single CPU
#SBATCH --mem=64gb           # Job memory request
#SBATCH --time=60-00:00:00        # Time limit days-hrs:min:sec
#SBATCH --output=MummerSR2_%j.log  # Standard output and error log
#SBATCH --chdir=/work/unckless/a948g501/SR2_Affinis/MaskerOutput.DaffSR2

module load mummer

echo "STARTING"

mummer Daffinis_STfemale_v5.1.masked.fasta SR2ONT.pilon3.masked.fasta 

echo "DONE"
```

## Nucmer


```python
#!/bin/bash
#SBATCH --job-name=NucmerSR2  # Job name
#SBATCH --partition=unckless      # Partition Name (Required)
#SBATCH --mail-type=END,FAIL     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=anjaligupta@ku.edu   # Where to send mail	
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1          # Run on a single CPU
#SBATCH --mem=64gb           # Job memory request
#SBATCH --time=60-00:00:00        # Time limit days-hrs:min:sec
#SBATCH --output=MummerSR2_%j.log  # Standard output and error log
#SBATCH --chdir=/work/unckless/a948g501/SR2_Affinis/MaskerOutput.DaffSR2

module load mummer

echo "STARTING"

nucmer Daffinis_STfemale_v5.1.masked.fasta SR2ONT.pilon3.masked.fasta -p SR2_masked

echo "DONE"
```


```python
mummerplot --prefix=SR2_masked --postscript --filter SR2_masked.delta
```

## Show-coords


```python
show-coords -rcl SR2_masked.delta > SR2_masked.coords
```


```python
show-coords -Br SR2_masked.delta > SR2_masked.btab
```
