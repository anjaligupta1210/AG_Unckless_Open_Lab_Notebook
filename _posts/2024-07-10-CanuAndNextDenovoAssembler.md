---
layout: post
title: Assembling nanopore long reads using Canu and NextDenovo
date: 10 July 2024  
category: [ Computational Pipelines ]
tags: [ NanoporeLongReadsAssembly, Canu, NextDenovo ]
---

### Assembling *D. affinis* SR2 genome using Canu

Make working directory


```python
mkdir /work/unckless/a948g501/Evol2024/AssemblySR2
```

Copy all my raw Nanopore long reads to the working directory


```python
scp a948g501@dtn.ku.edu:/resfs/GROUPS/MB/NIH0072614/Daffinis/RawReads/SR141ONT/allreads.fastq.gz .
```


```python
chmod +777 allreads.fastq.gz
```

Run Correction using CanuAssembler

`CanuCorrect.sh`


```python
#!/bin/bash
#SBATCH --job-name=CanuCorrect  # Job name
#SBATCH --partition=unckless      # Partition Name (Required)
#SBATCH --mail-type=END,FAIL     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=anjaligupta@ku.edu   # Where to send mail
#SBATCH --chdir=/work/unckless/a948g501/Evol2024/AssemblySR2/
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1          # Run on a single CPU
#SBATCH --mem=64gb           # Job memory request
#SBATCH --time=60-00:00:00        # Time limit days-hrs:min:sec
#SBATCH --output=CanuCorrect_%j.log  # Standard output and error log


module load conda
conda activate canu

canu -correct \
-p DaffSR2canucorrect -d canuCorrect \
genomeSize=185000000 \
-untrimmed -nanopore allreads.fastq.gz \
gridOptions="--partition=unckless"

conda deactivate
```

Assemble corrected reads


```python
#!/bin/bash
#SBATCH --job-name=CanuAssemble  # Job name
#SBATCH --partition=unckless      # Partition Name (Required)
#SBATCH --mail-type=END,FAIL     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=anjaligupta@ku.edu   # Where to send mail
#SBATCH --chdir=/work/unckless/a948g501/Evol2024/AssemblySR2/
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1          # Run on a single CPU
#SBATCH --mem=64gb           # Job memory request
#SBATCH --time=60-00:00:00        # Time limit days-hrs:min:sec
#SBATCH --output=CanuAssemble_%j.log  # Standard output and error log


module load conda
conda activate canu

canu -assemble \
-p DaffSR2canuassemble -d canuAssemble \
genomeSize=185000000 \
-corrected -nanopore canuCorrect/DaffSR2canucorrect.correctedReads.fasta.gz \
gridOptions="--partition=unckless"

conda deactivate
```

### Assembling *D. affinis* SR2 genome using NextDenovo

#### Installation


```python
wget https://github.com/Nextomics/NextDenovo/releases/latest/download/NextDenovo.tgz
tar -vxzf NextDenovo.tgz && cd NextDenovo
```


```python
cd /work/unckless/a948g501/Evol2024/AssemblySR2/NextDenovo/NextDenovo
```

#### Test if NextDenovo works


```python
module load python

python nextDenovo test_data/run.cfg
```

#### Running for SR2 genome assembly


```python
mkdir -p /work/unckless/a948g501/Evol2024/AssemblySR2/NextDenovo/NextDenovo/DaffSR2

#Prepare input.fofn
scp a948g501@dtn.ku.edu:/resfs/GROUPS/MB/NIH0072614/Daffinis/RawReads/SR141ONT/allreads.fastq.gz /work/unckless/a948g501/Evol2024/AssemblySR2/NextDenovo/NextDenovo/DaffSR2

ls /work/unckless/a948g501/Evol2024/AssemblySR2/NextDenovo/NextDenovo/DaffSR2/allreads.fastq.gz > input.fofn

#Create run.cfg

cp /work/unckless/a948g501/Evol2024/AssemblySR2/NextDenovo/NextDenovo/doc/run.cfg /work/unckless/a948g501/Evol2024/AssemblySR2/NextDenovo/NextDenovo/DaffSR2/
```

Edit `run.cfg` - 


```python
[General]
job_type = local # local, slurm, sge, pbs, lsf
job_prefix = SR2nextDenovo
task = all # all, correct, assemble
rewrite = yes # yes/no
deltmp = yes 
parallel_jobs = 8 # number of tasks used to run in parallel
input_type = raw # raw, corrected
read_type = ont # clr, ont, hifi
input_fofn = /work/unckless/a948g501/Evol2024/AssemblySR2/NextDenovo/NextDenovo/DaffSR2/input.fofn
workdir = /work/unckless/a948g501/Evol2024/AssemblySR2/NextDenovo/NextDenovo/Daff_SR2_NextDenovo_Again

[correct_option]
read_cutoff = 1k
genome_size = 185000000 # estimated genome size
sort_options = -m 20g -t 15
minimap2_options_raw = -t 8
pa_correction = 3 # number of corrected tasks used to run in parallel, each corrected task requires ~TOTAL_INPUT_BASES/4 bytes of memory usage.
correction_options = -p 15

[assemble_option]
minimap2_options_cns = -t 8 
nextgraph_options = -a 1

# see https://nextdenovo.readthedocs.io/en/latest/OPTION.html for a detailed introduction about all the parameters
```


```python
mkdir /work/unckless/a948g501/Evol2024/AssemblySR2/NextDenovo/NextDenovo/Daff_SR2_NextDenovo
```


```python
#!/bin/bash
#SBATCH --job-name=NextDenovo # Job name
#SBATCH --partition=unckless      # Partition Name (Required)
#SBATCH --mail-type=END,FAIL     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=anjaligupta@ku.edu   # Where to send mail
#SBATCH --chdir=/work/unckless/a948g501/Evol2024/AssemblySR2/NextDenovo/NextDenovo
#SBATCH --ntasks=8          # Number of tasks (adjust based on total CPUs and your needs)
#SBATCH --cpus-per-task=1     # CPUs per task
#SBATCH --mem=256GB    # Job memory request per CPU (adjust as needed)
#SBATCH --time=60-00:00:00    # Time limit days-hrs:min:sec
#SBATCH --output=NextDenovo_%j.log  # Standard output and error log


module python
module load conda
conda activate paralleltask

python nextDenovo /work/unckless/a948g501/Evol2024/AssemblySR2/NextDenovo/NextDenovo/DaffSR2/run.cfg

conda deactivate
```

### Running RepeatModeler and RepeatMasker on Canu Assembly SR2


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
#SBATCH --chdir=/work/unckless/a948g501/Evol2024/AssemblySR2/RM

module load repeatmodeler
module load repeatmasker/4.0.9

echo "STARTING"

BuildDatabase -name Daff_SR2_canu -engine ncbi /work/unckless/a948g501/Evol2024/AssemblySR2/canuAssemble/DaffSR2canuassemble.contigs.fasta

RepeatModeler -engine ncbi -pa 8 -database Daff_SR2_canu

RepeatMasker -pa 8 -gff -lib Daff_SR2_canu-families.fa -dir MaskerOutput.DaffSR2_canu /work/unckless/a948g501/Evol2024/AssemblySR2/canuAssemble/DaffSR2canuassemble.contigs.fasta

cp /work/unckless/a948g501/Evol2024/AssemblySR2/RM/MaskerOutput.DaffSR2_canu/DaffSR2canuassemble.contigs.fasta.masked DaffSR2canuassemble.contigs.masked.fasta

echo "DONE"

```
