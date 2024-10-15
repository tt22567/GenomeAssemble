#!/bin/bash

#SBATCH --job-name=sra_release
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --time=1-00:00:00

cd /mnt/storage/scratch/tt22567/Neuroptera_rawdata/Neuroptera_rawdata_AWS/rnaseq
module load apps/sratoolkit/3.0.0

mkdir releaseq

# Decompress rawdata into paired-end sequences (fasterq-dump does not have compression functionality, fastq-dump has compression functionality)
for i in *.rna; do fasterq-dump -e 6 -p --split-3 --skip-technical ./$i; done &&
# Compress the decompressed paired-end sequences into rna_1/2.fastq.gz and move them into the releaseq folder (the gzip command is built-in, not part of sratoolkit)
for i in *rna_?.fastq; do gzip ./$i; cp *rna_?.fastq.gz releaseq; done &&
# Sometimes there may be a folder, which can be deleted
for i in *rna_?.fastq; do rm *rna_?.fastq; done

