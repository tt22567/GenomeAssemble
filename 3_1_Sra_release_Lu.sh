#!/bin/bash

#SBATCH --job-name=sra release
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --time=1-00:00:00

cd /mnt/storage/scratch/tt22567/Neuroptera_rawdata/Neuroptera_rawdata_AWS/rnaseq

module load apps/sratoolkit/3.0.0

for i in *.rna; do fasterq-dump -e 6 -p --split-3 --skip-technical ./&i; done