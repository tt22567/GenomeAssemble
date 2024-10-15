#!/bin/bash

#SBATCH --job-name=Fastp
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --time=4-00:00:00

cd /mnt/storage/scratch/tt22567/Neuroptera_rawdata/Neuroptera_rawdata_GCP/Euclimacia

module load apps/sratoolkit/3.0.0

for i in *.gz; do fastp -i $i -o $i.clean; done 
