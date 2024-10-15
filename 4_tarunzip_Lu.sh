#!/bin/bash

#SBATCH --job-name=tarzip
#SBATCH --mem=100G
#SBATCH --partition=hmem
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=24

cd $PWD

#unzip and remove
tar -zxvf 1_rawdata_fastq.tar.gz
rm -rf *rna
