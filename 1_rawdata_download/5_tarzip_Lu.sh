#!/bin/bash
#SBATCH --job-name=tarzip
#SBATCH --mem=100G
#SBATCH --partition=hmem
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=24

cd $PWD

#Package and compress, and delete the original file
tar -zcvf 1_rawdata_fastq.tar.gz ./1_rawdata_fastq     
rm -rf 1_rawdata_fastq
