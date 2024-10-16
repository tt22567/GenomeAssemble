#!/bin/bash

#SBATCH --job-name=taxdump
#SBATCH --mem=50G
#SBATCH --partition=hmem
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=24

cd $PWD

wget -c -t 0 https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz 
tar zxf taxdump.tar.gz
