#!/bin/bash
#SBATCH --job-name=nr_db_index
#SBATCH --nodes=1
#SBATCH --partition=hmem
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=24Gb

cd $PWD

module load languages/anaconda3/2021-3.9-bioconda
conda activate Lacewings

# -c means resuming a download; -t specifies the number of seconds to retry after a connection is lost, 0 means retry indefinitely. Note: Make sure to execute this command in the original directory, otherwise it will restart the download if the file is not found.
# wget -c -t 0 https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz 

# Decompress nr.gz to obtain nr, the original file will be removed
# /usr/bin/gunzip nr.gz 

# Create a database to get nr.dmnd (the input file follows --in, which is the decompressed file, the output file follows -d, and the output format is .dmnd)
diamond makedb --in nr -d nr
