#!/bin/bash

#SBATCH --job-name=UniVec
#SBATCH --mem=50G
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=24

source /sw/languages/anaconda/anaconda.3.9-2021.12-bioconda/bin/activate   /user/home/tt22567/.conda/envs/Lacewings  #path to the env where blast installed

cd $PWD
UniVec=/mnt/storage/scratch/tt22567/Database/UniVec_db/UniVec
makeblastdb_PATH=/user/home/tt22567/.conda/envs/Lacewings/bin

#wget -c -t 0 https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec 

# Build the index
$makeblastdb_PATH/makeblastdb -in $UniVec -input_type fasta -dbtype nucl -title UniVectitle -out UniVec.dmnd -parse_seqids

# contigs.fasta is the database, CL is the database name, and out will be used as the -db parameter in future BLAST+ searches. You can give it a meaningful name.
# -in is followed by the input file, which is the FASTA sequence you want to format.
# -dbtype is followed by the sequence type: 'nucl' for nucleotides, 'prot' for proteins.
# -title gives the database a name (just for display; this cannot be used as the -db parameter in searches).
# -parse_seqids is recommended to add, although the exact reason is unclear for now.
# -out is followed by the database name, which you can choose to be meaningful. It will be used as the -db parameter in future BLAST+ searches.
