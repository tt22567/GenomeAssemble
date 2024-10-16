#!/bin/bash

#SBATCH --job-name=taxonomap
#SBATCH --mem=80G
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=24

module load languages/anaconda3/2021-3.9-bioconda
conda activate Lacewings

cd /mnt/storage/scratch/tt22567/Database/taxonomap_db

#wget -c -t 0 https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz

nr_db_path=/mnt/storage/scratch/tt22567/Database/nr_db
nr_path=$nr_db_path/nr
nr_dmnd_path=$nr_db_path/nr.dmnd
taxonomap_path=./prot.accession2taxid.gz
taxonnodes_path=/mnt/storage/scratch/tt22567/Database/taxdump_db/taxdump/nodes.dmp
#taxonnames_path=/mnt/storage/scratch/tt22567/Database/taxdump_db/taxdump/names.dmp

diamond makedb --in $nr_path --db $nr_dmnd_path --taxonmap $taxonomap_path --taxonnodes $taxonnodes_path 
#--taxonnames $taxonnames_path 

