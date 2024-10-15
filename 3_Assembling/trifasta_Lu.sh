#!/bin/bash

#SBATCH --job-name=inr_trifasta
#SBATCH --mem=180G
#SBATCH --partition=hmem
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=24

cd $PWD
output=/mnt/storage/scratch/tt22567/Neuroptera/4_trifasta/1_inRnaTri

module load languages/anaconda3/2021-3.9-bioconda
conda activate Lacewings

for i in *_1.*
do id=${i%%_1.*}
if [ ! -f $id_trinity.Trinity.fasta ]
Trinity --normalize_reads --seqType fq --max_memory 100G --left $id\_1.clean.gz --right $id\_2.clean.gz --output $id\_trinity --CPU 8 --full_cleanup --inchworm_cpu 8
fi
done  
