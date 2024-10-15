#!/bin/bash

#SBATCH --job-name=bbmap
#SBATCH --mem=280G
#SBATCH --partition=hmem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --time=5-00:00:00

cd $PWD

# For paired-end sequences with large size that cannot be assembled, remove duplicates first, then perform quality control and assembly.
source /sw/languages/anaconda/anaconda.3.9-2021.12-bioconda/bin/activate   /user/home/tt22567/.conda/envs/Lacewings

BBMP_PATH=/user/home/tt22567/.conda/envs/Lacewings/bin
FASTP_PATH=/user/home/tt22567/.conda/envs/Lacewings/bin

for i in *_1.fastq.gz
do
id=${i%%_1.fastq.gz}
fq1=$id\_1.fastq.gz
fq2=$id\_2.fastq.gz
outb1=$id\_1.ddp.fastq.gz
outb2=$id\_2.ddp.fastq.gz
outp1=$id\_1.ddp_clean.fastq.gz
outp2=$id\_2.ddp_clean.fastq.gz
$BBMP_PATH/clumpify.sh in1=$fq1 in2=$fq2 outb1=$outb1 outb2=$outb2 dedupe &&
$FASTP_PATH/fastp -i $outb1 -o $outp1 -I $outb2 -O $outp2 -c -W 4 -M 20 -3 -q 15 -u 40 -n 10 -l 36 -h ./html_json/$id.html -j ./html_json/$id.json &&
Trinity --normalize_reads --seqType fq --max_memory 120G --left $outp1 --right $outp2 --output $id\_trinity --CPU 8 --full_cleanup --inchworm_cpu 8
done
