#!/bin/bash

#SBATCH --job-name=Spades—4
#SBATCH --mem=300G
#SBATCH --partition=hmem
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=24

cd $PWD

module load languages/anaconda3/2021-3.9-bioconda
module load apps/spades/3.15.2
STOR_PATH=/mnt/storage/scratch/tt22567/Neuroptera/4_trifasta/4_outnDaTri #path used for transfer the assembly

##spades can only recognize sequences with specific endings, such as fq.gz, so the name needs to be changed
#for name in `ls *.clean.gz`
#do
#mv $name ${name%%.gz}.fq.gz
#done 

# check if there is scaffolds.fasta in the folder
#for document in *_spades
#do
#docuid=${document%%_spades}
#docunm=$docuid\_spades
  #if [ ! -f $docunm/scaffolds.fasta ]
  #then
  #rm -rf $docunm
  #fi
#done

#start the assembling
for i in *_1.clean.fq.gz
do 
id=${i%%_1.clean.fq.gz}
fq1=$id\_1.clean.fq.gz
fq2=$id\_2.clean.fq.gz
docunm=$id\_spades
 if [ ! -f  $docunm/$id\_scaffolds.fasta ] 
 then
 echo "Assembling $id starts"
 python3 /sw/apps/SPAdes-3.15.2-Linux/bin/spades.py -1 $fq1 -2 $fq2 -o $id\_spades -t 22 --isolate >>spades.log 2>&1 &&
 mv $docunm/scaffolds.fasta $docunm/$id\_scaffolds.fasta && cp $docunm/$id\_scaffolds.fasta $PATH
 echo "Assembling $id is done" >>spades.log
 fi
done

#>>spades.log 2>&1 this step could be ignored, as log file will automatically exist


