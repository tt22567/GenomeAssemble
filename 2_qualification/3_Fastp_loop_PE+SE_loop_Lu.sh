#!/bin/bash

#SBATCH --job-name=Fastp
#SBATCH --mem=80G
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=24

#SBATCH --time=1-00:00:00

cd $PWD
output=/mnt/storage/scratch/tt22567/Neuroptera/3_cleanfastq/1_inrna_clean

#For paired-end sequences
for i in *_1.*
do 
id=${i%%_1.*}

  if [ ! -f ./html_json/$id.html ] 
   then 
    echo "$id starts" >>inrnafastp.log
    /user/home/tt22567/.conda/envs/Lacewings/bin/fastp -i $id\_1.* -o $id\_1.clean.gz -I $id\_2.* -O $id\_2.clean.gz -c -W 4 -M 20 -3 -q 15 -u 40 -n 10 -l 36 -h ./html_json/$id.html -j ./html_json/$id.json;mv $id\_1.clean.gz $id\_2.clean.gz $output &&
    echo "$id is done" >>inrnafastp.log
  fi
  done 


#For single-end sequences
for i in *rna.fastq.gz
do 
id=${i%%.fastq.gz}

  if [ ! -f ./html_json/$id.html ] 
   then
    echo "$id starts" >>inrnafastp.log
    /user/home/tt22567/.conda/envs/Lacewings/bin/fastp -i $id.fastq.gz -o $id\_clean.gz -W 4 -M 20 -3 -q 15 -u 40 -n 10 -l 36 -h ./html_json/$id.html -j ./html_json/$id.json &&
	mv $id\_clean.gz $output &&
    echo "$id is done" >>inrnafastp.log
  fi
  done 