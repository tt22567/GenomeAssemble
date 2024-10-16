#!/bin/bash

#SBATCH --job-name=unite
#SBATCH --output=./unite_%A.out
#SBATCH --mem=80G
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=24

source /user/work/tt22567/Miniconda3/ENTER/bin/activate   /user/work/tt22567/Miniconda3/ENTER/envs/Lacewings

cd $PWD

PWD=/user/work/tt22567/Neuroptera/5_buscogene
UniteBUSCOs=/user/work/tt22567/Neuroptera/5_buscogene/UniteBUSCOs.py

target_file=$PWD/DIR_gene
logs_file=$PWD/DIR_log

unite_file=DIR_Unite
mkdir -p $unite_file

chmod -R 777 $UniteBUSCOs 

echo -e "#######################\nStarting to unite the BUSCOs of same species.\n$(date)\n#######################\n\n" >> $logs_file/UniteEcho.log 
$UniteBUSCOs --indir $target_file --outdir $unite_file --sortmissing >> $logs_file/UniteBUSCOs.log 2>&1 && 
echo -e "#######################\nBUSCOs united.\n$(date)\n#######################\n\n" >> $logs_file/UniteEcho.log

