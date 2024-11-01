#!/bin/bash

#SBATCH --job-name=1_mafft
#SBATCH --account=gely018542
#SBATCH --partition=cpu
#SBATCH --mem=50GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --time=1-00:00:00
#SBATCH --output=./1_mafft_%A.out

source /user/work/tt22567/Miniconda3/ENTER/bin/activate   /user/work/tt22567/Miniconda3/ENTER/envs/twist

cd $PWD

PWD=/user/work/tt22567/Neuroptera/6_alignment

target_file=$PWD/DIR_Unite      
mafft_path=/user/work/tt22567/Miniconda3/ENTER/envs/twist/bin
seqkit_path=/user/work/tt22567/Miniconda3/ENTER/envs/twist/bin

mafft_file=./DIR_mafft
logs_file=./DIR_log
csv_file=./DIR_csv

mkdir -p $mafft_file
mkdir -p $logs_file
mkdir -p $csv_file

#calculating the number of genes before alignment
touch $csv_file/calculate.txt  &&
num_gens_befcoverage=$(ls -l $target_file/*faa | wc -l) 
echo -e "A total of $num_gens_befcoverage in total from BUSCO" >>$csv_file/1_calculate.txt
  
#aligning the genes
for i in $target_file/BSC_*.faa;do
  id=$(echo $i | sed "s#$target_file\/##g")
  nm=$(echo $i | sed -e "s#$target_file\/##g;s#.faa##g")
  od=$nm\_mft.faa
  
  if [ ! -f $mafft_file/$od ];then
    echo "maffting $i starts" >>$logs_file/1_mafft.log &&
    $mafft_path/mafft-linsi --thread 4 --leavegappyregion $target_file/$id > $mafft_file/$od && >>$logs_file/1_mafft.log 2>&1 && 
    echo "maffting $i is finished" >>$logs_file/1_mafft.log
  else 
    echo "maffting $i is already finished, no need to do it again" >>$logs_file/1_mafft.log
  fi  &&
done

#calculating the length of alignments after maffting
echo "calculating the length of alignments (including gaps) after maffting starts" >>$logs_file/1_mafft.log &&
touch $csv_file/2_seqkit_stats_maffted_align.txt &&
$seqkit_path/seqkit stats $mafft_file/BSC_*_mft.faa 2>&1 >> $csv_file/2_seqkit_stats_maffted_align.txt &&
echo "calculating the length of alignments (including gaps) after maffting is finished" >>$logs_file/1_mafft.log


