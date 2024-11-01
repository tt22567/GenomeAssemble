#!/bin/bash
#SBATCH --job-name=2_BMGE
#SBATCH --account=gely018542
#SBATCH --partition=cpu
#SBATCH --mem=80GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --time=14-00:00:00
#SBATCH --output=./2_BMGE_%A.out

source /user/work/tt22567/Miniconda3/ENTER/bin/activate   /user/work/tt22567/Miniconda3/ENTER/envs/twist

cd $PWD

PWD=/user/work/tt22567/Neuroptera/6_alignment

target_file=$PWD/DIR_mafft   
bmge_path=/user/work/tt22567/Miniconda3/ENTER/envs/twist/share/bmge-1.12-1
seqkit_path=/user/work/tt22567/Miniconda3/ENTER/envs/twist/bin

mafft_file=./DIR_mafft
bmge_file=./DIR_bmge
logs_file=./7_DIR_log

mkdir -p $mafft_file
mkdir -p $bmge_file
mkdir -p $logs_file

#trimming alignments with BMGE
for i in $mafft_file/BSC_*_mft.faa;do
  id=$(echo $i | sed "s#$mafft_file\/##g")
  nm=$(echo $i | sed -e "s#$mafft_file\/##g;s#.faa##g")
  od=$nm\.bmg.faa

  if [ ! -f $bmge_file/$od ];then
    echo "trimming $mafft_file/$id with bmge start" >>$logs_file/2_bmge.log  &&
	java -Xmx90G -jar $bmge_path/BMGE.jar -i $mafft_file/$id -t AA -m BLOSUM90 -h 0.5 -s NO -of $logs_file/$od >>$logs_file/2_bmge.log 2>&1 &&
	echo "trimming $mafft_file/$id with bmge are finished" >>$logs_file/2_bmge.log
  else
    echo "trimming $mafft_file/$id with BMGE is already done, no need to do it again" >>$logs_file/2_bmge.log
  fi
done


