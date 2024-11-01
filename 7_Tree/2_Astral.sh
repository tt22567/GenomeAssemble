#!/bin/bash
#SBATCH --job-name=7_Astral
#SBATCH --account=gely018542
#SBATCH --partition=cpu
#SBATCH --mem=80GB
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --time=14-00:00:00
#SBATCH --output=./7_Astral_%A.out

NT=$[$SLURM_NNODES*$SLURM_NTASKS_PER_NODE]
FNP=$[$SLURM_NTASKS-2]

source /user/work/tt22567/Miniconda3/ENTER/bin/activate   /user/work/tt22567/Miniconda3/ENTER/envs/twist

cd $PWD

PWD=/user/work/tt22567/Neuroptera/6_alignment

long_file=./DIR_long
Astral_file=./DIR_Astral
logs_file=./DIR_log

mkdir -p $long_file
mkdir -p $Astral_file
mkdir -p $logs_file

module add apps/astral/5.7.8

if [ ! -f $long_file/Astral_input.treefile ];then
  echo "constructing the single species-tree with Astral starts" >>$logs_file/5.2_Astral_tree.log &&
  cat $Astral_file/*treefile > $Astral_file/Astral_input.treefile &&
  java -jar /mnt/storage/software/apps/ASTRAL/astral.5.7.8.jar -i $Astral_file/Astral_input.treefile -o $Astral_file/Astral_output.treefile >>$logs_file/5.1_IQtree_Astral.log 2>&1 &&
  echo "constructing the single species-tree with Astral is finished" >>$logs_file/5.2_Astral_tree.log
fi

