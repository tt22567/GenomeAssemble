#! /bin/bash

#SBATCH --job-name=MFP1_BMGE
#SBATCH --account=gely018542
#SBATCH --partition=cpu
#SBATCH --mem=80GB
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=10
#SBATCH --time=14-00:00:00
#SBATCH --output=./MFP1_BMGE_%A.out

NT=$[$SLURM_NNODES*$SLURM_NTASKS_PER_NODE]
FNP=$[$SLURM_NTASKS-2]

source /user/work/tt22567/Miniconda3/ENTER/bin/activate   /user/work/tt22567/Miniconda3/ENTER/envs/twist

cd $PWD
PWD=/user/work/tt22567/Neuroptera/6_treefile

target_file=./DIR_bmge
storage_path=./MFP_BMGE
log_file=./DIR_log

mkdir -p $storage_path
mkdir -p $log_file

# Multithread preparation
tmp_fifofile="/tmp/$$.fifo" 
trap "exec 9>&-;exec 9<&-;exit 0" 2 
mkfifo $tmp_fifofile
exec 9<>$tmp_fifofile 
rm $tmp_fifofile

for ((i=1;i<=$FNP;i++))
do
    echo >&9
done

for i in DIR_bmge/*mft.bmg.faa;do 
    read -u9
    {
  id=$(echo $i | sed 's#DIR_bmge\/##g')
  if [[ $(ls -l $i | awk '{ print $5 }') -gt 0  ]] && [ ! -f $storage_path/$id.treefile ];then
    sleep 5s
	if [ ! -s $storage_path/$id.log ];then
      echo "constructing phylogeny for $id with IQTREE starts" >>$log_file/MFP_BMGE.log
      iqtree -s $i -m MFP -T AUTO -bb 1000 -pre MFP_BMGE/$id -alninfo &&  
      rm -rf $storage_path/$id.bionj $storage_path/$id.ckp.gz $storage_path/$id.log $storage_path/$id.mldist &&
      echo "constructing phylogeny for $id with IQTREE is done" >>$log_file/MFP_BMGE.log
	fi
  fi 
    echo >&9 
    } & 
done

