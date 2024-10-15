#! /bin/bash
#SBATCH --job-name=in_megahit
#SBATCH --mem=250G
#SBATCH --partition=hmem
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --time=10-00:00:00

source /user/work/tt22567/Miniconda3/ENTER/bin/activate   /user/work/tt22567/Miniconda3/ENTER/envs/Lacewings

DIR_Megahit=/user/work/tt22567/Miniconda3/ENTER/envs/Lacewings/bin

for i in *_1.clean.fq.gz;do 
  id=${i%%_1.clean.fq.gz}
  fq1=$id\_1.clean.fq.gz
  fq2=$id\_2.clean.fq.gz
  docunm=$id\_spades
  documn1=$id\_Megahit
  if [ ! -f  $docunm/$id\_scaffolds.fasta ];then
    echo "Assembling $id with megahit starts" >>spades.log
    $DIR_Megahit/megahit -1 $fq1 -2 $fq2 -t 24 -m 0.95 --min-contig-len 500 --out-dir $documn1 --out-prefix $id >>spades.log 2>&1 &&
	echo "Assembling $id with megahit finished" >>spades.log
  fi
done
