#!/bin/bash
#SBATCH --job-name=GenDown
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=24Gb

DIR_current=$(pwd)
DIR_csv=$DIR_current/$1
NLine=$(cat $DIR_csv |wc -l)

if [ -f $DIR_current/CLine.log ]; then 
    st=$(sed -n 1p $DIR_current/CLine.log)
else
    st=2
fi

mkdir -p $DIR_current/assembly
mkdir -p $DIR_current/rnaseq
mkdir -p $DIR_current/wgs

####Multithread preparation
tmp_fifofile="/tmp/$$.fifo"
trap "exec 9>&-;exec 9<&-;exit 0" 2
mkfifo $tmp_fifofile
exec 9<>$tmp_fifofile 
rm $tmp_fifofile

for ((i=1;i<=20;i++))
do
    echo >&9
done

for ((n=$st;n<=$NLine;n++))
do
    line=$(sed -n "$n"p $DIR_csv | sed -e 's/\r//g;s/ftp:/https:/g') # remove ending \r and replace ftp with https
    orgn=$(echo $line | awk -F ',' '{print$1}')
    accn=$(echo $line | awk -F ',' '{print$2}')
    type=$(echo $line | awk -F ',' '{print$3}')
    furl=$(echo $line | awk -F ',' '{print$4}')

    opfx="$orgn"_"$accn"

    echo "*****************"
    echo "Requesting $opfx"

        case $type in 
            "Assembly_Chromosome")
                wget -q -O ./assembly/$opfx.fna.gz $furl && 
                # curl -s -o ./assembly/$opfx.fna.gz $furl && 
                gzip -d ./assembly/$opfx.fna.gz 
                # wget -q -O ./assembly/$opfx.faa.gz ${furl%/*}/*_protein.faa.gz && 
                # if [ -s "./assembly/$opfx.faa.gz" ]; then
                #     gzip -d ./assembly/$opfx.faa.gz
                # else
                #     rm ./assembly/$opfx.faa.gz
                # fi
            ;;

            "Assembly_Contig" | "Assembly_Scaffold")
                wget -q -O ./assembly/$opfx.fna.gz $furl && 
                gzip -d ./assembly/$opfx.fna.gz
            ;;

            "RNA-Seq-TRANSCRIPTOMIC")
                wget -q -O ./rnaseq/$opfx.rna $furl 
                /user/home/tt22567/.conda/envs/Lacewings/bin/fastq-dump -e 6 -p --split-3 --skip-technical ./assembly/$opfx.rna
            ;;

            "WGS-GENOMIC")
                wget -q -O ./wgs/$opfx.dna $furl 
                /user/home/tt22567/.conda/envs/Lacewings/bin/fastq-dump -e 6 -p --split-3 --skip-technical ./assembly/$opfx.rna
            ;;
        esac
    wait 
    echo $n > CLine.log
done

#Xiumei Lu modified from LL
#csv should be modified with four column: Organism  Run_num	Library_source	Download_address
#i.e. Neoneuromus_ignobilis	SRR19893053	RNA-Seq-TRANSCRIPTOMIC	gs://sra-pub-run-5/SRR19893053/SRR19893053.1
#i.e. Chrysoperla_carnea	SRR10418476	WGS-GENOMIC	gs://sra-pub-run-12/SRR10418476/SRR10418476.1

