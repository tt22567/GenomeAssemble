#! /bin/bash
#SBATCH --job-name=outRNABUSCO
#SBATCH --mem=60G
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=24

cd $PWD

PWD=/user/work/tt22567/Neuroptera/5_buscogene
DIR_trans=/user/work/tt22567/Neuroptera/4_trifasta/out_Faa #in_FAA,out_Faa
DIR_db=/user/work/tt22567/Database/Busco_db/endopterygota_odb10
DIR_busco=/user/work/tt22567/Miniconda3/ENTER/envs/Megaloptera/lib/python3.9/site-packages/busco
DIR_logfile=./DIR_log
DIR_csvfile=./DIR_Busco_csv

db_name=endopterygota_odb10
lengths_cutoff_eff=2

#intrans,outrans
DIR_output=./DIR_Busco_outrans
log_file=busco_outrans.log 
csv_file=summarizeBusco_outrans.csv

mkdir -p $DIR_output
mkdir -p $DIR_logfile
mkdir -p $DIR_csvfile

source /user/work/tt22567/Miniconda3/ENTER/bin/activate   /user/work/tt22567/Miniconda3/ENTER/envs/Megaloptera

#a way to include more frament gene into complete gene
if [ ! -f $DIR_db/lengths_cutoff.ori ];then 
  echo "preparing lengths_cutoff by doubling its $3 starts" >> $DIR_logfile/$log_file
  cp $DIR_db/lengths_cutoff $DIR_db/lengths_cutoff.ori &&
  awk -v x=$lengths_cutoff_eff 'BEGIN{OFS=FS="\t"}{gen = $1; B = $2; delta = $3; len = $4; ndelta = delta * x; print gen, B, ndelta, len }' $DIR_db/lengths_cutoff.ori > $DIR_db/lengths_cutoff &&
  echo "lengths_cutoff is prepared with doubled $3" >> $DIR_logfile/$log_file
else 
  echo "lengths_cutoff is prepared, no need to do it again" >> $DIR_logfile/$log_file
fi

#Busco assessment and get single-copy orthology gene
for i in $DIR_trans/*.Trans.DX.L.pep.cdhit.faa;do
  id=$(echo $i | sed -e "s#$DIR_trans\/##g")
  od=$id\_$db_name
  txt=short_summary.specific.$db_name.$od.txt
  if [ ! -d $DIR_output/$od ] || [ ! -f $DIR_output/$od/$txt ];then
    rm -rf $DIR_output/$od
    echo "Busco assessment and geting single-copy othology genes of $i start" >> $DIR_logfile/$log_file
    python $DIR_busco/run_BUSCO.py -m proteins -i $i -l $DIR_db -o $DIR_output/$od --cpu 10 &&
    echo "Busco assessment and geting single-copy othology genes of $i is finished" >> $DIR_logfile/$log_file
  else 
    echo "Busco assessment and geting single-copy othology genes of $i is already done, no need to do it again" >> $DIR_logfile/$log_file
  fi
done

#--config $DIR_config/config.ini --augustus_species aedes 
#DIR_config=/user/work/tt22567/Miniconda3/ENTER/envs/Megaloptera/config
#-o output_file_name

#summarizeBUSCO_assessment
if [ ! -f $DIR_csvfile/$csv_file ];then
  touch $DIR_csvfile/$csv_file
  echo -e "id,c,s,d,f,m,t,cp,sp,dp,fp,mp,md,bv,hs,ld,cd,ng,tm" >>$DIR_csvfile/$csv_file

  for document in $DIR_output/*$db_name;do 
    dm=$(echo $document | sed "s#$DIR_output\/##g")
    id=$(echo $document | sed "s#$DIR_output\/##g;s#\_$db_name##g")
    txt=short_summary.specific.$db_name.$dm.txt
    for tt in $document/$txt;do
      c=$(cat $tt | grep "Complete BUSCOs" | awk '{print $1}')
      s=$(cat $tt | grep "Complete and single-copy BUSCOs" | awk '{print $1}')
      d=$(cat $tt | grep "Complete and duplicated BUSCOs" | awk '{print $1}')
      f=$(cat $tt | grep "Fragmented BUSCOs" | awk '{print $1}')
      m=$(cat $tt | grep "Missing BUSCOs" | awk '{print $1}')
	  t=$(cat $tt | grep "Total BUSCO groups searched" | awk '{print $1}')
      #p(percentage)
	  cp=$(cat $tt | grep "C:" | awk -F "%" '{print $1}' | awk -F ":" '{print $2}')
	  sp=$(cat $tt | grep "C:" | awk -F "%" '{print $2}' | awk -F ":" '{print $2}')
	  dp=$(cat $tt | grep "C:" | awk -F "%" '{print $3}' | awk -F ":" '{print $2}')
	  fp=$(cat $tt | grep "C:" | awk -F "%" '{print $4}' | awk -F ":" '{print $2}')
	  mp=$(cat $tt | grep "C:" | awk -F "%" '{print $5}' | awk -F ":" '{print $2}')
	  #related infor-softwares (Busco and hmmsearch),mode,linage dataset,creation date,number of genomes
	  md=$(cat $tt | grep "# BUSCO was run in mode" | awk -F ": " '{print $2}')
	  bv=$(cat $tt | grep "# BUSCO version" | awk -F ": " '{print $2}')
	  hs=$(cat $tt | grep "hmmsearch" | awk -F ": " '{print $2}')
	  ld=$(cat $tt | grep "# The lineage dataset" | awk -F ": " '{print $2}' | awk -F "(" '{print $1}')
	  cd=$(cat $tt | grep "Creation date" | awk -F "Creation date: " '{print $2}' | awk -F "," '{print $1}')
	  ng=$(cat $tt | grep "number of genomes: " | awk -F "number of genomes: " '{print $2}' | awk -F "," '{print $1}')
	  echo -e "$id,$c,$s,$d,$f,$m,$t,$cp,$sp,$dp,$fp,$mp,$md,$bv,$hs,$ld,$cd,$ng,$(date)" >>$DIR_csvfile/$csv_file
    done
  done
fi

#gather single_copy_busco_sequences together
gene_file=DIR_gene 
mkdir -p $gene_file

for document in $DIR_output/*$db_name;do 
  id=$(echo $document | sed "s#$DIR_output\/##g")
  gen=$document/run_$db_name/busco_sequences/single_copy_busco_sequences
  od=out_$id
  if [ ! -d $gene_file/$od ];then
    echo "copy the single-copy-busco genes of $id starts" >> $DIR_logfile/$log_file &&
	cp -r $gen  $gene_file/$od  2>&1 && rm -rf $gene_file/$od/*fna && #outgroup exits both faa and fna fasta files.
	echo echo "copy the single-copy-busco genes of $id is finished, no need to do it again" >> $DIR_logfile/$log_file
  fi
done

