#! /bin/bash
#SBATCH --job-name=inDnBuFly
#SBATCH --mem=80G
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=24

cd $PWD

PWD=/user/work/tt22567/Neuroptera/5_buscogene
DIR_genomes=/user/work/tt22567/Neuroptera/4_trifasta/2_inDnaTri  #2_inDnaTri,4_outDnaTri
DIR_db=/user/work/tt22567/Database/Busco_db/endopterygota_odb10
DIR_busco=/user/work/tt22567/Miniconda3/ENTER/envs/Megaloptera/lib/python3.9/site-packages/busco
DIR_logfile=./DIR_log
DIR_csvfile=./DIR_Busco_csv

lengths_cutoff_eff=2
db_name=endopterygota_odb10
ref_species=aedes
#aedes伊蚊；culex库蚊；Anopheles_gambiae冈比亚蚊；fly苍蝇；
#bombus_impatiens1美洲东部熊蜂；bombus_terrestris2熊蜂；nasonia金小蜂；honeybee1蜜蜂；
#heliconius_melpomene1袖蝶；
#tribolium2012赤拟谷盗

#ingenome,outgenome
DIR_output=./DIR_Busco_ingenome
DIR_csvfile=./DIR_Busco_csv
log_file=busco_ingenome_$ref_species.log
csv_file=summarizeBusco_ingenome_$ref_species.csv

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

#Busco assessment and get single-copy othology gene
for i in $DIR_genomes/*_scaffolds.fasta;do
  id=$(echo $i | sed -e "s#$DIR_genomes\/##g")
  od=$id\_$ref_species\_$db_name
  txt=short_summary.specific.$db_name.$od.txt
  if [ ! -d $DIR_output/$od ] || [ ! -f $DIR_output/$od/$txt ]; then
    rm -rf $DIR_output/$od
	echo "Busco assessment and getting single-copy othology genes of $i start" >> $DIR_logfile/$log_file
    python $DIR_busco/run_BUSCO.py -m genome -i $i -l $DIR_db -o $DIR_output/$od --cpu 10 --augustus_species $ref_species &&
	echo "Busco assessment and getting single-copy othology genes of $i finished" >> $DIR_logfile/$log_file
  fi
done

#--config $DIR_config/config.ini --augustus_species $species_name
#DIR_config=/user/work/tt22567/Miniconda3/ENTER/envs/Megaloptera/config
#-o output_file_name

#summarizeBUSCO_assessment
if [ ! -f $DIR_csvfile/$csv_file ];then
  touch $DIR_csvfile/$csv_file
  echo -e "id,rf,c,s,d,f,m,t,cp,sp,dp,fp,mp,ns,nc,tl,pgp,sn50,cn50,bv,hs,bb,mb,tb,ag,md,ld,cd,ng,tm" >>$DIR_csvfile/$csv_file

  for document in $DIR_output/*_$ref_species\_$db_name;do 
    dm=$(echo $document | sed "s#$DIR_output\/##g")
    id=$(echo $document | sed "s#$DIR_output\/##g;s#_$ref_species\_$db_name##g")
    txt=short_summary.specific.$db_name.$dm.txt
    for tt in $document/$txt;do
      rf=$ref_species
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
	  #Assembly Statistics
	  ns=$(cat $tt | grep "Number of scaffolds" | awk '{print $1}')
	  nc=$(cat $tt | grep "Number of contigs" | awk '{print $1}')
	  tl=$(cat $tt | grep "Total length" | awk '{print $1}')
	  pgp=$(cat $tt | grep "Percent gaps" | awk -F "%" '{print $1}')
	  sn50=$(cat $tt | grep "Scaffold N50" | awk '{print $1}')
	  cn50=$(cat $tt | grep "Contigs N50" | awk '{print $1}')
	  #softwares (Busco,hmmsearch,bbtools,makeblastdb,tblastn,agustus)
	  bv=$(cat $tt | grep "# BUSCO version" | awk -F ": " '{print $2}')
	  hs=$(cat $tt | grep "hmmsearch" | awk -F ": " '{print $2}')
	  bb=$(cat $tt | grep "bbtools" | awk -F ": " '{print $2}')
	  mb=$(cat $tt | grep "makeblastdb" | awk -F ": " '{print $2}')
	  tb=$(cat $tt | grep "tblastn" | awk -F ": " '{print $2}')
	  ag=$(cat $tt | grep "augustus: " | awk -F ": " '{print $2}')
	  #related-info(mode,linage dataset,creation date,number of genomes)
	  md=$(cat $tt | grep "# BUSCO was run in mode" | awk -F ": " '{print $2}')
	  ld=$(cat $tt | grep "# The lineage dataset" | awk -F ": " '{print $2}' | awk -F "(" '{print $1}')
	  cd=$(cat $tt | grep "Creation date" | awk -F "Creation date: " '{print $2}' | awk -F "," '{print $1}')
	  ng=$(cat $tt | grep "number of genomes: " | awk -F "number of genomes: " '{print $2}' | awk -F "," '{print $1}')
	  echo -e "$id,$rf,$c,$s,$d,$f,$m,$t,$cp,$sp,$dp,$fp,$mp,$ns,$nc,$tl,$pgp,$sn50,$cn50,$bv,$hs,$bb,$mb,$tb,$ag,$md,$ld,$cd,$ng,$(date)" >>$DIR_csvfile/$csv_file
    done
  done
fi

#gather single_copy_busco_sequences together and remove files with less than three faa
gene_file=DIR_gene 
mkdir -p $gene_file


for document in $DIR_output/*$db_name;do 
  id=$(echo $document | sed "s#$DIR_output\/##g")
  gen=$document/run_$db_name/busco_sequences/single_copy_busco_sequences
  od=in_$id  #in;out
  nm_faa=$(ls $gene_file/$od/*faa | wc -l)
  
  #gather single_copy_busco_sequences together
  if [ ! -d $gene_file/$od ];then
    echo "copy the single-copy-busco genes of $id starts" >> $DIR_logfile/$log_file &&
	cp -r $gen  $gene_file/$od  2>&1 && rm -rf $gene_file/$od/*fna && #outgroup exits both faa and fna fasta files.
	echo echo "copy the single-copy-busco genes of $id is finished, no need to do it again" >> $DIR_logfile/$log_file
  fi
  
  #remove files with less than three faa   
  for file in $gene_file/$od;do
    if [ $nm_faa <= 3 ];then 
	  echo "$gene_file/$od has less than 3 faa, and therefore should be deleted" >> $DIR_logfile/$log_file &&
	  rm -rf $gene_file/$od  &&
      echo "the delation of $gene_file/$od has finished, no need to do it again" >> $DIR_logfile/$log_file
	fi
  done
	  
done



