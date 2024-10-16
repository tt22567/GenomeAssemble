#! /bin/bash
#SBATCH --job-name=inRnaBUSCO
#SBATCH --mem=60G
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=24

cd $PWD

PWD=/user/work/tt22567/Neuroptera/5_BuscoNonDecont
DIR_trans=/user/work/tt22567/Neuroptera/4_trifasta/3_outRnaTri  #1_inRnaTri, 3_outRnaTri
DIR_db=/user/work/tt22567/Database/Busco_db/endopterygota_odb10
DIR_busco=/user/work/tt22567/Miniconda3/ENTER/envs/Megaloptera/lib/python3.9/site-packages/busco
DIR_Longest_isoform=/user/work/tt22567/Miniconda3/ENTER/envs/Lacewings/opt/trinity-2.8.5/util/misc
DIR_cdhit=/user/work/tt22567/Miniconda3/ENTER/envs/Lacewings/bin

#from_Tri means we do Busco without Decontmination
DIR_output1=./DIR_outfaa  #DIR_infaa_intrans, DIR_outfaa_intrans
DIR_output2=./DIR_Busco_outrans  #DIR_Busco_intrans, DIR_Busco_outrans
DIR_logfile=./DIR_log

db_name=endopterygota_odb10
lengths_cutoff_eff=2

#Intrans,Outrans
log_file1=Outrans_L.pep.cdhit.log
log_file2=Outrans_Busco.log
csv_file=summarizeBusco_outrans.csv

mkdir -p $DIR_output1
mkdir -p $DIR_output2
mkdir -p $DIR_logfile

source /user/work/tt22567/Miniconda3/ENTER/bin/activate   /user/work/tt22567/Miniconda3/ENTER/envs/Lacewings

# Transdecode and cdhit trinity assemblies
for rnatri in $DIR_trans/*_trinity.Trinity.fasta;do
   id=$(echo $rnatri | sed -e "s#$DIR_trans\/##g")
   od=$(echo $rnatri | sed -e "s#$DIR_trans\/##g;s#\_trinity.Trinity.fasta##g")

#transdecoder find LongOrfs and predict protein
if [[ ! -f $DIR_output1/$od.L.trans.pep.faa ]]; then
    echo "\n========================\nSelecting the longest isoform.\n$(date)\n========================\n" >> $DIR_logfile/$log_file1 && 
    perl $DIR_Longest_isoform/get_longest_isoform_seq_per_trinity_gene.pl $rnatri >> $DIR_output1/$od.L.fasta && 
        
    echo "\n========================\nTransDecoding $rnatri LongestORF starts.\n$(date)\n========================\n" >> $DIR_logfile/$log_file1 && 
    TransDecoder.LongOrfs -t $DIR_output1/$od.L.fasta >> $DIR_logfile/$log_file1 2>&1 && 
    TransDecoder.Predict -t $DIR_output1/$od.L.fasta >> $DIR_logfile/$log_file1 2>&1 && 
    
	mv $PWD/"$od"*.transdecoder.* $DIR_output1 && 
    mv $DIR_output1/"$od"*.transdecoder.pep $DIR_output1/$od.L.trans.pep.faa && 
    rm -rf $PWD/$od*.transdecoder_dir* &&  
	rm -rf $PWD/pipeliner.*.cmds &&
    echo "\n========================\nTransDecoding $rnatri LongestORF finished.\n$(date)\n========================\n" >> $DIR_logfile/$log_file1 
else
    echo -e "\n########################\nTransdecoding has already been accomplished, no need to do it again.\n$(date)\n########################\n" >> $DIR_logfile/$log_file1 
fi

#Cdhit removing duplications
if [[ ! -f $DIR_output1/$od.L.trans.pep.cdhit.faa ]]; then
    echo "\n========================\nremoving the duplications of $DIR_output1/$od.L.trans.pep.faa starts.\n$(date)\n========================\n" >> $DIR_logfile/$log_file1 && 
    $DIR_cdhit/cd-hit-est -i $DIR_output1/$od.L.trans.pep.faa -o $DIR_output1/$od.L.trans.pep.cdhit.faa -c 0.98 -n 10 &&
	rm -rf $DIR_output1/$od.L.trans.pep.cdhit.faa.clstr &&
	echo "\n========================\nremoving the duplications of $DIR_output1/$od.L.trans.pep.faa is done.\n$(date)\n========================\n" >> $DIR_logfile/$log_file1
else
    echo "\n========================\nremoving the duplications of $DIR_output1/$od.L.trans.pep.faa is finished, no need to do it again.\n$(date)\n========================\n" >> $DIR_logfile/$log_file1
fi
done

source /user/work/tt22567/Miniconda3/ENTER/bin/deactivate   /user/work/tt22567/Miniconda3/ENTER/envs/Lacewings
source /user/work/tt22567/Miniconda3/ENTER/bin/activate   /user/work/tt22567/Miniconda3/ENTER/envs/Megaloptera

#a way to include more frament gene into complete gene
if [ ! -f $DIR_db/lengths_cutoff.ori ];then 
    echo "preparing lengths_cutoff by doubling its $3 starts" >> $DIR_logfile/$log_file2
    cp $DIR_db/lengths_cutoff $DIR_db/lengths_cutoff.ori &&
    awk -v x=$lengths_cutoff_eff 'BEGIN{OFS=FS=","}{gen = $1; B = $2; delta = $3; len = $4; ndelta = delta * x; print gen, B, ndelta, len }' $DIR_db/lengths_cutoff.ori > $DIR_db/lengths_cutoff &&
    echo "lengths_cutoff is prepared with doubled $3" >> $DIR_logfile/$log_file2
else 
    echo "lengths_cutoff is prepared, no need to do it again" >> $DIR_logfile/$log_file2
fi

#Busco assessment and get single-copy othology gene
for i in $DIR_output1/*.L.trans.pep.cdhit.faa;do
    id=$(echo $i | sed -e "s#$DIR_output1\/##g")
    od=$id\_$db_name
	txt=short_summary.specific.$db_name.$od.txt
if [ ! -d $DIR_output2/$od ] || [ ! -f $DIR_output2/$od/$txt ];then
    rm -rf $DIR_output2/$od
    echo "Busco assessment and geting single-copy othology genes of $i start" >> $DIR_logfile/$log_file2
    python $DIR_busco/run_BUSCO.py -m proteins -i $i -l $DIR_db -o $DIR_output2/$od --cpu 10 &&
    echo "Busco assessment and geting single-copy othology genes of $i is finished" >> $DIR_logfile/$log_file2
else 
    echo "Busco assessment and geting single-copy othology genes of $i is already done, no need to do it again" >> $DIR_logfile/$log_file2
fi
done

#--config $DIR_config/config.ini
#DIR_config=/user/work/tt22567/Miniconda3/ENTER/envs/Megaloptera/config
#-o output_file_name

#summarizeBUSCO_assessment
touch $csv_file
echo -e "id,c,s,d,f,m,t,cp,sp,dp,fp,mp,md,bv,hs,ld,cd,ng,tm" >>$csv_file

for document in $DIR_output2/*$db_name;do 
  dm=$(echo $document | sed "s#$DIR_output2\/##g") 
  id=$(echo $document | sed "s#$DIR_output2\/##g;s#\_$db_name##g") 
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
	echo -e "$id,$c,$s,$d,$f,$m,$t,$cp,$sp,$dp,$fp,$mp,$md,$bv,$hs,$ld,$cd,$ng,$(date)" >>$csv_file
  done
done

#unitBusco
