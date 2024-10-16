#! /bin/bash

cd $PWD

#this script is to assess the assembly quality of fasta

PWD=/user/work/tt22567/Neuroptera/3_tridecont/in_Fasta #before assembly:1_inRnaTri,2_inDnaTri,3_outRnaTri,4_outDnaTri;after assembly:in_Fasta,out_Fasta
DIR_seqkit=/user/work/tt22567/Miniconda3/ENTER/envs/Lacewings/bin

target_group=inRnaDXTri  #inRnaTri,outRnaTri,inDnaTri,outDnaTri,inRnaDXTri;outRnaDXTri
target_seq=trinity.Trinity.Sifted.DeContX.fasta  #trinity.Trinity.fasta；scaffolds.fasta;trinity.Trinity.Sifted.DeContX.fasta

DIR_logfile=DIR_N50
log_file=echo.log
txt_file=N50_$target_group.txt
txt_rev=N50_$target_group\_rev.txt
csv_file=$target_group\_N50.csv

mkdir -p $DIR_logfile

source /user/work/tt22567/Miniconda3/ENTER/bin/activate   /user/work/tt22567/Miniconda3/ENTER/envs/Lacewings

#N50 assessment
for i in $PWD/*_$target_seq;do
    echo "assessment the assembly of $i starts" >>$DIR_logfile/$log_file &&
    $DIR_seqkit/seqkit stat -a $i 2>&1 >> $DIR_logfile/$txt_file  &&
    echo "assessment the assembly of $i is done" >>$DIR_logfile/$log_file
done

#prepare for txt_rev: Here we replace comma with full stop, so we can manage the csv file well, but remember to change them back in the final csv file.
if [ ! -f $DIR_logfile/$txt_rev ];then 
cat $DIR_logfile/$txt_file | sed "s#,#.#g;s#%#p#g" >> $DIR_logfile/$txt_rev
fi

#transfer txt to csv file
touch $DIR_logfile/$csv_file
echo -e "id,format,type,num_seqs,sum_len,min_len,avg_len,max_len,Q1,Q2,Q3,sum_gap,N50,Q20P,Q30P,GCP" >>$DIR_logfile/$csv_file

for i in $PWD/*_$target_seq;do
  
  id=$(echo $i | sed "s#$PWD\/##g")
  format=$(cat $DIR_logfile/$txt_rev | grep "$i" | awk -F "[ ]*" '{print $2}')
  type=$(cat $DIR_logfile/$txt_rev | grep "$i" | awk -F "[ ]*" '{print $3}')
  num_seqs=$(cat $DIR_logfile/$txt_rev | grep "$i" | awk -F "[ ]*" '{print $4}')
  sum_len=$(cat $DIR_logfile/$txt_rev | grep "$i" | awk -F "[ ]*" '{print $5}')
  min_len=$(cat $DIR_logfile/$txt_rev | grep "$i" | awk -F "[ ]*" '{print $6}')
  avg_len=$(cat $DIR_logfile/$txt_rev | grep "$i" | awk -F "[ ]*" '{print $7}')
  max_len=$(cat $DIR_logfile/$txt_rev | grep "$i" | awk -F "[ ]*" '{print $8}')
  Q1=$(cat $DIR_logfile/$txt_rev | grep "$i" | awk -F "[ ]*" '{print $9}')
  Q2=$(cat $DIR_logfile/$txt_rev | grep "$i" | awk -F "[ ]*" '{print $10}')
  Q3=$(cat $DIR_logfile/$txt_rev | grep "$i" | awk -F "[ ]*" '{print $11}')
  sum_gap=$(cat $DIR_logfile/$txt_rev | grep "$i" | awk -F "[ ]*" '{print $12}')
  N50=$(cat $DIR_logfile/$txt_rev | grep "$i" | awk -F "[ ]*" '{print $13}')
  Q20P=$(cat $DIR_logfile/$txt_rev | grep "$i" | awk -F "[ ]*" '{print $14}')
  Q30P=$(cat $DIR_logfile/$txt_rev | grep "$i" | awk -F "[ ]*" '{print $15}')
  GCP=$(cat $DIR_logfile/$txt_rev | grep "$i" | awk -F "[ ]*" '{print $16}')
  
  echo -e "$id,$format,$type,$num_seqs,$sum_len,$min_len,$avg_len,$max_len,$Q1,$Q2,$Q3,$sum_gap,$N50,$Q20P,$Q30P,$GCP" >>$DIR_logfile/$csv_file
  #reminder: change the comma back
done 



  
