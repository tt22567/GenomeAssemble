#! /bin/bash

#This script is used for the adding of order and family name for busco gene files ahead of the species name

cd $PWD

PWD=/user/work/tt22567/Neuroptera/5_buscogene
target_file=$PWD/DIR_gene
abb_file=$PWD/addfamilyname.csv
log_file=$PWD/DIR_log

ct=out  #in,out

for file in $target_file/$ct\_*ERR*_endpg10;do    #SRR OR ERR                                  #in_Agulla_sp._AD-2015_SRR1811836.rna.Trans.DX.L.pep.cdhit.faa_en
    id=$(echo $file | sed "s#$target_file\/out_##g")    #in or out                             #Agulla_sp._AD-2015_SRR1811836.rna.Trans.DX.L.pep.cdhit.faa_endpg10
	srr=ERR$(echo $id | awk -F "ERR" '{print $2}' | awk -F "." '{print $1}')  #SRR OR ERR      #SRR1811836   
    c1=$(cat $abb_file | grep "$srr" |  awk -F "," '{print $1}')                               #RAP_RAP
    new=$ct\_$c1\_$id                                                                          #in_RAP_RAP_Agulla_sp._AD-2015_SRR1811836.rna.Trans.DX.L.pep.cdhit.faa_endpg10
    
    echo "move $id to $new starts"  >>$log_file/addfamilyname.log
    mv $target_file\/$ct\_$id  $target_file\/$new &&
    echo "move $id to $new is finished"  >>$log_file/addfamilyname.log
done

#csv format
#SHORT	Order	Family	ScientificName	Run_num
#NEU_CHR	Neuroptera	Chrysopidae	Chrysopa pallens	SRR4181653


