#!/bin/sh
#SBATCH --job-name=out_Decont
#SBATCH --partition=cpu
#SBATCH --output=./outrnatri_%A.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=10
#SBATCH --mem=80GB
#SBATCH --time=10-00:00:00

DNP=$SLURM_CPUS_PER_TASK  # Number of CPUs per task
NT=$[$SLURM_NNODES*$SLURM_NTASKS_PER_NODE]  # Total number of tasks = number of nodes * number of tasks per node
FNP=$[$SLURM_NTASKS]  # Total number of tasks - 2  # If the number of nodes is changed to 1, remove the '-2' inside the square brackets

cd /user/work/tt22567/Neuroptera/4_trifasta
#storage file for *_trinity.Trinity.fasta

DIR_cur=$(pwd)
#DIR_cnd=/user/home/vn21703/anaconda3  #trinity,blast,transdecoder,diamond should be installed in one conda env, within the same virtural env.
UniVecDB=/user/work/tt22567/Database/UniVec_db/UniVec.dmnd # path to UniVec_db
NRDiaDB=/user/work/tt22567/Database/nr_db/nr.dmnd #path to nr_db
SeqOrinPiePlot=user/work/tt22567/Neuroptera/4_trifasta/SeqOriginPiePlot.py #path to script for drawing
TaxDumpDB=/user/work/tt22567/Database/taxdump_db/taxdump #path to taxdump
taxonomap_path=/user/work/tt22567/Database/taxonomap_db/prot.accession2taxid.gz  
DIR_asm=$DIR_cur/3_outRnaTri  #path to storage file of tri.fasta
DIR_sif=$DIR_cur/out_Fasta 
DIR_dec=$DIR_cur/out_Outfmt 
DIR_tdc=$DIR_cur/out_Faa #path to Trans.DX.L.pep.faa
DIR_Trc=$DIR_cur/out_TranCont #path to drawings

#-p if there is already a file that has the same name, then it won't creat a same one.
mkdir -p $DIR_sif
mkdir -p $DIR_dec
mkdir -p $DIR_tdc
mkdir -p $DIR_Trc

InputFasta=$(ls ./3_outRnaTri/*_trinity.Trinity.fasta | sed -e "s/.fasta//g;s/.\/3_outRnaTri\///g") # all input fasta files  #output=speciesname_rna.trinity.Trinity, such as Xanthostigma_gobicola_SRR1811918.rna_trinity.Trinity
NNN=$(echo $InputFasta | wc -w)  #calculate the number

EcdysFasta=$(ls ./3_outRnaTri/*_trinity.Trinity.fasta | sed -e "s/.fasta//g;s/.\/3_outRnaTri\///g" | grep -E "^Xenos|^Mengenilla|^Stylops") # exceptions: Species of Strepsiptera should be decontaminated by searching Holometabola genes, such as Mengenilla_moldrzyki_SRR619393.rna_trinity.Trinity

### Functions
runDeVect(){

    local rnatri=$(echo $1 | sed -e "s#$DIR_asm\/##g;s#.fasta##g")

    if [[ $EcdysFasta =~ $rnatri ]]; then 
        local DeConProg=/user/work/tt22567/Neuroptera/4_trifasta/DeContaminRNA-tx.py
    else
        local DeConProg=/user/work/tt22567/Neuroptera/4_trifasta/DeContaminRNA-rk.py
    fi
    
    echo -e "\n########################\nRemoval of vector sequences...\n$(date)\n########################\n" 
    echo -e "\n========================\nScreening $rnatri.fasta (round 1)...\n========================\n" 

    cp $1 $DIR_sif/$rnatri.Sifted.fasta
        
    echo -e "\nStart screening $rnatri.fasta.\n" 
    echo -e "\n$(blastn -version)\n" 
    blastn -db $UniVecDB -query $DIR_sif/$rnatri.Sifted.fasta -out $DIR_dec/$rnatri.decont -outfmt "6 qseqid sseqid qstart qend  score evalue pident ppos staxids" -task blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true -evalue 1e-5 -searchsp 1750000000000 && 
    echo -e "\nFinished screening $rnatri.fasta.\n$(date)\n" 
    
    local counter=1

    while [ $(ls -l $DIR_dec/$rnatri.decont | awk '{ print $5 }') -gt 0 -a $counter -le 10 ] 
    do
        echo -e "\n========================\nScreening $rnatri.fasta (round $[$counter+1])...\n$(date)\n========================\n" 
        echo -e "\nRemoving vector sequences from $rnatri.fasta...\n" 
        $DeConProg --fasta $DIR_sif/$rnatri.Sifted.fasta --csv $DIR_dec/$rnatri.decont --out $DIR_sif --keep && 
        mv -f $DIR_sif/$rnatri.Sifted.Sifted.fasta $DIR_sif/$rnatri.Sifted.fasta &&         
        echo -e "\nFinished. Clean assembly saved to $DIR_sif/$rnatri.Sifted.fasta...\n$(date)\n" 

        echo -e "\n========================\nScreening $DIR_sif/$rnatri.Sifted.fasta to confirm no vector sequences...\n$(date)\n========================\n" 
        rm $DIR_dec/$rnatri.decont && 
        blastn -db $UniVecDB -query $DIR_sif/$rnatri.Sifted.fasta -out $DIR_dec/$rnatri.decont -outfmt "6 qseqid sseqid qstart qend score evalue pident ppos staxids" -task blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true -evalue 1e-5 -searchsp 1750000000000 && 
        echo -e "\nFinished screening $DIR_sif/$rnatri.Sifted.fasta.\n$(date)\n" 
        
        ((counter++))
    done
    
    echo -e "\n########################\nCleaning Vectors from $rnatri is accomplished after $counter rounds.\n$(date)\n########################\n" 
}

runFindAlien(){

    local rnatri=$(echo $1 | sed -e "s#$DIR_sif\/##g;s#.Sifted.fasta##g")

    echo -e "\n########################\nSearching contaminated sequences...\n$(date)\n########################\n" 

    if [[ ! -f $DIR_dec/$rnatri.decontX ]] || [[ $(ls -l $DIR_dec/$rnatri.decontX | awk '{ print $5 }') -eq 0 ]]; then 
        echo -e "\n========================\nSearching the unwanted transcripts by Diamond and non-redundant protein database...\n$(date)\n========================\n" && 
        diamond blastx --taxonmap $taxonomap_path --threads $DNP --db $NRDiaDB --query $1 --out $DIR_dec/$rnatri.decontX --outfmt 6 qseqid sseqid qstart qend score evalue pident ppos staxids --max-target-seqs 10 --evalue 1e-5 && 
        echo -e "\nFinished screening $DIR_sif/$rnatri.Sifted.fasta.\n$(date)\n" 
    else
        echo -e "\n########################\nScreening $rnatri.Sifted.fasta has already been accomplished, no need to do it again.\n$(date)\n########################\n" 
    fi
}

runDeCont(){

    local rnatri=$(echo $1 | sed -e "s#$DIR_sif\/##g;s#.Sifted.fasta##g")

    if [[ $EcdysFasta =~ $rnatri ]]; then 
        local DeConProg=/user/work/tt22567/Neuroptera/4_trifasta/DeContaminRNA-tx.py
        local DeConVad=Holometabola #this is for the specific groups, strepsiptera, will decont in a broader taxonomic level.
    else
        local DeConProg=/user/work/tt22567/Neuroptera/4_trifasta/DeContaminRNA-rk.py
        local DeConVad=order   #lower case; ingroup Neuropterida--superorder; outgroup most Holometabolous orders--order, except Strepsiptera--Holometabola; how to decide the taxonomy level depends on the conditions, if the group you want to decontaminate has enough genes on the NCBI, then you should just search and decont it with a strict taxonomic level. However, if there are no much genes on NCBI, then its better to decont the quarry sequence in a broader taxonomic level. 
    fi
    
    echo -e "\n########################\nRemoval of contaminated sequences...\n$(date)\n########################\n" 

    echo -e "\n========================\nRemoving alien transcripts from the assembly $rnatri ...\n$(date)\n========================\n" && 
    $DeConProg --fasta $DIR_sif/$rnatri.Sifted.fasta --csv $DIR_dec/$rnatri.decontX --out $DIR_sif --sensitive 0 --valid $DeConVad --db $TaxDumpDB &&
    echo -e "\nFinished. Clean assembly saved to $DIR_sif/$rnatri.Sifted.DeContX.fasta...\n$(date)\n" 

    echo -e "\n########################\nCleaning alien transcripts from $rnatri is accomplished...\n$(date)\n########################\n" 
}

chmod -R 777 /user/work/tt22567/Neuroptera/4_trifasta/*.py 

#source /sw/languages/anaconda/anaconda.3.9-2021.12-bioconda/bin/activate   /user/home/tt22567/.conda/envs/Lacewings  #path to the env where blast installed
source /user/work/tt22567/Miniconda3/ENTER/bin/activate      /user/work/tt22567/Miniconda3/ENTER/envs/Lacewings

### Remove vector fragments
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

# vector-removing processes
for rnatri in $InputFasta
do
    read -u9
    {

    # remove vector fragments
    if [[ ! -f $DIR_dec/$rnatri.decont ]] || [[ $(ls -l $DIR_dec/$rnatri.decont | awk '{ print $5 }') -gt 0 ]]; then 
        runDeVect $DIR_asm/$rnatri.fasta >> $DIR_dec/$rnatri.log 2>&1 
    else
        echo -e "\n########################\nCleaning Vectors from $rnatri has already been accomplished, no need to do it again.\n$(date)\n########################\n" >> $DIR_dec/$rnatri.log 
    fi 

    echo >&9 
    } & 
done
wait

### Finding decontamination
# Multithread preparation
tmp_fifofile="/tmp/$$.fifo"
trap "exec 9>&-;exec 9<&-;exit 0" 2
mkfifo $tmp_fifofile
exec 9<>$tmp_fifofile 
rm $tmp_fifofile

for ((i=1;i<=$NT;i++))
do
    echo >&9
done

# Finding decontamination process
for rnatri in $InputFasta
do
    read -u9
    {

    # remove transcripts that are likely to be from unwanted taxa
    while [[ ! -f $DIR_dec/$rnatri.decontX ]] || [[ $(ls -l $DIR_dec/$rnatri.decontX | awk '{ print $5 }') -eq 0 ]] 
    do
        runFindAlien $DIR_sif/$rnatri.Sifted.fasta >> $DIR_dec/$rnatri.log 2>&1 
    done

    echo >&9 
    } & 
done
wait

#source /sw/languages/anaconda/anaconda.3.9-2021.12-bioconda/bin/deactivate  /user/home/tt22567/.conda/envs/Lacewings

#source /sw/languages/anaconda/anaconda.3.9-2021.12-bioconda/bin/activate    /user/home/tt22567/.conda/envs/Lacewings #transdecoder is installed in Lacewings, together with trinity

### Decontamination, Transdecode decontaminated assemblies, and Display transcript origins
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

# 
SiftFasta=$(ls $DIR_dec/*_trinity.Trinity.decontX | sed -e "s#.decontX##g;s#$DIR_dec\/##g")

for rnatri in $SiftFasta
do
    read -u9
    {

    # remove transcripts that are likely to be from unwanted taxa
    while [[ ! -f $DIR_sif/$rnatri.Sifted.DeContX.fasta ]] || [[ $(ls -l $DIR_sif/$rnatri.Sifted.DeContX.fasta | awk '{ print $5 }') -lt 100 ]] 
    do
        runDeCont $DIR_sif/$rnatri.Sifted.fasta >> $DIR_dec/$rnatri.log 2>&1 
    done

    # Display transcript origins
    if [[ ! -f $DIR_Trc/$rnatri.png ]]; then
        echo -e "\n========================\nDrawing $nrna ...\n$(date)\n========================\n" >> $DIR_dec/$rnatri.log && 
        $SeqOrinPiePlot --infile $DIR_dec/$rnatri.decontX --out $DIR_Trc --db $TaxDumpDB >> $DIR_dec/$rnatri.log 2>&1 && 
        echo -e "\n========================\nFinished $nrna ...\n$(date)\n========================\n" >> $DIR_dec/$rnatri.log         
    else
        echo -e "\n########################\nDrawing has already been accomplished, no need to do it again.\n$(date)\n########################\n" >> $DIR_dec/$rnatri.log 
    fi
    
    # Transdecode decontaminated assemblies
    if [[ ! -f $DIR_tdc/$(echo $rnatri | sed -e "s/_trinity.Trinity//g").Trans.DX.L.pep.faa ]]; then
        echo "\n========================\nSelecting the longest isoform.\n$(date)\n========================\n" >> $DIR_dec/$rnatri.log && 
        perl /user/work/tt22567/Miniconda3/ENTER/envs/Lacewings/opt/trinity-2.8.5/util/misc/get_longest_isoform_seq_per_trinity_gene.pl $DIR_sif/$rnatri.Sifted.DeContX.fasta > $DIR_sif/$rnatri.Sifted.DeContX.L.fasta && 
        
        echo "\n========================\nTransDecoding $rnatri (Sifted+DeContX) LongestORF starts.\n$(date)\n========================\n" >> $DIR_dec/$rnatri.log && 
        TransDecoder.LongOrfs -t $DIR_sif/$rnatri.Sifted.DeContX.L.fasta >> $DIR_dec/$rnatri.log 2>&1 && 
        TransDecoder.Predict -t $DIR_sif/$rnatri.Sifted.DeContX.L.fasta >> $DIR_dec/$rnatri.log 2>&1 && 
        mv $DIR_cur/"$rnatri"*.transdecoder.* $DIR_tdc && 
        mv $DIR_tdc/"$rnatri"*.transdecoder.pep $DIR_tdc/$(echo $rnatri | sed -e "s/_trinity.Trinity//g").Trans.DX.L.pep.faa && 
        rm -rf $DIR_cur/$rnatri.Sifted.DeContX.L.fasta.transdecoder_dir* pipeliner.*.cmds && 
        echo "\n========================\nTransDecoding $rnatri (Sifted+DeContX) LongestORF finished.\n$(date)\n========================\n" >> $DIR_dec/$rnatri.log 
    else
        echo -e "\n########################\nTransdecoding has already been accomplished, no need to do it again.\n$(date)\n########################\n" >> $DIR_dec/$rnatri.log 
    fi

    #Cdhit removing duplications
	if [[ ! -f $DIR_tdc/$(echo $rnatri | sed -e "s/_trinity.Trinity//g").Trans.DX.L.pep.cdhit.faa ]]; then
        echo "\n========================\nremoving the duplications of $DIR_tdc/$(echo $rnatri | sed -e "s/_trinity.Trinity//g").Trans.DX.L.pep.faa starts.\n$(date)\n========================\n" >> $DIR_dec/$rnatri.log && 
        /user/work/tt22567/Miniconda3/ENTER/envs/Lacewings/bin/cd-hit-est -i $DIR_tdc/$(echo $rnatri | sed -e "s/_trinity.Trinity//g").Trans.DX.L.pep.faa -o $DIR_tdc/$(echo $rnatri | sed -e "s/_trinity.Trinity//g").Trans.DX.L.pep.cdhit.faa -c 0.98 -n 10 &&
		rm $DIR_tdc/$(echo $rnatri | sed -e "s/_trinity.Trinity//g").Trans.DX.L.pep.cdhit.faa.clstr &&
		echo "\n========================\nremoving the duplications of $DIR_tdc/$(echo $rnatri | sed -e "s/_trinity.Trinity//g").Trans.DX.L.pep.faa is done.\n$(date)\n========================\n" >> $DIR_dec/$rnatri.log
	else
        echo "\n========================\nremoving the duplications of $DIR_tdc/$(echo $rnatri | sed -e "s/_trinity.Trinity//g").Trans.DX.L.pep.faa is finished, no need to do it again.\n$(date)\n========================\n" >> $DIR_dec/$rnatri.log
	fi	
	
	echo -e "\n########################\nAll data processing for $nrna (including Decontamination and removing duplications) has been done.\n$(date)\n########################\n" >> $DIR_dec/$rnatri.log 
    
	echo >&9 
    } & 
done
wait


