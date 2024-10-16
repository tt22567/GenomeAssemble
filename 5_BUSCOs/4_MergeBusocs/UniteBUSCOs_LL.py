#!/user/work/tt22567/Miniconda3/ENTER/envs/Lacewings/bin/python

import os, re
from Bio import SeqIO
import matplotlib.pyplot as plt
from matplotlib import gridspec as gs
import numpy as np
import pandas as pd
import argparse

def TaxaDirMap(dir):  #DirList
    #taxa = dir.str.split(r'(_SRS|_ERP|_ERS|_GCA|_GCF|_UP|.MG|_PRJEB|_NC|_NEB|_N1|_N2)', expand=True, n=1).iloc[:,0].drop_duplicates()  
    taxa = dir.str.split(r'(_SRR|_GCF|_ERR|_GCA)', expand=True, n=1).iloc[:,0].drop_duplicates()  #only left species names
    taxa.reset_index(drop=True, inplace=True)
    tdlist = taxa.apply(lambda x: dir[dir.str.find(x)==0])  #map of taxa and DirList；
    return tdlist  #TaxaDir

def DirFastaMap(dir, prefix='.'):  #DirList, prefix=InDir
    dflist = dir.apply(lambda x: pd.DataFrame([[x, i] for i in os.listdir(prefix+'/'+x) if i[-4:]=='.faa']))  # For each faa file, build a dataframe
    dflist.index = dir.values
    return dflist  #DirFasta

def TaxaFastaTable(dirfasta, prefix='.'):  #DirFasta, prefix=InDir
    df = pd.concat(dirfasta.tolist(), ignore_index=True)
    #taxa = df.iloc[:,0].str.split(r'(_SRS|_ERP|_ERS|_GCA|_GCF|_UP|.MG|_PRJEB|_NC|_NEB|_N1|_N2)', expand=True, n=1).iloc[:,0]
    taxa = df.iloc[:,0].str.split(r'(_SRR|_GCF|_ERR|_GCA)', expand=True, n=1).iloc[:,0]
    gene = df.iloc[:,1]
    fastapath = df.apply(lambda x: prefix + '/' + x[0] + '/' + x[1], axis=1)
    return pd.DataFrame({'Taxon':taxa, 'Gene':gene, 'Path':fastapath}, index=range(len(taxa)))  #TaxaFasta

def selectLongest(file, id):
    seqlen = 0
    for i in file:
        fa = SeqIO.read(i, "fasta")
        if len(fa) > seqlen:
            fasta = fa
            fasta.id = id
            fasta.name = ''
            fasta.description = ''
            seqlen = len(fasta)
    return fasta

def writeSeq(gene):
    fp = TaxaGeneUni.loc[TaxaGeneUni.loc[:,'Gene']==gene,:]
    genefasta = fp.apply(lambda x: selectLongest(file=TaxaFasta.loc[(TaxaFasta.loc[:,'Taxon']==x['Taxon']) & (TaxaFasta.loc[:,'Gene']==gene),'Path'].tolist(), id=x['Taxon']), axis=1)
    #SeqIO.write(genefasta, OutDir+'/BSC_'+gene.replace('at6231',''), "fasta")
    SeqIO.write(genefasta, OutDir+'/BSC_'+gene.replace('at33392',''), "fasta")
    res = fp.copy()
    res.loc[:,'Max.Length'] = genefasta.apply(len)
    return res

def TaxaGeneMatrix(genestattable):
    genestattable.index=genestattable.loc[:,'Taxon']
    genuni = genestattable.loc[:,'Gene'].drop_duplicates()
    tgmat = list(map(lambda x: genestattable.loc[genestattable.loc[:,'Gene']==x, 'Max.Length'], genuni))
    tgmat = pd.concat(tgmat, sort=True, axis=1)
    tgmat.columns = genuni
    tgmat.sort_index(axis=1, inplace=True)
    totlen = tgmat.apply(lambda x: sum(x.dropna()), axis=1)
    totlen.name='Total.Length'
    tgmat = pd.concat([totlen, tgmat], sort=True, axis=1)
    tgmat.reset_index(inplace=True)    
    return tgmat

def HeatMap(taxastat, outdir, dpi=300, nfsize=50, tickstep=10):
    ts_log = taxastat.iloc[:-2,2:].apply(np.log).fillna(0)
    ts_bin = taxastat.iloc[:-2,2:].isnull()
    xlabpos = np.arange(len(ts_log.columns), step=tickstep*3)
    ylabpos = np.arange(len(ts_log.index), step=tickstep)
    # plt.rc('xtick',labelsize=3)
    # plt.rc('ytick',labelsize=3)
    fig = plt.figure(figsize=(3*nfsize,2*nfsize)) # width, height
    gsp = gs.GridSpec(3, 1, height_ratios=[1, 1, 1])
    axes0 = plt.subplot(gsp[0])
    htp0 = axes0.imshow(ts_log, cmap=plt.cm.jet, aspect='auto')
    axes0.set_title('Log-length of gene vs. taxon', fontdict={'fontsize':25})
    axes0.set_xticks([])
    axes0.set_yticks(ylabpos)
    axes0.set_yticklabels(labels=ts_log.index[ylabpos].tolist(), fontdict={'fontsize':3})
    axes0.set_ylabel('Taxon', fontdict={'fontsize':18})
    axes1 = plt.subplot(gsp[1])
    htp1 = axes1.imshow(ts_bin, cmap=plt.cm.gray, aspect='auto')
    axes1.set_title('Presence of gene vs. taxon', fontdict={'fontsize':25})
    axes1.set_xticks([])
    axes1.set_yticks(ylabpos)
    axes1.set_yticklabels(labels=ts_bin.index[ylabpos].tolist(), fontdict={'fontsize':3})
    axes1.set_ylabel('Taxon', fontdict={'fontsize':18})
    axes2 = plt.subplot(gsp[2])
    axes2.bar(taxastat.iloc[-1,2:].index, taxastat.iloc[-1,2:], width=1, color='orange')
    axes2.margins(x=0)
    axes2.set_xlabel('BUSCO gene', fontdict={'fontsize':18})
    axes2.set_ylabel('Missing', fontdict={'fontsize':18})
    axes2.set_title('Count of missing genes', fontdict={'fontsize':25})
    axes2.set_xticks(xlabpos)
    axes2.set_xticklabels(labels=ts_bin.columns[xlabpos].tolist(), fontdict={'fontsize':3})
    plt.setp(axes2.get_xticklabels(), rotation=270)
    fig.tight_layout()
    plt.subplots_adjust(hspace=0.2)
    fig.savefig(outdir+'/'+'TaxaGeneCount.png', dpi=dpi)
    # fig.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', '-i',
                        type=str,
                        help='Sequences requiring masking')
    parser.add_argument('--outdir', '-o',
                        type=str,
                        help='The output path of the masked sequence',
                        default=".")
    parser.add_argument('--sortmissing', '-s',
                        action='store_true',
                        help='The final output TaxaStat table and the TaxaGeneCount figure are sorted alphabetically by indeces (taxon and gene names), but if this is specified, they will be sorted by the number of missing genes ascendingly.')
    args = parser.parse_args()

    # initializing and collecting information
    InDir = args.indir.replace('\\','/')    #InDir = "/user/work/tt22567/Neuroptera/5_deconfasta/DIR_gene_test"
    OutDir = args.outdir.replace('\\','/')  #OutDir = "/user/work/tt22567/Neuroptera/5_deconfasta/DIR_Unite"
    CurDir = os.getcwd().replace('\\','/')  #CurDir = "/user/work/tt22567/Neuroptera/5_deconfasta"
    if not os.path.isabs(InDir):  #to detect whether InDir is an absolute path or not
        InDir = CurDir + '/' + InDir.replace('./','')
    if not os.path.isabs(OutDir):
        OutDir = CurDir + '/' + OutDir.replace('./','')
    if not os.path.exists(OutDir):
        os.makedirs(OutDir, mode=0o777, exist_ok=True)    
    DirList = os.listdir(InDir)   #os.listdir(InDir) will return the list of files in the file of InDir
    DirList = pd.Series([i for i in DirList if os.path.isdir(InDir+'/'+i)])  #os.path.isdir () is a built-in Python function that is used to check whether the specified path is an existing directory or not
    TaxaDir = TaxaDirMap(DirList)  #TaxaDirMap()  dir=DirList  TaxaDir=tdlist
    DirFasta = DirFastaMap(DirList, prefix=InDir)  #DirFastaMap() dir=DirList, prefix='.'; DirFasta=dflist
    TaxaFasta = TaxaFastaTable(DirFasta, prefix=InDir) #dataframe including three columes:taxon,faa,path
    GeneCount = TaxaFasta.groupby(['Taxon','Gene']).count().reset_index()  #发生了什么？taxon少了
    TaxaGeneUni = TaxaFasta[['Taxon','Gene']].drop_duplicates(inplace=False, ignore_index=True)  #
    GeneUni = TaxaGeneUni['Gene'].drop_duplicates(inplace=False)

    # generating gene-fasta files and summary dataframe    
    GeneLenght = list(map(writeSeq, GeneUni))
    GeneLenght = pd.concat(GeneLenght, ignore_index=True)
    GeneStat = pd.merge(GeneCount, GeneLenght, on=['Taxon', 'Gene'], sort=True)
    #GeneStat['Gene'] = GeneStat['Gene'].apply(lambda x: 'BSC_'+x.replace('at6231.faa',''))
    GeneStat['Gene'] = GeneStat['Gene'].apply(lambda x: 'BSC_'+x.replace('at33392.faa',''))
    GeneStat.columns = ['Taxon', 'Gene', 'Source', 'Max.Length']
    del(GeneLenght)

    # generating taxa-gene summary
    TaxaStat1 = TaxaGeneUni.groupby(['Taxon']).count().reset_index()
    TaxaStat1.sort_values(by='Taxon', ascending=True, inplace=True)
    TaxaStat1.columns = ['Taxon', 'Total.Genes']
    TaxaStat2 = TaxaGeneMatrix(GeneStat)
    TaxaStat = pd.merge(TaxaStat1, TaxaStat2, on=['Taxon'], sort=True)
    TaxaStat.index = TaxaStat.loc[:, 'Taxon']
    TaxaStat = TaxaStat.iloc[:, 1:]
    TSavg = pd.Series(TaxaStat.median(axis=0), name='MEDIAN').astype(int)
    TSmissing = pd.Series(TaxaStat.isnull().sum(axis=0), name='MISSING').astype(int)
    TSmissing[0] = sum(TSmissing[2:])
    TaxaStat.loc['AVERAGE',:] = TSavg
    TaxaStat.loc['MISSING',:] = TSmissing
    if args.sortmissing:
        colindex = TaxaStat.columns[:2].tolist() + TaxaStat.iloc[-1, 2:].sort_values(ascending=True).index.tolist()
        rowindex = TaxaStat.iloc[:-2, 0].sort_values(ascending=False).index.tolist() + TaxaStat.index[-2:].tolist()
        TaxaStat = TaxaStat.loc[rowindex, colindex]
    del(TaxaStat1, TaxaStat2, TSavg, TSmissing)

    # writing uniting summaries
    GeneStat.to_excel(OutDir+'/'+'UnitedGeneStat.xlsx')
    TaxaStat.to_excel(OutDir+'/'+'TaxaGeneCount.xlsx')
    print('In summary, the UniteBUSCOs processing has treated {DL} protein prediction files and extracted {NG} genes for {NT} species. In the output {NT} * {NG} data matrix, there are {NM} missing loci ({PM}%).\n'.format(DL=str(len(DirList)), NG=str(len(GeneUni)), NT=str(len(TaxaStat.index)-2), NM=str(TaxaStat.isna().sum().sum()), PM=str(round(TaxaStat.isna().sum().sum()*100/(len(TaxaStat.index)-2)/len(GeneUni),2))))
    # TaxaStat = pd.read_excel('E:/Nema_backup/data/busco/united/TaxaGeneCount.xlsx', index_col=0, header=0)

    # plotting heat map
    HeatMap(taxastat=TaxaStat, outdir=OutDir, nfsize=30, tickstep=1)
