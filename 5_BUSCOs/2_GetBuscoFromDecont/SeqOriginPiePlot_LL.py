#!/sw/languages/anaconda/anaconda.3.9-2021.12-bioconda/bin/python3.9

## Pie charts of phylum based transcripts classification
## libs
import matplotlib.pyplot as plt
import pandas as pd
import argparse, os
# import seaborn as sb
from math import ceil
from Bio import Entrez
Entrez.email ='lianges.luex@gmail.com' # please feel free to get in touch with me, when you have any questions or problems. 

def qsmaker(i):
    qs = pd.Series(str(Dec.iloc[i]['staxids']).split(';')).apply(lambda x: x if x.isdigit() else None).dropna()
    qs.index = [Dec.index[i]] * len(qs)
    return qs

def idcorr(oldid):
    newid = TaxMerge.loc[oldid]['validid'] if oldid in TaxMerge.index else oldid
    return newid

def getHigherLineage(taxid, rank=''):
    LineageID = [taxid]
    while LineageID[-1]!=1:
        LineageID.append(TaxNodes.loc[LineageID[-1]]['parent'])
    LineageID=pd.Series(LineageID, index=TaxNodes.loc[LineageID]['rank'])
    if rank=='':
        return TaxNames.loc[LineageID]['nametxt']
    elif rank in LineageID.index:
        return TaxNames.loc[LineageID.iloc[LineageID.index.tolist().index(rank):]]['nametxt']
    else:
        return None

def qclass(GenName, rank, remote=False):
    if remote:
        qtaxrefId = Entrez.read(Entrez.esearch(db='taxonomy', term=GenName))['IdList'][0]
        qtaxrefLin = Entrez.read(Entrez.efetch(db='taxonomy', id=qtaxrefId))[0]['LineageEx']
        qreftax = [rk['ScientificName'] for rk in qtaxrefLin if rk['Rank'].lower() == rank][0] 
    else:
        qtaxrefId = TaxNames[TaxNames['nametxt']==GenName].index.values[0] 
        qreftax = getHigherLineage(qtaxrefId, rank=rank).iloc[0]
    return qreftax

def freq(dec, target, stepmax=100, remote=False):
    qs = list(map(qsmaker, range(len(dec.index))))
    qsid = pd.concat(qs)
    qsid.dropna(inplace=True)
    if remote:
        SidUni = qsid.drop_duplicates().astype(str)
        stp1 = [x*stepmax for x in range(ceil(len(SidUni)/stepmax))]
        stp2 = stp1[1:] + [len(SidUni)]
        rk = []
        for i in range(len(stp1)):
            handle = Entrez.efetch(db="taxonomy", id=SidUni[stp1[i]:stp2[i]].tolist())
            TaxInfo = Entrez.read(handle)
            handle.close()
            rk.extend([[t['AkaTaxIds'][0] if 'AkaTaxIds' in t.keys() else t['TaxId'], t['Lineage']] for t in TaxInfo])
        rk = pd.DataFrame(rk).drop_duplicates() # 2-column dataframe with no index
        rk.index = rk.iloc[:,0]
        rk = pd.Series(rk.iloc[:,1].apply(lambda x: x.split('; '))) # series with index
    else:
        qsid = qsid.astype(int)
        qsid = qsid.apply(idcorr)
        qsid = qsid[~qsid.isin(TaxDel)]
        SidUni = qsid.drop_duplicates()
        rk = SidUni.apply(lambda x: getHigherLineage(x).tolist()) # series with index
        rk.index = SidUni.values 
        rk.dropna(inplace=True)
        rk.drop_duplicates(inplace=True)    
    qsid = qsid[qsid.isin(rk.index)].astype(int)
    TAGclade=pd.Series([target,'Holometabola','Insecta','Arthropoda', 'Ecdysozoa', 'Spiralia', 'Protostomia', 'Deuterostomia', 'Xenacoelomorpha', 'Bilateria', 'Cnidaria', 'Ctenophora', 'Placozoa', 'Porifera', 'Viridiplantae', 'Fungi', 'Bacteria', 'Archaea', 'Viruses'], index=[target,'Holom+','Insec+','Arthr+','Ecdys+','Spira','Proto+','Deute','Xenac','Bilat+','Cnida','Cteno','Placo','Porif','Virid','Fungi','Bacte','Archa','Virus'])
    def tagrk(x):
        tagin = TAGclade.isin(rk[x])
        if any(tagin):
            return TAGclade.index[tagin][0]
        else:
            return 'others' 
    qstag = qsid.apply(tagrk)
    mx = lambda x: max(x, key=x.tolist().count)
    QSidvote = qstag.groupby([qstag.index]).apply(mx) # assign the Qid to the taxon that most of its Sids belong to.
    xx = QSidvote.groupby(QSidvote).count() # summarize how many Qids of each taxon.
    xx.dropna(inplace=True)
    res = pd.DataFrame(columns=['Count','Percentage','Explode'])
    res['Count'] = xx
    res.sort_values(['Count'], ascending=False, inplace=True)
    sumCount=res['Count'].sum()
    res['Percentage'] = res['Count'].div(sumCount) * 100
    res['Explode'] = [0.1 if i == target else 0 for i in res.index]
    return res

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--infile', '-i',
                        type=str,
                        help='DecContamination file ends with ".DeConX"')
    parser.add_argument('--out', '-o',
                        type=str,
                        help='Output directory')
    parser.add_argument('--db', '-d',
                        type=str,
                        help='Designate the database location of NCBI Taxonomy: "remote", search online Taxonomy Database; or provide the directory that contains all the uncompressed files of standlone version of Taxonomy Database (e.g. /YOUR/OWN/PATH/TO/taxdump).',
                        default="remote")
    args = parser.parse_args()
    # reading data
    decontfile = args.infile.replace('\\','/')
    #decontheader = ['', 'sseqid', 'qstart', 'qend', 'staxids']
    decontheader = ['', 'sseqid', 'qstart', 'qend', 'score', 'evalue', 'pident', 'ppos', 'staxids']
    # preprocessing
    Dec = pd.read_csv(decontfile, index_col=0, names=decontheader, sep="\t")
    gnm = decontfile.split('/')[-1].split('_')
    gnsp = (gnm[0].isupper() and gnm[1:4] or gnm[0:3])
    if args.db.lower()=="remote": # search online NCBI Taxonomy database
        ReMoTe = True
    elif all([os.path.exists(args.db + '/nodes.dmp'), os.path.exists(args.db + '/names.dmp'), os.path.exists(args.db + '/merged.dmp'), os.path.exists(args.db + '/delnodes.dmp')]): # search the local NCBI database, or custom database of same format
        ReMoTe = False
        TaxNodes = pd.read_csv(args.db + '/nodes.dmp', names=['taxid','parent','rank'], index_col=0, header=None, usecols=[0,2,4], sep="\t")
        TaxNames = pd.read_csv(args.db + '/names.dmp', names=['taxid','nametxt','uniquename','nameclass'], index_col=0, header=None, usecols=[0,2,4,6], sep="\t")
        TaxNames = TaxNames[TaxNames['nameclass']=='scientific name']
        TaxMerge = pd.read_csv(args.db + '/merged.dmp', names=['oldid','validid'], index_col=0, header=None, usecols=[0,2], sep="\t")
        TaxDel = pd.read_csv(args.db + '/delnodes.dmp', index_col=None, header=None, usecols=[0], sep="\t").iloc[:,0]
    else:
        print('Wrong input for the argument of "--db".\n')
        exit()
    # counting
    RankRef = qclass(gnsp[0], rank='phylum', remote=ReMoTe)
    RankFreq = freq(Dec, RankRef, stepmax=50, remote=ReMoTe)
    # drawing  # n126-141
    tit = RankRef + '+' + ' '.join(gnsp)
    print('Drawing' + tit + '\n') 
    taxcol = pd.Series(['darkgray','mediumturquoise','lavender','cyan','royalblue','brown','coral','red','darkgreen','hotpink','purple','orange','yellow','fuchsia','goldenrod','khaki','olive','orchid','salmon','wheat'], index=[RankRef,'Holom+','Insec+','Arthr','Ecdys+','Spira','Proto+','Deute','Xenac','Bilat+','Cnida','Cteno','Placo','Porif','Virid','Fungi','Bacte','Archa','Virus','others'])
    # RankFreq = RankFreq[:-1]
    col = taxcol[RankFreq.index].tolist()
    plt.figure(facecolor='snow')
    plt.axes(aspect="equal") # circle
    plt.xlim(0,8)
    plt.ylim(0,8)
    # plt.subplot(3, 4, fl, frameon = False)
    plt.pie(x=RankFreq['Count'], colors=col, labels=RankFreq.index, labeldistance=1.1, radius=1, explode=RankFreq['Explode'], frame=0)
    plt.title(tit.replace('+','\n'))
    plt.xticks(())
    plt.yticks(())
    plt.savefig(args.out + '/' + decontfile.split('/')[-1].replace('.decontX','.png'))
    RankFreq = RankFreq.append(pd.Series(RankFreq.sum(axis=0), name='Sum'))
    RankFreq['Explode'] = RankFreq['Explode'].apply(lambda x: '*' if x > 0 else None)
    RankFreq.at['Sum','Explode'] = None
    RankFreq.to_excel(args.out + '/' + decontfile.split('/')[-1].replace('.decontX','.xlsx'))
    