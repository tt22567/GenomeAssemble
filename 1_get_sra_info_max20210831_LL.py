from Bio import Entrez
from datetime import datetime
from bs4 import BeautifulSoup
from pandas import DataFrame
from math import ceil
import pandas as pd
import argparse
Entrez.email ='A.N.Other@example.com'


def ProcBar(percent, StartStr='', EndStr='', TotalLength=50):
    bar = ''.join(["\033[0;42m%s\033[0m"%' '] * int(percent * TotalLength)) + ''
    bar = '\r' + StartStr + bar.ljust(TotalLength) + ' {:0>4.1f}%'.format(percent*100) + EndStr
    print(bar, end='', flush=True)

def get_SRAinfo_df(id):
    Toatal_length=len(id)
    stp = [x * 10000 for x in range(ceil(Toatal_length / 10000))]
    df_SRAinfo = DataFrame(columns=['Accession','SpeciesTaxid','Organism', 'Run_num', 'Library_strategy',
                                    'Library_source', 'Total_bases', 'Download_address', 'Org'])
    for i in stp:
        print(' Get data from NCBI start:'+str(i)+" - "+str(i+10000))
        raw_info_id = Entrez.efetch(db="sra",retstart=i, id=id,retmode="xml")
        print(' Load data information...')
        record_id = raw_info_id.read()
        list_info = BeautifulSoup(record_id, "lxml")
        print(' Analyze data information')
        num = 1
        for id_info in list_info.find_all("experiment_package"):
                length=len(list_info.find_all("experiment_package"))
                ProcBar((num) / length, StartStr='>>>', EndStr='|' + str(length) + ' records')
                num=num+1
                #library_layout
                library_layout=str(id_info.find_all("library_layout")[0].contents[0])
                if library_layout.endswith('</paired>'):
                    # download_address
                    taxid = int(id_info.find_all('taxon_id')[0].string)
                    for srainfo in id_info.find_all("alternatives"):
                        if srainfo['org'] == 'NCBI':
                            org = srainfo['org']
                            download_address = srainfo['url']
                            #accession
                            accession=str(id_info.find_all("member")[0]['accession'])
                            #organism
                            organism=id_info.find_all("member")[0]['organism']
                            #data_type
                            library_strategy = id_info.find_all("library_strategy")[0].string
                            #data_source
                            library_source = id_info.find_all("library_source")[0].string
                            #run_num
                            run_num=id_info.find_all("run")[0]['accession']
                            #total_bases
                            total_bases=int(id_info.find_all("run")[0]['total_bases'])
                    id_sum = [accession,taxid,organism,run_num,library_strategy,library_source,total_bases,download_address,org]
                    #return {'accession':accession,'organism':organism,'run_num':run_num,'szie(k)':size,'download_address':download_address}
                    df_SRAinfo.loc[accession] = id_sum
        #df_SRAinfo.to_excel('0-rawdata_{0}_{1}.xlsx'.format('SRA',str(i)))
        print(' ')
    print('Finally, the amount of paired data (NCBI) is '+str(len(df_SRAinfo)))
    #df_SRAinfo.to_excel('0-rawdata_SRA.xlsx')
    return df_SRAinfo

def get_taxinfo_df(taxid_list):
    taxinfo_tab=DataFrame(columns=['tax_id','Phylum', 'SubPhylum', 'Class', 'Order', 'SubOrder', 'InfraOrder',
                                   'SuperFamily', 'Family', 'SubFamily','ScientificName'])
    Toatal_length=len(taxid_list)
    stp = [x * 10000 for x in range(ceil(Toatal_length / 10000))]
    for i in stp:
        print(' Get data from NCBI: ' + str(i) + " - " + str(i + 10000))
        TaxInfo=Entrez.read(Entrez.efetch(db="taxonomy",retstart=i, id=taxid_list))
        print(' Analyze data information...')
        for num_id in range(len(TaxInfo)):
            ProcBar((num_id) / len(TaxInfo), StartStr='>>>', EndStr='|' + str(len(TaxInfo)) + ' records')
            taxinfo_dict = {}
            tax_id=int(TaxInfo[num_id]['TaxId'])
            taxinfo_dict['tax_id'] = tax_id
            taxinfo_dict['ScientificName'] = TaxInfo[num_id]['ScientificName']
            for rank in ['Phylum', 'SubPhylum', 'Class', 'Order', 'SubOrder', 'InfraOrder', 'SuperFamily', 'Family', 'SubFamily']:
                try:
                    info_rank = [i['ScientificName'] for i in TaxInfo[num_id]['LineageEx'] if i['Rank'].lower() == rank.lower()][0]
                except IndexError:
                    info_rank = None
                taxinfo_dict[rank] = info_rank
            taxinfo_tab.loc[tax_id] = taxinfo_dict
        print(' ')
    print('Finally, the amount of taxinfo is ' + str(len(taxinfo_tab)))
    return taxinfo_tab


def esearch_workfolw(Find_What):
    print("0-Start esearch :"+Find_What)
    FullRecords_id = Entrez.read(Entrez.esearch(db="sra", RetMax=10000000, term=Find_What+"[Organism]"))
    idlist=FullRecords_id['IdList']
    print("1-The total data :"+str(len(idlist)))
    print("2-Get SRA information :")
    df_SRAinfo=get_SRAinfo_df(idlist)
    print("3-Get taxon info :")
    taxid_list=list(set(list(df_SRAinfo["SpeciesTaxid"])))
    df_taxinfo=get_taxinfo_df(taxid_list)
    print("4-Concat df_taxinfo and df_SRAinfo")
    #Summarydf = concat([df_taxinfo, df_SRAinfo], axis=1)
    Summarydf = pd.merge(df_taxinfo, df_SRAinfo, how='right', right_on="SpeciesTaxid", left_index=True)

    print("5-Output Summary Table")
    # RNA-seq to excel
    print("     >Delete duplicate species data and only keep the maximum length. ")
    RNA_Sort_table = Summarydf.loc[Summarydf["Library_source"] == "TRANSCRIPTOMIC"].sort_values(by=['Total_bases'], ascending=False)

    raw_rna_table_rows=len(RNA_Sort_table.index)
    print(" The number of RNA-Seq original data items: "+str(raw_rna_table_rows))

    RNA_final_table = RNA_Sort_table.drop_duplicates(subset=['Organism'], keep='first')

    filtered_table_rows=len(RNA_final_table.index)
    print(" The number of RNA-Seq filtered data is: "+str(filtered_table_rows))
    # WGS_info to excel
    WGS_Sort_table = Summarydf.loc[Summarydf["Library_source"] == 'GENOMIC'].sort_values(by=['Total_bases'], ascending=False)

    raw_wgs_table_rows=len(WGS_Sort_table.index)
    print(" The number of WGS original data items: "+str(raw_wgs_table_rows))

    WGS_filter_table=WGS_Sort_table[WGS_Sort_table['Total_bases'] > 5000000000]
    WGS_final_table = WGS_filter_table.drop_duplicates(subset=['Organism'], keep='first')
    print("     >Delete item that Total_bases < 5000000000")
    filtered_table_rows=len(WGS_final_table.index)
    print(" The number of WGS filtered data is: "+str(filtered_table_rows))

    #output_rawdata
    data = datetime.now().strftime('%Y%m%d')
    prefix=Find_What
    Summarydf.to_excel('0-rawdata_{0}_{1}.xlsx'.format(prefix, data))
    WGS_Sort_table.to_excel('1-WGS_all_{0}_{1}.xlsx'.format(prefix, data))
    RNA_Sort_table.to_excel('1-RNA_all_{0}_{1}.xlsx'.format(prefix, data))
    RNA_final_table.to_excel("2-{0}_RNA-Seq_{1}.xlsx".format(prefix, data))
    WGS_final_table.to_excel("2-{0}_WGS_{1}.xlsx".format(prefix, data))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', '-i',
                        type=str,
                        help='Search content or can be a file containing search content line by line.'
                             'Format requirements according to NCBI')
    args = parser.parse_args()
    esearch_workfolw(args.input)
    #esearch_workfolw('Mecoptera')
