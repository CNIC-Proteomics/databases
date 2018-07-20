import sys, os, logging
import urllib.request
import datetime
import re
from Bio import SwissProt

__author__ = 'jmrodriguezc'

class creator:
    URL_UNIPROT = 'https://www.uniprot.org/uniprot/?'
    SPECIES = {
        'human': {
            'scientific': 'Homo sapiens',
            'proteome': 'UP000005640'
        },
        'mouse': {
            'scientific': 'Mus musculus',
            'proteome': 'UP000000589'
        },
        'pig': {
            'scientific': 'Sus scrofa',
            'proteome': 'UP000008227'
        },
        'rabbit': {
            'scientific': 'Oryctolagus cuniculus',
            'proteome': 'UP000001811'
        }
    }    
    TMP_DIR = 'tmp'
    LIST_TERMS = ['GO', ]
    HEADER = ['Category', 'Fasta']

    '''
    Creates the databases
    '''
    def __init__(self, s):
        # species variables
        self.species = s.lower()
        if self.species in self.SPECIES:
            self.proteome_id = self.SPECIES[ self.species ]['proteome']
        else:
            sys.exit("ERROR: Species parameter has been not found")
        # date time
        self.datetime = datetime.datetime.now().strftime("%Y%m%d")
        # delete any temporal file
        self._delete_tmp_dir(self.TMP_DIR)

    def _delete_tmp_dir(self, dir):
        files = [ f for f in os.listdir(dir) ]
        for f in files:
            try:
                os.remove(os.path.join(dir, f))
            except Exception as e:
                logging.error(e)

    def download_raw_db(self):
        '''
        Download the raw database file
        '''
        url = self.URL_UNIPROT
        if self.proteome_id:
            url += "query=proteome:"+ self.proteome_id
            url_txt = url+"&format=txt"
            logging.debug('get '+url_txt)
            self.db_dat = self.TMP_DIR + '/' + self.species + '_' + self.proteome_id + '_' + self.datetime + '.dat'
            urllib.request.urlretrieve(url_txt, self.db_dat)
            url_txt = url+"&format=fasta"
            logging.debug('get '+url_txt)
            self.db_fasta = self.TMP_DIR + '/' + self.species + '_' + self.proteome_id + '_' + self.datetime + '.fasta'
            urllib.request.urlretrieve(url_txt, self.db_fasta)

    def parse_raw_db(self):
        '''
        Parse the raw database file
        '''
        if self.db_dat:
            for record in SwissProt.parse( open(self.db_dat) ):
# PTHR22443~FAMILY NOT NAMED	>sp|A0AUZ9|KAL1L_HUMAN KAT8 regulatory NSL complex subunit 1-like protein OS=Homo sapiens GN=KANSL1L PE=1 SV=2
# metastatic colorectal cancer	>sp|A0AV02|S12A8_HUMAN Solute carrier family 12 member 8 OS=Homo sapiens GN=SLC12A8 PE=2 SV=3
# 778063:centrosomal protein 250kDa	>tr|Q15136|Q15136_HUMAN Protein kinase A-alpha (Fragment) OS=Homo sapiens GN=KIN27 PE=2 SV=1
# 778063:centrosomal protein 250kDa	>tr|Q5JPD5|Q5JPD5_HUMAN Centromere protein J OS=Homo sapiens GN=DKFZp667E025 PE=2 SV=1
                dsc = record.description
                org = re.search(r'OS=(\w* \w*)', record.organism, re.I | re.M)
                gen = re.search(r'Name=(\w*)', record.gene_name, re.I | re.M)
                comm = ">"+ record.accessions[0] +"|"+ record.entry_name +"|"+ dsc +" OS="+ org +" GN="+gen

                print("*****")
                print( comm )
                print( "--" )
                print(record.entry_name)
                print(",".join(record.accessions))
                print(record.keywords)
                print(repr(record.organism))
                print(record.sequence[:20] + "...")
                print("----")
                
                for reference in record.cross_references:
                    print(reference)
                # for reference in record.references:
                #     for ref in reference.references:
                #         print(ref)
                print( record['DR'] )
                # for feature in record.features:
                #     print(feature)

    def to_file(self, outdir):
        '''
        Print to file
        '''
        # create output directory if does not exist
        if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
        outfile = self.TMP_DIR + '/' + self.species + '_' + self.proteome_id + '_' + self.datetime + '.txt'


