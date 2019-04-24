import sys, os, logging
import urllib.request
import datetime
import re
import zipfile
import json
from Bio import SwissProt
from Bio import SeqIO
from Bio.KEGG import REST
from Bio.KEGG import Enzyme


class creator:
    URL_UNIPROT = 'https://www.uniprot.org/uniprot/?'
    URL_CORUM   = 'http://mips.helmholtz-muenchen.de/corum/download/allComplexes.json.zip'
    URL_PANTHER = 'ftp://ftp.pantherdb.org/sequence_classifications/current_release/PANTHER_Sequence_Classification_files/'
    SPECIES_LIST = {
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
    LIST_TERMS = ['GO','KEGG','PANTHER','Reactome','CORUM','DrugBank']
    HEADER = ['Category','Hit']
    TIME = datetime.datetime.now().strftime("%Y%m")

    '''
    Creates the databases
    '''
    def __init__(self, s, o, i=None, f=None):
        # assign species
        species = s.lower()
        if species in self.SPECIES_LIST:
            self.species = species
            self.proteome_id = self.SPECIES_LIST[self.species]['proteome']
        else:
            sys.exit( "ERROR: Species parameter has been not found. Try with: "+", ".join(self.SPECIES_LIST.keys()) )
        
        # create output directory if does not exist
        self.outdir = o
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir, exist_ok=True)
        # create temporal file
        self.TMP_DIR = os.path.dirname(os.path.abspath(__file__)) + '/../tmp/'+ self.species
        os.makedirs(self.TMP_DIR, exist_ok=True)
        logging.debug(self.TMP_DIR)
        # get/create the fasta sequences
        if i:
            self.db_fasta = i
            self.outfname = ".".join(os.path.basename(i).split(".")[:-1])
        else:
            self.outfname = species +'_'+ self.proteome_id +'_'+ self.TIME
            self.db_fasta = self.outdir +'/'+ self.outfname +'.fasta'
            self.download_fasta_db(self.db_fasta, f)
        # create data files
        self.db_uniprot = self.TMP_DIR +'/'+ self.outfname +'.uniprot.dat'
        self.db_corum   = self.TMP_DIR +'/'+ ".".join(os.path.basename( self.URL_CORUM ).split(".")[:-1]) # get the filename from the URL (without 'zip' extension)
        self.db_panther = self.TMP_DIR +'/'+ self.outfname +'.panther.dat'
        # create output file
        self.outfile = self.outdir +'/'+ self.outfname +'.tsv'
        # delete any temporal file
        self._delete_tmp_dir(self.TMP_DIR)

    def _delete_tmp_dir(self, dir):
        files = [ f for f in os.listdir(dir) ]
        for f in files:
            try:
                os.remove(os.path.join(dir, f))
            except Exception as e:
                logging.error(e)

    def download_raw_dbs(self, filt):
        '''
        Download the raw databases
        '''
        # UniProt
        # filter by SwissProt (Reviewd) if apply
        url = self.URL_UNIPROT +'query=proteome:'+ self.proteome_id
        url = url +'&format=txt'
        if filt and filt == "sw":
            url += '%20reviewed:yes'
        url += '&format=fasta'
        logging.debug("get "+url)
        urllib.request.urlretrieve(url, self.db_uniprot)
        
        # CORUM
        # download all complexes file (using the same name)
        # unzip the file
        url = self.URL_CORUM
        db_dat = self.TMP_DIR +'/'+ os.path.basename(url)
        logging.debug("get "+url)
        urllib.request.urlretrieve(url, db_dat)
        zip_ref = zipfile.ZipFile(db_dat, 'r')
        zip_ref.extractall(self.TMP_DIR)
        zip_ref.close()
        
        # PANTHER
        # get the list of species and extract the file name
        url = self.URL_PANTHER
        result = urllib.request.urlopen(url).read().decode('utf-8')
        if result:
            pattern = re.search(r'\s*(PTHR[^\_]*\_'+self.species+'\_)', result, re.I | re.M)
            if pattern:
                url = self.URL_PANTHER + pattern[1]
                logging.debug("get "+url)
                urllib.request.urlretrieve(url, self.db_panther)
        

    def download_fasta_db(self, outfile, filt):
        '''
        Download the fasta database file
        '''
        url = self.URL_UNIPROT +'query=proteome:'+ self.proteome_id        
        if filt and filt == "sw": # filter by SwissProt
            url += '%20reviewed:yes'
        url += '&format=fasta'
        logging.debug('get '+url)
        urllib.request.urlretrieve(url, outfile)

    def extract_identifiers(self, regex):
        '''
        Extract the identifiers from the FASTA file
        '''
        ids = []
        pattern = regex if regex else '[^\|]*\|([^\|]*)\|' # by default is UniProt
        records = SeqIO.parse(self.db_fasta, "fasta")
        for record in records:
            match = re.search(pattern, record.description, re.I | re.M)
            if match:
                ids.append(match[1])
        return ids
        
    def extract_categories(self, ids, type):
        '''
        Parse the raw database file
        '''
        output = ''
        if self.db_uniprot:
            # create reports from external data
            logging.info('create reports from external data...')
            corum_json = None
            panther_txt = None
            if os.path.isfile(self.db_corum):
                with open(self.db_corum, 'r') as f:
                    corum_json = json.load(f)
            logging.debug('corum done')
            if os.path.isfile(self.db_panther):
                with open(self.db_panther, 'r') as f:
                    panther_txt = f.read()
            logging.debug('panther done')
            # Extract the info from the main database (UniProt), if apply
            logging.info('extract the info from the main database...')
            for record in SwissProt.parse( open(self.db_uniprot) ):
                prot_acc = record.accessions[0]
                if prot_acc in ids:
                    prot_accs = ";".join(record.accessions[1:])
                    pattern = re.search(r'Name=(\w*)', record.gene_name, re.I | re.M)
                    gene_name = pattern[1] if pattern else record.gene_name         
                    pattern = re.search(r'[RecName|SubName]: Full=([^\;|\{]*)', record.description, re.I | re.M)
                    dsc = pattern[1] if pattern else record.description
                    dclass = record.data_class
                    # identification column
                    if type == "gene":
                        id = gene_name
                    else:
                        id = prot_acc
                    main_columns = id
                    # save cross references
                    for reference in record.cross_references:
                        extdb = reference[0]
                        extdesc = ''
                        if extdb in self.LIST_TERMS:
                            extid = reference[1]
                            extdesc = "|".join(reference[1:])
                            if extdb == "KEGG":
                                extdesc = self._extract_cat_kegg(extid)
                            elif extdb == "CORUM":
                                extdesc = self._extract_cat_corum(id, corum_json)
                            elif extdb == "PANTHER":
                                extdesc = self._extract_cat_panther(id, panther_txt)
                            if extdesc != '':
                                output += extdesc +"\t"+ main_columns +"\n"
        return output

    def _extract_cat_kegg(self, id):
        '''
        Parse the raw database file
        '''            
        out = id+'|'
        try:
            record = REST.kegg_get(id).read()
            if record:
                pattern = re.search(r'DEFINITION\s*([^\n]*)', record, re.I | re.M)
                out += pattern[1] if pattern else ''
                pattern = re.search(r'PATHWAY\s*([^\n]*)', record, re.I | re.M)
                out += ";"+pattern[1] if pattern else ''
            pass
        except:
            pass
        return out

    def _extract_cat_corum(self, id, allComp):
        '''
        Parse the raw database file
        '''
        out = ''
        if allComp:
            comps = list(filter(lambda person: id in person['subunits(UniProt IDs)'], allComp))
            if comps:
                out += ";".join([ str(comp['ComplexID'])+'|'+comp['ComplexName'] for comp in comps if 'ComplexID' in comp and 'ComplexName' in comp ])
        return out

    def _extract_cat_panther(self, id, allFam):
        '''
        Parse the raw database file
        '''
        out = ''
        if allFam:
            pattern = re.search(r'UniProtKB='+id+'\t*([^\t]*)\t*([^\t]*)', allFam, re.I | re.M)
            out += pattern[1]+'|'+pattern[2] if pattern else ''
        return out

    def to_file(self, output):
        '''
        Print to file
        '''
        f = open(self.outfile, "w")
        f.write("\t".join(self.HEADER)+"\n")
        f.write(output)
        f.close()



