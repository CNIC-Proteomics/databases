import sys, os, logging
import urllib.request
import datetime
import re
import zipfile
import json
from Bio import SwissProt
from Bio import SeqIO
from Bio.KEGG import REST


class creator:
    # https://www.uniprot.org/uniprot/?query=proteome:up000005640&format=fasta&include=yes&fil=reviewed:yes
    URL_UNIPROT = 'https://www.uniprot.org/uniprot/?'
    URL_UNIPROT += 'include=yes&' # include all isoforms
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
        'rat': {
            'scientific': 'Rattus norvegicus',
            'proteome': 'UP000002494'
        },
        'pig': {
            'scientific': 'Sus scrofa',
            'proteome': 'UP000008227'
        },
        'rabbit': {
            'scientific': 'Oryctolagus cuniculus',
            'proteome': 'UP000001811'
        },
        'zebrafish': {
            'scientific': 'Danio rerio',
            'proteome': 'UP000000437'
        }
    }
    LIST_IDS = ['Name','IsoIDs','Accession','Accessions','Gene','Class','Species','Description']
    LIST_TERMS = ['Ensembl','RefSeq','CCDS','GO','KEGG','PANTHER','Reactome','CORUM','DrugBank']
    HEADER = ['Category','Hit']
    TIME = datetime.datetime.now().strftime("%Y%m")

    '''
    Creates the databases
    '''
    def __init__(self, s, o, f=None, d=False):
        
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
        self.TMP_DIR = os.path.dirname(os.path.abspath(__file__)) +'/../tmp/'+ self.TIME +'/'+ self.species
        os.makedirs(self.TMP_DIR, exist_ok=True)
        logging.debug(self.TMP_DIR)
        
        # download sequences from UniProt
        self.outfname = species +'_'+ self.proteome_id +'_'+ self.TIME +'_'+ f if f else ''
        self.db_fasta = self.outdir +'/'+ self.outfname +'.fasta'        
        self._download_fasta_db(self.db_fasta, f)
        
        # remove duplicate sequences
        if d:
            logging.debug("remove duplicate sequences")
            self._remove_duplicates(self.db_fasta, self.db_fasta)
            
        # create data files
        self.db_uniprot = self.TMP_DIR +'/'+ self.outfname +'.uniprot.dat'
        self.db_corum   = self.TMP_DIR +'/'+ ".".join(os.path.basename( self.URL_CORUM ).split(".")[:-1]) # get the filename from the URL (without 'zip' extension)
        self.db_panther = self.TMP_DIR +'/'+ self.outfname +'.panther.dat'
        
        # create output files
        self.outfile = self.outdir +'/'+ self.outfname +'.categories.tsv'


    def _delete_tmp_dir(self, dir):
        files = [ f for f in os.listdir(dir) ]
        for f in files:
            try:
                os.remove(os.path.join(dir, f))
            except Exception as e:
                logging.error(e)


    def _download_fasta_db(self, outfile, filt):
        '''
        Download the fasta database file
        '''
        url = self.URL_UNIPROT +'query=proteome:'+ self.proteome_id        
        if filt and filt == "sw": # filter by SwissProt
            url += '&fil=reviewed:yes'            
        elif filt and filt == "tr": # filter by TrEMBL
            url += '&fil=reviewed:no'
        url += '&format=fasta'
        logging.debug('get '+url)
        urllib.request.urlretrieve(url, outfile)


    def _remove_duplicates(self, infile, outfile=None):
        '''
        Remove duplicated sequences
        '''
        seqs = dict()
        precords = SeqIO.parse(infile, "fasta")
        records = SeqIO.to_dict(precords)
        # read all sequences
        # delete the duplicated sequences based on the sorted id's
        for i in list(records):
            record = records[i]
            s = str(record.seq)
            if s in seqs:
                l = [seqs[s], i]
                logging.warning( "duplicated sequences: {}".format(",".join(l)) , exc_info=False)
                if outfile:
                    l.sort()
                    logging.warning( "deleting the sequences {}".format(l[1]) , exc_info=False)
                    del records[l[1]]
                    seqs[s] = l[0]
            else:
                seqs[s] = i
        # write file if apply
        if outfile:
            with open(outfile, 'w') as handle:
                SeqIO.write(records.values(), handle, 'fasta')
        return None


    def download_raw_dbs(self, filt):
        '''
        Download the raw databases
        '''
        # delete any temporal file
        # self._delete_tmp_dir(self.TMP_DIR)
        
        # UniProt
        # filter by SwissProt (Reviewd) if apply
        if not os.path.isfile(self.db_uniprot):
            url = self.URL_UNIPROT +'query=proteome:'+ self.proteome_id
            if filt and filt == "sw":
                url += '&fil=reviewed:yes'
            url += '&format=txt'
            logging.debug("get "+url)
            urllib.request.urlretrieve(url, self.db_uniprot)
        else:
            logging.debug('cached uniprot')
        
        # CORUM
        # download all complexes file (using the same name)
        # unzip the file
        if not os.path.isfile(self.db_corum):
            url = self.URL_CORUM
            db_dat = self.TMP_DIR +'/'+ os.path.basename(url)
            logging.debug("get "+url)
            urllib.request.urlretrieve(url, db_dat)
            zip_ref = zipfile.ZipFile(db_dat, 'r')
            zip_ref.extractall(self.TMP_DIR)
            zip_ref.close()
        else:
            logging.debug('cached corum')
        
        # PANTHER
        # get the list of species and extract the file name
        if not os.path.isfile(self.db_panther):
            url = self.URL_PANTHER
            result = urllib.request.urlopen(url).read().decode('utf-8')
            if result:
                pattern = re.search(r'\s*(PTHR[^\_]*\_'+self.species+')', result, re.I | re.M)
                if pattern:
                    url = self.URL_PANTHER + pattern[1]
                    logging.debug("get "+url)
                    urllib.request.urlretrieve(url, self.db_panther)
        else:
            logging.debug('cached panther')
    
    
    def create_qreport(self):
        '''
        Create protein report
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
            # create cross-references data
            logging.info('create cross-references data from UniProtKB database...')
            for record in SwissProt.parse( open(self.db_uniprot) ):
                # local variable
                outs = dict()
                
                # extract main info ---
                name = record.entry_name
                acc = record.accessions[0]
                accs = ";".join(record.accessions[1:])
                pattern = re.search(r'Name=(\w*)', record.gene_name, re.I | re.M)
                gene = pattern[1] if pattern else record.gene_name  
                pattern = re.search(r'[RecName|SubName]: Full=([^\;|\{]*)', record.description, re.I | re.M)
                dsc = pattern[1] if pattern else record.description
                dclass = record.data_class
                pattern = re.search(r'([\w|\s]*)\s+\(\w*\)', record.organism, re.I | re.M)
                species = pattern[1] if pattern else record.organism
                # extract isoforms IDs
                comm = [c for c in record.comments if 'ALTERNATIVE PRODUCTS:' in c]
                if comm:
                    IsoIds = re.findall(r'IsoId=([^\;]*)\;', comm[0], re.I | re.M | re.DOTALL)
                    IsoIds = ";".join(IsoIds)
                    # delete *-1 prefix from isoform Ids
                    IsoIds = re.sub(rf"{acc}-1",f"{acc}", IsoIds)
                else:
                    IsoIds = acc
                # save to out report
                outs['Name']        = name
                outs['IsoIDs']      = IsoIds
                outs['Accession']   = acc
                outs['Accessions']  = accs
                outs['Gene'] = gene
                outs['Class'] = dclass
                outs['Species'] = species
                outs['Description'] = dsc

                # create cross-references data ---
                # filter by given list of terms
                xs = [x for x in record.cross_references if x[0] in self.LIST_TERMS]
                # create dictionary with the common Xreferences
                xrefs = dict()
                for xref in xs:
                    xrefs.setdefault(xref[0], []).append(xref[1:])
                # convert the xref data to plain text in one line
                for extdb,xref in xrefs.items():
                    extdesc = ''
                    if extdb == "Ensembl" or extdb == "RefSeq" or extdb == "CCDS":
                        extdesc = self._extract_xref_ids(xref, acc)
                    elif extdb == "GO":
                        extdesc = self._extract_cat_go(xref)
                    elif extdb == "KEGG": # remote access
                        extdesc = self._extract_cat_kegg(xref)
                    elif extdb == "PANTHER":
                        extdesc = self._extract_cat_panther(panther_txt, xref, acc)
                    elif extdb == "Reactome":
                        extdesc = self._extract_cat_reactome(xref)
                    elif extdb == "CORUM":
                        extdesc = self._extract_cat_corum(corum_json, acc)
                    elif extdb == "DrugBank":
                        extdesc = self._extract_cat_drugbank(xref)
                    if extdesc != '':
                        # replace bad characters
                        extdesc = extdesc.replace("\t"," ")
                        extdesc = extdesc.replace('–','-')
                        # delete *-1 prefix from isoform Ids
                        extdesc = re.sub(rf"\[{acc}-1\]",f"[{acc}]", extdesc)
                    # save to out report
                    outs[extdb] = extdesc
                
                # create the line of output text
                output += "\t".join([outs[o] if o in outs else '' for o in self.LIST_IDS+self.LIST_TERMS])
                output += "\n"
        return output
      
          
    def _extract_xref_ids(self, xref, acc):
        '''
        Parse the xref data
        '''
        out = ''
        for x in xref:
            id = x[0]
            dsc = ''
            for y in x[1:]:
                m = re.search(rf"\[({acc}[^\]]*)\]\s*$", y, re.I | re.M)
                if m:
                    id += f"[{m[1]}]"
                    y = re.sub(rf"\.\s*\[{acc}[^\]]*\]\s*",'', y)
                dsc += f"{y}|"
            dsc = re.sub(r'[-|\|]*\s*$','', dsc) # delete - or | at the end of string
            out += f"{id}>{dsc};"
        out = re.sub(r'[>|;]*$','', out)# delete ; or > at the end of string
        return out

    def _extract_cat_go(self, xref):
        '''
        Parse the xref data
        '''
        out = ''
        for x in xref:
            id = x[0]
            dsc = "|".join(x[1:])
            out += f"{id}>{dsc};"
        out = re.sub(r'[>|;]*$','', out)# delete ; or > at the end of string
        return out

    def _extract_cat_kegg(self, xref):
        '''
        Parse the raw database file
        '''
        out = ''
        for x in xref:
            id = x[0]
            dsc = ''
            try:
                record = REST.kegg_get(id).read()
                if record:
                    pattern = re.search(r'DEFINITION\s*([^\n]*)', record, re.I | re.M)
                    dsc += pattern[1] if pattern else ''
                    pattern = re.search(r'PATHWAY\s*([\w\W]*)MODULE', record, re.I | re.M)
                    dsc += "|"+re.sub(r'\s*\n\s*','|', pattern[1]) if pattern else ''
                pass
            except:
                pass
            dsc = re.sub(r'[-|\|]*\s*$','', dsc) # delete - or | at the end of string
            out += f"{id}>{dsc};"
        out = re.sub(r'[>|;]*$','', out)# delete ; or > at the end of string
        return out

    def _extract_cat_panther(self, datatxt, xref, acc):
        '''
        Parse the raw database file
        '''
        out = ''
        if datatxt:
            pattern = re.search(rf"UniProtKB={acc}\t*([^\t]*)\t*([^\t]*)", datatxt, re.I | re.M)
            out += pattern[1]+'|'+pattern[2] if pattern else ''
        if not datatxt or out == '':
            out = ";".join([x[0] for x in xref])
        out = re.sub(r'[>|;]*$','', out)# delete ; or > at the end of string
        return out

    def _extract_cat_reactome(self, xref):
        '''
        Parse the xref data
        '''
        out = ''
        for x in xref:
            id = x[0]
            dsc = "|".join(x[1:])
            out += f"{id}>{dsc};"
        out = re.sub(r'[>|;]*$','', out)# delete ; or > at the end of string
        return out

    def _extract_cat_corum(self, datatxt, acc):
        '''
        Parse the raw database file
        '''
        out = ''
        if datatxt:
            comps = list(filter(lambda person: acc in person['subunits(UniProt IDs)'], datatxt))
            if comps:
                out += ";".join([ 'compID_'+str(comp['ComplexID'])+'>'+comp['ComplexName'] for comp in comps if 'ComplexID' in comp and 'ComplexName' in comp ])
        return out

    def _extract_cat_drugbank(self, xref):
        '''
        Parse the xref data
        '''
        out = ''
        for x in xref:
            id = x[0]
            dsc = "|".join(x[1:])
            out += f"{id}>{dsc};"
        out = re.sub(r'[>|;]*$','', out)# delete ; or > at the end of string
        return out


    def to_file(self, output):
        '''
        Print to file
        '''
        f = open(self.outfile, "w")
        # create header of output
        header = "\t".join([o for o in self.LIST_IDS+self.LIST_TERMS])
        header += "\n"
        f.write(header)
        f.write(output)
        f.close()

