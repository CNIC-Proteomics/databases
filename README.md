# SCRIPTS FOR THE CREATION OF DATABASES IN THE PROTEOMIC UNIT OF CNIC

### Requirements:
Python3 and the following Python packages:
- sys
- os
- logging
- urllib.request
- datetime
- re
- zipfile
- json
- from Bio import SwissProt
- from Bio import SeqIO
- from Bio.KEGG import REST
- from Bio.KEGG import Enzyme


## Executions

The following scripts, download the FASTA sequences from UniProt proteomes and create the system biology database for the given species:
```bash
python src/create_db_sb.py -s human   -o databases -vv  &> logs/create_db_sb.human.log
python src/create_db_sb.py -s mouse   -o databases -vv  &> logs/create_db_sb.mouse.log
python src/create_db_sb.py -s pig     -o databases -vv  &> logs/create_db_sb.pig.log
python src/create_db_sb.py -s rabbit  -o databases -vv  &> logs/create_db_sb.rabbit.log
```

Create the system biology database from the given FASTA sequence:
```bash
python src/create_db_sb.py -s human   -i test/test.fa  -o test/human -vv 
python src/create_db_sb.py -s human   -i test/test.fa  -r '[^\|]*\|([^\|]*)\|' -o test/human -vv 
```

# DecoyPYrat.v2 - Fast Hybrid Decoy Sequence Database Creation for Proteomic Mass Spectromtery Analyses

Improvement on DecoyPYrat script from Sanger (https://www.sanger.ac.uk/science/tools/decoypyrat)

## Execution

```bash
python src/decoyPYrat.v2.py  --output_fasta test/uniprot_MusMusculus_dic2016.decoy.fasta  test/uniprot_MusMusculus_dic2016.fasta

cat uniprot_MusMusculus_dic2016.target.fa uniprot_MusMusculus_dic2016.decoy.fasta > uniprot_MusMusculus_dic2016.target-decoy.fa
```



