# SCRIPTS FOR THE CREATION OF DATABASES IN THE PROTEOMIC UNIT OF CNIC

Execute the script which create the database for the system biology

```bash
.venv_win/Scripts/activate && python src/create_db_sb.py -s pig -o test
```

# DecoyPYrat.v2 - Fast Hybrid Decoy Sequence Database Creation for Proteomic Mass Spectromtery Analyses

Improvement on DecoyPYrat script from Sanger (https://www.sanger.ac.uk/science/tools/decoypyrat)

## Execution
1. Execute the decoyPYrat.v2.bat script.
2. Add the Target FASTA file (with *.fasta or *.fa extension)
2. Add the Decoy FASTA file (with *.fasta or *.fa extension)


```bash
/cygdrive/d/programs/Python3x/python src/decoyPYrat.v2.py  --output_fasta test/uniprot_MusMusculus_dic2016.decoy.fasta  test/uniprot_MusMusculus_dic2016.fasta

cat uniprot_MusMusculus_dic2016.target.fa uniprot_MusMusculus_dic2016.decoy.fasta > uniprot_MusMusculus_dic2016.target-decoy.fa
```

