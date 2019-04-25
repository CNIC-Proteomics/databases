#!/usr/bin/bash

# Declare variables
CODEDIR="d:/projects/databases"
OUTDIR="D:/Dropbox/databases/$(date +"%Y%m")" # with date folder
WSDIR="D:/Dropbox/databases/current_release"

# Function that executes the input command
run_cmd () {
  echo "-- $1"  
  eval $1
}

# for the following species...
# create the System biology database
SPECIES_LIST=(human mouse pig rabbit)
for SPECIES in "${SPECIES_LIST[@]}"
do
    # execute commands
    CMD="time python '${CODEDIR}/src/create_db_sb.py' -s ${SPECIES} -o '${OUTDIR}' -vv  &> '${CODEDIR}/logs/create_db_sb.${SPECIES}.log' "
    run_cmd "${CMD}"
done

# for the following species...
# create the Target/Decoy database
for FASTA in $(ls -1 "${OUTDIR}"/*.fa)
do
    # get local variables
    filename=$(basename "${FASTA}")
    filename="${filename%.*}"
    OUTFILE_dc="${OUTDIR}/${filename}.decoy.fa"
    OUTFILE_tg="${OUTDIR}/${filename}.target.fa"
    OUTFILE="${OUTDIR}/${filename}.target-decoy.fa"
    LOGFILE="${CODEDIR}/logs/decoyPYrat.${filename}.log"
    # execute commands
    CMD="time python '${CODEDIR}/src/decoyPYrat.v2.py' --output_fasta '${OUTFILE_dc}' --decoy_prefix=DECOY '${FASTA}' &> '${LOGFILE}' && cat ${OUTFILE_tg} ${OUTFILE_dc} > ${OUTFILE}"
    run_cmd "${CMD}"
done

# Delete the last version and copy the new version to the folder
rm -rf "${WSDIR}/*"
cp -rp "${OUTDIR}/." "${WSDIR}/."

