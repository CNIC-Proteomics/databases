#!/usr/bin/bash

# Declare variables
DATE="$(date +"%Y%m")" # with date folder
CODEDIR="d:/projects/databases"
BASEDIR="//tierra.cnic.es/SC/U_Proteomica/UNIDAD/Databases/NextCloud"
OUTDIR="${BASEDIR}/${DATE}" # with date folder
WSDIR="${BASEDIR}/current_release"
LOGDIR="${CODEDIR}/logs/${DATE}"

# Function that executes the input command
run_cmd () {
  echo "-- $1"  
  eval $1
}

# prepare workspaces
mkdir "${LOGDIR}"
mkdir "${WSDIR}"

# for the following species...
# create the System biology database
SPECIES_LIST=(human mouse pig rabbit)
for SPECIES in "${SPECIES_LIST[@]}"
do
    # get local variables
    LOGFILE="${LOGDIR}/create_db_sb.${SPECIES}.log"
    # execute commands
    CMD="time python '${CODEDIR}/src/create_db_sb.py' -s ${SPECIES} -o '${OUTDIR}' -vv  &> '${LOGFILE}' "
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
    LOGFILE="${LOGDIR}/decoyPYrat.${filename}.log"
    # execute commands
    CMD="time python '${CODEDIR}/src/decoyPYrat.v2.py' --output_fasta '${OUTFILE_dc}' --decoy_prefix=DECOY '${FASTA}' &> '${LOGFILE}' && cat ${OUTFILE_tg} ${OUTFILE_dc} > ${OUTFILE}"
    run_cmd "${CMD}"
done

# Delete the last version and copy the new version to the folder
rm -rf "${WSDIR}/*"
cp -rp "${OUTDIR}/." "${WSDIR}/."

