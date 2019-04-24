#!/usr/bin/bash

# Declare variables
CODEDIR="d:/projects/databases"
OUTDIR="D:/Dropbox/databases/$(date +"%Y%m")" # with date folder
WSDIR="D:/Dropbox/databases/current_release"

# Create the System biology database for each given species
SPECIES_LIST=(human mouse pig rabbit)
for SPECIES in "${SPECIES_LIST[@]}"
do
    CMD="time python '${CODEDIR}/src/create_db_sb.py' -s ${SPECIES} -o '${OUTDIR}' -vv  &> '${CODEDIR}/logs/create_db_sb.${SPECIES}.log'"
    echo "-- ${CMD}"
    eval ${CMD}
done

# Delete the last version and copy the new version to the folder
rm -rf "${WSDIR}/*"
cp -rp "${OUTDIR}/." "${WSDIR}/."

