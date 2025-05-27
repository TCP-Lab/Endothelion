#!/bin/bash

# ICs: all
# AQPs: all
# SLCs: only inorganic solute carriers
# Pumps: all
# ABCs: none
# GPCR: none

db_file="MTPDB.sqlite"
archive_path="./data/in/MTP-DB/${db_file}.gz"
db_path="./data/${db_file}"

# Extract the archive
if [[ ! -f "$db_path" ]]; then
    if [[ -f "$archive_path" ]]; then
        printf "Extracting MTP-DB ... "
        gzip -dc "$archive_path" > "$db_path"
        printf "Done\n"
    else
        printf "\nMTP-DB not found locally!"
        printf "Run 'kerblam data fetch' to download fresh database.\n"
    fi
fi

Rscript --vanilla "./src/make_genesets.R" \
    "$db_path" \
    "./data/out"

echo "DONE!"
