#!/bin/bash

# ==============================================================================
#? Make the global set of ICTs of interest by querying the MTP-DB
#?
#? Transportome sub-setting is done according to the following criteria:
#? - ICs: all
#? - AQPs: all
#? - SLCs: only fully-inorganic solute carriers
#? - Pumps: all
#? - ABCs: none
#? - GPCRs/RTKs: 52 receptors of interest (hard-coded in 'make_geneset.R')
# ==============================================================================

db_filename="MTPDB.sqlite"
archive_path="./data/in/MTP-DB/${db_filename}.gz"
db_path="./data/${db_filename}"
out_dir="./data/in"

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

# Make the gene-set
Rscript --vanilla "./src/make_geneset.R" \
    "$db_path" \
    "$out_dir"

echo "DONE!"
