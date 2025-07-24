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

# --- General settings and variables -------------------------------------------
source ./src/bash_commons.sh

db_path="./data/MTPDB.sqlite"
out_dir="./data/in"

# --- Extract the archive ------------------------------------------------------
_extract_mtpdb "$db_path"

# --- Make the gene-set --------------------------------------------------------
printf "Making the geneset...\n"
Rscript --vanilla "./src/make_geneset.R" \
    "$db_path" \
    "$out_dir"

echo "DONE!"
