#!/bin/bash

# ==============================================================================
#? Makes the global set of ICTs of interest by querying the MTP-DB (v1.25.24)
#?
#? Transportome sub-setting is done according to the following criteria:
#? - ICs: all (436)
#? - AQPs: all (14)
#? - SLCs: only fully-inorganic solute carriers (81)
#? - Pumps: all (90)
#? - ABCs: none
#? - GPCRs/RTKs: 51 receptors of interest (hard-coded in 'make_geneset.R')
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
