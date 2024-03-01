#!/bin/bash

# ============================================================================ #
#  The Endothelion Pipeline
# ============================================================================ #

# --- General settings and variables -------------------------------------------

set -e # "exit-on-error" shell option
set -u # "no-unset" shell option

# --- Pipeline starts here------------------------------------------------------

in_path="./data/in"
count_type="TPM"
threshold=true
GOIs="./ICT_set.csv"

# Looping through files with spaces in their names or paths is not such a
# trivial thing...
OIFS="$IFS"
IFS=$'\n'
for count_matrix in `find "${in_path}" -maxdepth 4 \
	-type f -iname "*.tsv" | sort`
do
	echo "Analyzing file: ${count_matrix}"
	Rscript endo_profiler.R \
		"$count_matrix" \
		"$count_type" \
		"$threshold" \
		"$GOIs"
	echo
done
IFS="$OIFS"
echo "DONE"
exit 0 # Success exit status

# GSE76528_TPM
# GSE138309_TPM
# GSE139133_TPM
# GSE195781_TPM
# GSE205739_TPM
