#!/bin/bash

# ==============================================================================
#? The Endothelion Pipeline
#?
#? Add the long description here as a multi-line comment.
#? It also supports some *markdown*!
# ==============================================================================

# --- General settings and variables -------------------------------------------

set -e # "exit-on-error" shell option
set -u # "no-unset" shell option

# For a friendlier use of colors in Bash
red=$'\e[1;31m'
grn=$'\e[1;32m'
yel=$'\e[1;33m'
end=$'\e[0m'

# Runtime option file
OPTS="./data/in/runtime_options.json"
# Set variables
in_path="$(cat $OPTS | jq -r ".in_path")"
count_type="$(cat $OPTS | jq -r ".count_type")"
threshold="$(cat $OPTS | jq -r ".threshold")"
GOIs="$(cat $OPTS | jq -r ".GOIs")"

# --- The pipeline starts here -------------------------------------------------

# Global count of TSV target files
files_found=$(find "${in_path}" -maxdepth 4 -type f -iname "*.tsv" | wc -l)
echo -e "\nFound ${files_found} TSV files to analyze.\n"

# Looping through files with spaces in their names or paths is not such a
# trivial thing...
OIFS="$IFS"
IFS=$'\n'
counter=1
for count_matrix in `find "${in_path}" -maxdepth 4 \
	-type f -iname "*.tsv" | sort`
do
	echo "----------------------------------------------------"
	echo "${yel}Analyzing file ${counter}/${files_found}:${end} ${count_matrix}"
	Rscript "./src/endo_profiler.R" \
		"$count_matrix" \
		"$count_type" \
		"$threshold" \
		"$GOIs" \
		"./data/out/$(dirname "${count_matrix#${in_path}}")"
	echo
	((counter++))
done
IFS="$OIFS"
echo -e "${grn}PIPELINE SUCCESSFULLY COMPLETED${end}\n"

exit 0 # Success exit status
