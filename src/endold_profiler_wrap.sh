#!/bin/bash

# ==============================================================================
# This Bash wrapper is invoked to parse the `runtime_options` JSON file and get
# the user-defined parameter values set therein. Then, it searches the input
# directory (including possible subfolders) for all the available TSV count
# matrices, which will be sequentially fed to the associated R script for actual
# data analysis. Note that the last parameter of `endold_profiler.R` (i.e., the
# output directory) is not retrieved from the JSON option file. On the contrary
# it is hard-coded by the wrapper so as to recreate within the `./data/out`
# directory (as per Kerblam! standard), a filesystem closely tracing the one
# found in `./data/in`, but containing just the results of the analysis.
# ==============================================================================

# --- General settings and variables -------------------------------------------
source ./src/bash_commons.sh

# Set the input path (based on which pipeline is running)
in_path="$1"
count_type="TPM"

# Set variables from runtime option JSON file
OPTS="./data/in/runtime_options.json"
GOIs="$(cat $OPTS | jq -r ".GOIs")"
threshold_adapt="$(cat $OPTS | jq -r ".threshold_adapt")"
threshold_value="$(cat $OPTS | jq -r ".threshold_value")"

# --- Main program -------------------------------------------------------------
# Global count of TSV target files
files_found=$(find "$in_path" -maxdepth 4 -type f \
	-iname "*_CountMatrix_*.tsv" | wc -l)
echo -e "          found ${files_found} TSV count matrices\n"

# Looping through files with spaces in their names or paths is not such a
# trivial thing...
OIFS="$IFS"
IFS=$'\n'
counter=1
for count_matrix in $(find "$in_path" -maxdepth 4 -type f \
	-iname "*_CountMatrix_*.tsv")
do
	echo "----------------------------------------------------"
	echo "${yel}Analyzing file ${counter}/${files_found}:${end}"
	echo "$count_matrix"
	Rscript "./src/endold_profiler.R" \
		"$count_matrix" \
		"$count_type" \
		"$threshold_adapt" \
		"$threshold_value" \
		"$GOIs" \
		"${in_path/\/in\//\/out\/}"
	echo
	((counter++))
done
IFS="$OIFS"
