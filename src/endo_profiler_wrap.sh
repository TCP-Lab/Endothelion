#!/bin/bash

# ==============================================================================
# This endo_profiler Bash wrapper is invoked to parses the `runtime_options`
# JSON file and get the user-defined parameter values set therein. Then, it
# searches the input directory (including possible subfolders) for all the
# available TSV count matrices, which will be sequentially fed to the related R
# script for actual data analysis. The last parameter of `endo_profiler.R`
# (i.e., the output directory) is not retrieved from the JSON option file. On
# the contrary it is hard-coded by the wrapper so as to recreate within the
# `./data/out` directory (Kerblam! standard), a filesystem closely tracing the
# one found in `./data/in`, but containing only the results of the analysis.
# ==============================================================================

# --- General settings and variables -------------------------------------------

# The so-called strict mode (see https://mywiki.wooledge.org/BashFAQ/105)
set -e           # "exit-on-error" shell option
set -u           # "no-unset" shell option
set -o pipefail  # exit on within-pipe error

# For a friendlier use of colors in Bash
red=$'\e[1;31m' # Red
grn=$'\e[1;32m' # Green
yel=$'\e[1;33m' # Yellow
end=$'\e[0m'

# Set the input path, depending on the pipeline that is running
in_path="$1"
count_type="TPM"

# Set variables from runtime option JSON file
OPTS="./data/in/runtime_options.json"
GOIs="$(cat $OPTS | jq -r ".GOIs")"
threshold_adapt="$(cat $OPTS | jq -r ".threshold_adapt")"
threshold_value="$(cat $OPTS | jq -r ".threshold_value")"

# --- Main program -------------------------------------------------------------

# Global count of TSV target files
files_found=$(find "${in_path}" -maxdepth 4 -type f \
	-iname "*_CountMatrix_*.tsv" | wc -l)
echo -e "Found ${files_found} TSV files to analyze.\n"

# Looping through files with spaces in their names or paths is not such a
# trivial thing...
OIFS="$IFS"
IFS=$'\n'
counter=1
for count_matrix in $(find "${in_path}" -maxdepth 4 -type f \
	-iname "*_CountMatrix_*.tsv")
do
	echo "----------------------------------------------------"
	echo "${yel}Analyzing file ${counter}/${files_found}:${end}"
	echo "$count_matrix"
	Rscript "./src/endo_profiler.R" \
		"$count_matrix" \
		"$count_type" \
		"$threshold_adapt" \
		"$threshold_value" \
		"$GOIs" \
		"./data/out/$(dirname "${count_matrix#${in_path}}")"
	echo
	((counter++))
done
IFS="$OIFS"
