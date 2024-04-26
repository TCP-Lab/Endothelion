#!/bin/bash

# ==============================================================================
# Description here ...
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
target_path="$1"

# Set variables from runtime option JSON file
OPTS="./data/in/runtime_options.json"
GOIs="$(cat $OPTS | jq -r ".GOIs")"
subGOIs_prefix="$(cat $OPTS | jq -r ".subGOIs_prefix")"
central_tendency="$(cat $OPTS | jq -r ".central_tendency")"

# --- Main program -------------------------------------------------------------

# Search endo-models
models_found=$(find "$target_path" -type f -iname "*profileReport.csv" \
	| xargs -d "\n" -n 1 dirname | xargs -d "\n" -n 1 dirname \
	| sort | uniq | wc -l)
echo -e "          Found ${models_found} endo-model(s)\n"

# Looping through files with spaces in their names or paths is not such a
# trivial thing...
OIFS="$IFS"
IFS=$'\n'
counter=1
for endo_model in $(find "$target_path" -type f -iname "*profileReport.csv" \
	| xargs -d "\n" -n 1 dirname | xargs -d "\n" -n 1 dirname \
	| sort | uniq)
do
	echo "----------------------------------------------------"
	echo "${yel}Analyzing endo-model ${counter}/${models_found}:${end}"
	echo "$endo_model"
	Rscript "./src/endo_function.R" \
		"$endo_model" \
		"$central_tendency" \
		"$GOIs" \
		"$subGOIs_prefix" \
		"$endo_model"
	echo
	((counter++))
done
IFS="$OIFS"
