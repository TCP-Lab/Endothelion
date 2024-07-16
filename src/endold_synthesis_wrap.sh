#!/bin/bash

# ==============================================================================
# To find the different endothelial models (endo-models), output directory is
# first searched for all "*profileReport.csv" files as returned from the
# previous analysis step, then basenames and the lowest levels of their paths
# (i.e., the GEO series ID part) are removed by applying `dirname` twice.
# Unique results of such `find` statement represent the endo-models addressed by
# the pipeline, under the hypothesis that input filesystem was properly set.
# ==============================================================================

# --- General settings and variables -------------------------------------------
source ./src/bash_commons.sh

# Set the input path (based on which pipeline is running)
target_path="$1"

# Set variables from runtime option JSON file
OPTS="./data/in/runtime_options.json"
GOIs="$(cat $OPTS | jq -r ".GOIs")"
central_tendency="$(cat $OPTS | jq -r ".central_tendency")"

# --- Main program -------------------------------------------------------------
# Search endo-models
models_found=$(find "$target_path" -type f -iname "*profileReport.csv" \
	| xargs -d "\n" -n 1 dirname | xargs -d "\n" -n 1 dirname \
	| sort | uniq | wc -l)
echo -e "          found ${models_found} endo-model(s)\n"

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
	Rscript "./src/endold_synthesis.R" \
		"$endo_model" \
		"$central_tendency" \
		"$GOIs" \
		"ICT_subset_" \
		"$endo_model"
	echo
	((counter++))
done
IFS="$OIFS"
