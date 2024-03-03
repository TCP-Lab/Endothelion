#!/bin/bash

# ==============================================================================
#? The Endothelion Pipeline - alpha
#?
#? This pipeline is used to visualize the absolute expression of some genes of
#? interest (GOIs) out of a count matrix containing the whole genome and many
#? biological replicates. The `endo_profiler` R script is used for this task.
#? It sequentially takes as its first argument all the TSV files found in the
#? input folder (including possible subfolders). Most of the parameter values
#? are retrieved from the `runtime_options` JSON file. Namely,
#?  1. *in_path*: the input folder (`./data/in` in a *Kerblam!* project);
#?  2. *count_type*: count units (usually *TPMs* for absolute expression);
#?  3. *threshold_adapt*: should a GMM be used to find the expression threshold?
#?  4. *threshold_value*: threshold value or the number of Gaussian components;
#?  5. *GOIs*: path to the list of GOIs.
#? Specifically, when `threshold_adapt == true`, the integer entered as the
#? `threshold_value` indicates the number of Gaussian components to be used in
#? the mixture for the adaptive calculation of the expression threshold. In
#? contrast, when `threshold_adapt == false`, the `threshold_value` is the real
#? value to be used as constant threshold in all the experiments.
#? The last parameter of `endo_profiler.R` (i.e., the output directory) is not
#? retrieved from the JSON option file. On the contrary it is hard-coded so as
#? to recreate within the `./data/out` directory (*Kerblam!* standard), a
#? filesystem similar to that found in `./data/in`, but containing only the
#? results of the analysis.
#? More details on the actual analysis algorithm and its expected arguments can
#? be found in the header of the R script file `./src/endo_profiler.R`.
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
threshold_adapt="$(cat $OPTS | jq -r ".threshold_adapt")"
threshold_value="$(cat $OPTS | jq -r ".threshold_value")"
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
for count_matrix in `find "${in_path}" -maxdepth 4 -type f -iname "*.tsv"`
do
	echo "----------------------------------------------------"
	echo "${yel}Analyzing file ${counter}/${files_found}:${end} ${count_matrix}"
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
echo -e "${grn}PIPELINE SUCCESSFULLY COMPLETED${end}\n"

exit 0 # Success exit status
