#!/bin/bash

# ==============================================================================
#? ICT+ absolute expression analysis in hCMEC-D3 cell line using SeqLoader
#?
#? This pipeline performs the *Phase I* Endothelion task for ICT+ absolute
#? expression assessment applied to multiple hCMEC-D3 cell line datasets.
#? Unlike hCMEC, this new pipeline takes advantage of SeqLoader for defining the
#? the S3 xSeries and xModel classes in R.
# ==============================================================================

# For a friendlier use of colors in Bash
red=$'\e[1;31m' # Red
grn=$'\e[1;32m' # Green
yel=$'\e[1;33m' # Yellow
end=$'\e[0m'


# Set variables from runtime option JSON file
OPTS="./data/in/runtime_options.json"
GOIs="$(cat $OPTS | jq -r ".GOIs")"
threshold_adapt="$(cat $OPTS | jq -r ".threshold_adapt")"
threshold_value="$(cat $OPTS | jq -r ".threshold_value")"

echo -e "\n\e[1;35mSTARTING hCMEC-D3 PROFILING\e[0m"
echo -e "\e[1;35m===========================\e[0m"
# --- The pipeline starts here -------------------------------------------------

# Set the input path (based on which pipeline is running)
in_path="./data/in/Lines/hCMEC_D3/"


# ICT+ absolute expression profiling in hCMEC-D3 cell line
echo -e "\n\e[1;32mSTEP 1 :: absolute expression profiling\e[0m"

echo "----------------------------------------------------"
echo "${yel}Analyzing file ${counter}/${files_found}:${end}"
echo "$count_matrix"
Rscript "./src/endo_profiler_x.R" \
	"$in_path" \
	"$threshold_adapt" \
	"$threshold_value" \
	"$GOIs" \
	"${in_path/\/in\//\/out\/}"


# --- The pipeline ends here ---------------------------------------------------
if [[ $? -eq 0 ]]; then
	echo -e "\e[1;32mPIPELINE COMPLETED SUCCESSFULLY\e[0m\n"
else
	echo -e "\e[1;31m\nPIPELINE FAILED\e[0m\n"
fi
