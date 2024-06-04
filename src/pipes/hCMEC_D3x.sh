#!/bin/bash

# ==============================================================================
#? ICT+ absolute expression analysis in hCMEC-D3 cell line using SeqLoader
#?
#? This pipeline performs the *Phase I* Endothelion task for ICT+ absolute
#? expression assessment applied to multiple hCMEC-D3 cell line datasets.
#? Unlike hCMEC, this new pipeline takes advantage of SeqLoader for defining the
#? the S3 xSeries and xModel classes in R.
# ==============================================================================

# General settings and variables 
source ./src/commons.sh

# Set variables from runtime option JSON file
OPTS="./data/in/runtime_options.json"
GOIs="$(cat $OPTS | jq -r ".GOIs")"
central_tendency="$(cat $OPTS | jq -r ".central_tendency")"
threshold_adapt="$(cat $OPTS | jq -r ".threshold_adapt")"
threshold_value="$(cat $OPTS | jq -r ".threshold_value")"

echo -e "\n${mag}STARTING hCMEC-D3 PROFILING${end}"
echo -e "${mag}===========================${end}"
# --- The pipeline starts here -------------------------------------------------

# Set the input path (pipeline-specific)
in_path="./data/in/Lines/hCMEC_D3/"

# ICT+ absolute expression profiling in hCMEC-D3 cell line
echo -e "\n${grn}STEP 1 :: absolute expression profiling${end}"
Rscript "./src/endo_profiler_x.R" \
	"$in_path" \
	"$central_tendency" \
	"$threshold_adapt" \
	"$threshold_value" \
	"$GOIs" \
	"${in_path/\/in\//\/out\/}"

# --- The pipeline ends here ---------------------------------------------------
if [[ $? -eq 0 ]]; then
	echo -e "${grn}PIPELINE COMPLETED SUCCESSFULLY${end}\n"
else
	echo -e "${red}\nPIPELINE FAILED${end}\n"
fi
