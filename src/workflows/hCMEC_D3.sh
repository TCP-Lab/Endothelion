#!/bin/bash

# ==============================================================================
#? ICT+ absolute expression analysis in hCMEC/D3 cell line (using SeqLoader)
#?
#? This pipeline performs the *Phase I* task of the Endothelion project for the
#? analysis of ICT+ absolute expression in multiple hCMEC/D3 cell line datasets.
#? Unlike the first implementation (i.e., `hCMECold.sh`), this new pipeline uses
#? *SeqLoader* for the definition of S3 *xSeries* and *xModel* classes in R,
#? allowing, among other things, direct invocation of a single, simpler, and
#? more maintainable R script, without the need for additional wrappers.
# ==============================================================================
 
# --- General settings and variables -------------------------------------------
source ./src/bash_commons.sh

# --- Extract the archive ------------------------------------------------------
db_path="./data/MTPDB.sqlite"
_extract_mtpdb "$db_path"

# Set variables from runtime option JSON file
OPTS="./data/in/runtime_options.json"
GOIs="$(cat $OPTS | jq -r ".GOIs")"
threshold_adapt="$(cat $OPTS | jq -r ".threshold.adapt")"
threshold_value="$(cat $OPTS | jq -r ".threshold.value")"
descriptive="$(cat $OPTS | jq -r ".meta_analysis.central_tendency.descriptive")"
chart_type="$(cat $OPTS | jq -r ".meta_analysis.chart.type")"

echo -e "\n${mag}STARTING hCMEC/D3 PROFILING${end}"
echo -e "${mag}===========================${end}"
# --- The pipeline starts here -------------------------------------------------

# Set the input path (pipeline-specific)
in_path="./data/in/Lines/hCMEC_D3/"

# Sample Quality Control
echo -e "\n${grn}STEP 00 :: Quality Control${end}"
if [[ ! -f ./src/pca_hc.R ]]; then
	# Temporary solution while waiting for x.FASTQ to be put in pipeline
	echo -e "\n${yel}Downloading the R script for PCA (by x.FASTQ)${end}"
	wget -P "./src" \
		https://raw.githubusercontent.com/TCP-Lab/x.FASTQ/main/workers/pca_hc.R
fi
echo -e "\n${yel}Expression Boxes, PCA, and HC${end}"
Rscript --vanilla "./src/pca_hc.R" \
	'countmatrix.*\.tsv' \
	"${in_path/\/in\//\/out\/}" \
	"$in_path"

# ICT+ absolute expression profiling in hCMEC-D3 cell line
Rscript --vanilla "./src/endo_profiler.R" \
	"$in_path" \
	"$descriptive" \
	"$threshold_adapt" \
	"$threshold_value" \
	"$GOIs" \
	"./data/MTPDB.sqlite" \
	"$chart_type" \
	"${in_path/\/in\//\/out\/}"

# --- The pipeline ends here ---------------------------------------------------
if [[ $? -eq 0 ]]; then
	echo -e "${mag}===============================${end}"
	echo -e "${mag}PIPELINE COMPLETED SUCCESSFULLY\n${end}"
else
	echo -e "${red}\nPIPELINE FAILED${end}\n"
fi
