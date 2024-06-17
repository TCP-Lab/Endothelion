#!/bin/bash

# ==============================================================================
#? ICT+ absolute expression analysis in hCMEC/D3 cell line (without SeqLoader)
#?
#? This pipeline performs the *Phase I* task of the Endothelion project for the
#? analysis of ICT+ absolute expression in multiple hCMEC/D3 cell line datasets.
#? In this first implementation of the workflow (written before SeqLoader
#? definition) two Bash wrappers are invoked to parse the `runtime_options` JSON
#? file and get the user-defined parameter values set therein. After that, they
#? search the input directory for suitable expression data files, which will be
#? sequentially fed to the associated R scripts for actual data analysis. See
#? comments in the header of `endold_profiler.R` and `endold_synthesis.R` for
#? more info about the analysis algorithm and its expected arguments.
# ==============================================================================

# General settings and variables 
source ./src/bash_commons.sh

echo -e "\n${mag}STARTING hCMEC-D3 PROFILING${end}"
echo -e "${mag}===========================${end}"
# --- The pipeline starts here -------------------------------------------------

# ICT+ absolute expression profiling in hCMEC/D3 cell line
echo -e "\n${grn}STEP 1 :: absolute expression profiling${end}"
bash ./src/endold_profiler_wrap.sh "./data/in/Lines/hCMEC_D3/"

# Cross-study statistics from functional gene sets in hCMEC-D3 cell line
echo -e "\n${grn}STEP 2 :: cross-study statistics${end}"
bash ./src/endold_synthesis_wrap.sh "./data/out/Lines/hCMEC_D3/"

# --- The pipeline ends here ---------------------------------------------------
if [[ $? -eq 0 ]]; then
	echo -e "${grn}PIPELINE COMPLETED SUCCESSFULLY${end}\n"
else
	echo -e "${red}\nPIPELINE FAILED${end}\n"
fi
