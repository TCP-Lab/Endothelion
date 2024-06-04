#!/bin/bash

# ==============================================================================
#? ICT+ absolute expression analysis in all endo-models
#?
#? This pipeline performs the *Phase I* Endothelion task for ICT+ absolute
#? expression assessment applied to all the available datasets.
#? The *endo_profiler* Bash wrapper is invoked to parse the `runtime_options`
#? JSON file and get the user-defined parameter values set therein (see
#? `./README.md` for details). Then, it searches the input directory (including
#? possible subfolders) for all the available TSV count matrices, which will be
#? sequentially fed to the associated R script for actual data analysis.
#? Further details about the analysis algorithm and its expected arguments can
#? be found in the header of the R script file `./src/endo_profiler.R`.
# ==============================================================================

# General settings and variables 
source ./src/commons.sh

echo -e "\n${mag}STARTING MASTER PROFILING${end}\n"
echo -e "${mag}=========================${end}"
# --- The pipeline starts here -------------------------------------------------

# ICT+ absolute expression profiling in all endo-models
echo -e "\n${grn}STEP 1 :: absolute expression profiling${end}"
bash ./src/endo_profiler_wrap.sh "./data/in/"

# Cross-study statistics from functional gene sets in all endo-models
echo -e "\n${grn}STEP 2 :: cross-study statistics${end}"
bash ./src/endo_function_wrap.sh "./data/out/"

# --- The pipeline ends here ---------------------------------------------------
if [[ $? -eq 0 ]]; then
	echo -e "${grn}PIPELINE COMPLETED SUCCESSFULLY${end}\n"
else
	echo -e "${red}\nPIPELINE FAILED${end}\n"
fi
