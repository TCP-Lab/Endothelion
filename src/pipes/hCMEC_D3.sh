#!/bin/bash

# ==============================================================================
#? ICT+ absolute expression analysis in hCMEC-D3 cell line
#?
#? This pipeline performs the *Phase I* Endothelion task for ICT+ absolute
#? expression assessment applied to multiple hCMEC-D3 cell line datasets.
#? The *endo_profiler* Bash wrapper is invoked to parse the `runtime_options`
#? JSON file and get the user-defined parameter values set therein (see
#? `./README.md` for details). Then, it searches the input directory (including
#? possible subfolders) for all the available TSV count matrices, which will be
#? sequentially fed to the associated R script for actual data analysis.
#? Further details about the analysis algorithm and its expected arguments can
#? be found in the header of the R script file `./src/endo_profiler.R`.
# ==============================================================================

echo -e "\n\e[1;35mSTARTING hCMEC-D3 PROFILING\e[0m"
echo -e "\e[1;35m===========================\e[0m"
# --- The pipeline starts here -------------------------------------------------

# ICT absolute expression profiling in hCMEC-D3 cell line
echo -e "\n\e[1;32mSTEP 1 :: absolute expression profiling\e[0m"
bash ./src/endo_profiler_wrap.sh "./data/in/Lines/hCMEC_D3/"

# Cross-study statistics from functional gene sets in hCMEC-D3 cell line
echo -e "\n\e[1;32mSTEP 2 :: cross-study statistics\e[0m"
bash ./src/endo_function_wrap.sh "./data/out/Lines/hCMEC_D3/"

# --- The pipeline ends here ---------------------------------------------------
if [[ $? -eq 0 ]]; then
	echo -e "\e[1;32mPIPELINE COMPLETED SUCCESSFULLY\e[0m\n"
else
	echo -e "\e[1;31m\nPIPELINE FAILED\e[0m\n"
fi
