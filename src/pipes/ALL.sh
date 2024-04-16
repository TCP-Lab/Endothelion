#!/bin/bash

# ==============================================================================
#? Global ICT absolute expression profiling
#?
#? This pipeline performs the *Phase I* Endothelion task for ICT absolute
#? expression assessment applied to all the available datasets.
#? The *endo_profiler* Bash wrapper is invoked to parses the `runtime_options`
#? JSON file and get the user-defined parameter values set therein (see
#? `./README.md` for details). Then, it searches the input directory (including
#? possible subfolders) for all the available TSV count matrices, which will be
#? sequentially fed to the related R script for actual data analysis. Note that
#? the last parameter of `endo_profiler.R` (i.e., the output directory) is not
#? retrieved from the JSON option file. On the contrary it is *hard-coded* by
#? the wrapper so as to recreate within the `./data/out` directory (*Kerblam*!
#? standard), a filesystem closely tracing the one found in `./data/in`, but
#? containing only the results of the analysis.
#? Further details about the analysis algorithm and its expected arguments can
#? be found in the header of the R script file `./src/endo_profiler.R`.
# ==============================================================================

# --- The pipeline starts here -------------------------------------------------
echo -e "\n\e[1;32mSTARTING GLOBAL PROFILING\e[0m\n"

# Global ICT absolute expression profiling
bash ./src/endo_profiler_wrap.sh "./data/in/"

if [[ $? -eq 0 ]]; then
	echo -e "\e[1;32mPIPELINE COMPLETED SUCCESSFULLY\e[0m\n"
else
	echo -e "\e[1;31m\nPIPELINE FAILED\e[0m\n"
fi
