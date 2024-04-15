#!/bin/bash

# ==============================================================================
#? Global ICT absolute expression
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

bash ./src/endo_profiler_wrap.sh "./data/in/"
