# Endothelion
Large-scale comparison of ICT expression in healthy endothelium.

## Dependencies
### Bash
- jq (>= 1.7.1)
- kerblam (>= 1.0.0-rc.1) [optional]

### R
- ggplot2 (>= 3.5.0)
- r4tcpl (>= 1.5.1)

## Phase I: Absolute Assessment
The first phase of __Endothelion__ was designed to profile (i.e., to provide a
graphical representation in terms of bar charts of) the absolute expression of
some given genes of interest (GOIs) for different endothelial models (both cell
lines and tissues) starting from available count matrices containing the whole
genome and many biological replicates.

For the purpose of the Endothelion project,
GOI list is the list of Ion Channels and Transporters (ICTs) genes as defined
in the `ICT_set.csv` file (which actually also contains some GPCR receptors),
but in principle it can be *any* set of genes. Most importantly, counts
within the expression matrices are assumed to be already in TPM units to
actually be suitable for absolute expression evaluation.
The `endo_profiler` R script and the related Bash wrapper are used for this
task. The shell wrapper first parses the the `runtime_options` JSON file to
get the user-defined parameter values set therein. Then, it searches the
input directory (including possible subfolders) for all the available count
matrices, which will be sequentially fed to the R script for actual data
analysis. Currently, the runtime options are just three:
1. *GOIs*: path to the list of GOIs;
2. *threshold_adapt*: should a GMM be used to find the expression threshold?
3. *threshold_value*: threshold value or the number of Gaussian components.
Specifically, when `threshold_adapt == true`, the integer entered as the
`threshold_value` indicates the number of Gaussian components to be used in
the mixture for the adaptive calculation of the expression threshold. In
contrast, when `threshold_adapt == false`, the `threshold_value` is the real
value to be used as constant threshold in all the experiments.
The last parameter of `endo_profiler.R` (i.e., the output directory) is not
retrieved from the JSON option file. On the contrary it is *hard-coded* by
the wrapper so as to recreate within the `./data/out` directory (*Kerblam!*
standard), a filesystem tracing the one found in `./data/in`, but containing
only the results of the analysis.
More details about the analysis algorithm and its expected arguments can be
found in the header of the R script file `./src/endo_profiler.R`.