# Endothelion
Large-scale comparison of ICT expression in healthy endothelium.

## Dependencies
### Bash
- jq (>= 1.7.1)
- kerblam (>= 1.0.0-rc.1) [optional]

### R
- ggplot2 (>= 3.5.0)
- tidyr (>= 1.3.1)
- dplyr (>= 1.1.4)
- httr (>= 1.4.7)
- AnnotationDbi (>= 1.60.2)
- org.Hs.eg.db (>= 3.16.0)
- PCAtools (>= 2.10.0)
- r4tcpl (>= 1.5.1)

### SeqLoader
- rlang (>= 1.1.3)
- magrittr (>= 2.0.3)

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


## hCMEC/D3 (Blood-Brain Barrier Cell Line)
### Query
- DB: [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra)
- Search: `hCMEC D3`
- Results: 12 _BioProjects_
- Zenodo: https://zenodo.org/records/12729454

### Included Studies
| ENA BioProject ID | Study Alias   |
|:-----------------:|:-------------:|
| PRJNA307652       | GSE76528      |
| PRJNA575504       | GSE138309     |
| PRJNA578611       | GSE139133     |
| PRJNA777606       | GSE187565     |
| PRJNA847413       | GSE205739     |
| PRJEB48614        | E-MTAB-11129  |
| PRJNA667281       | --            |
| PRJNA896725       | --            |

### Excluded Studies
| ENA BioProject ID | Study Alias   | Reason for Exclusion       |
|:-----------------:|:-------------:|:--------------------------:|
| PRJNA802135       | GSE195781     | duplicated/bad runs        |
| PRJNA607654       | GSE145581     | miRNA-Seq                  |
| PRJNA1073892      | GSE255171     | regulatory T cells (Tregs) |
| PRJNA307651       | GSE76530      | miRNA-Seq                  |

#### The PRJNA802135/GSE195781 Case
Out of the 5 control Runs of PRJNA802135/GSE195781 (2022), 3 turned out to be
___identical___ to other runs already published in the previous study of the
same group of authors (2016), under the different ID PRJNA307652/GSE76528.
Namely:

| PRJNA802135 (2022) | PRJNA307652 (2016) |
|:------------------:|:------------------:|
| SRR17833475        | SRR3085451         |
| SRR17833476        | SRR3085449         |
| SRR17833478        | SRR3085446         |

When comparing the corresponding FASTQ files between the two studies, the single
reads actually appear to be identical, and the different file size only depends
on the different pattern used for the heading line (_Field 1_) of each read.
```bash
# Check it out

zcat SRR17833478_1.fastq.gz | wc -l
zcat SRR3085446_1.fastq.gz | wc -l

zcat SRR17833478_1.fastq.gz | head
zcat SRR3085446_1.fastq.gz | head

zcat SRR17833478_1.fastq.gz | tail
zcat SRR3085446_1.fastq.gz | tail
```

This explains the batch effect that affects PRJNA802135 control samples (see PCA
and hierarchical clustering), the abnormal standard deviation levels that are
produced when all the samples are pooled together, and the _unlikely_ levels of
correlation between the 2 datasets (especially when removing the 2
different-looking runs from PRJNA802135).

As for the 2 non-duplicated samples of PRJNA802135, they showed a very poor
quality profile, such as severe problems of adapter contamination, eventually
featuring percentages of alignment of just 20%.

__For all these reasons, PRJNA802135 study was completely purged from the global
cohort.__
