# Endothelion
Large-scale comparison of ICT expression in healthy endothelium.

## Project Structure
___Endothelion___ is a project dedicated to exploring a special subset of the transportome--the ensemble of ion channels, pumps, and solute carriers (SLCs) that control inorganic ion transport and homeostasis--in human healthy endothelium.
### Phase I: Absolute Assessment
__Phase One__ focuses on characterizing the transportome within individual endothelial models, including both established cell lines and tissue samples.
By analyzing publicly available RNA-Seq datasets, each gene of interest (GOI) is assessed for presence or absence and visualized with intuitive charts.
This phase provides a clear snapshot of which transportome elements are expressed in each model, forming a foundation for further functional studies.
### Phase II: Differential Expression
__Phase Two__ moves beyond single models to examine differential expression across endothelial cells from different anatomical regions of healthy tissues.
Quantitative comparison of RNA-Seq data highlights genes that are enriched or depleted in specific vascular beds, uncovering regional specialization within the endothelium.
This analysis aims at revealing patterns in transportome composition that may underlie functional differences between endothelial populations.

## Transportome Analysis
### Making the Gene Set
For the purposes of the Endothelion project, the list of GOIs is compiled by querying our Membrane Transport Protein Database ([MTP-DB](https://github.com/TCP-Lab/MTP-DB)).
This gene set can be recreated on the fly by running the _Kerblam!_ workflow:
```bash
kerblam run make_geneset
```
The resulting list comprises 672 transportome elements selected for their relevance to inorganic ion transport and homeostasis, and includes the following categories:
- all known human ion channels (436 genes, i.e., the complete channelome, including possible auxiliary or modulatory subunits),
- the entire set of aquaporins (14 genes),
- all ATPase pumps (90 genes),
- solute carriers (SLCs) specific for inorganic solutes (81 genes),
- a set of 51 receptors (GPCRs and RTKs) with roles in inorganic ion dynamics.

> [!NOTE]
> Pay attention to protein vs. gene nomenclature for VEGFRs.
> Gene Symbols are very confusing:
> - _FLT1_ -_is for_-> VEGFR-1
> - _FLT2_ -_is for_-> this symbol doesn't exist any more (but it was the old name of FGFR1)
> - _FLT3_ -_is for_-> CD135 (i.e., the RTK receptor for the cytokine Flt3 ligand, _FLT3LG_)
> - _FLT4_ -_is for_-> VEGFR-3
> - _KDR_  -_is for_-> VEGFR-2

### Running the Analysis
The entire analysis workflow is implemented within the [_Kerblam!_](https://www.kerblam.dev/) project management open-source framework [(Visentin _et al._, 2025)](https://f1000research.com/articles/14-88), to ensure full transparency and result reproducibility.

You can retrieve the processed expression tables (i.e., outputs of the x.FASTQ pipeline) for all included studies directly from the corresponding Zenodo repository by:
```bash
kerblam data fetch
```
You can replicate all the steps of the Endothelion pipeline for transportome profiler in this way:
```bash
kerblam run hCMEC_D3
```

### Dependencies
#### Bash
- jq (>= 1.7.1)
- kerblam (>= 1.0.0-rc.1) [optional]
#### R
- limma (>= 3.58.1)
- sva (>= 3.50.0)
- ggplot2 (>= 3.5.0)
- Hmisc (>= 5.2-3)
- tidyr (>= 1.3.1)
- dplyr (>= 1.1.4)
- purrr (>= 1.0.2)
- httr (>= 1.4.7)
- DBI (>= 1.2.3)
- RSQLite (>= 2.3.5)
- AnnotationDbi (>= 1.60.2)
- org.Hs.eg.db (>= 3.16.0)
- PCAtools (>= 2.10.0)
- r4tcpl (>= 1.5.1)
### SeqLoader
- rlang (>= 1.1.3)
- magrittr (>= 2.0.3)

## hCMEC/D3 (Blood-Brain Barrier Cell Line)
### Query
- DB: [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra)
- Search: `hCMEC D3`
- Results: 12 _BioProjects_
- Zenodo: https://zenodo.org/records/12729454
- _Kerblam!_ workflow: `hCMEC_D3`

### Included Studies
| ENA BioProject ID | Study Alias | Ctrl Runs | Library | Median Read Length | Average Depth | Uniquely Mapped Reads | Platform | Reference |
|:-----------------:|:-----------:|:---------:|:-------:|:------------------:|:-------------:|:---------------------:|:--------:|:---------:|
| [PRJNA307652](https://www.ebi.ac.uk/ena/browser/view/PRJNA307652) | [GSE76528](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76528)   | 8 | PE | 2 × 51 bp  | 57.1 M | 78.9 % | Illumina HiSeq 2000   | [PMID: 26973449](https://pubmed.ncbi.nlm.nih.gov/26973449/) |
| [PRJNA575504](https://www.ebi.ac.uk/ena/browser/view/PRJNA575504) | [GSE138309](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138309) | 3 | PE | 2 × 78 bp  | 22.9 M | 91.3 % | Illumina NextSeq 550  | [PMID: 32757312](https://pubmed.ncbi.nlm.nih.gov/32757312/) |
| [PRJNA578611](https://www.ebi.ac.uk/ena/browser/view/PRJNA578611) | [GSE139133](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139133) | 2 | PE | 2 × 150 bp | 24.3 M | 95.3 % | Illumina NovaSeq 6000 | [PMID: 32985481](https://pubmed.ncbi.nlm.nih.gov/32985481/) |
| [PRJNA777606](https://www.ebi.ac.uk/ena/browser/view/PRJNA777606) | [GSE187565](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE187565) | 2 | PE | 2 × 150 bp | 27.8 M | 94.2 % | Illumina NovaSeq 6000 | [PMID: 40097733](https://pubmed.ncbi.nlm.nih.gov/40097733/) |
| [PRJNA847413](https://www.ebi.ac.uk/ena/browser/view/PRJNA847413) | [GSE205739](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205739) | 4 | PE | 2 × 150 bp | 23.3 M | 60.2 % | Illumina NovaSeq 6000 | _NA_ |
| [PRJEB48614](https://www.ebi.ac.uk/ena/browser/view/PRJEB48614)   | [E-MTAB-11129](https://www.ebi.ac.uk/biostudies/ArrayExpress/studies/E-MTAB-11129?query=E-MTAB-11129) | 3 | PE | 2 × 41 bp | 23.1 M | 85.9 % | Illumina NextSeq 500 | [PMID: 35967327](https://pubmed.ncbi.nlm.nih.gov/35967327/) |
| [PRJNA667281](https://www.ebi.ac.uk/ena/browser/view/PRJNA667281) | -- | 3 | PE | 2 × 150 bp | 22.1 M | 96.2 % | Illumina NovaSeq 6000 | [PMID: 33631268](https://pubmed.ncbi.nlm.nih.gov/33631268/) |
| [PRJNA896725](https://www.ebi.ac.uk/ena/browser/view/PRJNA896725) | -- | 5 | PE | 2 × 150 bp | 25.6 M | 94.1 % | Illumina NovaSeq 6000 | [PMID: 38638822](https://pubmed.ncbi.nlm.nih.gov/38638822/) |

### Excluded Studies
| ENA BioProject ID | Study Alias   | Reason for Exclusion       |
|:-----------------:|:-------------:|:--------------------------:|
| PRJNA802135       | GSE195781     | duplicated/bad runs        |
| PRJNA607654       | GSE145581     | miRNA-Seq                  |
| PRJNA1073892      | GSE255171     | regulatory T cells (Tregs) |
| PRJNA307651       | GSE76530      | miRNA-Seq                  |

#### The PRJNA802135/GSE195781 Case
Out of the 5 control Runs of PRJNA802135/GSE195781 (2022), 3 turned out to be ___identical___ to other runs already published in the previous study of the same group of authors (2016), under the different ID PRJNA307652/GSE76528.
Namely:

| PRJNA802135 (2022) | PRJNA307652 (2016) |
|:------------------:|:------------------:|
| SRR17833475        | SRR3085451         |
| SRR17833476        | SRR3085449         |
| SRR17833478        | SRR3085446         |

When comparing the corresponding FASTQ files between the two studies, the single reads actually appear to be identical, and the different file size only depends on the different pattern used for the heading line (_Field 1_) of each read.
```bash
# Check it out

zcat SRR17833478_1.fastq.gz | wc -l
zcat SRR3085446_1.fastq.gz | wc -l

zcat SRR17833478_1.fastq.gz | head
zcat SRR3085446_1.fastq.gz | head

zcat SRR17833478_1.fastq.gz | tail
zcat SRR3085446_1.fastq.gz | tail
```

This explains the batch effect that affects PRJNA802135 control samples (see PCA and hierarchical clustering), the abnormal standard deviation levels that are produced when all the samples are pooled together, and the _unlikely_ levels of correlation between the 2 datasets (especially when removing the 2 different-looking runs from PRJNA802135).

As for the 2 non-duplicated samples of PRJNA802135, they showed a very poor quality profile, such as severe problems of adapter contamination, eventually featuring percentages of alignment of just 20%.

__For all these reasons, PRJNA802135 study was completely purged from the global cohort.__
