# Utility Scripts and Functions

## The ___prefixer___ patch
The pre-release version of x.FASTQ (until _x.FASTQ_ Ver.Sum x.13.62.0, April 11,
2024) didn't use any suffix ID for STAR/RSEM (i.e., _anqFASTQ_ module) output
files, thus preventing _MultiQC_ from including them in the final report.
`_prefixer` function can be used as a patch for such _x.FASTQ_ old output by
prepending the unique ENA-SRR run ID to every file and folder found within the
target directory, including possible subdirectories. Specifically, from each GEO
project directory, e.g.,
```bash
cd ~/WORKS/Endothelion/Lines/hCMEC_D3/GSE205739
```
1. I removed the whole _MultiQC_ output folder and the related log file, as well
	as the assembled expression matrix `Count_Matrix_genes_TPM.tsv` along with
	its log file.
	```bash
	rm -rf MultiQC_out/ Z_QC_MultiQC_GSE* Count_Matrix_genes_TPM.tsv Z_Counts_GSE*
	```
1. Then I moved to the `Count` subdirectory and run the `_prefixer` function
	over each ENA-SRR sub-directory (the function is supposed to have been
	_sourced_ beforehand).
	```bash
	 cd Counts
	 
	 _prefixer SRR19592963
	 _prefixer SRR19592964
	 _prefixer SRR19592965
	 _prefixer SRR19592966

	 cd ..
	```
	I also tried to automate the process using `find`
	```bash 
	find ./Counts/ -maxdepth 1 -type d -name "SRR*" -exec bash -c '_prefixer "$1"' _ {} \;
	```
	but for some reason it didn't work...
1. Having assigned a different name to each different alignment, I could rerun
	_MultiQC_ to include also _RSEM_ and _STAR_ output logs in the final report.
	```bash
	qcfastq --tool=MultiQC .
	```
1. Then I regenerated the assembled expression matrix by running the latest
	`assembler.R` version to have GSE-prefixed count matrix (i.e.,
	`GSE205739_CountMatrix_genes_TPM.tsv`).
	```bash
	countfastq -n .
	```
1. Also I added the metadata, as retrieved from ENA database, with the `extra`
	column filled with `1` by default, which I replaced with `0` for those
	non-control samples not to be included in Endothelion reanalysis.
	```bash
	metaharvest -x=1 -e GSE205739 > GSE205739_meta.csv
	```
1. After having checked that everything was OK, I could remove the (trimmed)
	FASTQ still present in the folder to free disk space.
	```bash
	rm *.fastq.gz
	```
1. Then, I moved one level up and made a _gzip_ compressed archive of the
	whole GEO project directory that I placed inside the project directory
	itself.
	```bash
	cd ..
	tar -czvf ./GSE205739/GSE205739.tar.gz GSE205739
	```
1. Finally, I downloaded the _tar.gz_ archive, the count matrix, and metadata on
	my local machine for subsequent uploading to _Zenodo_.
	```bash
	# From my local machine (actually this command downloaded also the
	# `GSE205739_wgets.sh` that I deleted before sending stuff to Zenodo).
	scp fear@130.192.101.241:/home/fear/WORKS/Endothelion/Lines/hCMEC_D3/GSE205739/GSE205739* .
	```

## Generate a test dataset
`make_test_dataset.sh` Bash script can be used to generate a lightweight dataset
suitable for testing changes or new features of _Endothelion_ pipelines. From
the project directory, just run
```bash
./src/utils/make_test_dataset.sh
```
and the script will
1. make the target directory `./data/in/Lines/hCMEC_D` (if it doesn't exist yet)
	or, upon user confirmation, clean it by removing possible files already
	present therein;
1. download 2 _Endothelion_ count matrices (and related metadata) from _Zenodo_,
	namely GEO series ___GSE138309___ and ___GSE139133___ featuring just 3 and 2
	control samples, respectively;
1. further resize them row-wise by keeping only the GOI entries and 5000 more
	random genes (actually, this last step is performed by the R script of the
	same name).
1. save both these lightweight count matrices and a copy of the related metadata
	after prepending a `test_` string to their names (according to the
	`[data.profiles.test]` settings in `kerblam.toml`).

__NOTE 1:__
Since some of the 5000 random genes could be already present among the GOIs, the
final length of the count matrices cannot be predicted exactly. However, a seed
is set at the beginning of the R script (`set.seed(7)`) to ensure
reproducibility.

__NOTE 2:__
When running the _hCMEC_D3_ (or _hCMEC_D3x_) pipeline using the test data set,
namely
```bash
kerblam run --profile test hCMEC_D3
```
the warnings _Bad TPM normalization in series..._ will be thrown. Clearly this
is expected and perfectly normal, given the row-wise reduction of the count
matrices.

__NOTE 3:__
For both count matrices, a further warning will be issued during the pipeline
run, namely
```
WARNING:
 Can't find these Genes of Interest in the Count Matrix:
  CACNA1C
  ANO5
```
Again, this is expected and perfectly normal, since count matrices were
intersected with the list of GOIS during the row-reduction step _without_
exploding their _Gene Symbol_ entries (as the `endo_profiler.R` script does),
thus loosing the `ANO5,LOC102723370` and `CACNA1C,CACNA1C-IT2` genes. I kept
this fortuitous feature to check `endo_profiler.R` capability of detecting
possible missing GOIs among the genes included in the count matrix. 
