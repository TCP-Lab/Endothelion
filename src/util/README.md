# Utility Scripts and Functions

## The prefixer patch

From each GEO project directory, e.g.,
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
	over each ENA-SRR sub-directory
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
	`assembler.R` (v.1.6.0) to have GSE-prefixed count matrix (i.e.,
	`GSE205739_CountMatrix_genes_TPM.tsv`).
	```bash
	countfastq -n .
	```
1. Also I added the metadata, as retrieved from ENA database, with the `extra`
	column filled with `1` by default, which I replaced with `0` for those
	non-control samples not included in the Endothelion reanalysis.
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
