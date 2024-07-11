# Validation workflow (checklist)

1. check prefixes in folder _Counts_: `tree Counts/`
1. check prefix in _CountMatrix_genes_TPM.tsv_ file
1. check prefixes in folder _MultiQC_out_: `tree MultiQC_out/`
1. add PCA report: `qcfastq --tool=PCA .`
1. for quality control, download:
	- MultiQC report
		`scp fear@130.192.101.241:/home/fear/WORKS/Endothelion/Lines/hCMEC_D3/GSE138309/MultiQC_out/GSExxxxxxx_multiqc_report.html .`
	- PCA report
		`scp fear@130.192.101.241:/home/fear/WORKS/Endothelion/Lines/hCMEC_D3/GSE138309/PCA_out/GSExxxxxxx_CountMatrix_genes_TPM/* .`
1. download metadata again and manually edit the `extra` field using _nano_
	(set 0 to exclude a Run)
	`metaharvest -x=1 -e GSExxxxxxx > GSExxxxxxx_meta.csv`
1. delete FASTQ files
	`rm *.fastq.gz`
1. make the TAR GZ compressed archive from the parent directory
	`tar -czvf ./GSExxxxxxx/GSExxxxxxx.tar.gz GSExxxxxxx`
1. locally download files to upload to Zenodo later
	`scp fear@130.192.101.241:/home/fear/WORKS/Endothelion/Lines/hCMEC_D3/GSExxxxxxx/GSExxxxxxx* .`
