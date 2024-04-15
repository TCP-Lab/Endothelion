#!/bin/bash

utils_path="$(dirname "$(realpath "$0")")"
test_path="${utils_path}/../../data/in/test/"

mkdir -p "$test_path"

wget -nc -P "$test_path" https://zenodo.org/records/10960146/files/GSE138309_CountMatrix_genes_TPM.tsv
wget -nc -P "$test_path" https://zenodo.org/records/10960146/files/GSE138309_meta.csv
wget -nc -P "$test_path" https://zenodo.org/records/10960146/files/GSE139133_CountMatrix_genes_TPM.tsv
wget -nc -P "$test_path" https://zenodo.org/records/10960146/files/GSE139133_meta.csv

head -n 1000 GSE138309_CountMatrix_genes_TPM.tsv > test_GSE138309_CountMatrix_genes_TPM.tsv
head -n 1000 GSE139133_CountMatrix_genes_TPM.tsv > test_GSE139133_CountMatrix_genes_TPM.tsv

head -n 10 ICT_set.csv > test_ICT_set.csv 

