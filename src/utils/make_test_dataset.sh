#!/bin/bash

utils_path="$(dirname "$(realpath "$0")")"
endo_path="$(echo "$utils_path" | sed 's/\(Endothelion\).*/\1/')"
target_path="${endo_path}/data/in/Lines/hCMEC_D3"
gois_path="${endo_path}/data/in/ICT_set.csv"

mkdir -p "$target_path"

projs=("GSE138309" "GSE139133")
for proj in "${projs[@]}"; do
	wget -O "${target_path}/${proj}_CountMatrix_genes_TPM.tsv" \
		https://zenodo.org/records/10960146/files/${proj}_CountMatrix_genes_TPM.tsv
	wget -O "${target_path}/${proj}_meta.csv" \
		https://zenodo.org/records/10960146/files/${proj}_meta.csv

	Rscript "${endo_path}/src/utils/make_test_dataset.R" \
		"$target_path" "$gois_path" $proj
done
