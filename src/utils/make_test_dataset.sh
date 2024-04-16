#!/bin/bash

utils_path="$(dirname "$(realpath "$0")")"
endo_path="$(echo "$utils_path" | sed 's/\(Endothelion\).*/\1/')"
target_path="${endo_path}/data/in/Lines/hCMEC_D3"
gois_path="${endo_path}/data/in/ICT_set.csv"

printf "\nThis procedure will delete all files possibly present in"
printf "\n  .${target_path#${endo_path}}\n\n"
valid_ans=false
while ! $valid_ans; do
    read -ep "Proceed anyway (Y/n)? " ans
    no_rgx="^[Nn][Oo]?$"
    yes_rgx="^[Yy](ES|es)?$"
    if [[ $ans =~ $no_rgx ]]; then
        printf "  Aborting test set creation... nothing changed.\n"
        exit 1
    elif [[ $ans =~ $yes_rgx || -z "$ans" ]]; then
        printf "\n"
        valid_ans=true
    else
        printf "  Invalid answer '$ans'\n"
    fi
done

mkdir -p "$target_path"
rm "${target_path}"/*

projs=("GSE138309" "GSE139133")
for proj in "${projs[@]}"; do
	wget -O "${target_path}/${proj}_CountMatrix_genes_TPM.tsv" \
		https://zenodo.org/records/10960146/files/${proj}_CountMatrix_genes_TPM.tsv
	wget -O "${target_path}/${proj}_meta.csv" \
		https://zenodo.org/records/10960146/files/${proj}_meta.csv

	Rscript "${endo_path}/src/utils/make_test_dataset.R" \
		"$target_path" "$gois_path" $proj

	cp "${target_path}/${proj}_meta.csv" "${target_path}/test_${proj}_meta.csv"
done
