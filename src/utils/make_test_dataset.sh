#!/bin/bash

# This Bash script can be used to generate a lightweight data set suitable for
# testing changes or new features of Endothelion pipelines.

# Retrieve paths
utils_path="$(dirname "$(realpath "$0")")"
endo_path="$(echo "$utils_path" | sed 's/\(Endothelion\).*/\1/')"
target_path="${endo_path}/data/in/Lines/test_hCMEC_D3"
gois_path="${endo_path}/data/in/ICT_set.csv"

# Ask for user permission to clean the target directory when already present
if [[ -d "$target_path" ]]; then
    printf "\nTarget directory .${target_path#${endo_path}} already exists!"
    printf "\nThis procedure will delete all files possibly present in it.\n\n"
    valid_ans=false
    while ! $valid_ans; do
        read -ep "Proceed anyway (Y/n)? " ans
        no_rgx="^[Nn][Oo]?$"
        yes_rgx="^[Yy](ES|es)?$"
        if [[ $ans =~ $no_rgx ]]; then
            printf "Aborting test set creation... nothing changed.\n"
            exit 1
        elif [[ $ans =~ $yes_rgx || -z "$ans" ]]; then
            printf "\n"
            valid_ans=true
        else
            printf "Invalid answer '$ans'\n"
        fi
    done
fi

# Make the target directory (if it doesn't exist yet) or clean it if exists.
mkdir -p "$target_path"
rm "${target_path}"/*

# Download two Endothelion count matrices from Zenodo and make them small
projs=("GSE138309" "GSE139133")
for proj in "${projs[@]}"; do
	wget -O "${target_path}/${proj}_CountMatrix_genes_TPM.tsv" \
		https://zenodo.org/records/10960146/files/${proj}_CountMatrix_genes_TPM.tsv
	wget -O "${target_path}/${proj}_meta.csv" \
		https://zenodo.org/records/10960146/files/${proj}_meta.csv

	Rscript "${endo_path}/src/utils/make_test_dataset.R" \
		"$target_path" "$gois_path" "${proj}_CountMatrix_genes_TPM.tsv"
    rm "${target_path}/${proj}_CountMatrix_genes_TPM.tsv"

	mv "${target_path}/${proj}_meta.csv" "${target_path}/mini${proj}_meta.csv"
done
