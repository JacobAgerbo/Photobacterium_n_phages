#!/bin/sh
for file in $1
do
    echo "generating gene call file to anvi'o for: "$file""
    file1=""$file".fixed_provirus.tsv"
    file2=""$file"_gene_calls.txt"
    echo -e "gene_callers_id\tsource\taccession\tfunction\te_value" > split_info_geNomad_functions_"$file"
    while IFS=$'\t' read -r contig start stop marker annotation e_value; do
    awk -F'\t' -v contig="$contig" -v start="$start" -v stop="$stop" -v marker="$marker" -v annotation="$annotation" -v e_value="$e_value" \
        '$2 == contig && $3 > start && $4 <= stop {print $1 "\t" "geNomad" "\t" "geNomad" "\t" "prophage" "\t" e_value}' \
        "$file2" >> geNomad_prophage_"$file"
    done < "$file1"
done