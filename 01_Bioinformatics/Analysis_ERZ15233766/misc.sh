


awk -F'|' '/reverse:False/ {split($2,a,":"); printf(">%s_%d\n", a[2], NR)}' test.fa > test.fa
sed -i.bak 's/^>\(.*sample_id:\)\([^|]*\)\(_out.*\)$/> \2_'$(echo $((LINENO-1)))'/' test.fa

awk 'BEGIN {OFS="\t"} {sub("_out", "", $1); $2 = "/Volumes/Jacob_Agerbo_HD/NWS/00_DATA/"$1"_1.fastq.gz"; $3 = "/Volumes/Jacob_Agerbo_HD/NWS/00_DATA/"$1"_2.fastq.gz"; print}' soi.txt > samples.txt

anvi-search-primers --primers-txt primers.txt \
                               --samples samples.txt \
                               --only-report-primer-match \
                               -o PRIMER-MATCHES


cat PRIMER-MATCHES/*_VR1-PRIMER-MATCHES.fa > LysB_VR1-PRIMER-MATCHES.fa


oligotype LysB_VR1-PRIMER-MATCHES.fa \
          LysB_VR1-PRIMER-MATCHES.fa-ENTROPY \
          --number-of-auto-components 10 \
          --min-substantive-abundance 2 \
          --output-directory LysB_OLIGOTYPING








anvi-gen-variability-profile -c ../HoloFood_ERR4918746_bin.2_Photobacterium.db \
                             -p ../PROFILE.db \
                             --gene-caller-ids 3503 \
                             -o LysB_3503_AA.txt \
                             --engine AA \
                             --kiefl-mode




#

for i in D100 D10 D11 D22 D34 D45 D51 D73 D87 D8 D90 MG98
    do
    cp /Volumes/Jacob_Agerbo_HD/NWS/00_DATA/"$i"_* /Users/bfg522/Dropbox/Arbejde/PostDoc/PHOTOBACTERIUM_2/NWS_Photobacterium/HYPERVARIABILITY/
done    






anvi-gen-variability-profile -c HoloFood_ERR4918746_bin.2_Photobacterium.db \
                             -p PROFILE.db \
                             -C DEFAULT \
                             -b CORE \
                             --min-coverage-in-each-sample 10 \
                             --min-occurrence 3 \
                             --include-split-names \
                             --quince-mode \
                             -o PHOTO_NWS-SNVs.txt



anvi-gen-variability-profile -c HoloFood_ERR4918746_bin.2_Photobacterium.db \
                             -p PROFILE.db \
                             --genes-of-interest geNomad_prophage_gene_call.txt \
                             --min-coverage-in-each-sample 10 \
                             --min-occurrence 3 \
                             --include-split-names \
                             --quince-mode \
                             -o PHOTO_NWS-geNomad-SNVs.txt
