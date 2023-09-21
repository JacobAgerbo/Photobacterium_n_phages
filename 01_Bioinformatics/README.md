# Comparative Analysis of salmonid-related Photobacterium and its relation to phages

## Description of bioninformatics

### Publicly available genomes of Photobacterium
We recruited all publicly available genomes with the bacterial genera Photobacterium from the National Center for Biotechnology Information (NCBI), using the command line-based NCBI Datasets, resulting in 361 publicly available Photobacterium genomes. Furthermore, were metagenomic assembled genomes from recent investigation (REF -Holofish).

```{bash}
# Set the search term to "Photobacterium[Organism]"
search_term="Photobacterium[Organism]"

# Use the datasets command to search for genomes matching the search term
datasets search genome "$search_term" --output-file genome_summary.json

# Use the datasets command to download the genomes in FASTA format
datasets download genome --input-file genome_summary.json --exclude-rna --exclude-protein

# Extract the downloaded files to a directory
mkdir photobacterium_genomes
unzip "*.fasta.zip" -d photobacterium_genomes

# Clean up the downloaded zip files
rm *.fasta.zip
```

### Phylogenomic Photobacterium genomes
All genomes were compared with salmonid-related MAGs using anvi’o/v7.1 [1,2] for phylogenomic and comparative analysis. The phylogenomic analysis was carried out based on bacterial single-copy core genes (SCGs) using an anvi’o bacterial database with 71 bacterial SCGs. Amino acid sequences were extracted from Hidden Markov Models (HMMs) of SCGs and concatenated into aligned amino acid sequences. Concatenated amino acid sequences were used to generate a Newick-based maximum-likelihood phylogeny using FastTree2 [3].

So lets start making a anvi'o database for each genome and annotate with HMMs and functional annotations, like **COG**, **KEGG**, and **PFAM**.

```{bash}
THREADS=4
for i in `ls *fa | awk 'BEGIN{FS=".fa"}{print $1}'`
do
    anvi-gen-contigs-database -f $i.fa -o $i.db -T $THREADS
    anvi-run-hmms -c $i.db
    anvi-run-ncbi-cogs -c $i.db
    anvi-run-kegg-kofams -c $i.db
    anvi-run-pfams -c $i.db
done
```
Now that we have the anvi'o databases, lets do some bash-hacking and make the external-genomes.txt file, which are needed for the smooth anvi'o.

```{bash}
# Create an empty file
> external-genomes.txt

# Loop through all .db files in the current directory
for i in *.fa; do
  # Extract filename without extension
  filename=$(basename "$i" .db)
  # Get absolute path of the file
  path=$(realpath "$i")
  # Append filename and path to external-genomes.txt
  echo -e "$filename\t$path" >> external-genomes.txt
done
```

Now we generate a concatenated fasta-file with protein calls from our HMMs. We use all genes in the Bacterial Single-copy core gene (SCGs) database and make a FastTree based on the concatenated fasta file.

```{bash}
anvi-get-sequences-for-hmm-hits --external-genomes external-genomes.txt \
                                -o concatenated-proteins.fa \
                                --hmm-source Bacteria_71 \
                                --return-best-hit \
                                --get-aa-sequences \
                                --concatenate

anvi-gen-phylogenomic-tree -f concatenated-proteins.fa \
                           -o phylogenomic-tree.txt                                
```     




### Comparative genomics of Photobacterium genomes
The comparative analysis was carried out as in the previous studies [4,5], where similarities of each amino acid sequence in every genome were calculated against every other amino acid sequence across all genomes using BLASTp. We implemented Minbit heuristics of 0.5 to eliminate weak matches between two amino acid sequences [6] and an MCL inflation of 2. We used the MCL algorithms to identify gene clusters in amino acid sequence similarity [7]. We calculated ANI using PyANI [8]. Euclidean distance and ward linkage were used to organise gene clusters and genomes. A summary of the pan-genome generated for this study is available [here]().

### Metabolic reconstruction of Photobacterium
Metabolic reconstruction of compared MAGs was based on KOfams and was carried out using the anvi’o platform. We calculated the level of completeness for a given KEGG module [9,10] in our genomes using the programme anvi-estimate-metabolism, which leveraged the previous annotation of genes with KEGG orthologs (KOs). The URL https://merenlab.org/m/anvi-estimate-metabolism serves as a tutorial for this program which details the modes of usage and output file formats. The Heatmap of completion scores was illustrated using the ComplexHeatmap [11] package for R. Genomes were clustered based on similarity across the completion of pathways, as previously done for other bacterial genera [4].

### Enrichment analysis of KOfams 
The statistical approach for enrichment analysis is previously defined [12]. Briefly, the programme anvi-compute-functional-enrichment determined enrichment scores for KOfams genomes of salmonid-related Photobacterium and non-salmonid-related Photobacterium by fitting a binomial generalised linear model (GLM) to the occurrence of each KOfam in each group and then computing a Rao test statistic. We considered any KOfam with a q-value less than 0.05 to be ‘enriched’ in its associated group. The volcano plot was visualised using the EnhancedVolcano package for R.

**1.** 	Eren AM, Kiefl E, Shaiber A, Veseli I, Miller SE, Schechter MS, et al. Community-led, integrated, reproducible multi-omics with anvi’o. Nat Microbiol. 2021;6: 3–6.

**2.** 	Murat Eren A, Esen ÖC, Quince C, Vineis JH, Morrison HG, Sogin ML, et al. Anvi’o: an advanced analysis and visualization platform for ‘omics data. PeerJ. 2015;3: e1319.

**3.** 	Price MN, Dehal PS, Arkin AP. FastTree 2--approximately maximum-likelihood trees for large alignments. PLoS One. 2010;5: e9490.

**4.** 	Rasmussen JA, Kiilerich P, Madhun AS, Waagbø R, Lock E-JR, Madsen L, et al. Co-diversification of an intestinal Mycoplasma and its salmonid host. ISME J. 2023. doi:10.1038/s41396-023-01379-z

**5.**  Rasmussen JA, Villumsen KR, Duchêne DA, Puetz LC, Delmont TO, Sveier H, et al. Genome-resolved metagenomics suggests a mutualistic relationship between Mycoplasma and salmonid hosts. Commun Biol. 2021;4: 579.

**6.** 	Benedict MN, Henriksen JR, Metcalf WW, Whitaker RJ, Price ND. ITEP: an integrated toolkit for exploration of microbial pan-genomes. BMC Genomics. 2014;15: 8.

**7.** 	van Dongen S, Abreu-Goodger C. Using MCL to extract clusters from networks. Methods Mol Biol. 2012;804: 281–295.

**8.** 	Pritchard L, Glover RH, Humphris S, Elphinstone JG, Toth IK. Genomics and taxonomy in diagnostics for food security: soft-rotting enterobacterial plant pathogens. Anal Methods. 2015;8: 12–24.

**9.** 	Kanehisa M, Goto S, Sato Y, Kawashima M, Furumichi M, Tanabe M. Data, information, knowledge and principle: back to metabolism in KEGG. Nucleic Acids Res. 2014;42: D199–205.

**10.** 	Kanehisa M, Furumichi M, Tanabe M, Sato Y, Morishima K. KEGG: new perspectives on genomes, pathways, diseases and drugs. Nucleic Acids Res. 2017;45: D353–D361.

**11.** 	Gu Z, Eils R, Schlesner M. Complex heatmaps reveal patterns and correlations in multidimensional genomic data. Bioinformatics. 2016;32: 2847–2849.

**12.** 	Shaiber A, Willis AD, Delmont TO, Roux S, Chen L-X, Schmid AC, et al. Functional and genetic markers of niche partitioning among enigmatic members of the human oral microbiome. Genome Biol. 2020;21: 292.