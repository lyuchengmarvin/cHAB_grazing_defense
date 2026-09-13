#########################################
## Author: Yu-Cheng Lin
## Date: 2026.02.08
## This is the workflow to create input files for the estimation of BGC effect sizes on defense traits using multiple linear regression. This workflow curated reference BGC sequences for a Microcystis strain and calculated tbe sequence similarity and the presence of referernce genes in the rest of the tested strains. 
#### BGC gene-trait correlation analysis

**Predict biosynthetic gene cluster genes with antiSMASH:**

1. Annotate your reference genomes with [antiSMASH](https://antismash.secondarymetabolites.org/#!/start)

2. Download GeneBank files and extract fasta gene sequences gene-by-gene with custome script:
- script: `get_antiSMASH_fasta.sh`
- output: an umbrella fasta file with gene sequences from all predicted BGCs `LE19-84_antiSMASH_BGC_genes.fasta `
```{command-line}
## Change parameters in the script: `antiSMASH_gbk_dir` & `antiSMASH_fasta_output`
bash ../../scripts/get_antiSMASH_fasta.sh
```

**Annotate reference genotypes with a curated non-redundant BGC database (Hart et al., 2026):**

To get the reference BGC sequences, We annotated the gene sequences predicted by PRODRIGAL from the genomes of the reference strain. Because there are some discrepancy between antiSMASH predicted result and Hart's BGC database. We combined the references from both databases.

1. Annotate the genes with both databases
```{command-line}
# This script iterates through each BLAST record and finds the best hit based on identity score and query coverage.
## Annotate with Hart's database
python ../../scripts/blast_targets.py -m blastn -d BGCdb_seed5.fa -q ~/Desktop/LE_annotation/bakta/LE19-84.ffn -c 70 -p blast_result/LE19-84_lauren
## Annotate the antiSMASH predicted result pulled from the genebank files using the script `get_antiSMASH_fasta.sh`
python ../../scripts/blast_targets.py -m blastn -d antiSMASH_results/extracted_fna/LE19-84_antiSMASH_BGC_genes.fasta -q  ~/Desktop/LE_annotation/bakta/LE19-84.ffn -c 70 -p blast_result/LE19-84_antiSMASH
```
Remember to manually convert the tab-delimited tables to csv so that the R script can process them. Outputing to csv file often mess up the columns because there are many other symbols in the annotation field.

2. Find consistency
- Combine the results into one table `LE1984_anno_combined.txt` with the R scrip: `curate_BGC_annotation.Rmd`.
- output the combined annotation: `LE1984_bgc_annotation_combined.csv`. We used Hart's database as the basis to organize antiSMASH predicted result because antiSMASH sometimes separate genes from the same clusters due to contig breaks. We used use the annotation from Lauren_db as the core genes for a given BGC, and throw in gene sequences predicted by antiSMASH only as additional. BGC clusters that are too fragmented to judge if they come from the same famliy will be categorized as miscallaneous. For example, there are cases where a few genes were blasted to microcystin reference but lacking too many genes to be considered a full cluster. These genes are also not in consequitive numbers, indicating scattering around in the genome, or the assembly is simply too fragemented to capture.
- Note: We observed that many BGCs were on the edge of a contig and can be annotated to the same cluster in the database. One of the caveat of this study is that some genes at the edge of the cutoff will not be annotated. Take Aeruginosin aerB for example in LE19-84, its on the edge of ctg225 and ctg 266.

**BLAST reference BGC genes to the LE genomes:**

1. Place all the manually curated BGC fasta files in each of its BGC folder under the reference strain directory. These are the reference BGC genotypes. Perform BLAST search for them in LE culture collection PRODRIGAL gene annotation .ffn files.

```command-line
## Before running, change the file paths and output names in the script. "$strain_name", "$path_to_reference" and "$path_to_db_fasta_dir". The script will iterate over all the BGC folders under the strain directory.
bash ../scripts/blast_ref_bgcgeno.sh
```

2. For gene-trait correlation analysis, move all `targets_summary.txt` files under each `blast_results/`.

**Highly correlated genes in the pangenome:**

To check if evolutionary history mediates the BGC-trait correlation, this section check gene-trait correlation and genetic distance based on the core genes.

1. The Microcystis pangenome (Kiledal et al., 2024) was built by Panaroo and annotated by Bakta, which can be accessed through this [link](https://zenodo.org/records/7847321). Anders used the criteria of at least presenting in 90% of the genomes to account for the reference incompletness. However, I want to get the most represetative sequences so I used 98%. Core gene seqeunces from the reference strains of four phenotypically representative groups will be pulled separately. 
- script: `pull_genes_from_pangenomes.Rmd`

2. Move the reference core gene sequences to the `reference_core_genes/` directory and run blast search for these core genes in the LE culture collection.
```command-line
bash ../scripts/blast_ref_coregenes.sh
```
3. For the core gene summary table files
- Remove '|strain_name' in the database column.
- Format: Reference genotype,Gene ID,gene name
- core gene summary tables `LE19-84_core_gene_table.csv` is made from copying from the blast result files.

**Estimate the effect size for each BGCs**

See `figS9_effect_size.Rmd`.