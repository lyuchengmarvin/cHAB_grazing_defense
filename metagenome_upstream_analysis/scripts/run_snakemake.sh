#!/bin/bash
#SBATCH --job-name=inStrain_pangenome
#SBATCH --mail-user=lyucheng@umich.edu
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=03-00:00:00
#SBATCH --account=vdenef0
#SBATCH --partition=standard
#SBATCH --export=ALL
#SBATCH --error=error-%x-%j.log
#SBATCH --output=output-%x-%j.log

# Initialize conda for shell interaction
source ~/miniconda3/etc/profile.d/conda.sh

# Load conda environment
conda activate run_metagenome

# Run Snakemake under our working directory
cd /scratch/vdenef_root/vdenef0/lyucheng/ND_feed1819_metagenome

# Make sure that snakemake is not locked
snakemake -s snakefile --unlock

# Submit Snakemake job to SLURM

## Run the initial fastqc
# snakemake --profile ./config/greatlakes \
#     --conda-frontend conda \
#     -s snakefile run_fastqc

## De-duplicate reads and trim adapters
# snakemake --profile ./config/greatlakes \
#     --conda-frontend conda \
#     -s snakefile run_fastp

## Remove contaminants
# snakemake --profile ./config/greatlakes \
#     --conda-frontend conda \
#     -s snakefile run_human_read_removal

## Run the post-processing fastqc
# snakemake --profile ./config/greatlakes \
#     --conda-frontend conda \
#     -s snakefile run_post_fastqc


## Run the recruitment mapping script
# snakemake --profile ./config/greatlakes \
#     --conda-frontend conda \
#     -s downstream_snakefile run_recruitment_mapping

## Filter bam files 
# snakemake --profile ./config/greatlakes \
#     --conda-frontend conda \
#     --rerun-incomplete \
#     -s downstream_snakefile run_filterBam

## Get reference allele frequency
# snakemake --profile ./config/greatlakes \
#     --conda-frontend conda \
#     --rerun-incomplete \
#     -s downstream_snakefile run_get_ref_freq

## Variant calling on the recruitment mapping results
# snakemake --profile ./config/greatlakes \
#     --conda-frontend conda \
#     --rerun-incomplete \
#     -s downstream_snakefile run_call_variants

## Run inStrain profiling on the recruitment mapping results
snakemake --profile ./config/greatlakes \
    --conda-frontend conda \
    --rerun-incomplete \
    -s downstream_snakefile run_inStrain_profile

# Contruct biosynthetic gene clusters for the WLECC culture metegenome samples
## Run biosyntheticSPAdes assembly
# snakemake --profile ./config/greatlakes \
#     --conda-frontend conda \
#     -s snakefile run_assemble_biosyntheticSPAdes