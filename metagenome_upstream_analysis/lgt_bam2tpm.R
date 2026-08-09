## This script is used to parse bam files and calculate TPM values for each gene.
library('paletteer')
library('stringr')
library('tidyverse')
library('purrr')
library('Rsamtools')

## Define function to read bam file and calculate TPM values
bam2TPM = function(metadata) {
  bam_count_std = data.frame()
  for (sampleid in 1:nrow(metadata)) {
    bam.counts = purrr::map_dfr(metadata[sampleid, 'Dir'], ~ Rsamtools::idxstatsBam(.x))
    std.tbl = bam.counts %>%
      mutate(
        per_million = sum(metadata[sampleid, 'total_reads']) / 1000000,
        rpm = mapped / per_million,
        rpkm = rpm / seqlength,
        SampleID = as.character(metadata[sampleid, 'sample']),
        lake = as.character(metadata[sampleid, 'lake']),
        Date = as.character(metadata[sampleid, 'exp_date']),
        Treatment = as.character(metadata[sampleid, 'treatment'])
      ) %>%
      mutate(
        total_rpkm = sum(rpkm, na.rm = T),
        rpkm_relative_abundance = rpkm / total_rpkm * 100
      ) %>%
      mutate(length_kb = seqlength / 1000, RPK = mapped / length_kb) %>%
      group_by(SampleID) %>%
      mutate(sum_RPK_per_million = sum(RPK, na.rm = T) / 1e6) %>%
      ungroup() %>%
      mutate(TPM = RPK / sum_RPK_per_million) %>%
      mutate(TPM_relative_abundance = TPM / sum(TPM, na.rm = T) * 100) %>%
      select(
        SampleID,
        Date,
        Treatment,
        lake,
        seqnames,
        seqlength,
        mapped,
        rpkm,
        rpkm_relative_abundance,
        TPM,
        TPM_relative_abundance
      )
    
    bam_count_std = rbind(bam_count_std, std.tbl)
  }
  return(bam_count_std)
}

## Main
setwd('/nfs/turbo/lsa-vdeneflab/vdenef1-mistorage/Shared/lyucheng/nikesh_metG_assembly/graz_resist_genes')
# load metadata
metadata = read.table('data/ND_samples_metadata',header = T,sep = '\t') %>%
  filter(!is.na(total_reads)) %>%
  mutate(bam_dir_prefix = paste0("3003-", sample, "_"))
# define bam directory
bam_dir = '/nfs/turbo/lsa-vdeneflab/vdenef1-mistorage/Shared/lyucheng/nikesh_metG_assembly/graz_resist_genes/data/recruitment_map/just_genes/'
# list bam files that match sample IDs in metadata and add to metadata
bam_files = list.files(bam_dir, pattern = "lgt_full_recruitment.id99.bam$", full.names = TRUE)
# Match each prefix in metadata to the corresponding file in the file list
metadata = metadata %>%
  mutate(Dir = map_chr(bam_dir_prefix, ~ bam_files[str_detect(bam_files, .x)]))
# make TPM table
filt99cov50 <- bam2TPM(metadata)
write.csv(filt99cov50, "data/ND_samples_lgt_full_id99.mapped.counts.csv",quote = F, row.names = F)
