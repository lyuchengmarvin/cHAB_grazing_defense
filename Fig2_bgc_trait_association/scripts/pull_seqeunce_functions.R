# Title: Retrieve sequences from a pangenome table

## The function to retrieve one gene
retreive_sequence <- function(fasta, gene){
  seqs <- Biostrings::readDNAStringSet(fasta) %>% 
    data.frame(seq = ., seq_name = names(.)) %>% 
    unglue::unglue_unnest(seq_name, patterns = "{gene_name} {gene_annotation}") %>% 
    filter(gene_name == gene)
}

## The function to retrieve all genes from the pangenome
# inputs: 
# pangenome_tab: a pangenome table with the first column as gene names, and the rest of the columns as the accession numbers of the genes in each genome.
# path_gene_anno: a path to the directory with the annotated gene fasta files.
pull_seq_from_pangenome <- function(pangenome_tab, path_gene_anno = '~/Desktop/LE_annotation/bakta/'){
  # Compile fasta file paths
  gene_table <- pangenome_tab  %>%
    pivot_longer(-Gene,names_to = "genome",values_to = "gene_name") %>%
    filter(!is.na(gene_name)) %>%
    separate_longer_delim(gene_name, delim = ";") %>% 
    mutate(genes_fasta = stringr::str_glue(paste0(path_gene_anno,"{genome}.ffn")))
  # Check if the fasta files exist
  gene_table <- gene_table %>%
    mutate(annotation_files_exist = fs::file_exists(genes_fasta))
  # Retrieve sequences
  seqs <- purrr::map2_dfr(gene_table$genes_fasta,gene_table$gene_name,retreive_sequence,.progress=TRUE)
  # Add new fasta headers
  seqs_w_new_headers <- seqs %>% 
    left_join(gene_table) %>% 
    mutate(new_fasta_header = stringr::str_glue("{Gene}|{genome}"))
  
  return(seqs_w_new_headers)
}

## The function to write the retrieved sequences
# inputs:
# seqs_w_new_headers: a data frame with the retrieved sequences and new fasta headers.
# output_dir: a directory to write the fasta files.
write_sequences <- function(seqs_w_new_headers, output_dir = "~/Desktop/LE_annotation/pangenome_sequences/"){
  # Create output directory if it does not exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Write sequences to fasta files
  gene_names <- levels(as.factor(seqs_w_new_headers$Gene))
  
  for (g in gene_names){
    sub_gene_table <- subset(seqs_w_new_headers, Gene == g)
    sub_gene_table$seq %>% 
      `names<-`(sub_gene_table$new_fasta_header) %>% 
      Biostrings::DNAStringSet() %>% 
      Biostrings::writeXStringSet(file.path(output_dir, paste0(g,".fasta")))
  }
}