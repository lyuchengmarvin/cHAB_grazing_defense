# Functions for the purpose of pulling target gene sequences from PRODRIGAL gene fastas

## The function to retrieve one gene
retreive_sequence <- function(fasta, gene){
  seqs <- Biostrings::readDNAStringSet(fasta) %>% 
    data.frame(seq = ., seq_name = names(.)) %>% 
    unglue::unglue_unnest(seq_name, patterns = "{gene_name} {gene_annotation}") %>% 
    filter(gene_name == gene)
  return(seqs)
}
## The function to retrieve all genes from the gene table
pull_sequences <-function(gene_table){
  # gene_table must have columns: genes_fasta, Query_gene_name, strain
  # genes_fasta is the path to the fasta file containing the gene sequences
  # Query_gene_name is the name of the gene to be pulled from the fasta file
  # strain is the name of the strain from which the gene is pulled
  # The function returns a data frame with columns: gene_name, gene_annotation, seq, new_fasta_header
  # new_fasta_header is a string in the format "gene_name|strain"
  # Pull sequences
  seqs <- purrr::map2_dfr(gene_table$genes_fasta,
                  gene_table$Query_gene_name,
                  retreive_sequence,.progress=TRUE) %>% 
    # Add new fasta headers
    left_join(gene_table, by = c('gene_name'='Query_gene_name')) %>%
    mutate(new_fasta_header = stringr::str_glue("{gene_name}|{strain}"))
  return(seqs)
}
  
## The function to write the retrieved sequences
# inputs:
# seqs_w_new_headers: a data frame with the retrieved sequences and new fasta headers.
# output_dir: a directory to write the fasta files.
write_sequences <- function(seqs_w_new_headers, 
                            output_dir = "curated_BGC_sequences/",
                            output_file = "all_retrieved_sequences.fasta",
                            write_individual_fastas = TRUE){
  # Create output directory if it does not exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Write individual fasta files for each gene
  if (write_individual_fastas){
    gene_names <- levels(as.factor(seqs_w_new_headers$gene_name))
    # Write each gene to an individual fasta file
    for (g in gene_names){
      sub_gene_table <- subset(seqs_w_new_headers, gene_name == g)
      sub_gene_table$seq %>% 
        `names<-`(sub_gene_table$new_fasta_header) %>% 
        Biostrings::DNAStringSet() %>% 
        Biostrings::writeXStringSet(file.path(output_dir, paste0(g,".fasta")))
    }
  }else{
    # write all sequences to a single fasta file
    all_fasta_path <- file.path(output_dir, output_file)
    seqs_w_new_headers$seq %>% 
      `names<-`(seqs_w_new_headers$new_fasta_header) %>% 
      Biostrings::DNAStringSet() %>%
      Biostrings::writeXStringSet(.,filepath = all_fasta_path)

  }
}
  
