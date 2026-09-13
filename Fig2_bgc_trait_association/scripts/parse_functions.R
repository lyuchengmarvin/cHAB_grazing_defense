## Author: Yu-Cheng (Marvin) Lin
########################################
# Functions to process annotation files
parse_single_anno_result <- function(annotation_dir, annotation_key){
  # Test if in csv format
  if (!grepl('\\.csv$', annotation_dir)){
    stop("The annotation file must be in csv format.")
  }
  # Find strain name
  strain_name <- basename(annotation_dir) %>%
    sub('_lauren_blastn_table.csv$', '', .)
  # Load annotation results and add filter
  bgc_lauren <- read.csv(annotation_dir, header= T) %>%
    # because there will be multiple hits due to fragmented genes, caused by contig breaks, we need to get the ones that has the highest query coverage of the same hit.
    group_by(Reference_Name) %>%
    slice_max(order_by = Query_Coverage, n = 1, with_ties = FALSE) %>%
    filter(Reference_Name!='Unmatched') 
  # Combine annotation with BGC key
  bgc_lauren_w_key <- annotation_key %>%
    left_join(bgc_lauren, by = c('Header' = 'Reference_Name')) %>%
    mutate(Strain = strain_name) %>%
    unglue::unglue_unnest(Query_Name, patterns = "{Query_gene_name} {Query_gene_annotation}") %>%
    ungroup() %>%
    select(Strain, Query_gene_name, Query_gene_annotation, New_BGC_name, Gene, Identity_Score, Query_Coverage, gene_type, Main_Product, Biosynthetic_Class) 
  
  # Return the combined annotation
  return(bgc_lauren_w_key)
}

parse_multiple_anno_results <- function(annotation_dirs, annotation_key){
  # annotation_dirs: a vector of annotation file paths
  # annotation_key: a data frame containing the BGC key
  
  # Initialize an empty list to store results
  annot_list <- list()
  
  # Loop through each annotation file and parse
  for (dir in annotation_dirs){
    annot_result <- parse_single_anno_result(dir, annotation_key)
    annot_list[[dir]] <- annot_result
  }
  
  # Combine all annotation results into a single data frame
  combined_annot <- bind_rows(annot_list)

  # Return the combined annotation
  return(combined_annot)
}

# Function to read blast result from targeted analysis
read_seq_identity <- function(file_path){
  seq_idt_df <- read.table(file_path, header = T) %>%
    group_by(Strain_Name) %>%
    mutate(pres_absence = ifelse(Reference_Name == 'Unmatched', 0, 1)) %>%
    mutate(perc_present = mean(pres_absence)) %>%
    # replace unmatched with NA for summarizing
    mutate(Reference_Name = ifelse(Reference_Name == 'Unmatched', NA, Reference_Name)) %>%
    mutate(Identity_Score = ifelse(Reference_Name == 'Unmatched', NA, Identity_Score),
           Query_Coverage = ifelse(Reference_Name == 'Unmatched', NA, Query_Coverage)) %>%
    # dont consider NA when summarizing
    summarise(ANI = mean(Identity_Score, na.rm = TRUE)/100,
              ACOV = mean(Query_Coverage, na.rm = TRUE)/100,
              Percent_presence = mean(perc_present, na.rm = TRUE)) %>%
    # Give 0 values to the NAs
    replace_na(list(ANI = 0, ACOV = 0, Percent_presence = 0)) %>%
    # Rename some strain names for consistency
    mutate(Strain = str_replace(Strain_Name, 'LE18-131', 'LE19-131')) %>%
    mutate(Strain = str_replace(Strain_Name, 'Microcystis_aeruginosa_LE3', 'LE3'))
  return(seq_idt_df)
}

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
  # gene_table must have columns: genes_fasta, Query_gene_name, Strain
  # genes_fasta is the path to the fasta file containing the gene sequences
  # Query_gene_name is the name of the gene to be pulled from the fasta file
  # Strain is the name of the strain from which the gene is pulled
  # The function returns a data frame with columns: gene_name, gene_annotation, seq, new_fasta_header
  # new_fasta_header is a string in the format "gene_name|Strain"
  # Pull sequences
  seqs <- purrr::map2_dfr(gene_table$genes_fasta,
                  gene_table$Query_gene_name,
                  retreive_sequence,.progress=TRUE) %>% 
    # Add new fasta headers
    left_join(gene_table, by = c('gene_name'='Query_gene_name')) %>%
    mutate(new_fasta_header = stringr::str_glue("{gene_name}|{Strain}|{Gene}"))
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

# Functions to facilitate plotting
## Function to plot identity score and query coverage variation
plot_variation <- function(data, title_prefix, if_save=TRUE){
  data_num <- data %>%
    drop_na() %>%
    group_by(Gene) %>%
    summarise(n = n())
  
  p_idt <- ggplot(data, aes(x=Gene, y=Identity_Score)) +
    geom_boxplot(box.color='grey70', fill='lavenderblush3', alpha=0.7, outlier.shape = NA) +
    geom_jitter(width=0.1, alpha=0.3) +
    labs(title=paste(title_prefix, 'Gene Identity Score Variation'),
         x='Gene',
         y='Identity Score (%)') +
    geom_text(data=data_num, aes(x=Gene, y=105, label=n), inherit.aes = FALSE) +
    theme_pubclean() +
    theme(plot.title= element_text(size=12,color='black',hjust=0.5,face = "bold"),
          axis.line = element_line(size=0, colour = "black"),
          axis.ticks = element_line(size=1, colour = "black"),
          axis.text.y = element_text(hjust = 0,size=6,color='black',face = "bold"),
          axis.text.x = element_text(hjust = 0.5,size=6,color='black',face = "bold",angle = 90),
          axis.title = element_text(size=15,color='black',face = "bold"),
          panel.background = element_rect(fill = "grey100"),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
          legend.title = element_text(size=10,color='black',hjust=0.5, face = "bold"),
          legend.position = 'right',
          panel.grid.major = element_blank())
  
  
  p_cov <- ggplot(data, aes(x=Gene, y=Query_Coverage)) +
    geom_boxplot(box.color='grey70', fill='lavenderblush3', alpha=0.7, outlier.shape = NA) +
    geom_jitter(width=0.1, alpha=0.3) +
    labs(title=paste(title_prefix, 'Gene Coverage Variation'),
         x='Gene',
         y='Query Coverage (%)') +
    geom_text(data=data_num, aes(x=Gene, y=105, label=n), inherit.aes = FALSE) +
    theme_pubclean() +
    theme(plot.title= element_text(size=12,color='black',hjust=0.5,face = "bold"),
          axis.line = element_line(size=0, colour = "black"),
          axis.ticks = element_line(size=1, colour = "black"),
          axis.text.y = element_text(hjust = 0,size=6,color='black',face = "bold"),
          axis.text.x = element_text(hjust = 0.5,size=6,color='black',face = "bold",angle = 90),
          axis.title = element_text(size=15,color='black',face = "bold"),
          panel.background = element_rect(fill = "grey100"),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
          legend.title = element_text(size=10,color='black',hjust=0.5, face = "bold"),
          legend.position = 'right',
          panel.grid.major = element_blank())
  if (if_save){
    ggsave(filename = paste0(title_prefix, '_Identity_Score_Variation.pdf'),
           plot = p_idt,
           device = 'pdf',
           width = max(6, n_distinct(data$Gene)*0.6),
           height = 5)
    ggsave(filename = paste0(title_prefix, '_Query_Coverage_Variation.pdf'),
                  plot = p_cov,
                  device = 'pdf',
                  width = max(6, n_distinct(data$Gene)*0.6),
                  height = 5)
  }
  
  print(p_idt)
  print(p_cov)
}
