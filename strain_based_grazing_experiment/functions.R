library(tidyverse)
library(ggpubr)
library(glue)
library(patchwork)
#Importing data from excel function
read_data <- function(file_list) {
  read_excel(file_list, 
             col_types = c("date", "numeric", "numeric", 
                           "numeric", "numeric", "numeric", 
                           "text","text","text","text")) |>
    na.omit()
}
#Importing data from excel function
cc_import <- function(file_list) {
  
  # Assuming metadata is in the first 36 lines of the file
  metadata_raw <- read_csv(file_list,n_max= 37, col_names = c("property","value"))
  
  wide_metadata <- metadata_raw |>
    mutate(property = str_to_lower(property),
           property = str_replace_all(property,":", ""),
           property = str_replace_all(property," +", "_")) |>
    as_tibble() |>
    pivot_wider(names_from = property, values_from = value)  
  
  # Read the rest of the data minus the header rows
  headers_raw <- read_csv(file_list, skip = 81,n_max= 3, col_names = FALSE)
  
  # Transpose the header rows to process them by columns
  headers_transposed <- t(headers_raw) |>
    tolower() |>
    as_tibble() |>
    mutate(V1 = case_when(V1 == "diff." ~ '', 
                          V1 == "bin number" ~ "bin_no",
                          str_detect(V1, "lower") ~ "diam_lb",
                          str_detect(V1, "upper") ~ "diam_ub",
                          TRUE ~ V1)) |>
    # Combine the header rows, and process them column-wise
    pmap_chr(
      ~{
        # Remove whitespace from each piece
        elements <- trimws(c(...))
        # Remove NA values and trim whitespace
        elements <- trimws(elements[!is.na(elements) & elements != ""])
        # Concatenate the non-empty elements, convert to lowercase, and replace spaces with underscores
        combined <- tolower(paste(elements, collapse = "_"))
        # Replace any leftover whitespace or multiple underscores with a single underscore
        combined <- gsub(" +", "_", combined)
        combined <- gsub("__+", "_", combined)
        combined
      }
    )
  
  
  cc_data <- read_csv(file_list, skip = 85, col_names = FALSE)
  colnames(cc_data) <- headers_transposed
  
  all_data <- cbind(wide_metadata ,cc_data)
  return(all_data)
}

cc_import2 <- function(file_list) {
  
  # Assuming metadata is in the first 36 lines of the file
  metadata_raw <- read_csv(file_list,n_max= 35, col_names = c("property","value"))
  
  wide_metadata <- metadata_raw |>
    mutate(property = str_to_lower(property),
           property = str_replace_all(property,":", ""),
           property = str_replace_all(property," +", "_")) |>
    as_tibble() |>
    pivot_wider(names_from = property, values_from = value)  
  
  # Read the rest of the data minus the header rows
  headers_raw <- read_csv(file_list, skip = 35,n_max= 3, col_names = FALSE)
  
  # Transpose the header rows to process them by columns
  headers_transposed <- t(headers_raw) |>
    tolower() |>
    as_tibble() |>
    mutate(V1 = case_when(V1 == "diff." ~ '', 
                          V1 == "bin number" ~ "bin_no",
                          str_detect(V1, "lower") ~ "diam_lb",
                          str_detect(V1, "upper") ~ "diam_ub",
                          TRUE ~ V1)) |>
    # Combine the header rows, and process them column-wise
    pmap_chr(
      ~{
        # Remove whitespace from each piece
        elements <- trimws(c(...))
        # Remove NA values and trim whitespace
        elements <- trimws(elements[!is.na(elements) & elements != ""])
        # Concatenate the non-empty elements, convert to lowercase, and replace spaces with underscores
        combined <- tolower(paste(elements, collapse = "_"))
        # Replace any leftover whitespace or multiple underscores with a single underscore
        combined <- gsub(" +", "_", combined)
        combined <- gsub("__+", "_", combined)
        combined
      }
    )
  
  
  cc_data <- read_csv(file_list, skip = 39, col_names = FALSE)
  colnames(cc_data) <- headers_transposed
  
  all_data <- cbind(wide_metadata ,cc_data)
  return(all_data)
}

#Plotting fluorescence change during a mussel experiment and outputting the file as a png
plot_fp_data <- function(split_df) {
  for (i in 1:length(split_df)) {
    print( ggplot(data = split_df[[i]], aes(x = hours, y = avg, shape = source, color = fluor)) +
             geom_point(data = filter(split_df[[i]], Timepoint != "TF"),size = 2, 
                        position = position_dodge(width = 0.025)) +
             geom_point(data = filter(split_df[[i]], Timepoint == "TF"),size = 2,
                        position = position_dodge(width = 0.025)) + 
             stat_smooth(data = filter(split_df[[i]], Timepoint != "TF"), method = "lm", se = F) +
             scale_shape_manual("",
                                values = c("Mixed" = 17,
                                           "Water Column" = 16),
                                labels = c("blugrn_avg" = "Bluegreen",
                                           "grn_avg" = "Green Algae")) +
             scale_color_manual("", values = c("blugrn" = "royalblue1",
                                               "grn" = "lightgreen"),
                                labels = c("blugrn" = "Bluegreen",
                                           "grn" = "Green Algae")) +
             geom_errorbar(mapping = aes(ymin = avg - sd,
                                         ymax = avg + sd), width = 0.025,
                           position = position_dodge(width = 0.025)) +
             scale_x_continuous(breaks = c(0,0.25,0.50,0.75,1,1.5,2)) +
             xlab("Time Point (Hours)") +
             ylab("Chlorophyll (μg/L)") +
             facet_wrap(facets = "Treatment") +
             theme_pubclean() +
             theme(legend.position = "right") )
    
    ggsave(glue("./figures/mf_results/raw_fp_results/","{names(split_df[i])}",".png"), width = 180,
           height = 200, units = 'mm', dpi = 450)
  }
}
## Function for plotting trait data
plot_mussel_trait <- function(plot_data, xname, yname, std_fr, trait){
  xname <- sym(xname)
  yname <- sym(yname)
  # Scatterplot of green and blugreen
  t_space <- ggplot(plot_data, aes(x = !!xname, y = !!yname)) +
    geom_point(aes(fill = Clade), shape = 21, size = 3,color = "grey50") +
    # fill color by clade
    # customize colors
    scale_fill_manual(values = c('tan2', 'lightskyblue','gold','cadetblue',
                                 'tomato1','plum1','darkseagreen2','grey30')) +
    # plot 1:1 line
    geom_abline(intercept = 0, slope = 1, color = 'grey50',linetype = 'dashed') +
    # plot the standard feeding rate (green-only control)
    geom_vline(xintercept = std_fr, color = 'grey50', linetype = 'dotted')+
    geom_vline(xintercept = std_fr, color = 'grey', linetype = 'dotdash')+
    # put label on point and separate from each other
    ggrepel::geom_text_repel(aes(label = Culture), 
                             size = 2, alpha = 0.5,fontface = 'bold') +
    #geom_text(aes(label = Culture), size = 1.5, hjust = 1, vjust = -0.5) +
    scale_x_continuous(limits = c(-15, 65)) +
    scale_y_continuous(limits = c(-15, 40)) +
    labs(title = 'Trait Space for IDM resistance',
         x = 'Clearance Rate on Green',
         y = 'Clearance Rate on Blue-Green',
         caption = paste0('Calculation: ',trait)) +
    theme_bw()+
    theme(plot.title= element_text(size=15,color='black',hjust=0.6),
          legend.title = element_text(size=10,color='black',hjust=0.5),
          legend.position = 'left',
          panel.grid.major = element_blank())
  # add marginal density
  t_plot_den <- ggExtra::ggMarginal(t_space, type = "density", fill = "grey", size = 5)
  t_plot_den
}
#Plotting fluorescence change during a mussel experiment and outputting the file as a png
plot_cc_data <- function(split_df) {
  # for (i in 1:length(split_df)) {
  #   print(ggplot(split_df[[i]], aes(x = mean_size,y = log10(number_per_ml))) +
  #            geom_col(position= "stack") +
  #            scale_x_continuous(breaks = c(2,4,6,8,10,15,20,30,40,50,60)) +
  #            facet_wrap(~media) +
  #            labs(title = glue("Size Distribution by Media Type ", "{names(split_df[i])}"),
  #                 x = "Size (μm)",
  #                 y = "Log 10 Counts per mL") +
  #            theme_pubclean(base_size = 12) +
  #            theme(axis.title.x = ggtext::element_markdown(),
  #                  legend.position = 'bottom',
  #                  axis.text.x = element_text(margin = margin(r = 0)),
  #                  plot.title = ggtext::element_markdown(),
  #                  legend.text = ggtext::element_markdown(),
  #                  panel.border = element_rect(fill = NA, colour = "black", 
  #                                              linewidth = 0.7),
  #                  strip.background = element_rect(fill = 'grey90', colour = "black", 
  #                                                  linewidth = 0.7),
  #                  axis.line = element_blank()) +
  #            rotate_x_text(angle = 90))
  #   
  #   ggsave(glue("./figures/size_results/","{names(split_df[i])}",".png"), width = 180,
  #          height = 200, units = 'mm', dpi = 450)
  # }
  for (i in 1:length(split_df)) {
    
    count_plot <- ggplot(split_df[[i]], aes(x = mean_size,y = number_percent)) +
            geom_col(position= "identity", fill = "coral") +
            scale_x_continuous(breaks = c(2,4,6,8,10,15,20,30,40,50,60)) +
            facet_grid(~media) +
            labs(title = glue("Size Distribution by Media Type: ", "{names(split_df[i])}"),
                 subtitle = "Number",
                 x = "Particle Diameter (μm)",
                 y = "Number Percent") +
            ylim(c(0,25)) +
            theme_pubclean(base_size = 12) +
            theme(axis.title.x = ggtext::element_markdown(),
                  legend.position = 'bottom',
                  axis.text.x = element_text(margin = margin(r = 0)),
                  plot.title = ggtext::element_markdown(),
                  legend.text = ggtext::element_markdown(),
                  panel.border = element_rect(fill = NA, colour = "black", 
                                              linewidth = 0.7),
                  strip.background = element_rect(fill = 'grey90', colour = "black", 
                                                  linewidth = 0.7),
                  axis.line = element_blank()) +
            rotate_x_text(angle = 90)

    surfacearea_plot <- ggplot(split_df[[i]], aes(x = mean_size,y = surface_area_percent)) +
      geom_col(position= "identity", fill = "limegreen") +
      scale_x_continuous(breaks = c(2,4,6,8,10,15,20,30,40,50,60)) +
      facet_grid(~media) +
      labs(subtitle = "Surface Area",
           x = "Particle Diameter (μm)",
           y = "Surface Area Percent") +
      ylim(c(0,10)) +
      theme_pubclean(base_size = 12) +
      theme(axis.title.x = ggtext::element_markdown(),
            legend.position = 'bottom',
            axis.text.x = element_text(margin = margin(r = 0)),
            plot.title = ggtext::element_markdown(),
            legend.text = ggtext::element_markdown(),
            panel.border = element_rect(fill = NA, colour = "black", 
                                        linewidth = 0.7),
            strip.background = element_rect(fill = 'grey90', colour = "black", 
                                            linewidth = 0.7),
            axis.line = element_blank()) +
            rotate_x_text(angle = 90)

    volume_plot <- ggplot(split_df[[i]], aes(x = mean_size,y = volume_percent)) +
      geom_col(position= "identity", fill = "plum") +
      scale_x_continuous(breaks = c(2,4,6,8,10,15,20,30,40,50,60)) +
      facet_grid(~media) +
      labs(subtitle = "Biovolume",
           x = "Particle Diameter (μm)",
           y = "Volume Percent") +
      ylim(c(0,15)) +
      theme_pubclean(base_size = 12) +
      theme(axis.title.x = ggtext::element_markdown(),
            legend.position = 'bottom',
            axis.text.x = element_text(margin = margin(r = 0)),
            plot.title = ggtext::element_markdown(),
            legend.text = ggtext::element_markdown(),
            panel.border = element_rect(fill = NA, colour = "black", 
                                        linewidth = 0.7),
            strip.background = element_rect(fill = 'grey90', colour = "black", 
                                            linewidth = 0.7),
            axis.line = element_blank()) +
            rotate_x_text(angle = 90)

    plot <- (count_plot / surfacearea_plot / volume_plot)
    
    print(plot)
    ggsave(glue("./figures/media_results/","{names(split_df[i])}_percent",".png"), width = 180,
           height = 200, units = 'mm', dpi = 450)
  }
  
}

#Plotting fluorescence change during a mussel experiment and outputting the file as a png
plot_cc_exp <- function(split_df) {
  
  for (i in 1:length(split_df)) {
    
    count_plot <- ggplot(split_df[[i]], aes(x = mean_size,y = number_percent)) +
      geom_col(position= "identity", fill = "coral") +
      scale_x_continuous(breaks = c(1,2,4,6,8,10,15,20)) +
      ylim(c(0,15)) +
      labs(title = glue("Size Characteristics of ", "{names(split_df[i])}"),
           subtitle = "Number",
           x = "Particle Diameter (μm)",
           y = "Number (%)") +
      theme_pubclean(base_size = 12) +
      theme(axis.title.x = ggtext::element_markdown(),
            legend.position = 'bottom',
            axis.text.x = element_text(margin = margin(r = 0)),
            plot.title = ggtext::element_markdown(),
            legend.text = ggtext::element_markdown(),
            panel.border = element_rect(fill = NA, colour = "black", 
                                        linewidth = 0.7),
            strip.background = element_rect(fill = 'grey90', colour = "black", 
                                            linewidth = 0.7),
            axis.line = element_blank()) +
      rotate_x_text(angle = 90)
    
    surfacearea_plot <- ggplot(split_df[[i]], aes(x = mean_size,y = surface_area_percent)) +
      geom_col(position= "identity", fill = "limegreen") +
      scale_x_continuous(breaks = c(1,2,4,6,8,10,15,20)) +
      ylim(c(0,15)) +
      labs(subtitle = "Surface Area",
           x = "Particle Diameter (μm)",
           y = "Surface Area (%)") +
      theme_pubclean(base_size = 12) +
      theme(axis.title.x = ggtext::element_markdown(),
            legend.position = 'bottom',
            axis.text.x = element_text(margin = margin(r = 0)),
            plot.title = ggtext::element_markdown(),
            legend.text = ggtext::element_markdown(),
            panel.border = element_rect(fill = NA, colour = "black", 
                                        linewidth = 0.7),
            strip.background = element_rect(fill = 'grey90', colour = "black", 
                                            linewidth = 0.7),
            axis.line = element_blank()) +
      rotate_x_text(angle = 90)
    
    volume_plot <- ggplot(split_df[[i]], aes(x = mean_size,y = volume_percent)) +
      geom_col(position= "identity", fill = "plum") +
      scale_x_continuous(breaks = c(1,2,4,6,8,10,15,20)) +
      ylim(c(0,15)) +
      labs(subtitle = "Biovolume",
           x = "Particle Diameter (μm)",
           y = "Biovolume (%)") +
      theme_pubclean(base_size = 12) +
      theme(axis.title.x = ggtext::element_markdown(),
            legend.position = 'bottom',
            axis.text.x = element_text(margin = margin(r = 0)),
            plot.title = ggtext::element_markdown(),
            legend.text = ggtext::element_markdown(),
            panel.border = element_rect(fill = NA, colour = "black", 
                                        linewidth = 0.7),
            strip.background = element_rect(fill = 'grey90', colour = "black", 
                                            linewidth = 0.7),
            axis.line = element_blank()) +
      rotate_x_text(angle = 90)
    
    plot <- (count_plot / surfacearea_plot / volume_plot)
    
    print(plot)
    ggsave(glue("./figures/size_results/","{names(split_df[i])}_percent",".png"), width = 180,
           height = 200, units = 'mm', dpi = 450)
  }
  
}

