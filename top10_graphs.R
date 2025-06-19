#Rscript to plot top 5 / 10 most abundant taxa for each sampler

library(ggplot2)    
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(scales)  

# Set plotting theme
theme_set(theme_bw())

#sample meta data (this hasn't changed)
meta <- read.csv("old_parameters_MARTI_samp_comp_read_data/Phyloseq_data/Sample_table.csv")

# Marti summarised data ---------------------------------------------------

# Same MARTi data but exported June 2025 (just summed counts)
marti <- read.delim("samp_comp_0624_marti/marti_assignments_lca_0.1_all_levels_2025-JUN-18_9-27-47.tsv")

# Clean up column headers
n_cols <- ncol(marti)

# Clean only from column 4 onwards
if (n_cols > 3) {
  new_colnames <- colnames(marti)
  new_colnames[4:n_cols] <- sub("\\.\\.samp_comp.*", "", new_colnames[4:n_cols])
  colnames(marti) <- new_colnames
}
marti_sum <- read.csv("../samp_comp_0624_marti/samp_comp_summed_0624.csv")


marti_long <- pivot_longer(marti, cols = starts_with(c("CF_", "NHM_")),
                           names_to = "Sample_ID",
                           values_to = "Count")

marti_long$Count <- as.numeric(marti_long$Count)

marti_meta <- merge(marti_long, meta)
marti_meta$NumReads <- as.numeric(gsub(",","",marti_meta$NumReads)) #Remove commas from the NumRead col
marti_meta <- marti_meta %>% 
  mutate(
    percent_classified_read = Count / NumReads * 100,
    HPM = Count / NumReads * 1e6 #hits per million
  )

# Plotting aesthetics ----
#Colour palette - want consistency between graphs
Tol_muted <- c('#88CCEE', '#44AA99', '#117733', '#332288', '#DDCC77', '#999933','#CC6677', '#882255', '#AA4499', '#DDDDDD')

custom_theme <- theme_minimal(base_size = 12) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.3),
    panel.grid.minor = element_blank(),
  )

sampler_labels <- c(
  "Compact" = "Coriolis\nCompact",
  "Micro" = "Coriolis\nMicro",
  "Bobcat" = "InnovaPrep\nBobcat",
  "Cub" = "InnovaPrep\nCub",
  "Sass" = "SASS\n3100"
)

sampler_levels <- c(
  "Compact",
  "Micro",
  "Bobcat",
  "Cub",
  "Sass"
)

location_labels <- c(
  "NHM" = "Natural History Museum",
  "Cfarm" = "Church Farm"
)


duration_labels <- c(
  "25" = "25 minutes",
  "50" = "50 minutes"
)


# Functions ---------
# Function to extract top n taxa
get_top_taxa <- function(data, rank, top_n = 5) {
  data %>%
    filter(NCBI.Rank == rank) %>%
    group_by(Sample_ID) %>%
    arrange(desc(HPM), .by_group = TRUE) %>%
    slice_head(n = top_n) %>%
    filter(Count != 0) %>%
    ungroup()
}

# Function to plot all samplers on one graph 
plot_grouped_top_taxa <- function(data, rank, color_mapping, output_path){
  
  # Ensure proper sampler ordering
  data$Sampler <- factor(data$Sampler, levels = sampler_levels)
  
  # Capitalise rank for title and filename
  rank_clean <- paste0(toupper(substring(rank, 1, 1)), substring(rank, 2))
  
  # Plot
  p <- ggplot(data, aes(x = Sampler, y = HPM, fill = Name)) +
    geom_bar(stat = "identity") +
    facet_grid(rows = vars(Location), cols = vars(Sample_length),
               labeller = labeller(
                 Location = location_labels,
                 Sample_length = duration_labels
               )) +
    labs(
      x = "Sample",
      y = "Hits per Million"
    ) +
    scale_fill_manual(values = color_mapping, na.value = "grey80") +
    scale_y_continuous(
      labels = label_number(scale_cut = cut_short_scale(), accuracy = 1)
    ) +
    scale_x_discrete(labels = sampler_labels) +
    custom_theme
  
  # Save the plot
  ggsave(output_path, plot = p, width = 12, height = 8, dpi = 300)
}


# To plot samplers individually ----
plot_top_taxa_by_sampler <- function(data, samplers, rank, color_mapping, output_dir = ".", top_n = 5) {
  
  for (sampler in samplers) {
    
    # Filter and process data
    sampler_data <- data %>%
      filter(Sampler == sampler, NCBI.Rank == rank) %>%
      group_by(Sample_ID) %>%
      arrange(desc(HPM), .by_group = TRUE) %>%
      slice_head(n = top_n) %>%
      filter(Count != 0) %>%
      ungroup()
    
    # Plot
    p <- ggplot(sampler_data, aes(x = Sample_ID, y = HPM, fill = Name)) +
      geom_bar(stat = "identity") +
      facet_grid(~Location, scales = "free_x",
                 labeller = labeller(Location = location_labels)) +
      labs(
        x = "Sample",
        y = "Hits per Million"
      ) +
      scale_y_continuous(
        labels = label_number(scale_cut = cut_short_scale(), accuracy = 1)) +  # Option 1: short format with whole numbers
      custom_theme +
      theme(
        legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1)
      ) +
      scale_fill_manual(values = color_mapping)
    
    # Capitalise first letter of rank for filename
    rank_clean <- paste0(toupper(substring(rank, 1, 1)), substring(rank, 2))
    
    # Save plot
    file_name <- file.path(output_dir, paste0(rank_clean, "_Top", top_n, "_", sampler, ".png"))
    ggsave(file_name, plot = p, width = 12, height = 8)
  }
}


#Phylum ---------------
marti_phylum <- get_top_taxa(marti_meta, "phylum", top_n = 5)

# Extract unique Phylum values
unique_phylums <- unique(marti_phylum$Name) 

# Create a named vector
phy_color_mapping <- setNames(Tol_muted[1:length(unique_phylums)], unique_phylums)

# Call the function to plot all samplers
plot_grouped_top_taxa(
  data = marti_phylum,
  rank = "phylum",
  color_mapping = phy_color_mapping,
  output_path = "Images/graphs_marti_0624/top5/Phylum_Top5_AllSamp.png"
)


# Phylum plot per sampler 
plot_top_taxa_by_sampler(
  data = marti_meta,
  samplers = c("Compact", "Micro", "Bobcat", "Cub", "Sass"),
  rank = "phylum",
  color_mapping = phy_color_mapping,
  output_dir = "Images/graphs_marti_0624/top5"
)

# Genus -------------------------------------------------------------------
marti_genus <- get_top_taxa(marti_meta, "genus", top_n = 5)

# 29 Unique genera
unique_genus <- unique(marti_genus$Name)
uniq_genus_df <- as.data.frame(unique_genus)

# Create a color palette with 29 colors
gen_color_mapping <- colorRampPalette(brewer.pal(12, "Paired"))(29)

# Create a named colour vector
color_mapping_gen <- setNames(colour_palette_29[1:length(unique_genus)], unique_genus)

plot_grouped_top_taxa(
  data = marti_genus,
  rank = "genus",
  color_mapping = gen_color_mapping,
  output_path = "Images/graphs_marti_0624/top5/Genus_Top5_AllSamp.png"
)

plot_top_taxa_by_sampler(
  data = marti_meta,
  samplers = c("Compact", "Micro", "Bobcat", "Cub", "Sass"),
  rank = "genus",
  color_mapping = gen_color_mapping,
  output_dir = "Images/graphs_marti_0624/top5"
) 


#Species ----------------------------------

marti_species <- get_top_taxa(marti_meta, "species", top_n = 5)

# 33 unique species - not that diff to genus level 
unique_species <- unique(marti_species$Name)
unique_species_df <- as.data.frame(unique_species)

# Create a color palette with 33 colors
colour_palette_33 <- colorRampPalette(brewer.pal(12, "Paired"))(33)

# Create a named vector
spe_color_mapping <- setNames(colour_palette_33[1:length(unique_species)], unique_species)

plot_grouped_top_taxa(
  data = marti_species,
  rank = "species",
  color_mapping = spe_color_mapping,
  output_path = "Images/graphs_marti_0624/top5/Species_Top5_AllSamp.png"
)

plot_top_taxa_by_sampler(
  data = marti_meta,
  samplers = c("Compact", "Micro", "Bobcat", "Cub", "Sass"),
  rank = "species",
  color_mapping = spe_color_mapping,
  output_dir = "Images/graphs_marti_0624/top5"
) 







# Older code --- PHYLOSEQ data -----------------------------------------------------------
#These datasets are created in the Updated_phyloseq_analysis.R script

#Reading in data ----------
top_species <- read.csv("MARTI_samp_comp_updated_parameters_0823/top10_species_2508.csv")%>% 
  select(-X)

top_phylum <- read.csv("MARTI_samp_comp_updated_parameters_0823/top10_phylum_2508.csv") %>% 
  select(-X)

#sample meta data (this hasn't changed)
meta <- read.csv("old_parameters_MARTI_samp_comp_read_data/Phyloseq_data/Sample_table.csv")

#Combining
top_phylum_meta <- merge(top_phylum, meta)

#Top 10 Phylum graphs from phyloseq --------------------------------
#With faceting 
ggplot(top_phylum_meta, aes(x=Sample_ID, y = abundance, fill = Phylum)) +
  geom_bar(stat="identity") +
  facet_grid(~Sampler, scales ="free_x") +
  labs(title = "Top 10 Phylum per sample",
       x = "Sample",
       y = "Relative Abundance") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1))

#Bobcat
top_phylum_bobcat <- top_phylum_meta %>% 
  filter(Sampler == "Bobcat")

ggplot(top_phylum_bobcat, aes(x=Sample_ID, y = abundance, fill = Phylum)) + 
  geom_bar(stat="identity") +
  facet_grid(~Location, scales ="free_x") +
  labs(title = "Top 10 Phylum Bobcat",
       x = "Sample",
       y = "Relative Abundance") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1))

#Compact
top_phylum_Compact <- top_phylum_meta %>% 
  filter(Sampler == "Compact")

ggplot(top_phylum_Compact, aes(x=Sample_ID, y = abundance, fill = Phylum)) + 
  geom_bar(stat="identity") +
  facet_grid(~Location, scales ="free_x") +
  labs(title = "Top 10 Phylum Compact",
       x = "Sample",
       y = "Relative Abundance") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1))

#Micro
top_phylum_Micro <- top_phylum_meta %>% 
  filter(Sampler == "Micro")

ggplot(top_phylum_Micro, aes(x=Sample_ID, y = abundance, fill = Phylum)) + 
  geom_bar(stat="identity") +
  facet_grid(~Location, scales ="free_x") +
  labs(title = "Top 10 Phylum Micro",
       x = "Sample",
       y = "Relative Abundance") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1))

#Cub
top_phylum_Cub <- top_phylum_meta %>% 
  filter(Sampler == "Cub")

ggplot(top_phylum_Cub, aes(x=Sample_ID, y = abundance, fill = Phylum)) + 
  geom_bar(stat="identity") +
  facet_grid(~Location, scales ="free_x") +
  labs(title = "Top 10 Phylum Cub",
       x = "Sample",
       y = "Relative Abundance") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1))

#SASS
top_phylum_Sass <- top_phylum_meta %>% 
  filter(Sampler == "Sass")

ggplot(top_phylum_Sass, aes(x=Sample_ID, y = abundance, fill = Phylum)) + 
  geom_bar(stat="identity") +
  facet_grid(~Location, scales ="free_x") +
  labs(title = "Top 10 Phylum Sass",
       x = "Sample",
       y = "Relative Abundance") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1))







