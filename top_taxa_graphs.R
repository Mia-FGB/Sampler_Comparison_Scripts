#Rscript to plot top 5 / 10 most abundant taxa for each sampler

library(ggplot2)    
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(scales)  

# Load in data -------
#sample meta data (this hasn't changed)
meta <- read.csv("metadata/Sample_table.csv")

# Marti summarised data - Exported Jul 25, MinReadLength 150, Min Identity 85

marti <- read.delim("marti_outputs/marti_Jul_150_v2/summed_marti_assignments_lca_0.1_all_levels_2025-JUL-9_11-52-34.tsv")

# Clean up column headers
n_cols <- ncol(marti)

# Clean only from column 4 onwards
if (n_cols > 3) {
  new_colnames <- colnames(marti)
  new_colnames[4:n_cols] <- sub("\\.\\.samp_comp.*", "", new_colnames[4:n_cols])
  colnames(marti) <- new_colnames
}

marti_long <- pivot_longer(marti, cols = starts_with(c("CF_", "NHM_")),
                           names_to = "Sample_ID",
                           values_to = "Count")

marti_long$Count <- as.numeric(marti_long$Count)


## Merge on metadata and calculate ------

marti_meta <- merge(marti_long, meta)
marti_meta$NumReads <- as.numeric(gsub(",","",marti_meta$NumReads)) #Remove commas from the NumRead col
marti_meta <- marti_meta %>% 
  mutate(
    percent_classified_read = Count / NumReads * 100,
    HPM = Count / NumReads * 1e6,  #hits per million
    HP100k = Count / NumReads * 100000
  )

#rename sampler
marti_meta$Sampler <- gsub("Micro", "μ", marti_meta$Sampler)
marti_meta$Sample_ID <- gsub("Micro", "μ", marti_meta$Sample_ID)

# Filter the data on read count and HPM
filtered_marti <- marti_meta  %>% 
  filter(HPM >= 100, Count > 5)


# Plotting aesthetics ----
#Colour palette - want consistency between graphs where possible
palette <- c('#8B4513', '#9e0142', '#d53e4f', '#f46d43',
                       '#fdae61','#fee08b', '#e6f598', '#abdda4',
                       '#66c2a5','#3288bd', '#5e4fa2', '#A672A7')

custom_theme <- theme_minimal(base_size = 12) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.3),
    panel.grid.minor = element_blank(),
  )

sampler_labels <- c(
  "Compact" = "Coriolis\nCompact",
  "μ" = "Coriolis μ",
  "Bobcat" = "InnovaPrep\nBobcat",
  "Cub" = "InnovaPrep\nCub",
  "Sass" = "SASS\n3100"
)

sampler_levels <- c(
  "Compact",
  "μ",
  "Bobcat",
  "Cub",
  "Sass"
)

sampler_filename_map <- c(
  "Bobcat" = "Bobcat",
  "Compact" = "Compact",
  "μ" = "Mu",       # Use "Mu" in filenames
  "Cub" = "Cub",
  "Sass" = "Sass"
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
  p <- ggplot(data, aes(x = Sampler, y = HP100k, fill = Name)) +
    geom_bar(stat = "identity") +
    facet_grid(rows = vars(Location), cols = vars(Sample_length),
               labeller = labeller(
                 Location = location_labels,
                 Sample_length = duration_labels
               )) +
    labs(
      x = "Sample",
      y = "Hits per 100,000"
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
    p <- ggplot(sampler_data, aes(x = Sample_ID, y = HP100k, fill = Name)) +
      geom_bar(stat = "identity") +
      facet_grid(~Location, scales = "free_x",
                 labeller = labeller(Location = location_labels)) +
      labs(
        x = "Sample",
        y = "Hits per 100,000"
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
    
    sampler_safe <- sampler_filename_map[[sampler]]  # gets "Mu" if sampler is "μ"
    
    # Save plot
    file_name <- file.path(output_dir, paste0(rank_clean, "_Top", top_n, "_", sampler_safe, ".png"))
    ggsave(file_name, plot = p, width = 12, height = 8)
  }
}

# Plot with replicates not combined
plot_grouped_top_taxa_replicates <- function(data, rank, color_mapping, output_path) {
  
  # Ensure proper sampler order
  data$Sampler <- factor(data$Sampler, levels = sampler_levels)
  
  # Create replicate-aware x-axis label
  data <- data %>%
    mutate(
      Sample_Label = paste(Sampler, Repeat, sep = "_"),
      Sample_Label = factor(Sample_Label, levels = unique(Sample_Label[order(Sampler, Repeat)]))
    )
  
  # Capitalise rank for title/filename
  rank_clean <- paste0(toupper(substring(rank, 1, 1)), substring(rank, 2))
  
  # Plot
  p <- ggplot(data, aes(x = Sample_Label, y = HP100k, fill = Name)) +
    geom_bar(stat = "identity") +
    facet_grid(rows = vars(Location), cols = vars(Sample_length),
               labeller = labeller(
                 Location = location_labels,
                 Sample_length = duration_labels
               )) +
    labs(
      x = "Sample replicate",
      y = "Hits per 100,000"
    ) +
    scale_fill_manual(values = color_mapping, na.value = "grey80") +
    scale_y_continuous(
      labels = label_number(scale_cut = cut_short_scale(), accuracy = 1)
    ) +
    custom_theme +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.title.x = element_blank(),
      plot.margin = margin(10, 10, 30, 10)
    ) 
  
  # Save the plot
  ggsave(output_path, plot = p, width = 12, height = 8, dpi = 300)
}


# Calling the functions to plot top taxa ------

##Phylum ---------------
marti_phylum <- get_top_taxa(filtered_marti, "phylum", top_n = 5)

# Extract unique Phylum values - 11
unique_phylums <- unique(marti_phylum$Name) 

# Create a named vector
phy_color_mapping <- c('#8B4513', '#9e0142', '#d53e4f', '#f46d43',
                       '#fdae61','#fee08b', '#abdda4',
                       '#66c2a5','#3288bd', '#A672A7')

# Call the function to plot all samplers
plot_grouped_top_taxa(
  data = marti_phylum,
  rank = "phylum",
  color_mapping = phy_color_mapping,
  output_path = "Images/top5/Phylum_Top5_AllSamp.png"
)



# Phylum plot per sampler - calculates top 5 indpeendetly for each sampler
plot_top_taxa_by_sampler(
  data = filtered_marti,
  samplers = c("Compact", "μ", "Bobcat", "Cub", "Sass"),
  rank = "phylum",
  color_mapping = phy_color_mapping,
  output_dir = "Images/top5"
)


## Genus ---------------------------------------------------
marti_genus <- get_top_taxa(filtered_marti, "genus", top_n = 5)

# 32 Unique genera
unique_genus <- unique(marti_genus$Name)
uniq_genus_df <- as.data.frame(unique_genus)

# Create a color palette with 32 colors
# gen_color_mapping <- colorRampPalette(brewer.pal(12, "Paired"))(32)
# 
# # Create a named colour vector
# color_mapping_gen <- setNames(colour_palette_29[1:length(unique_genus)], unique_genus)

# Specific colours
gen_colors <- c(
  "#d96567", "#658eaf", "#6ea66d", "#f2a961", "#f2f27e", "#9e7056", "#eba7cb", "#84b8a8",
  "#efb098", "#9ea8c1", "#dba6c7", "#b1cd82", "#dac7ab", "#aaaaaa", "#4b9680", "#ce8954",
  "#db6fa6", "#799e50", "#dab959", "#9e8250", "#85633b", "#cdb88c", "#eae2cd", "#cadedb",
  "#78aba6", "#644a81", "#a19bb9", "#e5b880", "#ebebeb", "#97c3bc", "#d4c39c", "#337e73"
)

gen_color_mapping <- setNames(gen_colors[1:length(unique_genus)], unique_genus)


plot_grouped_top_taxa(
  data = marti_genus,
  rank = "genus",
  color_mapping = gen_color_mapping,
  output_path = "Images/top5/Genus_Top5_AllSamp.png"
)

plot_top_taxa_by_sampler(
  data = filtered_marti,
  samplers = c("Compact", "μ", "Bobcat", "Cub", "Sass"),
  rank = "genus",
  color_mapping = gen_color_mapping,
  output_dir = "Images/top5/"
) 

##Species ----------------------------------

marti_species <- get_top_taxa(filtered_marti, "species", top_n = 5)

# 34 unique species - not that diff to genus level 
unique_species <- unique(marti_species$Name)
unique_species_df <- as.data.frame(unique_species)

# Create a color palette with 39 colors
spe_colors <- c(
  "#d96567", "#658eaf", "#6ea66d", "#956a9b", "#f2a961", "#f2f27e", "#9e7056", "#eba7cb",
  "#919191", "#84b8a8", "#efb098", "#9ea8c1", "#dba6c7", "#b1cd82", "#f2dd7c", "#dac7ab",
  "#aaaaaa", "#4b9680", "#ce8954", "#8784aa", "#db6fa6", "#799e50", "#dab959", "#9e8250",
  "#616161", "#85633b", "#cdb88c", "#eae2cd", "#cadedb", "#78aba6", "#27615c", "#644a81",
  "#a19bb9", "#e5b880", "#ebebeb", "#b59262", "#97c3bc", "#d4c39c", "#337e73"
)

# Map to your species (or genera)
spe_color_mapping <- setNames(spe_colors[1:length(unique_species)], unique_species)

plot_grouped_top_taxa(
  data = marti_species,
  rank = "species",
  color_mapping = spe_color_mapping,
  output_path = "Images/top5/Species_Top5_AllSamp.png"
)

plot_top_taxa_by_sampler(
  data = filtered_marti,
  samplers = c("Compact", "μ", "Bobcat", "Cub", "Sass"),
  rank = "species",
  color_mapping = spe_color_mapping,
  output_dir = "Images/top5"
) 

plot_grouped_top_taxa_replicates(
  data = marti_species,
  rank = "species",
  color_mapping = spe_color_mapping,
  output_path = "Images/top5/Species_Top5_Repeats.png"
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







