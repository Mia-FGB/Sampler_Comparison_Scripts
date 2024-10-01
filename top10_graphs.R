library(ggplot2)    
library(dplyr)
library(tidyr)
library(RColorBrewer)

# Set plotting theme
theme_set(theme_bw())

#sample meta data (this hasn't changed)
meta <- read.csv("../old_parameters_MARTI_samp_comp_read_data/Phyloseq_data/Sample_table.csv")

# Marti summarised data ---------------------------------------------------

#Data manipulation
#This is data from a previous marti run
#marti_sum <- read.csv("../Older MARTi Runs/MARTI_samp_comp_updated_parameters_0823/marti_output_august23_summarised_taxaname.csv")

#Now running the script with newer MARTi run data
#I have slightly modified the MARTi output csv as before. Splitting it into assigned & summed files also changed the headers

marti_sum <- read.csv("../samp_comp_0624_marti/samp_comp_summed_0624.csv")


marti_long <- pivot_longer(marti_sum, cols = starts_with(c("CF_", "NHM_")),
                           names_to = "Sample_ID",
                           values_to = "Count")

marti_long$Count <- as.numeric(marti_long$Count)

marti_meta <- merge(marti_long, meta)
marti_meta$NumReads <- as.numeric(gsub(",","",marti_meta$NumReads)) #Remove commas from the NumRead col
marti_meta <- marti_meta %>% 
  mutate(percent_classified_read = Count/NumReads * 100)

#Phylum ---------------
marti_phylum <- marti_meta %>% 
  filter(Rank %in% c("phylum")) %>%
  group_by(Sample_ID) %>%
  arrange(desc(percent_classified_read)) %>% 
  slice_head(n = 5) %>% # keeps top 5 rows per group 
  filter(Count != 0)

#Colour palette - want consistency between graphs
Tol_muted <- c('#88CCEE', '#44AA99', '#117733', '#332288', '#DDCC77', '#999933','#CC6677', '#882255', '#AA4499', '#DDDDDD')

# Extract unique Phylum values
unique_phylums <- unique(marti_phylum$Taxon)

# Create a named vector
color_mapping <- setNames(Tol_muted[1:length(unique_phylums)], unique_phylums)

#Plotting Top 5 Phylum ------

#All samplers - too busy
all_phylum <- ggplot(marti_phylum, aes(x=Sample_ID, y = percent_classified_read, fill = Taxon)) +
  geom_bar(stat = "identity") +
  facet_grid(~Sampler, scales ="free_x") + #Can change this to  
  labs(title = "Top 5 Phylum per sample",
       x = "Sample",
       y = "Percent of Classified Reads") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1)) +
  scale_fill_manual(values = color_mapping)

ggsave("../Images/graphs_marti_0624/top5/Phylum_Top5_AllSamp.svg", plot = all_phylum , device = "svg", width = 10, height = 8)

#Bobcat
marti_phylum_bobcat <-  marti_phylum %>% 
  filter(Sampler == "Bobcat")

bobcat_phy_5 <- ggplot(marti_phylum_bobcat, aes(x=Sample_ID, y = percent_classified_read, fill = Taxon)) +
  geom_bar(stat = "identity") +
  facet_grid(~Location, scales ="free_x") +
  labs(title = "Top 5 Phylum Bobcat",
       x = "Sample",
       y = "Percent of Classified Reads") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1)) +
  scale_fill_manual(values = color_mapping)

# Save the plot as an SVG
ggsave("../Images/graphs_marti_0624/top5/Phylum_Top5_Bobcat.svg", plot = bobcat_phy_5 , device = "svg", width = 10, height = 8)

#Compact
marti_phylum_Compact <-  marti_phylum %>% 
  filter(Sampler == "Compact")

compact_phy_5 <- ggplot(marti_phylum_Compact, aes(x=Sample_ID, y = percent_classified_read, fill = Taxon)) +
  geom_bar(stat = "identity") +
  facet_grid(~Location, scales ="free_x") +
  labs(title = "Top 5 Phylum Compact",
       x = "Sample",
       y = "Percent of Classified Reads") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1)) +
  scale_fill_manual(values = color_mapping)

ggsave("../Images/graphs_marti_0624/top5/Phylum_Top5_Compact.svg", plot = compact_phy_5 , device = "svg", width = 10, height = 8)

#Micro
marti_phylum_Micro <-  marti_phylum %>% 
  filter(Sampler == "Micro")

micro_phy_5 <- ggplot(marti_phylum_Micro, aes(x=Sample_ID, y = percent_classified_read, fill = Taxon)) +
  geom_bar(stat = "identity") +
  facet_grid(~Location, scales ="free_x") +
  labs(title = "Top 5 Phylum Micro",
       x = "Sample",
       y = "Percent of Classified Reads") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1))  +
  scale_fill_manual(values = color_mapping)

ggsave("../Images/graphs_marti_0624/top5/Phylum_Top5_Micro.svg", plot = micro_phy_5 , device = "svg", width = 10, height = 8)

#Cub
marti_phylum_Cub <-  marti_phylum %>% 
  filter(Sampler == "Cub")

cub_phy_5 <- ggplot(marti_phylum_Cub, aes(x=Sample_ID, y = percent_classified_read, fill = Taxon)) +
  geom_bar(stat = "identity") +
  facet_grid(~Location, scales ="free_x") +
  labs(title = "Top 5 Phylum Cub",
       x = "Sample",
       y = "Percent of Classified Reads") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1)) +
  scale_fill_manual(values = color_mapping)

ggsave("../Images/graphs_marti_0624/top5/Phylum_Top5_Cub.svg", plot = cub_phy_5 , device = "svg", width = 10, height = 8)

#Sass
marti_phylum_Sass <-  marti_phylum %>% 
  filter(Sampler == "Sass")

sass_phy_5 <- ggplot(marti_phylum_Sass, aes(x=Sample_ID, y = percent_classified_read, fill = Taxon)) +
  geom_bar(stat = "identity") +
  facet_grid(~Location, scales ="free_x") +
  labs(title = "Top 5 Phylum Sass",
       x = "Sample",
       y = "Percent of Classified Reads") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1)) +
  scale_fill_manual(values = color_mapping)

ggsave("../Images/graphs_marti_0624/top5/Phylum_Top5_Sass.svg", plot = sass_phy_5 , device = "svg", width = 10, height = 8)


# Genus -------------------------------------------------------------------
#Genus
marti_genus <- marti_meta %>% 
  filter(Rank %in% c("genus")) %>%
  group_by(Sample_ID) %>%
  arrange(desc(percent_classified_read)) %>% 
  slice_head(n = 5) %>% # keeps top 5 rows per group 
  filter(Count != 0)

# 29 Unique genera
unique_genus <- unique(marti_genus$Taxon)
uniq_genus_df <- as.data.frame(unique_genus)

# Create a color palette with 29 colors
colour_palette_29 <- colorRampPalette(brewer.pal(12, "Paired"))(29)

# Create a named vector
color_mapping_gen <- setNames(colour_palette_29[1:length(unique_genus)], unique_genus)

#All samplers together
all_gen_5 <- ggplot(marti_genus, aes(x=Sample_ID, y = percent_classified_read, fill = Taxon)) +
  geom_bar(stat = "identity") +
  facet_grid(~Sampler, scales ="free_x") +
  labs(title = "Top 5 Genus per sample",
       x = "Sample",
       y = "Percent of Classified Reads") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1)) +
  scale_fill_manual(values = color_mapping_gen)

ggsave("../Images/graphs_marti_0624/top5/Genus_Top5_All.svg", plot = all_gen_5 , device = "svg", width = 10, height = 8)

#Bobcat
marti_genus_bobcat <-  marti_genus %>% 
  filter(Sampler == "Bobcat")

bobcat_gen_5 <- ggplot(marti_genus_bobcat, aes(x=Sample_ID, y = percent_classified_read, fill = Taxon)) +
  geom_bar(stat = "identity") +
  facet_grid(~Location, scales ="free_x") +
  labs(title = "Top 5 Genus Bobcat",
       x = "Sample",
       y = "Percent of Classified Reads") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1)) +
  scale_fill_manual(values = color_mapping_gen)

ggsave("../Images/graphs_marti_0624/top5/Genus_Top5_Bobcat.svg", plot = bobcat_gen_5 , device = "svg", width = 10, height = 8)

#Compact
marti_genus_Compact <-  marti_genus %>% 
  filter(Sampler == "Compact")

compact_gen_5 <- ggplot(marti_genus_Compact, aes(x=Sample_ID, y = percent_classified_read, fill = Taxon)) +
  geom_bar(stat = "identity") +
  facet_grid(~Location, scales ="free_x") +
  labs(title = "Top 5 Genus Compact",
       x = "Sample",
       y = "Percent of Classified Reads") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1)) +
  scale_fill_manual(values = color_mapping_gen)

ggsave("../Images/graphs_marti_0624/top5/Genus_Top5_Compact.svg", plot = compact_gen_5 , device = "svg", width = 10, height = 8)

#Cub
marti_genus_Cub <-  marti_genus %>% 
  filter(Sampler == "Cub")

cub_gen_5 <- ggplot(marti_genus_Cub, aes(x=Sample_ID, y = percent_classified_read, fill = Taxon)) +
  geom_bar(stat = "identity") +
  facet_grid(~Location, scales ="free_x") +
  labs(title = "Top 5 Genus Cub",
       x = "Sample",
       y = "Percent of Classified Reads") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1)) +
  scale_fill_manual(values = color_mapping_gen)

ggsave("../Images/graphs_marti_0624/top5/Genus_Top5_Cub.svg", plot = cub_gen_5 , device = "svg", width = 10, height = 8)

#Micro
marti_genus_Micro <-  marti_genus %>% 
  filter(Sampler == "Micro")

micro_gen_5 <- ggplot(marti_genus_Micro, aes(x=Sample_ID, y = percent_classified_read, fill = Taxon)) +
  geom_bar(stat = "identity") +
  facet_grid(~Location, scales ="free_x") +
  labs(title = "Top 5 Genus Micro",
       x = "Sample",
       y = "Percent of Classified Reads") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1)) +
  scale_fill_manual(values = color_mapping_gen)

ggsave("../Images/graphs_marti_0624/top5/Genus_Top5_Micro.svg", plot = micro_gen_5 , device = "svg", width = 10, height = 8)


#Sass
marti_genus_Sass <-  marti_genus %>% 
  filter(Sampler == "Sass")

sass_gen_5 <-ggplot(marti_genus_Sass, aes(x=Sample_ID, y = percent_classified_read, fill = Taxon)) +
  geom_bar(stat = "identity") +
  facet_grid(~Location, scales ="free_x") +
  labs(title = "Top 5 Genus Sass",
       x = "Sample",
       y = "Percent of Classified Reads") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1)) +
  scale_fill_manual(values = color_mapping_gen)

ggsave("../Images/graphs_marti_0624/top5/Genus_Top5_Sass.svg", plot = sass_gen_5 , device = "svg", width = 10, height = 8)

#Species top 5 graphs -------------------------

#Species
marti_species <- marti_meta %>% 
  filter(Rank %in% c("species")) %>%
  group_by(Sample_ID) %>%
  arrange(desc(percent_classified_read)) %>% 
  slice_head(n = 5) %>% # keeps top 5 rows per group 
  filter(Count != 0)

# 33 unique species - not that diff to genus level 
unique_species <- unique(marti_species$Taxon)
unique_species_df <- as.data.frame(unique_species)

# Create a color palette with 33 colors
colour_palette_33 <- colorRampPalette(brewer.pal(12, "Paired"))(33)

# Create a named vector
color_mapping_spec <- setNames(colour_palette_33[1:length(unique_species)], unique_species)

#All samplers together
all_spe_5 <- ggplot(marti_species, aes(x=Sample_ID, y = percent_classified_read, fill = Taxon)) +
  geom_bar(stat = "identity") +
  facet_grid(~Sampler, scales ="free_x") +
  labs(title = "Top 5 Species per sample",
       x = "Sample",
       y = "Percent of Classified Reads") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1)) +
  scale_fill_manual(values = color_mapping_spec)

ggsave("../Images/graphs_marti_0624/top5/Species_Top5_All.svg", plot = all_spe_5 , device = "svg", width = 10, height = 8)

#Bobcat
 marti_species_bobcat <-  marti_species %>% 
   filter(Sampler == "Bobcat")

bobcat_spe_5 <- ggplot(marti_species_bobcat, aes(x=Sample_ID, y = percent_classified_read, fill = Taxon)) +
   geom_bar(stat = "identity") +
   facet_grid(~Location, scales ="free_x") +
   labs(title = "Top 5 Species Bobcat",
        x = "Sample",
        y = "Percent of Classified Reads") +
   theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1)) +
   scale_fill_manual(values = color_mapping_spec)

ggsave("../Images/graphs_marti_0624/top5/Species_Top5_Bobcat.svg", plot = bobcat_spe_5 , device = "svg", width = 10, height = 8)
 
#Compact
marti_species_Compact <-  marti_species %>% 
   filter(Sampler == "Compact")
 
compact_spe_5 <- ggplot(marti_species_Compact, aes(x=Sample_ID, y = percent_classified_read, fill = Taxon)) +
   geom_bar(stat = "identity") +
   facet_grid(~Location, scales ="free_x") +
   labs(title = "Top 5 Species Compact",
        x = "Sample",
        y = "Percent of Classified Reads") +
   theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1)) +
  scale_fill_manual(values = color_mapping_spec)

ggsave("../Images/graphs_marti_0624/top5/Species_Top5_Compact.svg", plot = compact_spe_5 , device = "svg", width = 10, height = 8)

#Cub
marti_species_Cub <-  marti_species %>% 
  filter(Sampler == "Cub")

cub_spe_5 <- ggplot(marti_species_Cub, aes(x=Sample_ID, y = percent_classified_read, fill = Taxon)) +
  geom_bar(stat = "identity") +
  facet_grid(~Location, scales ="free_x") +
  labs(title = "Top 5 Species Cub",
       x = "Sample",
       y = "Percent of Classified Reads") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1)) +
  scale_fill_manual(values = color_mapping_spec)

ggsave("../Images/graphs_marti_0624/top5/Species_Top5_Cub.svg", plot = cub_spe_5 , device = "svg", width = 10, height = 8)

#Micro
marti_species_Micro <-  marti_species %>% 
  filter(Sampler == "Micro")

micro_spe_5 <- ggplot(marti_species_Micro, aes(x=Sample_ID, y = percent_classified_read, fill = Taxon)) +
  geom_bar(stat = "identity") +
  facet_grid(~Location, scales ="free_x") +
  labs(title = "Top 5 Species Micro",
       x = "Sample",
       y = "Percent of Classified Reads") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1)) +
  scale_fill_manual(values = color_mapping_spec)

ggsave("../Images/graphs_marti_0624/top5/Species_Top5_Micro.svg", plot = micro_spe_5 , device = "svg", width = 10, height = 8)

 
#Sass
marti_species_Sass <-  marti_species %>% 
  filter(Sampler == "Sass")

sass_spe_5 <-ggplot(marti_species_Sass, aes(x=Sample_ID, y = percent_classified_read, fill = Taxon)) +
  geom_bar(stat = "identity") +
  facet_grid(~Location, scales ="free_x") +
  labs(title = "Top 5 Species Sass",
       x = "Sample",
       y = "Percent of Classified Reads") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1)) +
  scale_fill_manual(values = color_mapping_spec)

ggsave("../Images/graphs_marti_0624/top5/Species_Top5_Sass.svg", plot = sass_spe_5 , device = "svg", width = 10, height = 8)

# PHYLOSEQ data -----------------------------------------------------------
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







