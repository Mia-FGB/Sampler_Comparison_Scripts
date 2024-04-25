library(ggplot2)    
library(dplyr)
library(tidyr)

# Set plotting theme
theme_set(theme_bw())

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





# Marti summarised data ---------------------------------------------------

#Data manipulation
marti_sum <- read.csv("MARTI_samp_comp_updated_parameters_0823/marti_output_august23_summarised_taxaname.csv")

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
  slice_head(n = 5) # keeps top 5 rows per group 

#Plotting phylum all on one 
ggplot(marti_phylum, aes(x=Sample_ID, y = percent_classified_read, fill = Taxon)) +
  geom_bar(stat = "identity") +
  facet_grid(~Sampler, scales ="free_x") +
  labs(title = "Top 5 Phylum per sample",
       x = "Sample",
       y = "Percent of Classified Reads") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1))

#Bobcat
marti_phylum_bobcat <-  marti_phylum %>% 
  filter(Sampler == "Bobcat")

ggplot(marti_phylum_bobcat, aes(x=Sample_ID, y = percent_classified_read, fill = Taxon)) +
  geom_bar(stat = "identity") +
  facet_grid(~Location, scales ="free_x") +
  labs(title = "Top 5 Phylum Bobcat",
       x = "Sample",
       y = "Percent of Classified Reads") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1))

#Compact
marti_phylum_Compact <-  marti_phylum %>% 
  filter(Sampler == "Compact")

ggplot(marti_phylum_Compact, aes(x=Sample_ID, y = percent_classified_read, fill = Taxon)) +
  geom_bar(stat = "identity") +
  facet_grid(~Location, scales ="free_x") +
  labs(title = "Top 5 Phylum Compact",
       x = "Sample",
       y = "Percent of Classified Reads") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1))

#Micro
marti_phylum_Micro <-  marti_phylum %>% 
  filter(Sampler == "Micro")

ggplot(marti_phylum_Micro, aes(x=Sample_ID, y = percent_classified_read, fill = Taxon)) +
  geom_bar(stat = "identity") +
  facet_grid(~Location, scales ="free_x") +
  labs(title = "Top 5 Phylum Micro",
       x = "Sample",
       y = "Percent of Classified Reads") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1))

#Cub
marti_phylum_Cub <-  marti_phylum %>% 
  filter(Sampler == "Cub")

ggplot(marti_phylum_Cub, aes(x=Sample_ID, y = percent_classified_read, fill = Taxon)) +
  geom_bar(stat = "identity") +
  facet_grid(~Location, scales ="free_x") +
  labs(title = "Top 5 Phylum Cub",
       x = "Sample",
       y = "Percent of Classified Reads") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1))

#Sass
marti_phylum_Sass <-  marti_phylum %>% 
  filter(Sampler == "Sass")

ggplot(marti_phylum_Sass, aes(x=Sample_ID, y = percent_classified_read, fill = Taxon)) +
  geom_bar(stat = "identity") +
  facet_grid(~Location, scales ="free_x") +
  labs(title = "Top 5 Phylum Sass",
       x = "Sample",
       y = "Percent of Classified Reads") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1))

#Species top 5 graphs -------------------------
#Species
marti_species <- marti_meta %>% 
  filter(Rank %in% c("species")) %>%
  group_by(Sample_ID) %>%
  arrange(desc(percent_classified_read)) %>% 
  slice_head(n = 5)

#Too many to plot at once 
ggplot(marti_species, aes(x=Sample_ID, y = percent_classified_read, fill = Taxon)) +
  geom_bar(stat = "identity") +
  facet_grid(~Sampler, scales ="free_x") +
  labs(title = "Top 5 Species per sample",
       x = "Sample",
       y = "Percent of Classified Reads") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1))

#Bobcat
 marti_species_bobcat <-  marti_species %>% 
   filter(Sampler == "Bobcat")

 ggplot(marti_species_bobcat, aes(x=Sample_ID, y = percent_classified_read, fill = Taxon)) +
   geom_bar(stat = "identity") +
   facet_grid(~Location, scales ="free_x") +
   labs(title = "Top 5 Species Bobcat",
        x = "Sample",
        y = "Percent of Classified Reads") +
   theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1))
 
#Compact
marti_species_Compact <-  marti_species %>% 
   filter(Sampler == "Compact")
 
ggplot(marti_species_Compact, aes(x=Sample_ID, y = percent_classified_read, fill = Taxon)) +
   geom_bar(stat = "identity") +
   facet_grid(~Location, scales ="free_x") +
   labs(title = "Top 5 Species Compact",
        x = "Sample",
        y = "Percent of Classified Reads") +
   theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1))

#Cub
marti_species_Cub <-  marti_species %>% 
  filter(Sampler == "Cub")

ggplot(marti_species_Cub, aes(x=Sample_ID, y = percent_classified_read, fill = Taxon)) +
  geom_bar(stat = "identity") +
  facet_grid(~Location, scales ="free_x") +
  labs(title = "Top 5 Species Cub",
       x = "Sample",
       y = "Percent of Classified Reads") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1))

#Micro
marti_species_Micro <-  marti_species %>% 
  filter(Sampler == "Micro")

ggplot(marti_species_Micro, aes(x=Sample_ID, y = percent_classified_read, fill = Taxon)) +
  geom_bar(stat = "identity") +
  facet_grid(~Location, scales ="free_x") +
  labs(title = "Top 5 Species Micro",
       x = "Sample",
       y = "Percent of Classified Reads") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1))
 
#Sass
marti_species_Sass <-  marti_species %>% 
  filter(Sampler == "Sass")

ggplot(marti_species_Sass, aes(x=Sample_ID, y = percent_classified_read, fill = Taxon)) +
  geom_bar(stat = "identity") +
  facet_grid(~Location, scales ="free_x") +
  labs(title = "Top 5 Species Sass",
       x = "Sample",
       y = "Percent of Classified Reads") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1))
