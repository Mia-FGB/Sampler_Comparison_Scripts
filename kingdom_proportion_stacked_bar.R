library(ggplot2)    
library(dplyr)
library(data.table)
library(tidyr)

#Reading in data ----------
#sample meta data (this hasn't changed)
meta <- read.csv("old_parameters_MARTI_samp_comp_read_data/Phyloseq_data/Sample_table.csv")

#Kingdom reads - This data was calculated in updated_phyloseq_mia_analysis.R script
kingdom <- read.csv("MARTI_samp_comp_updated_parameters_0823/kingdom_proportions_3008.csv")  


# Percent graphs unclassified 30.08 ---------------------------------------

#Create dataframe
kingdom_meta <- merge(kingdom, meta)

kingdom_meta <- kingdom_meta %>% 
  pivot_longer(cols = c(unclassified_reads, bacteria_reads, eukaryote_reads, archaea_reads, virus_reads),
               names_to = "Kingdom",
               values_to = "Abundance"
  )

kingdom_meta <- kingdom_meta %>% 
  group_by(Sample_ID) %>% 
  mutate(percent_classified_reads = Abundance / sum(Abundance) * 100) %>% 
  ungroup() %>% 
  mutate(percent_tot_reads = Abundance / no_reads * 100)

#Plot
kingdom_colours <-  c( "#1E88E5", "#D81B60","#FFC107", "#0A794D", "grey")

#renaming the legend
kingdom_meta <- kingdom_meta %>% 
  mutate(Kingdom = recode(Kingdom,
                          "unclassified_reads" = "Unclassified",
                          "bacteria_reads" = "Bacteria",
                          "eukaryote_reads" = "Eukaryote",
                          "archaea_reads" = "Archaea",
                          "virus_reads" = "Virus"))

#reordering the legend
kingdom_meta$Kingdom <- factor(kingdom_meta$Kingdom ,
                               levels=c('Archaea', 'Bacteria', 
                                        'Eukaryote', 'Virus', 'Unclassified'))


#Facet Location
ggplot(kingdom_meta, aes(x = Sample_ID, y = percent_classified_reads, fill = Kingdom)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(~Location, scales ="free_x") +
  scale_fill_manual(values = kingdom_colours) +
  labs(title = "Percent classified reads by Kingdom",
       x = "Sample",
       y = "Percentage Classified Reads") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1))

#Facet sampler
ggplot(kingdom_meta, aes(x = Sample_ID, y = percent_classified_reads, fill = Kingdom)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(~Sampler, scales ="free_x") +
  scale_fill_manual(values = kingdom_colours) +
  labs(title = "Percent classified reads by Kingdom",
       x = "Sample",
       y = "Percentage Classified Reads") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1))

#Facet 
ggplot(kingdom_meta, aes(x = Sample_ID, y = percent_classified_reads, fill = Kingdom)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(~Air_volume, scales ="free_x") +
  scale_fill_manual(values = kingdom_colours) +
  labs(title = "Percent classified reads by Kingdom",
       x = "Sample",
       y = "Percentage Classified Reads") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1))

