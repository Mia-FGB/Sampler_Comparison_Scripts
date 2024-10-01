#setting up working directory and packages
setwd("~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Sampler comparison")
library(dplyr)
library(ggplot2)
library(lemon)
library(ggh4x)
library(forcats)
library(tidyverse)

#read in data
taxarank_data <- read.csv("../PHIbase_read_data/sampler_comparison_PHIbase_alldetails.csv")


#Need to normalise the reads for each location, e.g. hits per million reads
#Therefore we need the total number of reads for each sample, not just the classified reads
normalised_data <- 
  taxarank_data %>% 
  mutate(hits_per_100.000 = (read_count * 100000)/(total_read)) %>% 
  #log scale-Mutate the data to add a new row, making R see the barcodes as factors not numbers 
  mutate(log_hits = log10(hits_per_100.000), barcode = factor(barcode)) %>% 
  #adding a column for collection volume
  mutate(Collection_volume = (Collection_rate * Collection_time))
  
#adding a column for identity to the dataframe
normalised_data$identity <- paste( normalised_data$Air_Sampler,
                                  normalised_data$Collection_volume)

#----------------------------
#Graph-----------------------
#----------------------------

#top 10 taxa species------------

#going to start just with species
#each bar is a single barcode
#each bar is made up of the top 10 taxa for that barcode

#need to mutate the data to only include top 10 normalised read count for each barcode
#used this stack overflow https://stackoverflow.com/questions/31939643/make-a-table-showing-the-10-largest-values-of-a-variable-in-r

top10_taxa_species <- 
  normalised_data %>% 
  #first grouping them by barcode, as we want top 10 for each
  group_by(barcode) %>% 
  #arranging from biggest to smallest by normalised read count (arrange is from dplyr)
  arrange(desc(hits_per_100.000)) %>% 
  #slice selects the first 10 observations
  slice_head(n = 10)
  #noticed a problem that there are lots with just a 2 and so the 10th is being selected based on alphabetical order
  #going to filter out any reads <4

table_top10_species 


#Percentage stacked barchart for species --------------------
#This graph splits them all out individually which means there is less opportunity to compare
ggplot(top10_taxa_species, aes(fill = species, x = barcode, y = hits_per_100.000))+
  geom_bar(position = "fill", stat = "identity") +
  theme_bw()+ ggtitle("Species abundance (Top 10 taxa)") + ylab("Normalised read count (hits per 100,000 reads)") +
  labs(fill = "Species") + theme(legend.position="bottom") +
  facet_grid(~Air_Sampler + location, scales = "free_x")


#March 2023 - I prefer this, need to do it for my other variables though 
#This graph splits into location but need to also have air sampler marked
top10_taxa_species_graph <- ggplot(top10_taxa_species, aes(fill = species, x =identity, y = hits_per_100.000))+
  geom_bar(position = "fill", stat = "identity") +
  theme_bw()+ ggtitle("Species abundance (Top 10 taxa)") + ylab("Normalised read count (hits per 100,000 reads)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(fill = "Species") + theme(legend.position="bottom") +
  facet_grid(rows = vars(location))
#just need to save it in a way that all visible 
ggsave("Species_abundance_top10_0303.png", plot = top10_taxa_species_graph, width = 10, height = 10)

#Now want it not percentage but actual values
top10_taxa_species_stack_graph <- ggplot(top10_taxa_species, aes(fill = species, x =identity, y = hits_per_100.000))+
  geom_bar(position = "stack", stat = "identity") +
  theme_bw()+ ggtitle("Species abundance (Top 10 taxa)") + ylab("Normalised read count (hits per 100,000 reads)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(fill = "Species") + theme(legend.position="bottom") +
  facet_grid(rows = vars(location))
#To save
ggsave("Species_abundance_stack_top10_0703.png", plot = top10_taxa_species_stack_graph, width = 10, height = 10)


#Including other species ----------------------------------------------------
#Tried to have top10 per barcode and other but none of the code I tried would do it

#Need a way to group for each barcode
#Trying https://stackoverflow.com/questions/70304967/how-to-select-top-n-values-and-group-the-rest-of-the-remaining-ones
species_sort <- 
  normalised_data %>% 
  #first grouping them by barcode, as we want top 10 for each
  group_by(barcode) %>% 
  #arranging from biggest to smallest by normalised read count (arrange is from dplyr)
  arrange(desc(hits_per_100.000)) 

top10_species <- species_sort %>%  slice(1:10)

#Note the genus & family columns aren't relevant as have been grouped together and each have diff
other_species <- species_sort %>% anti_join(top10_species, by = c("species", "barcode")) %>% 
  mutate(species = "0thers") %>%  mutate(hits_per_100.000 = sum(hits_per_100.000)) %>% 
  distinct(barcode, .keep_all = TRUE)

top10_species_final <- top10_species %>%  bind_rows(other_species)

#Plotting
#Current best version (2023)
top10_other_species_stack_graph <- ggplot(top10_species_final, aes(fill = species, x =identity, y = hits_per_100.000))+
  geom_bar(position = "stack", stat = "identity") +
  theme_bw()+ ggtitle("Species abundance (Top 10 taxa)") + ylab("Normalised read count (hits per 100,000 reads)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(fill = "Species") + theme(legend.position="bottom") +
  facet_grid(rows = vars(location))
#Want the others to be in grey
ggsave("Species_abundance_stack_top10_other_1403.pdf", plot = top10_other_species_stack_graph, width = 12, height = 12)

#-------------------------------------------------------------------------------
#Genus -------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#Including other Genus ----------------------------------------------------
#Need to group by Genus and total the read count
#can then normalise this against total reads
normalised_genus_data <-
  taxarank_data %>% 
  #need to group by barcode, also want to retain the other variables
  group_by(barcode, genus, location, total_read, Air_Sampler) %>% 
  #once grouped by genus, want to count how many reads there are for each genus 
  mutate(genus_count = sum(read_count)) %>% 
  #Only want one row per grouped genus
  distinct(genus, .keep_all = TRUE) %>% 
  select(barcode, location, genus, total_read, Air_Sampler,
         Collection_rate, Collection_time, DNA_yield, genus_count) %>% 
  #normalise the read count
  mutate(hits_per_100.000 = (genus_count * 100000)/(total_read)) %>% 
  mutate(Collection_volume = (Collection_rate * Collection_time))

#adding a column for identity to the dataframe
normalised_genus_data$identity <- paste( normalised_genus_data$Air_Sampler,
                                         normalised_genus_data$Collection_volume)
genus_sort <- 
  normalised_genus_data %>% 
  #first grouping them by barcode, as we want top 10 for each
  group_by(barcode) %>% 
  #arranging from biggest to smallest by normalised read count (arrange is from dplyr)
  arrange(desc(hits_per_100.000)) 

top10_genus <- genus_sort %>%  slice(1:10)

#All the not top 10 genus
other_genus <- genus_sort %>% anti_join(top10_genus, by = c("genus", "barcode")) %>% 
  mutate(genus = "0thers") %>%  mutate(hits_per_100.000 = sum(hits_per_100.000)) %>% 
  distinct(barcode, .keep_all = TRUE)

top10_genus_final <- top10_genus %>%  bind_rows(other_genus)

#Plotting
top10_taxa_genus_other_stack_graph <- ggplot(top10_genus_final, aes(fill = genus, x =identity, y = hits_per_100.000))+
  geom_bar(position = "stack", stat = "identity") +
  theme_bw()+ ggtitle("Genus abundance (Top 10 taxa)") + ylab("Normalised read count (hits per 100,000 reads)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(fill = "Genus") + theme(legend.position="bottom") +
  facet_grid(rows = vars(location))
#just need to save it in a way that all visible 
ggsave("Genus_abundance_stack_other_top10_1403.pdf", plot = top10_taxa_genus_other_stack_graph, width = 12, height = 12)

#Genus not including other--------------------


#then group into top 10
top10_taxagenus <- 
  normalised_genus_data %>% 
  #first grouping them by barcode, as we want top 10 for each
  group_by(barcode) %>% 
  #arranging from biggest to smallest by normalised read count (arrange is from dplyr)
  arrange(desc(hits_per_100.000)) %>% 
  #slice selects the first 10 observations
  slice(1:10)

#then plot a faceted graph with percentage abundance
ggplot(top10_taxagenus, aes(fill = genus, x = barcode, y = hits_per_100.000))+
  geom_bar(position = "fill", stat = "identity") +
  theme_bw()+ ggtitle("Genus abundance (Top 10 taxa)") + ylab("Normalised read count (hits per 100,000 reads)") +
  labs(fill = "genus") + theme(legend.position="bottom") +
  facet_grid(~Air_Sampler + location, scales = "free_x")


#This graph splits into location but need to also have air sampler marked
#In order to make the same graph I've done for species I need to have identity column which is in the normalised data
top10_taxa_genus_graph <- ggplot(top10_taxagenus, aes(fill = genus, x =identity, y = hits_per_100.000))+
  geom_bar(position = "fill", stat = "identity") +
  theme_bw()+ ggtitle("Genus abundance (Top 10 taxa)") + ylab("Normalised read count (hits per 100,000 reads)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(fill = "Genus") + theme(legend.position="bottom") +
  facet_grid(rows = vars(location))
#just need to save it in a way that all visible 
ggsave("Genus_abundance_top10_0303.png", plot = top10_taxa_genus_graph, width = 10, height = 10)

#Same but not percentage 
top10_taxa_genus_stack_graph <- ggplot(top10_taxagenus, aes(fill = genus, x =identity, y = hits_per_100.000))+
  geom_bar(position = "stack", stat = "identity") +
  theme_bw()+ ggtitle("Genus abundance (Top 10 taxa)") + ylab("Normalised read count (hits per 100,000 reads)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(fill = "Genus") + theme(legend.position="bottom") +
  facet_grid(rows = vars(location))
#just need to save it in a way that all visible 
ggsave("Genus_abundance_stack_top10_0703.png", plot = top10_taxa_genus_stack_graph, width = 10, height = 10)


#---------------------------------------------------------
#Family--------------------------------------------------
#---------------------------------------------------------

#Including other family ----------------------------------------------------
#Need to group by Genus and total the read count
#can then normalise this against total reads
normalised_family_data <-
  taxarank_data %>% 
  #need to group by barcode, also want to retain the other variables
  group_by(barcode, family, location, total_read, Air_Sampler) %>% 
  #once grouped by family, want to count how many reads there are for each family 
  mutate(family_count = sum(read_count)) %>% 
  #Only want one row per grouped genus
  distinct(family, .keep_all = TRUE) %>% 
  select(barcode, location, family, total_read, Air_Sampler,
         Collection_rate, Collection_time, DNA_yield, family_count) %>% 
  #normalise the read count
  mutate(hits_per_100.000 = (family_count * 100000)/(total_read)) %>% 
  mutate(Collection_volume = (Collection_rate * Collection_time))

#adding a column for identity to the dataframe
normalised_family_data$identity <- paste( normalised_family_data$Air_Sampler,
                                          normalised_family_data$Collection_volume)
family_sort <- 
  normalised_family_data %>% 
  #first grouping them by barcode, as we want top 10 for each
  group_by(barcode) %>% 
  #arranging from biggest to smallest by normalised read count (arrange is from dplyr)
  arrange(desc(hits_per_100.000)) 

top10_family <- family_sort %>%  slice(1:10)

#Taking all those that aren't top 10 and relabelling them as others
other_family <- family_sort %>% anti_join(top10_family, by = c("family", "barcode")) %>% 
  mutate(family = "0thers") %>%  mutate(hits_per_100.000 = sum(hits_per_100.000)) %>% 
  distinct(barcode, .keep_all = TRUE)

#Binding back the top 10 and the summed others
top10_family_final <- top10_family %>%  bind_rows(other_family)

#Plotting
top10_taxa_family_other_stack_graph <- ggplot(top10_family_final, aes(fill = family, x =identity, y = hits_per_100.000))+
  geom_bar(position = "stack", stat = "identity") +
  theme_bw()+ ggtitle("family abundance (Top 10 taxa)") + ylab("Normalised read count (hits per 100,000 reads)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(fill = "family") + theme(legend.position="bottom") +
  facet_grid(rows = vars(location))
#just need to save it in a way that all visible 
ggsave("family_abundance_stack_other_top10_1603.pdf", plot = top10_taxa_family_other_stack_graph, width = 12, height = 12)

#Just top 10 Family ---------

#then group into top 10
top10_taxa_family <- 
  normalised_family_data %>% 
  #first grouping them by barcode, as we want top 10 for each
  group_by(barcode) %>% 
  #arranging from biggest to smallest by normalised read count (arrange is from dplyr)
  arrange(desc(hits_per_100.000)) %>% 
  #slice selects the first 10 observations
  slice(1:10)

#older graph
ggplot(top10_taxa_family, aes(fill = family, x = barcode, y = hits_per_100.000))+
  geom_bar(position = "fill", stat = "identity") +
  theme_bw()+ ggtitle("Family abundance (Top 10 taxa)") + ylab("Normalised read count (hits per 100,000 reads)") +
  labs(fill = "family") + theme(legend.position="bottom") +
  facet_grid(~Air_Sampler + location, scales = "free_x")

#March 2023 - This is better 
top10_taxa_family_graph <- ggplot(top10_taxa_family, aes(fill = family, x =identity, y = hits_per_100.000))+
  geom_bar(position = "fill", stat = "identity") +
  theme_bw()+ ggtitle("Family abundance (Top 10 taxa)") + ylab("Normalised read count (hits per 100,000 reads)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(fill = "Family") + theme(legend.position="bottom") +
  facet_grid(rows = vars(location))
#just need to save it in a way that all visible 
ggsave("Family_abundance_top10_0303.png", plot = top10_taxa_family_graph, width = 10, height = 10)

#Not percent
top10_taxa_family_stack_graph <- ggplot(top10_taxa_family, aes(fill = family, x =identity, y = hits_per_100.000))+
  geom_bar(position = "stack", stat = "identity") +
  theme_bw()+ ggtitle("Family abundance (Top 10 taxa)") + ylab("Normalised read count (hits per 100,000 reads)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(fill = "Family") + theme(legend.position="bottom") +
  facet_grid(rows = vars(location))
#just need to save it in a way that all visible 
ggsave("Family_abundance_stack_top10_0703.png", plot = top10_taxa_family_stack_graph, width = 10, height = 10)
