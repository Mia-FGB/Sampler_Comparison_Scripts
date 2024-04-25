#Packages
library(dplyr)
library(ggplot2)
library(tidyr)
library(reshape2)
library(scales)

#*-------------
#Plot aesthetics 
sampler_colours <- c('#BBCC33', '#77AADD', '#EE8866', '#EEDD88',  '#99DDFF')
theme_set(theme_bw())

#Marti-------------
#Loading in data 
marti_sum <- read.csv("MARTI_samp_comp_updated_parameters_0823/marti_output_august23_summarised_taxaname.csv")
marti_species <- marti_sum %>% 
  filter(Rank == "species") %>% #filter to retain species only
  subset(select = -NCBI_ID)

marti_genus <- marti_sum %>% 
  filter(Rank == "genus") %>%
  subset(select = -NCBI_ID)

marti_phylum <- marti_sum %>% 
  filter(Rank == "phylum") %>%
  subset(select = -NCBI_ID)

#Across whole group analysis ----------------------------------------------------
# Convert non-zero counts to 1 for presence
#marti_presence_data <- marti_species 
#marti_presence_data[, -c(1,2)] <-
#  ifelse(marti_species[,-c(1,2)] > 0, 1, 0) #returns 1 if count > 0, else 0; ignores first 2 cols

#marti_presence_data <- marti_genus 
#marti_presence_data[, -c(1,2)] <-
#  ifelse(marti_genus[,-c(1,2)] > 0, 1, 0) 

marti_presence_data <- marti_phylum 
marti_presence_data[, -c(1,2)] <-
  ifelse(marti_phylum[,-c(1,2)] > 0, 1, 0) 

#Unique taxa 
unique_phy <- marti_presence_data$Taxon[rowSums(marti_presence_data[, -c(1,2)]) == 1] #Sum row 1 = unique taxa.

marti_unique_phy <- #subset df to only contain uniq genus
  marti_presence_data[marti_presence_data$Taxon %in% unique_phy, ]

marti_unique_sample_name <- 
  melt(marti_unique_phy, id.vars = c("Taxon", "Rank")) %>% # Melt the data to long
  filter(value == 1) #Only keep present taxa
         
#write.csv(marti_unique_sample_name, "MARTI_samp_comp_updated_parameters_0823/unique_genus_name_sample.csv")  
                          
num_uniq <- marti_unique_sample_name %>% 
  group_by(variable) %>% 
  summarise(count = n()) #no uniq species per sample

#write.csv(num_uniq, "MARTI_samp_comp_updated_parameters_0823/unique_species_num_sample.csv") 
#write.csv(num_uniq, "MARTI_samp_comp_updated_parameters_0823/unique_genus_num_sample.csv") 
#Edit in Excel to include samples with 0 & Sampler, Location 

#uniq_sp <- read.csv("MARTI_samp_comp_updated_parameters_0823/unique_species_num_sample_edit.csv")
uniq <-  read.csv("MARTI_samp_comp_updated_parameters_0823/unique_genus_num_sample_edit.csv")


#Plot - Group by location 
location_uniq_ge <- ggplot(uniq, aes(x=sample_ID, y = uniq_ge_count, fill = Sampler)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  facet_grid(~Location, scales ="free_x") +
  scale_fill_manual(values = sampler_colours) +
  labs( x = "Sample",
        y = "Number of unique genus") 


#Total species in each sample
numeric_col_presence_data <- marti_presence_data[, 3:ncol(marti_presence_data)]
total <- data.frame(sample_ID = colnames(numeric_col_presence_data), 
                            Total = colSums(numeric_col_presence_data))

#Combining 
taxa_counts <- inner_join(uniq, total, by = "sample_ID")
colnames(taxa_counts)[c(5)] <-"total_count" #renaming col names

#Shared 
taxa_counts <- taxa_counts %>% 
  mutate(shared_ge_count = total_count - uniq_ge_count)


#Reshape data
taxa_count_long <- taxa_counts %>%
  pivot_longer(
    cols = c(uniq_ge_count, shared_ge_count),
    names_to = "Number_ge",
    values_to = "Count"
  ) %>%
  mutate(
    Number_ge = case_when(
      Number_ge == "uniq_ge_count" ~ "Unique",
      Number_ge == "shared_ge_count" ~ "Shared"
    )
  )

#Plot proportion of uniq & shared 
all_uniq_prop_bar <- ggplot(taxa_count_long, aes(x = sample_ID, y = Count, fill = Number_ge)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(~Location, scales ="free_x") +
  scale_fill_manual(values = sampler_colours) +
  labs(x = "Sample",
       y = "Percentage of Classified Genus") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1))


# Looking at uniq within a location ---------------------------------------
#Church Farm 
CF_samples <- marti_presence_data %>% 
  select(1:2, starts_with("CF"))

CF_u_ge <- CF_samples$Taxon[rowSums(CF_samples[, -c(1,2)]) == 1] 

CF_uniq_ge <- #subset df to only contain uniq genus
  CF_samples[marti_presence_data$Taxon %in% CF_u_ge, ]

CF_uniq_ge_name <- 
  melt(CF_uniq_ge, id.vars = c("Taxon", "Rank")) %>% # Melt the data to long
  filter(value == 1) #Only keep present taxa

CF_uniq_ge_name_num <- CF_uniq_ge_name %>% 
  group_by(variable) %>% 
  summarise(count = n()) 

#NHM 
NHM_samples <- marti_presence_data %>% 
  select(1:2, starts_with("NHM"))

NHM_u_ge <- NHM_samples$Taxon[rowSums(NHM_samples[, -c(1,2)]) == 1] 

NHM_uniq_ge <- #subset df to only contain uniq genus
  NHM_samples[marti_presence_data$Taxon %in% NHM_u_ge, ]

NHM_uniq_ge_name <- 
  melt(NHM_uniq_ge, id.vars = c("Taxon", "Rank")) %>% # Melt the data to long
  filter(value == 1) #Only keep present taxa

NHM_uniq_ge_name_num <- NHM_uniq_ge_name %>% 
  group_by(variable) %>% 
  summarise(count = n()) #no uniq geecies per sample

#Combine
Uniq_per_location <- rbind(CF_uniq_ge_name_num, NHM_uniq_ge_name_num)

#Samples with 0 uniq species - need to check and change this with each analysis 
missing_samples <- data.frame(
  variable = c("CF_SASS_1","CF_SASS_2", "CF_SASS_3", "CF_SASS_4", "NHM_SASS_3"),
  count = c("0", "0", "0", "0", "0"))

Uniq_per_location <- rbind(Uniq_per_location, missing_samples)

colnames(Uniq_per_location)[c(1,2)] <- c("sample_ID","uniq_ge_count") #renaming col names
Uniq_per_location$uniq_ge_count <- as.numeric(Uniq_per_location$uniq_ge_count)
Uniq_per_location$sample_ID <- as.character(Uniq_per_location$sample_ID)

#Need to be in the same order before binding 
order_rows <- order(Uniq_per_location$sample_ID) #Getting alphabetical order
Uniq_per_location <- Uniq_per_location[order_rows, ] #rearranging by this order

columns_to_add <- c("Sampler", "Location", "total_count")
Uniq_per_location <- cbind(Uniq_per_location, taxa_counts[columns_to_add]) #adding cols of metadata 

#Shared 
Uniq_per_location <- Uniq_per_location %>% 
  mutate(shared_ge_count = total_count - uniq_ge_count)

#Reshape data
Uniq_per_location_long <- Uniq_per_location %>%
  pivot_longer(
    cols = c(uniq_ge_count, shared_ge_count),
    names_to = "Number_ge",
    values_to = "Count"
  ) %>%
  mutate(
    Number_ge = case_when(
      Number_ge == "uniq_ge_count" ~ "Unique",
      Number_ge == "shared_ge_count" ~ "Shared"
    )
  )

#Plot proportion of uniq & shared geecies 
ggplot(Uniq_per_location_long, aes(x = sample_ID, y = Count, fill = Number_ge)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(~Location, scales ="free_x") +
  scale_fill_manual(values = sampler_colours) +
  labs(x = "Sample",
       y = "Percentage of Classified Genera") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1))

#Unique per Sampler ----------------------------------------------------
#Bobcat
Bobcat_samples <- marti_presence_data %>% 
  select(1:2, contains("Bobcat"))

Bobcat_u_ge <- Bobcat_samples$Taxon[rowSums(Bobcat_samples[, -c(1,2)]) == 1] 

Bobcat_uniq_ge <- #subset df to only contain uniq genus
  Bobcat_samples[marti_presence_data$Taxon %in% Bobcat_u_ge, ]

Bobcat_uniq_ge_name <- 
  melt(Bobcat_uniq_ge, id.vars = c("Taxon", "Rank")) %>% # Melt the data to long
  filter(value == 1) #Only keep present taxa

Bobcat_uniq_ge_name_num <- Bobcat_uniq_ge_name %>% 
  group_by(variable) %>% 
  summarise(count = n()) #no uniq genus per sample

#Compact
Compact_samples <- marti_presence_data %>% 
  select(1:2, contains("Compact"))

Compact_u_ge <- Compact_samples$Taxon[rowSums(Compact_samples[, -c(1,2)]) == 1] 

Compact_uniq_ge <- #subset df to only contain uniq genus
  Compact_samples[marti_presence_data$Taxon %in% Compact_u_ge, ]

Compact_uniq_ge_name <- 
  melt(Compact_uniq_ge, id.vars = c("Taxon", "Rank")) %>% # Melt the data to long
  filter(value == 1) #Only keep present taxa

Compact_uniq_ge_name_num <- Compact_uniq_ge_name %>% 
  group_by(variable) %>% 
  summarise(count = n()) #no uniq genus per sample

#Cub
Cub_samples <- marti_presence_data %>% 
  select(1:2, contains("Cub"))

Cub_u_ge <- Cub_samples$Taxon[rowSums(Cub_samples[, -c(1,2)]) == 1] 

Cub_uniq_ge <- #subset df to only contain uniq genus
  Cub_samples[marti_presence_data$Taxon %in% Cub_u_ge, ]

Cub_uniq_ge_name <- 
  melt(Cub_uniq_ge, id.vars = c("Taxon", "Rank")) %>% # Melt the data to long
  filter(value == 1) #Only keep present taxa

Cub_uniq_ge_name_num <- Cub_uniq_ge_name %>% 
  group_by(variable) %>% 
  summarise(count = n()) #no uniq genus per sample

#Micro
Micro_samples <- marti_presence_data %>% 
  select(1:2, contains("Micro"))

Micro_u_ge <- Micro_samples$Taxon[rowSums(Micro_samples[, -c(1,2)]) == 1] 

Micro_uniq_ge <- #subset df to only contain uniq genus
  Micro_samples[marti_presence_data$Taxon %in% Micro_u_ge, ]

Micro_uniq_ge_name <- 
  melt(Micro_uniq_ge, id.vars = c("Taxon", "Rank")) %>% # Melt the data to long
  filter(value == 1) #Only keep present taxa

Micro_uniq_ge_name_num <- Micro_uniq_ge_name %>% 
  group_by(variable) %>% 
  summarise(count = n()) #no uniq genus per sample

#Sass
Sass_samples <- marti_presence_data %>% 
  select(1:2, contains("Sass"))

Sass_u_ge <- Sass_samples$Taxon[rowSums(Sass_samples[, -c(1,2)]) == 1] 

Sass_uniq_ge <- #subset df to only contain uniq genus
  Sass_samples[marti_presence_data$Taxon %in% Sass_u_ge, ]

Sass_uniq_ge_name <- 
  melt(Sass_uniq_ge, id.vars = c("Taxon", "Rank")) %>% # Melt the data to long
  filter(value == 1) #Only keep present taxa

Sass_uniq_ge_name_num <- Sass_uniq_ge_name %>% 
  group_by(variable) %>% 
  summarise(count = n()) #no uniq genus per sample

#Combine
Uniq_per_sampler <- rbind(Bobcat_uniq_ge_name_num, Compact_uniq_ge_name_num, 
                          Cub_uniq_ge_name_num, Micro_uniq_ge_name_num, 
                          Sass_uniq_ge_name_num)

#Samples with 0 uniq  
#missing_samp <- data.frame(
#  variable = c("CF_SASS_4"),
#  count = c("0"))

#Uniq_per_sampler <- rbind(Uniq_per_sampler, missing_samp)

colnames(Uniq_per_sampler)[c(1,2)] <- c("sample_ID","uniq_ge_count") #renaming col names
Uniq_per_sampler$uniq_ge_count <- as.numeric(Uniq_per_sampler$uniq_ge_count)
Uniq_per_sampler$sample_ID <- as.character(Uniq_per_sampler$sample_ID)

#Need to be in the same order before binding 
order_rows <- order(Uniq_per_sampler$sample_ID) #Getting alphabetical order
Uniq_per_sampler <- Uniq_per_sampler[order_rows, ] #rearranging by this order

columns_to_add <- c("Sampler", "Location", "total_count")
Uniq_per_sampler <- cbind(Uniq_per_sampler, taxa_counts[columns_to_add]) #adding cols of metadata 

#Shared 
Uniq_per_sampler <- Uniq_per_sampler %>% 
  mutate(shared_ge_count = total_count - uniq_ge_count)

#Reshape data
Uniq_per_sampler_long <- Uniq_per_sampler %>%
  pivot_longer(
    cols = c(uniq_ge_count, shared_ge_count),
    names_to = "Number_ge",
    values_to = "Count"
  ) %>%
  mutate(
    Number_ge = case_when(
      Number_ge == "uniq_ge_count" ~ "Unique",
      Number_ge == "shared_ge_count" ~ "Shared"
    )
  )

#Plot proportion of uniq & shared species 
sampler_uniq_prop_bar <- ggplot(Uniq_per_sampler_long, aes(x = sample_ID, y = Count, fill = Number_ge)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(~Sampler, scales ="free_x") +
  scale_fill_manual(values = sampler_colours) +
  labs(x = "Sample",
       y = "Percentage of Classified Genera") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1))



#Older work ----------
#MARTi data looking for unique samples in species data
marti_old <- read.csv("old_parameters_MARTI_samp_comp_read_data/taxa_assignments_lca_0.1_species_assi_Samp_Comp.csv")

#Marti Genus data 
marti_genus <- read.csv("read_data/MARTI_samp_comp_read_data/taxa_assignments_lca_0.1_genus_Samp_Comp.csv")

# Convert non-zero counts to 1 for presence
#checks each element, returns 1 if the count is greater than 0, else returns 0
#want to ignore the first 3 columns
marti_genus_presence_data <- marti_genus
marti_genus_presence_data[, -c(1,2,3)] <- ifelse(marti_genus[,-c(1,2,3)] > 0, 1, 0)

#find unique taxa
#calculating the row sum, excluding the Genus. Where it sums to 1 this is a unique taxa.
marti_genus_unique_taxa <- marti_genus_presence_data$Taxon[rowSums(marti_genus_presence_data[, -c(1,2,3)]) == 1]
# Find columns where unique taxa are present
col_marti_genus_unique_taxa <- colnames(marti_genus_presence_data)[
  colSums(marti_genus_presence_data[rowSums(marti_genus_presence_data[, -c(1,2,3)]) == 1, -c(1,2,3)]) > 0]

#MEGAN data ----------------

#MEGAN Genus -------
#loading in the data - 17.05
Megan_genus_old <- read.csv("read_data/MEGAN_samp_comp_data/MEGAN_taxonNametoCount_assi_SampComp.csv")
#reexported data 26.05
Megan_assi_Gen <- read.csv("read_data/MEGAN_samp_comp_data/MEGAN_genus_summarised_all_leaves_SampComp.csv")
#rename 1st col
colnames(Megan_assi_Gen)[1] <- "Genus"

#Creating a sample count table

# Convert non-zero counts to 1 for presence
#checks each element, returns 1 if the count is greater than 0, else returns 0
#want to ignore the first column [, -1]
meg_gen_presence_data <- Megan_assi_Gen
meg_gen_presence_data[, -1] <- ifelse(Megan_assi_Gen[, -1] > 0, 1, 0)

#find unique taxa 
#calculating the row sum, excluding the Genus. Where it sums to 1 this is a unique taxa.
meg_gen_unique_taxa <- meg_gen_presence_data$Genus[rowSums(meg_gen_presence_data[, -1]) == 1]

# Find columns where unique taxa are present
cols_meg_gen_uniq_taxa <- 
  colnames(meg_gen_presence_data)[colSums(meg_gen_presence_data[
    rowSums(meg_gen_presence_data[, -1]) == 1, -1]) > 0]
#use colnames to get Genera
#checking which col have a sum greater than 0

#Results:
#"megan_NHM_Compact_barcode19 - Luteolibacter" "megan_NHM_Micro_barcode15 - Thrips" 

#Taxa that are present in less than 2 samples
#calculating the row sum, excluding the Genus. Where it sums to 1 this is a unique taxa.
less_two_taxa <- presence_data$Genus[rowSums(presence_data[, -1]) <= 2]

# Find columns where unique taxa are present
columns_with_less_two_taxa <- colnames(presence_data)[colSums(presence_data[rowSums(presence_data[, -1]) <=2, -1]) > 0]

#Shared across all samples
shared_taxa <- presence_data$Genus[rowSums(presence_data[, -1]) == 40]
shared_taxa_38 <- presence_data$Genus[rowSums(presence_data[, -1]) >= 38]



#MEGAN species----

Meg_Spe_data <- read.csv("read_data/MEGAN_samp_comp_data/MEGAN_species_summarised_all_leaves_SampComp.csv")
#rename 1st col
colnames(Meg_Spe_data)[1] <- "Species"

#Creating a sample count table
meg_spe_presence_data <- Meg_Spe_data
meg_spe_presence_data[, -1] <- ifelse( Meg_Spe_data[, -1] > 0, 1, 0)

#find unique taxa 
#calculating the row sum, excluding the Species. Where it sums to 1 this is a unique taxa.
meg_spe_unique_taxa <- meg_spe_presence_data$Species[rowSums(meg_spe_presence_data[, -1]) == 1]

# Find columns where unique taxa are present
cols_meg_spe_uniq_taxa <- colnames(meg_spe_presence_data)[colSums(meg_spe_presence_data[rowSums(meg_spe_presence_data[, -1]) == 1, -1]) > 0]


#MEGAN kingdom method-----
#loading in the data
Megan_assi_Kin <- read.csv("read_data/MEGAN_Kingdon_taxonNametoCount_assi_SampComp.csv")
Megan_sum_Kin <- read.csv("read_data/MEGAN_Kingdon_taxonNametoCount_sum_SampComp.csv")
Megan_percent_kin <- read.csv("read_data/MEGAN_sum_Kingdom_percent_Samp_Comp.csv")

#transpose the dataframe 
Megan_assi_Kin <- Megan_assi_Kin %>% 
  data.table::transpose(make.names = 'X.Datasets', keep.names = 'X.Datasets')

Megan_percent_kin <- Megan_percent_kin  %>% 
  data.table::transpose(make.names = 'X.Datasets', keep.names = 'X.Datasets')

#rename 1st col
colnames(Megan_assi_Kin)[1] <- "Sample"
colnames(Megan_percent_kin)[1] <- "Sample"


#write to csv
write.csv(Megan_assi_Kin, "read_data/MEGAN_samp_comp_data/MEGAN_assi_Kingdom_counts_Samp_Comp.csv", row.names=FALSE)
write.csv(Megan_percent_kin, "read_data/MEGAN_percent_Kingdom_counts_Samp_Comp.csv", row.names=FALSE)

#comparing to single barcode data
barcode01_kingdom_percent <- read.csv("read_data/megan_NHM_Bobcat_barcode01_kingdon_percent.csv")
#transpose the dataframe 
barcode01_kingdom_percent <- barcode01_kingdom_percent %>% 
  data.table::transpose(make.names = 'NCBI', keep.names = 'NCBI')