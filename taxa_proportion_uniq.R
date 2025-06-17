#Packages
library(dplyr)
library(ggplot2)
library(tidyr)
library(reshape2)
library(scales)


#sample meta data (this hasn't changed)
meta <- read.csv("old_parameters_MARTI_samp_comp_read_data/Phyloseq_data/Sample_table.csv")

#Marti-------------
#Loading in data 
#Older data
#marti_sum <-read.csv("../Older MARTi Runs/MARTI_samp_comp_updated_parameters_0823/marti_output_august23_summarised_taxaname.csv")

#New MARTi analysis 2024
marti_sum <- read.csv("samp_comp_0624_marti/samp_comp_summed_0624.csv")

extract_rank_df <- function(data, rank) {
  data %>%
    filter(Rank == rank) %>%
    select(-NCBI_ID)
}

marti_species <- extract_rank_df(marti_sum, "species")
marti_genus   <- extract_rank_df(marti_sum, "genus")
marti_phylum  <- extract_rank_df(marti_sum, "phylum")

# Data wrangling --------------
# Option to do diff taxonomic levels 

make_presence_absence <- function(df) {
  df[, -c(1,2)] <- ifelse(df[, -c(1,2)] > 0, 1, 0)
  return(df)
}

marti_presence_species <- make_presence_absence(marti_species)
marti_presence_genus   <- make_presence_absence(marti_genus)
marti_presence_phylum  <- make_presence_absence(marti_phylum)

#Unique taxa - change this depending on taxonomic level of interest
marti_presence_data <- marti_presence_species

# To treat samples as individuals -------------
# Step 1: Identify unique species (present in only one sample)
unique_taxa <- marti_presence_data$Taxon[
  rowSums(marti_presence_data[, -c(1,2)]) == 1
]

# Step 2: Filter data to only those unique taxa
marti_unique <- marti_presence_data %>%
  filter(Taxon %in% unique_taxa)

# Step 3: Melt to long format and filter to present taxa
unique_long <- marti_unique %>%
  melt(id.vars = c("Taxon", "Rank")) %>%
  filter(value == 1)

# Step 4: Count unique species per sample
uniq_counts <- unique_long %>%
  group_by(variable) %>%
  summarise(uniq_count = n(), .groups = "drop") 

# rename variable column, was behaving weird so have this workaround
uniq_counts$Sample_ID <- as.character(uniq_counts$variable)
uniq_counts <- uniq_counts %>%
  select(Sample_ID, uniq_count)

# Get full list of samples from your metadata
all_samples <- unique(meta$Sample_ID)

# Identify missing samples (not in uniq_counts)
missing_samples_df <- data.frame(
  Sample_ID = setdiff(all_samples, uniq_counts$Sample_ID),
  uniq_count = 0
)

# Combine counts with missing = 0 and join with metadata
uniq <- bind_rows(uniq_counts, missing_samples_df) %>%
  left_join(meta, by = "Sample_ID")

# Group by sampler first ------------
# Updated code to group by sampler then find unique June 2025

# Step 1: Add sampler identity
meta_info <- meta %>% select(Sample_ID, Sampler)
presence_with_sampler <- marti_presence_data %>%
  pivot_longer(-c(Taxon, Rank), names_to = "Sample_ID", values_to = "presence") %>%
  left_join(meta_info, by = "Sample_ID")

# Step 2: Summarise presence by sampler (1 if found in any of that sampler's samples)
sampler_presence <- presence_with_sampler %>%
  group_by(Taxon, Sampler) %>%
  summarise(present = as.integer(any(presence == 1)), .groups = "drop") %>%
  pivot_wider(names_from = Sampler, values_from = present, values_fill = 0)

# Step 3: Identify taxa unique to one sampler
sampler_presence$uniqueness <- rowSums(sampler_presence[ , -1]) == 1
unique_taxa <- sampler_presence %>%
  filter(uniqueness == TRUE)

# Step 4: Count unique taxa per sampler
unique_counts <- unique_taxa %>%
  pivot_longer(-c(Taxon, uniqueness), names_to = "Sampler", values_to = "presence") %>%
  filter(presence == 1) %>%
  group_by(Sampler) %>%
  summarise(n_unique_taxa = n(), .groups = "drop")

#Shared taxa per sampler ----
shared_taxa_all <- sampler_presence %>%
  filter(uniqueness == FALSE)

# Unique and shared taxa counts -----
long_taxa <- sampler_presence %>%
  pivot_longer(cols = -c(Taxon, uniqueness), names_to = "Sampler", values_to = "Present") %>%
  filter(Present == 1)  %>%
  mutate(Type = ifelse(uniqueness, "Unique", "Shared"))

taxa_counts <- long_taxa %>%
  group_by(Sampler, Type) %>%
  summarise(n_taxa = n(), .groups = "drop")

# Total species per sampler
species_per_sampler <- presence_with_sampler %>%
  filter(presence == 1) %>%
  group_by(Sampler) %>%
  summarise(n_unique_species = n_distinct(Taxon), .groups = "drop")


#Plotting------------------------

#Plot aesthetics 
# sampler_colours <- c('#BBCC33', '#77AADD', '#EE8866', '#EEDD88',  '#99DDFF')

sampler_colours <- setNames(brewer.pal(5, "Set2"), 
                            c("Bobcat", "Compact", 
                                "Micro", "Cub", 
                                "Sass"))
# Theme for plots ------
custom_theme <- theme_minimal(base_size = 12) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.3),
    panel.grid.minor = element_blank(),
  )

location_labels <- c(
  "NHM" = "Natural History Museum",
  "Cfarm" = "Church Farm"
)


duration_labels <- c(
  "25" = "25 minutes",
  "50" = "50 minutes"
)

sampler_levels <- c(
  "Compact",
  "Micro",
  "Bobcat",
  "Cub",
  "Sass"
)

# Same order as all other plots 
uniq$Sampler <- factor(uniq$Sampler, levels = sampler_levels)

#Plot - Group by location -----
# Get mean to plot 
uniq_summary_plot <- uniq %>%
  group_by(Location, Sample_length, Sampler) %>%
  summarise(
    uniq_mean = mean(uniq_count),
    uniq_se = sd(uniq_count) / sqrt(n()),
    .groups = "drop"
  )


location_uniq <- ggplot(uniq_summary_plot, aes(x=Sampler, y = uniq_mean, fill = Sampler)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = uniq_mean - uniq_se, ymax = uniq_mean + uniq_se),
                width = 0.2, colour = "black") +
  facet_grid(rows = vars(Location), cols = vars(Sample_length), 
             labeller = labeller(
               Location = location_labels,
               Sample_length = duration_labels
             )) +
  scale_fill_manual(values = sampler_colours) +
  labs( x = "Sample",
        y = "Mean number of unique genera ± SE") +
  scale_x_discrete(labels = c(
    "Compact" = "Coriolis\nCompact",
    "Micro" = "Coriolis\nMicro",
    "Bobcat" = "InnovaPrep\nBobcat",
    "Cub" = "InnovaPrep\nCub",
    "Sass" = "SASS\n3100"
  )) + 
  custom_theme 

location_uniq

ggsave("Images/graphs_marti_0624/uniq_shared/Unique_Species_per_Sample.pdf",
       plot = location_uniq , device = "svg", width = 12, height = 8)


#Shared & Unique plot 
# Working here - because I'm not sure if I trust it
taxa_counts$Sampler <- factor(taxa_counts$Sampler, levels = sampler_levels)
ggplot(taxa_counts, aes(x = Sampler, y = n_taxa, fill = Type)) +
  geom_bar(stat = "identity", position = "stack", colour = "black") +
  scale_fill_manual(values = c("Unique" = "#d53e4f", "Shared" = "#3288bd")) +
  labs(
    y = "Number of species detected",
    x = "Sampler",
    fill = "Taxon type"
  ) +  scale_x_discrete(labels = c(
    "Compact" = "Coriolis\nCompact",
    "Micro" = "Coriolis\nMicro",
    "Bobcat" = "InnovaPrep\nBobcat",
    "Cub" = "InnovaPrep\nCub",
    "Sass" = "SASS\n3100"
  )) + 
  custom_theme




#Plot - Group by Sample length 
length_uniq <- ggplot(uniq, aes(x=Sample_ID, y = uniq_count, fill = Sampler)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  facet_grid(~Sample_length, scales ="free_x") +
  scale_fill_manual(values = sampler_colours) +
  labs( x = "Sample",
        y = "Number of unique genera") 

ggsave("../Images/graphs_marti_0624/uniq_shared/Unique_Genus_Sample_Length.svg",
       plot = length_uniq , device = "svg", width = 10, height = 7)


#Total species in each sample ------------

numeric_col_presence_data <- marti_presence_data[, 3:ncol(marti_presence_data)]
total <- data.frame(Sample_ID = colnames(numeric_col_presence_data), 
                    total_count = colSums(numeric_col_presence_data))


#Combining the total number onto big uniq dataframe
taxa_counts <- inner_join(uniq, total, by = "Sample_ID")

#This order is important for later binding steps!
order_rows <- order(taxa_counts$Sample_ID) #Getting alphabetical order
taxa_counts <- taxa_counts[order_rows, ] #rearranging by this order


#Shared - Gettin a count for only the shared species (total includes the uniq)
taxa_counts <- taxa_counts %>% 
  mutate(shared_count = total_count - uniq_count)


#Reshape data - duplicating each row, so there is a seperate row for the unique and shared count
taxa_count_long <- taxa_counts %>%
  pivot_longer(
    cols = c(uniq_count, shared_count),
    names_to = "Number",
    values_to = "Count"
  ) %>%
  mutate(
    Number = case_when(
      Number == "uniq_count" ~ "Unique",
      Number == "shared_count" ~ "Shared"
    )
  )

#Plot proportion of uniq & shared 
all_uniq_prop_bar <- ggplot(taxa_count_long, aes(x = Sample_ID, y = Count, fill = Number)) +
  geom_bar(stat = "identity", 
          # position = "fill"
           ) +
  facet_grid(~Location, scales ="free_x") +
  scale_fill_manual(values = sampler_colours) +
  labs(x = "Sample",
       y = "Percentage of Classified Genera") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1))

ggsave("../Images/graphs_marti_0624/uniq_shared/Proportion_Unique_Genera_Location_count.svg",
       plot = all_uniq_prop_bar , device = "svg", width = 10, height = 7)


# Looking at uniq within a Location ---------------------------------------

#Church Farm, need to recalculate uniq counts ignoring those in NHM data ----
CF_samples <- marti_presence_data %>% 
  select(1:2, starts_with("CF"))

CF_u <- CF_samples$Taxon[rowSums(CF_samples[, -c(1,2)]) == 1] 

CF_uniq <- #subset df to only contain uniq taxon
  CF_samples[marti_presence_data$Taxon %in% CF_u, ]

CF_uniq_name <- 
  melt(CF_uniq, id.vars = c("Taxon", "Rank")) %>% # Melt the data to long
  filter(value == 1) #Only keep present taxa

CF_uniq_name_num <- CF_uniq_name %>% 
  group_by(variable) %>% 
  summarise(count = n()) 

#NHM same as with CF ----
NHM_samples <- marti_presence_data %>% 
  select(1:2, starts_with("NHM"))

NHM_u <- NHM_samples$Taxon[rowSums(NHM_samples[, -c(1,2)]) == 1] 

NHM_uniq <- #subset df to only contain uniq taxon
  NHM_samples[marti_presence_data$Taxon %in% NHM_u, ]

NHM_uniq_name <- 
  melt(NHM_uniq, id.vars = c("Taxon", "Rank")) %>% # Melt the data to long
  filter(value == 1) #Only keep present taxa

NHM_uniq_name_num <- NHM_uniq_name %>% 
  group_by(variable) %>% 
  summarise(count = n()) #no uniq species per sample

#Combine
Uniq_per_location <- rbind(CF_uniq_name_num, NHM_uniq_name_num)

#Samples with 0 uniq species - need to check and change this with each analysis 
missing_samples_2 <- data.frame(
  variable = c("CF_Micro_2","CF_SASS_2","CF_SASS_4", "NHM_Cub_3","NHM_SASS_2","NHM_SASS_3","NHM_SASS_4"),
  count = c("0", "0", "0", "0", "0", "0", "0"))

Uniq_per_location <- rbind(Uniq_per_location, missing_samples_2)

colnames(Uniq_per_location)[c(1,2)] <- c("Sample_ID","uniq_count") #renaming col names
Uniq_per_location$uniq_count <- as.numeric(Uniq_per_location$uniq_count)
Uniq_per_location$Sample_ID <- as.character(Uniq_per_location$Sample_ID)

#Need to be in the same order before binding 
order_rows <- order(Uniq_per_location$Sample_ID) #Getting alphabetical order
Uniq_per_location <- Uniq_per_location[order_rows, ] #rearranging by this order

#Add metadata 
Uniq_per_location <- Uniq_per_location %>%
  left_join(meta, by = "Sample_ID") 

#Add total taxon count calculated earlier in the analysis
Uniq_per_location <- inner_join(Uniq_per_location, total, by = "Sample_ID")

#Shared 
Uniq_per_location <- Uniq_per_location %>% 
  mutate(shared_count = total_count - uniq_count)

#Reshape data
Uniq_per_location <- Uniq_per_location %>%
  pivot_longer(
    cols = c(uniq_count, shared_count),
    names_to = "Number",
    values_to = "Count"
  ) %>%
  mutate(
    Number = case_when(
      Number == "uniq_count" ~ "Unique",
      Number == "shared_count" ~ "Shared"
    )
  )

#Plot proportion of uniq & shared Taxa, by location ------
prop_location <- ggplot(Uniq_per_location, aes(x = Sample_ID, y = Count, fill = Number)) +
  geom_bar(stat = "identity",
           #position = "fill"
           ) +
  facet_grid(~Location, scales ="free_x") +
  scale_fill_manual(values = sampler_colours) +
  labs(x = "Sample",
       y = "Percentage of Classified Genera") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1))

ggsave("../Images/graphs_marti_0624/uniq_shared/Proportion_Unique_Genera_Location_Seperate_count.svg",
       plot = prop_location , device = "svg", width = 10, height = 7)

#Unique per Sampler ----------------------------------------------------
#Bobcat
Bobcat_samples <- marti_presence_data %>% 
  select(1:2, contains("Bobcat"))

Bobcat_u <- Bobcat_samples$Taxon[rowSums(Bobcat_samples[, -c(1,2)]) == 1] 

Bobcat_uniq <- #subset df to only contain uniq genus
  Bobcat_samples[marti_presence_data$Taxon %in% Bobcat_u, ]

Bobcat_uniq_name <- 
  melt(Bobcat_uniq, id.vars = c("Taxon", "Rank")) %>% # Melt the data to long
  filter(value == 1) #Only keep present taxa

Bobcat_uniq_name_num <- Bobcat_uniq_name %>% 
  group_by(variable) %>% 
  summarise(count = n()) #no uniq genus per sample

#Compact
Compact_samples <- marti_presence_data %>% 
  select(1:2, contains("Compact"))

Compact_u <- Compact_samples$Taxon[rowSums(Compact_samples[, -c(1,2)]) == 1] 

Compact_uniq <- #subset df to only contain uniq genus
  Compact_samples[marti_presence_data$Taxon %in% Compact_u, ]

Compact_uniq_name <- 
  melt(Compact_uniq, id.vars = c("Taxon", "Rank")) %>% # Melt the data to long
  filter(value == 1) #Only keep present taxa

Compact_uniq_name_num <- Compact_uniq_name %>% 
  group_by(variable) %>% 
  summarise(count = n()) #no uniq genus per sample

#Cub
Cub_samples <- marti_presence_data %>% 
  select(1:2, contains("Cub"))

Cub_u <- Cub_samples$Taxon[rowSums(Cub_samples[, -c(1,2)]) == 1] 

Cub_uniq <- #subset df to only contain uniq genus
  Cub_samples[marti_presence_data$Taxon %in% Cub_u, ]

Cub_uniq_name <- 
  melt(Cub_uniq, id.vars = c("Taxon", "Rank")) %>% # Melt the data to long
  filter(value == 1) #Only keep present taxa

Cub_uniq_name_num <- Cub_uniq_name %>% 
  group_by(variable) %>% 
  summarise(count = n()) #no uniq genus per sample

#Micro
Micro_samples <- marti_presence_data %>% 
  select(1:2, contains("Micro"))

Micro_u <- Micro_samples$Taxon[rowSums(Micro_samples[, -c(1,2)]) == 1] 

Micro_uniq <- #subset df to only contain uniq genus
  Micro_samples[marti_presence_data$Taxon %in% Micro_u, ]

Micro_uniq_name <- 
  melt(Micro_uniq, id.vars = c("Taxon", "Rank")) %>% # Melt the data to long
  filter(value == 1) #Only keep present taxa

Micro_uniq_name_num <- Micro_uniq_name %>% 
  group_by(variable) %>% 
  summarise(count = n()) #no uniq genus per sample

#Sass
Sass_samples <- marti_presence_data %>% 
  select(1:2, contains("Sass"))

Sass_u <- Sass_samples$Taxon[rowSums(Sass_samples[, -c(1,2)]) == 1] 

Sass_uniq <- #subset df to only contain uniq genus
  Sass_samples[marti_presence_data$Taxon %in% Sass_u, ]

Sass_uniq_name <- 
  melt(Sass_uniq, id.vars = c("Taxon", "Rank")) %>% # Melt the data to long
  filter(value == 1) #Only keep present taxa

Sass_uniq_name_num <- Sass_uniq_name %>% 
  group_by(variable) %>% 
  summarise(count = n()) #no uniq genus per sample

#Combine
Uniq_per_sampler <- rbind(Bobcat_uniq_name_num, Compact_uniq_name_num, 
                          Cub_uniq_name_num, Micro_uniq_name_num, 
                          Sass_uniq_name_num)

#Samples with 0 uniq  
missing_samples_3 <- data.frame(
 variable = c("NHM_SASS_2", "CF_Micro_2"),
 count = c("0", "0"))

Uniq_per_sampler <- rbind(Uniq_per_sampler, missing_samples_3)

colnames(Uniq_per_sampler)[c(1,2)] <- c("Sample_ID","uniq_count") #renaming col names
Uniq_per_sampler$uniq_count <- as.numeric(Uniq_per_sampler$uniq_count)
Uniq_per_sampler$Sample_ID <- as.character(Uniq_per_sampler$Sample_ID)

#Need to be in the same order before binding 
order_rows <- order(Uniq_per_sampler$Sample_ID) #Getting alphabetical order
Uniq_per_sampler <- Uniq_per_sampler[order_rows, ] #rearranging by this order

columns_to_add <- c("Sampler", "Location", "total_count")
Uniq_per_sampler <- cbind(Uniq_per_sampler, taxa_counts[columns_to_add]) #adding cols of metadata 

#Shared 
Uniq_per_sampler <- Uniq_per_sampler %>% 
  mutate(shared_count = total_count - uniq_count)

#Reshape data
Uniq_per_sampler <- Uniq_per_sampler %>%
  pivot_longer(
    cols = c(uniq_count, shared_count),
    names_to = "Number",
    values_to = "Count"
  ) %>%
  mutate(
    Number = case_when(
      Number == "uniq_count" ~ "Unique",
      Number == "shared_count" ~ "Shared"
    )
  )

#Plot proportion of uniq & shared species 
sampler_uniq_prop_bar <- ggplot(Uniq_per_sampler, aes(x = Sample_ID, y = Count, fill = Number)) +
  geom_bar(stat = "identity",
           #position = "fill"
           ) +
  facet_grid(~Sampler, scales ="free_x") +
  scale_fill_manual(values = sampler_colours) +
  labs(x = "Sample",
       y = "Percentage of Classified Genera") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1))

ggsave("../Images/graphs_marti_0624/uniq_shared/Proportion_Unique_Genera_Sampler_Seperate_count.svg",
       plot = sampler_uniq_prop_bar , device = "svg", width = 10, height = 7)


##
#Older work ----------
##

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