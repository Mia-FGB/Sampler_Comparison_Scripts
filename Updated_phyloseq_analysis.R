

#Library -----------------------------------------------------------------------
library(phyloseq)   # Facilitate the import, storage, analysis, and graphical display of microbiome census data.
library(vegan)      # Analysis of variance using distance matrices using the adonis2 function
library(Maaslin2)   # Determine multivariable association between metadata and microbial meta-omics features; differential abundance analysis
library(ggplot2)    # Generate visualization plots 
library(ggsignif)   # Visualize comparisons and the significance between two groups
library(dplyr)
library(RColorBrewer)
library(BiocManager)
library(Biostrings)
library(knitr)
library(scater)
library(patchwork)
library(stringr)
library(fantaxtic)
library(data.table)
library(tidyr)
library(forcats)
library(Cairo)

# Theme for plots ------
custom_theme <- theme_minimal(base_size = 12) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.3),
    panel.grid.minor = element_blank(),
  )

#Loading in data -----------------------------------------------------------------------


#OTU table 
#Using data from MARTi run with updated parameters 

#otu <- read.csv("../Older MARTi Runs/MARTI_samp_comp_updated_parameters_0823/marti_output_august23_assigned.csv") #August 2023

otu <- read.csv("samp_comp_0624_marti/samp_comp_assigned_otu_0624.csv") #MARTi run from June 2024
rownames(otu) <- otu[,1] #Making the rownames the NCBI ID
otu <- otu[,-1]          #remove 1st first col

#taxa table - need to remake this
#taxa <- read.csv("../Older MARTi Runs/MARTI_samp_comp_updated_parameters_0823/Taxonomy/taxa_lineage_sep_na.csv") #August 2023
taxa <- read.csv("samp_comp_0624_marti/samp_comp_assigned_0624_taxaID_lineage.csv")

taxa <- taxa[-1, ] #Remove the first row

taxa<- taxa %>% 
  mutate_if(is.character, as.factor)

rownames(taxa) <- rownames(otu) #same rownames as the otu table
taxa <- taxa[,-1]               

#sample meta data (this hasn't changed between 2023 and 2024 MARTi parameters)
meta <- read.csv("old_parameters_MARTI_samp_comp_read_data/Phyloseq_data/Sample_table.csv")

rownames(meta) <- meta[,1]                              #1st col is rownames
meta$NumReads <- as.numeric(gsub(",","",meta$NumReads)) #Remove commas from the NumRead col

#tables need to be matrices
otu_mat<- as.matrix(otu)
tax_mat<- as.matrix(taxa)


# Phyloseq ----------------------------------------------------------------


#Make phyloseq object
phylo_OTU<- otu_table(otu_mat, taxa_are_rows = TRUE)
phylo_TAX<- tax_table(tax_mat)
phylo_samples <- sample_data(meta)
phylo_object<- phyloseq(phylo_OTU, phylo_TAX, phylo_samples) #Bring them together

#Read data visualisation https://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html ------

#Histogram - distribution of read data
sample_sum_df <- data.frame(sum = sample_sums(phylo_object)) #sample_sums = no assigned reads per sample

ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 15000) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") + ylab("Count") +
  theme(axis.title.y = element_blank())
#Very skewed towards lower read counts, but some pretty long reads also

smin  <- min(sample_sums(phylo_object)) #136
smean <- mean(sample_sums(phylo_object)) #13789.98
smax  <- max(sample_sums(phylo_object)) #97279

#Community composition ----
samp_phylum <- phylo_object %>%
  tax_glom(taxrank = "phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance - compared to pass reads
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa <2%
  arrange(phylum)                                      # Sort data frame alphabetically by phylum

# Plot ------------

# Setting up colour scheme 
phylum_colours <- c('#9e0142', '#d53e4f', '#f46d43', '#fdae61', '#fee08b', '#e6f598', '#abdda4', '#66c2a5',
                    '#3288bd', '#5e4fa2', '#A672A7', '#b5b3bd')

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

# Reorder phylum by total abundance, but keep "Unassigned" last
samp_phylum$phylum <- samp_phylum %>%
  mutate(phylum = fct_reorder(phylum, Abundance, .fun = sum)) %>%
  mutate(phylum = fct_relevel(phylum, "Unassigned", after = Inf)) %>%
  pull(phylum)

samp_phylum$Sampler <- factor(samp_phylum$Sampler, levels = sampler_levels)


# Phylum level plot
phylum_stacked_bar <- ggplot(samp_phylum, aes(x = Sampler, y = Abundance, fill = phylum)) + 
  # facet_grid(Location~., labeller = labeller(Location = location_labels)) +
  facet_grid(rows = vars(Location), cols = vars(Sample_length), 
             labeller = labeller(
               Location = location_labels,
               Sample_length = duration_labels
             )) +
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "Phylum", values = phylum_colours) +
  theme(axis.title.x = element_blank()) +   # Remove x axis title
  ylab("Relative Abundance (Phyla > 2%) \n") +
  scale_x_discrete(labels = c(
    "Compact" = "Coriolis\nCompact",
    "Micro" = "Coriolis\nMicro",
    "Bobcat" = "InnovaPrep\nBobcat",
    "Cub" = "InnovaPrep\nCub",
    "Sass" = "SASS\n3100"
  )) + 
  custom_theme 

  # + ggtitle("Phylum Composition of Different Air Samplers by Sampling Site") 

phylum_stacked_bar

ggsave("Images/graphs_marti_0624/taxonomic_abundance/Phylum_Stacked_Bar.png",
       plot = phylum_stacked_bar, width = 12, height = 8, dpi = 600)



#Taxonomy breakdown---------
df_phylo <- psmelt(phylo_object)

tot_sequences <-  meta %>% 
  group_by(Sampler, Location, Sample_length) %>% 
  summarise(avg_no_reads = mean(NumReads), avg_DNA_yield = mean(DNA_yield))

#Calculating reads per taxonomic group
tot_reads <-  meta %>% 
  group_by(Sample_ID) %>% 
  summarise(Sample = (Sample_ID), no_reads = (NumReads)) %>% 
  #Need to have same header to merge properly later
  select(-Sample_ID)

# Define a function for summarizing data by a specific column value (code from ChatGPT, made what I had a function)
summarise_by_column <- function(data, column_name, values, summary_column_name) {
  data %>%
    filter({{ column_name }} %in% values) %>% 
    group_by(Sample) %>% 
    summarise(!!summary_column_name := sum(Abundance))
}

#Run function for all Kingdoms
Unclassified_sample <- summarise_by_column(df_phylo, kingdom, NA, "unclassified_reads")
Eukaryote_sample    <- summarise_by_column(df_phylo, kingdom, "Eukaryota", "eukaryote_reads")
Bacteria_sample     <- summarise_by_column(df_phylo, kingdom, "Bacteria", "bacteria_reads")
Archaea_sample      <- summarise_by_column(df_phylo, kingdom, "Archaea", "archaea_reads")
Virus_sample        <- summarise_by_column(df_phylo, kingdom, "Viruses", "virus_reads")

fungi_phyla  <- c("Oomycota", "Mucoromycota", "Ascomycota", "Basidiomycota")
algae_phyla  <- c("Bacillariophyta", "Chlorophyta")
animal_phyla <- c("Nematoda", "Mollusca", "Arthropoda", "Chordata")

Plant_sample  <- summarise_by_column(df_phylo, phylum, "Streptophyta", "plant_reads")
Fungi_sample  <- summarise_by_column(df_phylo, phylum, fungi_phyla, "fungi_reads")
Algae_sample  <- summarise_by_column(df_phylo, phylum, algae_phyla, "algae_reads")
Animal_sample <- summarise_by_column(df_phylo, phylum, animal_phyla, "animal_reads")


#Bringing it all together
list_sample <- list(Unclassified_sample, Bacteria_sample, Eukaryote_sample, Archaea_sample,
             Virus_sample, Plant_sample, Fungi_sample, Algae_sample, Animal_sample, tot_reads)

all_reads <- Reduce(function(x, y) merge(x, y, all=TRUE), list_sample, accumulate=FALSE)

#Adding Read proportion-
all_reads <- all_reads %>% 
  mutate(percent_unclassified_reads = (unclassified_reads/no_reads)*100,
         percent_bacteria_reads     = (bacteria_reads/no_reads)*100,
         percent_eukaryote_reads    = (eukaryote_reads/no_reads)*100,
         percent_archaea_reads      = (archaea_reads/no_reads)*100,
         percent_virus_reads        = (virus_reads/no_reads)*100,
         percent_plant_reads        = (plant_reads/no_reads)*100,
         percent_fungi_reads        = (fungi_reads/no_reads)*100,
         percent_algae_reads        = (algae_reads/no_reads)*100,
         percent_animal_reads       = (animal_reads/no_reads)*100)

#Just Kingdom 
kingdom_list <- list(Unclassified_sample, Bacteria_sample, Eukaryote_sample, Archaea_sample,
                     Virus_sample, tot_reads)
kingdom_reads <- Reduce(function(x, y) merge(x, y, all=TRUE), kingdom_list, accumulate=FALSE)

#Eukaryote 
eukaryote_list <- list(Plant_sample, Fungi_sample, Algae_sample, Animal_sample, tot_reads)
eukaryote_reads <- Reduce(function(x, y) merge(x, y, all=TRUE), eukaryote_list, accumulate=FALSE)

#Older data now in Older Marti Runs
#write.csv(all_reads, "MARTI_samp_comp_updated_parameters_0823/taxa_proportions_2108.csv")  
#write.csv(kingdom_reads, "MARTI_samp_comp_updated_parameters_0823/kingdom_proportions_3008.csv")  
#write.csv(eukaryote_reads, "MARTI_samp_comp_updated_parameters_0823/eukaryoteproportions_3008.csv")  

# write.csv(all_reads, "../samp_comp_0624_marti/taxa_proportions_0624.csv")  
# write.csv(kingdom_reads, "../samp_comp_0624_marti/kingdom_proportions_0624.csv")  
# write.csv(eukaryote_reads, "../samp_comp_0624_marti/eukaryoteproportions_0624.csv")  


#Plot Eukaryote phylum proportions----------------------------------------------

phylo_eukaryote = subset_taxa(phylo_object, kingdom=="Eukaryota")

#Relative abundance graphs ----------------------------------------------
#Just animal/fungi/plant/algae

samp_phylum_euk <- phylo_eukaryote %>% 
  tax_glom(taxrank = "phylum") %>%   
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance - compared to pass reads
  psmelt()   %>%                                      
  #Adding new column to look at higher level
  mutate(Taxa = case_when(
    phylum == "Streptophyta" ~ "Plant",
    phylum == "Unassigned" ~ "Unassigned",
    phylum %in% algae_phyla ~ "Algae",
    phylum %in% animal_phyla ~ "Animal",
    phylum %in% fungi_phyla ~ "Fungi"))

taxa_colours <-  c( '#3288bd', '#d53e4f','#fee08b', '#abdda4', '#b5b3bd')

samp_phylum_euk$Sampler <- factor(samp_phylum_euk$Sampler, levels = sampler_levels)


# To also facet by time 
eukaryote_plot <- ggplot(samp_phylum_euk, aes(x = Sampler, y = Abundance, fill = Taxa)) + 
  facet_grid(rows = vars(Location), cols = vars(Sample_length), 
             labeller = labeller(
               Location = location_labels,
               Sample_length = duration_labels
             )) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = taxa_colours) +
  theme(axis.title.x = element_blank()) +
  ylab("Relative Abundance") +
  xlab("Sampler") +
  custom_theme +  
  scale_x_discrete(labels = c(
    "Compact" = "Coriolis\nCompact",
    "Micro" = "Coriolis\nMicro",
    "Bobcat" = "InnovaPrep\nBobcat",
    "Cub" = "InnovaPrep\nCub",
    "Sass" = "SASS\n3100"
  )) 


ggsave("Images/graphs_marti_0624/taxonomic_abundance/Eukaryote_Stacked_Bar.png",
       plot = eukaryote_plot, width = 12, height = 8, dpi = 600)

# To group time and sampler
samp_phylum_euk$Sampler_Length <- paste(samp_phylum_euk$Sampler,
                                        samp_phylum_euk$Sample_length)

eukaryote_taxa_graph <-  ggplot(samp_phylum_euk, aes(x = Sampler_Length, y = Abundance, fill = Taxa)) + 
  facet_grid(rows = vars(Location), labeller = labeller(Location = location_labels)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = taxa_colours) +
  theme(axis.title.x = element_blank()) + # Remove x axis title
  ylab("Relative Abundance") +
  xlab("Sampler and Duration") +
  scale_x_discrete(labels = c(
    "Compact 25" = "Coriolis\nCompact\n25",
    "Micro 25" = "Coriolis\nMicro\n25",
    "Bobcat 25" = "InnovaPrep\nBobcat\n25",
    "Cub 25" = "InnovaPrep\nCub\n25",
    "Sass 25" = "SASS\n3100\n25",
    "Compact 50" = "Coriolis\nCompact\n50",
    "Micro 50" = "Coriolis\nMicro\n50",
    "Bobcat 50" = "InnovaPrep\nBobcat\n50",
    "Cub 50" = "InnovaPrep\nCub\n50",
    "Sass 50" = "SASS\n3100\n50"
  )) + 
  custom_theme 

eukaryote_taxa_graph

ggsave("/Images/graphs_marti_0624/taxonomic_abundance/Eukaryote_Stacked_Bar.png",
       plot = eukaryote_taxa_graph, width = 6, height = 8, dpi = 600)


#All phyla in balck and white, prefer the graph from earlier in the code
eukaryote_phylum_graph <- ggplot(samp_phylum_euk, aes(x = Sampler_Length, y = Abundance, fill = phylum)) + 
  facet_grid(Location~.) +
  geom_bar(stat = "identity") +
  scale_fill_grey(start = 0, end = .9) +
  theme(axis.title.x = element_blank()) + # Remove x axis title
  ylab("Relative Abundance Phylum") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
  ggtitle("Phylum Composition of Eukaryote reads showing different Air Samplers by Sampling Site")

ggsave("../Images/graphs_marti_0624/taxonomic_abundance/Eukaryote_Phylum_Stacked_Bar.svg",
       plot = eukaryote_phylum_graph , device = "svg", width = 10, height = 8)

#Kingdom Proportions-------------
samp_kingdom_rel <- phylo_object %>% 
  tax_glom(taxrank = "kingdom") %>%   
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance - compared to pass reads
  psmelt() #Make long

samp_kingdom_rel$Sampler_Length_ID <- paste(samp_kingdom_rel$Sample_ID,
                                     samp_kingdom_rel$Sample_length)

samp_kingdom_rel$Sampler_Length <- paste(samp_kingdom_rel$Sampler,
                                            samp_kingdom_rel$Sample_length)

kingdom_rel_graph <- ggplot(samp_kingdom_rel, aes(x = Sampler_Length, y = Abundance, fill = kingdom)) + 
  facet_grid(Location~.) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = taxa_colours) +
  theme(axis.title.x = element_blank()) + # Remove x axis title
  ylab("Relative Abundance") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
  ggtitle("Kingdom Composition of different Air Samplers by Sampling Site")

ggsave("../Images/graphs_marti_0624/taxonomic_abundance/Kingdom_Stacked_Bar.svg",
       plot = kingdom_rel_graph , device = "svg", width = 10, height = 8)

#All samples  - Using phylo object but not faceting how I want and doesn't have unclassified reads
ggplot(samp_kingdom_rel, aes(x = Sample_ID, y = Abundance, fill = kingdom)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = taxa_colours) +
  theme(axis.title.x = element_blank()) + # Remove x axis title
  ylab("Relative Abundance") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
  ggtitle("Kingdom Composition of different Air Samplers by Sampling Site")




#Normalised reads graphs --------------------------------------------------------------
samp_phylum_euk_norm <- phylo_eukaryote %>% #Using the already subset data
  tax_glom(taxrank = "phylum") %>%   # agglomerate at phylum level
  psmelt()   %>%                                      # Melt to long format
  mutate(norm_read_count = (Abundance * 100000)/(NumReads)) %>%  #Normalise to read per 100,000 
  #Adding new column to look at higher level
  mutate(Taxa = case_when(
    phylum == "Streptophyta" ~ "Plant",
    phylum %in% algae_phyla ~ "Algae",
    phylum %in% animal_phyla ~ "Animal",
    phylum %in% fungi_phyla ~ "Fungi"))

samp_phylum_euk_norm$Sampler_Length <- paste(samp_phylum_euk_norm$Sampler,
                                        samp_phylum_euk_norm$Sample_length)

#animal/fungi/plant/algae
eukaryote_norm_taxa_graph <-  ggplot(samp_phylum_euk_norm, aes(x = Sampler_Length, y = norm_read_count, fill = Taxa)) + 
  facet_grid(rows = vars(Location)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = taxa_colours) +
  theme(axis.title.x = element_blank()) + # Remove x axis title
  ylab("Reads per 100,000") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
  custom_theme
  # ggtitle("Composition of Eukaryote reads showing different Air Samplers by Sampling Site") 

ggsave("../Images/graphs_marti_0624/taxonomic_abundance/Eukaryote_Count_Stacked_Bar.svg",
       plot = eukaryote_norm_taxa_graph , device = "svg", width = 10, height = 8)


#All phyla
eukaryote_norm_phylum_graph <- ggplot(samp_phylum_euk_norm, aes(x = Sampler_Length, y = Abundance, fill = phylum)) + 
  facet_grid(Location~.) +
  geom_bar(stat = "identity") +
  scale_fill_grey(start = 0, end = .9) +
  theme(axis.title.x = element_blank()) + # Remove x axis title
  ylab("Reads per 100,000") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
  ggtitle("Phylum Composition of Eukaryote reads showing different Air Samplers by Sampling Site")

ggsave("../Images/graphs_marti_0624/taxonomic_abundance/Eukaryote_Phylum_Count_Stacked_Bar.svg",
       plot = eukaryote_norm_phylum_graph  , device = "svg", width = 10, height = 8)

#Kingdom
samp_kingdom_norm <- phylo_object %>% 
  tax_glom(taxrank = "kingdom")   %>%   
  psmelt()                        %>% 
  mutate(norm_read_count = (Abundance * 100000)/(NumReads))

samp_kingdom_norm$Sampler_Length <- paste(samp_kingdom_norm$Sampler,
                                         samp_kingdom_norm$Sample_length)

kingdom_norm_graph <- ggplot(samp_kingdom_norm, aes(x = Sampler_Length, y = Abundance, fill = kingdom)) + 
  facet_grid(Location~.) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = taxa_colours) +
  theme(axis.title.x = element_blank()) + # Remove x axis title
  ylab("Reads per 100,000") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
  ggtitle("Kingdom Composition of different Air Samplers by Sampling Site")

ggsave("../Images/graphs_marti_0624/taxonomic_abundance/Kingdom_Count_Stacked_Bar.svg",
       plot = kingdom_norm_graph , device = "svg", width = 10, height = 8)


# Top 10 taxa per sample --------------------------------------------------
#Decided not to repeat this and instead did it in a sperate script using the summed data - 0624

#devtools::install_github("gmteunisse/fantaxtic")
top_sp <- top_taxa(phylo_object, 
                tax_level = "Species", 
                n_taxa = 10,
                grouping = "Sample_ID")

top10_species <- top_sp$top_taxa %>% 
  reframe(Sample_ID, tax_rank, abundance, Species)

write.csv(top10_species, "MARTI_samp_comp_updated_parameters_0823/top10_species_2508.csv") 

top_phy <- top_taxa(phylo_object, 
                tax_level = "Phylum", 
                n_taxa = 10,
                grouping = "Sample_ID")

top10_phylum <- top_phy$top_taxa %>% 
  reframe(Sample_ID, tax_rank, abundance, Phylum)

write.csv(top10_phylum, "MARTI_samp_comp_updated_parameters_0823/top10_phylum_2508.csv")  

#To plot something similar 
ggplot(top10_species, aes(x = Sample_ID, y = abundance, fill = Species)) +
         geom_col(color = "black") +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) + ylab("Average relative abundance") 

# Pulling out some stats -------

# Define phyla of interest
target_phyla <- c("Streptophyta", "Arthropoda", "Chordata")

# Summarise average relative abundance
summary_table <- samp_phylum_euk %>%
  filter(phylum %in% target_phyla) %>%
  group_by(Location, Sampler, Sample_length, phylum) %>%
  summarise(
    mean_abundance = mean(Abundance, na.rm = TRUE),
    sd_abundance = sd(Abundance, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

summary_table <- summary_table %>%
  mutate(across(contains("abundance"), ~ round(. * 100, 1))) %>%
  arrange(desc(mean_abundance))
