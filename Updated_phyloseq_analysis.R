#Using data from MARTi run with updated parameters 
#August 2023

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

# Set plotting theme
theme_set(theme_bw())

#Loading in data -----------------------------------------------------------------------


#OTU table 
otu <- read.csv("MARTI_samp_comp_updated_parameters_0823/marti_output_august23_assigned.csv")
rownames(otu) <- otu[,1] #Making the rownames the NCBI ID
otu <- otu[,-1]          #remove 1st first col

#taxa table 
taxa <- read.csv("MARTI_samp_comp_updated_parameters_0823/Taxonomy/taxa_lineage_sep_na.csv")
taxa<- taxa %>% 
  mutate_if(is.character, as.factor)

rownames(taxa) <- rownames(otu) #same rownames as the otu table
taxa <- taxa[,-1]               

#sample meta data (this hasn't changed)
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
  geom_histogram(color = "black", fill = "indianred", binwidth = 10000) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") + ylab("Count") +
  theme(axis.title.y = element_blank())
#Very skewed towards lower read counts, but some pretty long reads also

smin  <- min(sample_sums(phylo_object)) #136
smean <- mean(sample_sums(phylo_object)) #13789.98
smax  <- max(sample_sums(phylo_object)) #97279

#Community composition ----
samp_phylum <- phylo_object %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance - compared to pass reads
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa <2%
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum

# Plot 
phylum_colours <- c('#9e0142', '#d53e4f', '#f46d43', '#fdae61', '#fee08b', '#e6f598', '#abdda4', '#66c2a5', 
                    '#3288bd', '#5e4fa2', '#b5b3bd')

phylum_stacked_bar <- ggplot(samp_phylum, aes(x = Sampler, y = Abundance, fill = Phylum)) + 
  facet_grid(Location~.) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colours) +
  theme(axis.title.x = element_blank()) +   # Remove x axis title
  ylab("Relative Abundance (Phyla > 2%) \n") +
  ggtitle("Phylum Composition of Different Air Samplers by Sampling Site") 




#Taxonomy breakdown---------
df_phylo <- psmelt(phylo_object)

tot_sequences <-  meta %>% 
  group_by(Sampler, Location, Sample_length) %>% 
  summarise(avg_no_reads = mean(NumReads), avg_DNA_yield = mean(DNA_yield))

# Taxonomy proportion


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
Unclassified_sample <- summarise_by_column(df_phylo, Kingdom, NA, "unclassified_reads")
Eukaryote_sample    <- summarise_by_column(df_phylo, Kingdom, "Eukaryota", "eukaryote_reads")
Bacteria_sample     <- summarise_by_column(df_phylo, Kingdom, "Bacteria", "bacteria_reads")
Archaea_sample      <- summarise_by_column(df_phylo, Kingdom, "Archaea", "archaea_reads")
Virus_sample        <- summarise_by_column(df_phylo, Kingdom, "Viruses", "virus_reads")

fungi_phyla  <- c("Oomycota", "Mucoromycota", "Ascomycota", "Basidiomycota")
algae_phyla  <- c("Bacillariophyta", "Chlorophyta")
animal_phyla <- c("Nematoda", "Mollusca", "Arthropoda", "Chordata")

Plant_sample  <- summarise_by_column(df_phylo, Phylum, "Streptophyta", "plant_reads")
Fungi_sample  <- summarise_by_column(df_phylo, Phylum, fungi_phyla, "fungi_reads")
Algae_sample  <- summarise_by_column(df_phylo, Phylum, algae_phyla, "algae_reads")
Animal_sample <- summarise_by_column(df_phylo, Phylum, animal_phyla, "animal_reads")


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

write.csv(all_reads, "MARTI_samp_comp_updated_parameters_0823/taxa_proportions_2108.csv")  
write.csv(kingdom_reads, "MARTI_samp_comp_updated_parameters_0823/kingdom_proportions_3008.csv")  
write.csv(eukaryote_reads, "MARTI_samp_comp_updated_parameters_0823/eukaryoteproportions_3008.csv")  


#Plot Eukaryote phylum proportions----------------------------------------------

phylo_eukaryote = subset_taxa(phylo_object, Kingdom=="Eukaryota")

#Relative abundance graphs ----------------------------------------------

samp_phylum_euk <- phylo_eukaryote %>% 
  tax_glom(taxrank = "Phylum") %>%   
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance - compared to pass reads
  psmelt()   %>%                                      
  #Adding new column to look at higher level
  mutate(Taxa = case_when(
    Phylum == "Streptophyta" ~ "Plant",
    Phylum %in% algae_phyla ~ "Algae",
    Phylum %in% animal_phyla ~ "Animal",
    Phylum %in% fungi_phyla ~ "Fungi"))


samp_phylum_euk$Sampler_Length <- paste(samp_phylum_euk$Sampler,
                                        samp_phylum_euk$Sample_length)

taxa_colours <-  c( "#1E88E5", "#D81B60","#FFC107", "#0A794D")

#Just animal/fungi/plant/algae
eukaryote_taxa_graph <-  ggplot(samp_phylum_euk, aes(x = Sampler_Length, y = Abundance, fill = Taxa)) + 
  facet_grid(rows = vars(Location)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = taxa_colours) +
  theme(axis.title.x = element_blank()) + # Remove x axis title
  ylab("Relative Abundance") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
  ggtitle("Composition of Eukaryote reads showing different Air Samplers by Sampling Site") 

#All phyla
eukaryote_phylum_graph <- ggplot(samp_phylum_euk, aes(x = Sampler_Length, y = Abundance, fill = Phylum)) + 
  facet_grid(Location~.) +
  geom_bar(stat = "identity") +
  scale_fill_grey(start = 0, end = .9) +
  theme(axis.title.x = element_blank()) + # Remove x axis title
  ylab("Relative Abundance Phylum") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
  ggtitle("Phylum Composition of Eukaryote reads showing different Air Samplers by Sampling Site")

#Kingdom Proportions-------------
samp_kingdom_rel <- phylo_object %>% 
  tax_glom(taxrank = "Kingdom") %>%   
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance - compared to pass reads
  psmelt() #Make long

samp_kingdom_rel$Sampler_Length_ID <- paste(samp_kingdom_rel$Sample_ID,
                                     samp_kingdom_rel$Sample_length)

kingdom_rel_graph <- ggplot(samp_kingdom_rel, aes(x = Sampler_Length, y = Abundance, fill = Kingdom)) + 
  facet_grid(Location~.) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = taxa_colours) +
  theme(axis.title.x = element_blank()) + # Remove x axis title
  ylab("Relative Abundance") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
  ggtitle("Kingdom Composition of different Air Samplers by Sampling Site")

#All samples  - Using phylo object but not faceting how I want and doesn't have unclassified reads
ggplot(samp_kingdom_rel, aes(x = Sample_ID, y = Abundance, fill = Kingdom)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = taxa_colours) +
  theme(axis.title.x = element_blank()) + # Remove x axis title
  ylab("Relative Abundance") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
  ggtitle("Kingdom Composition of different Air Samplers by Sampling Site")




#Normalised reads graphs --------------------------------------------------------------
samp_phylum_euk_norm <- phylo_eukaryote %>% #Using the already subset data
  tax_glom(taxrank = "Phylum") %>%   # agglomerate at phylum level
  psmelt()   %>%                                      # Melt to long format
  mutate(norm_read_count = (Abundance * 100000)/(NumReads)) %>%  #Normalise to read per 100,000 
  #Adding new column to look at higher level
  mutate(Taxa = case_when(
    Phylum == "Streptophyta" ~ "Plant",
    Phylum %in% algae_phyla ~ "Algae",
    Phylum %in% animal_phyla ~ "Animal",
    Phylum %in% fungi_phyla ~ "Fungi"))

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
  ggtitle("Composition of Eukaryote reads showing different Air Samplers by Sampling Site") 


#All phyla
eukaryote_norm_phylum_graph <- ggplot(samp_phylum_euk_norm, aes(x = Sampler_Length, y = Abundance, fill = Phylum)) + 
  facet_grid(Location~.) +
  geom_bar(stat = "identity") +
  scale_fill_grey(start = 0, end = .9) +
  theme(axis.title.x = element_blank()) + # Remove x axis title
  ylab("Reads per 100,000") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
  ggtitle("Phylum Composition of Eukaryote reads showing different Air Samplers by Sampling Site")

#Kingdom
samp_kingdom_norm <- phylo_object %>% 
  tax_glom(taxrank = "Kingdom")   %>%   
  psmelt()                        %>% 
  mutate(norm_read_count = (Abundance * 100000)/(NumReads))

samp_kingdom_norm$Sampler_Length <- paste(samp_kingdom_norm$Sampler,
                                         samp_kingdom_norm$Sample_length)

kingdom_norm_graph <- ggplot(samp_kingdom_norm, aes(x = Sampler_Length, y = Abundance, fill = Kingdom)) + 
  facet_grid(Location~.) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = taxa_colours) +
  theme(axis.title.x = element_blank()) + # Remove x axis title
  ylab("Reads per 100,000") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
  ggtitle("Kingdom Composition of different Air Samplers by Sampling Site")



# Top 10 taxa per sample --------------------------------------------------

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
