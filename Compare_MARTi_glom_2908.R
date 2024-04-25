#Working on things Richard and I discussed 29/08

#Library -----------------------------------------------------------------------
library(phyloseq)
library(dplyr)
library(stringr)


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


#Checking phylum agglomeration ---------------------------------------
phylum_check <- phylo_object %>% 
  tax_glom(taxrank = "Phylum") %>%   # agglomerate at phylum level
  psmelt() 

phylum_check %>% 
  filter(Phylum =="Actinomycetota") %>% 
  summarise(tot_art = sum(Abundance))

#Checking kingdom agglomeration ---------------------------------------
kingdom_check <- phylo_object %>% 
  tax_glom(taxrank = "Kingdom") %>%   # agglomerate at phylum level
  psmelt() 

kingdom_check %>% 
  filter(Kingdom =="Viruses") %>% 
  summarise(tot = sum(Abundance))


# See if Taxa are in both data --------------------------------------------
#MARTi data 
marti <- read.csv("MARTI_samp_comp_updated_parameters_0823/compare_taxa_assignments_lca_0.1_all_levels_2023-AUG-10_16-8-58.csv")

marti_kingdom <- marti %>% 
  filter(NCBI.Rank == "kingdom")

