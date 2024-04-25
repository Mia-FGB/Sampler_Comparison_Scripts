library("phyloseq")
library("ggplot2")
library("dplyr")
library("ape")

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

# Tree --------------------------------------------------------------------
#https://joey711.github.io/phyloseq/import-data  
#All taxa, too many to be useful
random_tree <-  rtree(ntaxa(phylo_object), rooted=TRUE, tip.label=taxa_names(phylo_object))
plot(random_tree)

#Merge this into the phyloseq object
physeq_tree <-  merge_phyloseq(phylo_object, phylo_samples, random_tree)
plot_tree(physeq_tree, color="Location", label.tips="taxa_names", ladderize="left", plot.margin=0.3)

#https://www.bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-analysis.html#exploratory-tree-plots
SC.Euk <- subset_taxa(phylo_object, Kingdom=="Eukaryota")
random_euk_tree <-  rtree(ntaxa(SC.Euk), rooted=TRUE, tip.label=taxa_names(SC.Euk))
physeq_euk_tree <-  merge_phyloseq(SC.Euk, phylo_samples, random_euk_tree)

plot_tree(physeq_euk_tree, color="Sampler", shape="Location", label.tips="Phylum")


# Without phyloseq --------------------------------------------------------


