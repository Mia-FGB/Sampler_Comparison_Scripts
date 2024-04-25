#Using data from MARTi run with updated parameters 
#September 2023

#Library -----------------------------------------------------------------------
library(phyloseq)   # Facilitate the import, storage, analysis, and graphical display of microbiome census data.
library(ggplot2)
library(plyr)

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
GP <- phylo_object 

# PCoA  ---------------------------------------------------------
#following http://joey711.github.io/phyloseq/plot_ordination-examples 

#Remove taxa that do not appear more than 3 times in more than 50% of samples
wh0 = genefilter_sample(GP, filterfun_sample(function(x) x > 3), A=0.5*nsamples(GP))
GP1 = prune_taxa(wh0, GP)

#Transform to even sampling depth
GP1 = transform_sample_counts(GP1, function(x) 1E6 * x/sum(x))

#keep 10 most abundant phyla
phylum.sum = tapply(taxa_sums(GP1), tax_table(GP1)[, "Phylum"], sum, na.rm=TRUE)
top10phyla = names(sort(phylum.sum, TRUE))[1:10]
GP1 = prune_taxa((tax_table(GP1)[, "Phylum"] %in% top10phyla), GP1) 
ntaxa(GP1) #only leaves 53 

#Plot taxa and shade by phylum 
GP.ord <- ordinate(GP1, "NMDS", "bray")
p1 = plot_ordination(GP1, GP.ord, type="taxa", color="Phylum", title="taxa")
p1 + facet_wrap(~Phylum, 3) #stops the overplotting

#plot samples and shade by sampler - not splitting out
p2 = plot_ordination(GP1, GP.ord, type="samples", color="Sampler", shape="Location") 
p2 + geom_polygon(aes(fill=Sampler)) + geom_point(size=5) + ggtitle("samples")

p4 = plot_ordination(GP1, GP.ord, type="split", color="Phylum", shape="Sampler", label="Insect", title="split") 
p4
