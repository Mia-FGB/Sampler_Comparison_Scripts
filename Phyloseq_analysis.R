library(phyloseq)   # Facilitate the import, storage, analysis, and graphical display of microbiome census data.
library(vegan)      # Analysis of variance using distance matrices using the adonis2 function
library(Maaslin2)   # Determine multivariable association between metadata and microbial meta-omics features; differential abundance analysis
library(ggplot2)    # Generate visualization plots 
library(ggsignif)   # Visualize comparisons and the significance between two groups
library(microbiome) # Tools for microbiome analysis
library(dplyr)
library(RColorBrewer)

# Set plotting theme
theme_set(theme_bw())

#Creating a phyloseq object----------
#following this tutorial https://mvuko.github.io/meta_phyloseq/
#using summed species counts from MARTi

#OTU table 
otu <- read.csv("read_data/MARTI_samp_comp_read_data/Phyloseq_data/OTU_table.csv")
#was being weird so i removed my rownames to use these 
rownames(otu) <- paste0("OTU", 1:nrow(otu))

#taxa table 
taxa <- read.csv("read_data/MARTI_samp_comp_read_data/Phyloseq_data/taxa_table.csv")
#want taxa to be recognised as factors
taxa<- taxa %>% 
  mutate_if(is.character, as.factor)
#same rownames as the otu table
rownames(taxa) <- rownames(otu)

#sample data 
meta <- read.csv("read_data/MARTI_samp_comp_read_data/Phyloseq_data/Sample_table.csv")
#weird empty row on the bottom
meta <- meta %>%  na.omit
#think the 1st col need to be the rownames
rownames(meta) <- meta[,1]
#remove 1st col so not repeated
#meta[,1] <- NULL

#combining into phyloseq
#tables need to be matrices
otu_mat<- as.matrix(otu)
tax_mat<- as.matrix(taxa)

#transform data to phyloseq object
phylo_OTU<- otu_table(otu_mat, taxa_are_rows = TRUE)
phylo_TAX<- tax_table(tax_mat)
phylo_samples <- sample_data(meta)

#and put them in one object
phylo_object<- phyloseq(phylo_OTU, phylo_TAX, phylo_samples)

#checking it all works
plot_bar(phylo_object, fill = "Family")
plot_heatmap(phylo_object, taxa.label = "Kingdom")
plot_bar(phylo_object, "Kingdom", fill="Phylum", facet_grid=~Sampler)


#Diversity -----
#another tutorial https://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html 

#read data -Histogram to see distribution 
#sample_sums is a phyloseq function, counts number of assigned reads per sample
sample_sum_df <- data.frame(sum = sample_sums(phylo_object))
#Can then plot this dataframe 
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 1500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") + ylab("Count") +
  theme(axis.title.y = element_blank())

#this data does not have a normal distribution, skewed towards low read counts, only a few samples have high read counts

# mean, max and min of sample read counts
smin <- min(sample_sums(phylo_object)) #37
smean <- mean(sample_sums(phylo_object)) #2301.325
smax <- max(sample_sums(phylo_object)) #13549


#Stacked barplot of phyla to get a sense of community composition ----
#too many to reasonably distinguish between diff colours so only include those that contribute more than 2% relative abundance 
samp_phylum <- phylo_object %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance - compared to pass reads
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum

# Plot 
Phylum_colours <- c('#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#42d4f4', '#f032e6', '#fabed4', 
                    '#469990', '#dcbeff', '#9A6324', '#800000', '#aaffc3', '#000075')

ggplot(samp_phylum, aes(x = Sampler, y = Abundance, fill = Phylum)) + 
  facet_grid(Location~.) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  ylab("Relative Abundance (Phyla > 2%) \n") +
  ggtitle("Phylum Composition of Different Air Samplers by Sampling Site") 


### Unconstrained ordination------

#PCoA of communities
#Scale reads to even depth 

# Scales reads by 
# 1) taking proportions,
# 2) multiplying by a given library size of n
# 3) rounding down
scale_reads <- function(physeq, n) {
  physeq.scale <-
    transform_sample_counts(physeq, function(x) {
      (n * x/sum(x))
    })
  otu_table(physeq.scale) <- floor(otu_table(physeq.scale))
  physeq.scale <- prune_taxa(taxa_sums(physeq.scale) > 0, physeq.scale)
  return(physeq.scale)
}

#have chosen 500 completely randomly here - need to find a better way to normalise
#plot looks pretty similar with 2
samp_scale <- phylo_object %>%
  scale_reads(200) 

# Ordinate
samp_pcoa <- ordinate(
  physeq = samp_scale, 
  method = "PCoA", 
  distance = "bray"
)

#can play around with all the metadat when plotting this
# Plot 
plot_ordination(
  physeq = samp_scale,
  ordination = samp_pcoa,
  color = "Location",
  shape = "Sampler",
  title = "PCoA of Sampler Comparison (Bray Curtis)"
) + 
  scale_color_manual(values = c( '#e6194B', '#3cb44b')
  ) +
  geom_point(aes(color = Location, size = Air_volume), alpha = 0.9)

#Trying NMDS instead of PCoA
set.seed(1)

# Ordinate
samp_nmds <- ordinate(
  physeq = samp_scale, 
  method = "NMDS", 
  distance = "bray"
)
#Run 20 stress 0.1755471 
#generally anything below .2 is considered an acceptable stress for a NMDS plot

#PCoA plot--------
plot_ordination(
  physeq = samp_scale,
  ordination = samp_nmds,
  color = "Location",
  shape = "Sampler",
  title = "NMDS of Sampler Comparison (Bray Curtis)"
) + 
  scale_color_manual(values = c( '#e6194B', '#3cb44b')
  ) +
  geom_point(aes(color = Location, size = Air_volume))


#Permanova 
#Hyp: Different samplers have the same centroids (? - think it might be to do with clustering)
#Not sure hoe to interpret the results of this

set.seed(1)

# Calculate bray curtis distance matrix
samp_bray <- phyloseq::distance(samp_scale, method = "bray")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(phylo_object))

# Adonis test
adonis(samp_bray ~ Sampler, data = sampledf)

# Homogeneity of dispersion test
beta <- betadisper(samp_bray, sampledf$Sampler)
permutest(beta)
             
#Alpha Diversity ------
#This method involves subsampling the libraries with replacement to estimate the species abundance of the real population while standardising sampling effort 
min_lib <- min(sample_sums(phylo_object))

#We will subsample to 37 the minimum number of reads.
#We will repeat this 100 times and average the diversity estimates from each trial.

# Initialize matrices to store richness and evenness estimates
nsamp = nsamples(phylo_object)
trials = 100

richness <- matrix(nrow = nsamp, ncol = trials)
row.names(richness) <- sample_names(phylo_object)

evenness <- matrix(nrow = nsamp, ncol = trials)
row.names(evenness) <- sample_names(phylo_object)

# It is always important to set a seed when you subsample so your result is replicable 
set.seed(3)

for (i in 1:100) {
  # Subsample
  r <- rarefy_even_depth(phylo_object, sample.size = min_lib, verbose = FALSE, replace = TRUE)
  
  # Calculate richness
  rich <- as.numeric(as.matrix(estimate_richness(r, measures = "Observed")))
  richness[ ,i] <- rich
  
  # Calculate evenness
  even <- as.numeric(as.matrix(estimate_richness(r, measures = "InvSimpson")))
  evenness[ ,i] <- even
}

#Calculate the Mean and sd per sample for observed richness and inverse simpsons index to store in a df
# Create a new dataframe to hold the means and standard deviations of richness estimates
Sample <- row.names(richness)
mean <- apply(richness, 1, mean)
sd <- apply(richness, 1, sd)
measure <- rep("Richness", nsamp)
rich_stats <- data.frame(Sample, mean, sd, measure)

# Create a new dataframe to hold the means and standard deviations of evenness estimates
Sample <- row.names(evenness)
mean <- apply(evenness, 1, mean)
sd <- apply(evenness, 1, sd)
measure <- rep("Inverse Simpson", nsamp)
even_stats <- data.frame(Sample, mean, sd, measure)

#combine richness and eveness into one df
alpha <- rbind(rich_stats, even_stats)

#add sample metadata using merge 
s <- data.frame(sample_data(phylo_object))
alphadiv <- merge(alpha, s, by = "Sample") 

#Plot alpha diversity measures for different samplers in a facet
ggplot(alphadiv, aes(x = Sampler, y = mean, color = Location, group = Location, shape = Location)) +
  geom_point(size = 3) + 
  facet_wrap(~measure, ncol = 1, scales = "free") +
  scale_color_manual(values = c('#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231')) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

