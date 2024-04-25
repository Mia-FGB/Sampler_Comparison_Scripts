library(vegan)      # Analysis of variance using distance matrices using the adonis2 function
library(Maaslin2)   # Determine multivariable association between metadata and microbial meta-omics features; differential abundance analysis
library(ggplot2)    # Generate visualization plots 
library(ggsignif)   # Visualize comparisons and the significance between two groups
library(dplyr)
library(RColorBrewer)
library(mia)
library(miaViz)
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
counts <- otu            #to make the mia code work

#taxa table 
taxa <- read.csv("MARTI_samp_comp_updated_parameters_0823/Taxonomy/taxa_lineage_sep_na.csv")
taxa<- taxa %>% 
  mutate_if(is.character, as.factor)

rownames(taxa) <- rownames(otu) #same rownames as the otu table
taxa <- taxa[,-1]               
tax <-  taxa                    #to make the mia code work

#sample meta data (this hasn't changed)
meta <- read.csv("old_parameters_MARTI_samp_comp_read_data/Phyloseq_data/Sample_table.csv")

rownames(meta) <- meta[,1]                              #1st col is rownames
meta$NumReads <- as.numeric(gsub(",","",meta$NumReads)) #Remove commas from the NumRead col
samples <-  meta                                        #to make the mia code work


### MIA package ###-------------
counts<- as.matrix(counts)

#data checks
all(rownames(samples) == colnames(counts)) # TRUE
class(samples)                             # DataFrame
all(rownames(tax) == rownames(counts))     # TRUE
class(tax)                                 # Dataframe
class(counts)                              # numeric

# Create a TreeSE
tse_taxa <- TreeSummarizedExperiment(assays =  SimpleList(counts = counts),
                                     colData = DataFrame(samples),
                                     rowData = DataFrame(tax))
#Exploration and QC

#Calculate Frequencies
rowData(tse_taxa)$Kingdom %>% table()
rowData(tse_taxa)$Phylum %>% table()

tse_taxa <- transformCounts(tse_taxa, method = "relabundance") #Add relative Abundance

plotAbundanceDensity(tse_taxa, layout = "jitter", assay_name = "relabundance",
                     n = 40, #Top 40 plotted
                     point_size=1, point_shape=19, point_alpha=0.1) + 
  scale_x_log10(label=scales::percent)


plotAbundanceDensity(tse_taxa, layout = "density", assay_name = "relabundance",
                     n = 5, #Plot top 5
                     colour_by="Sampler", point_alpha=1/10) + scale_x_log10()

top_features <- getTopTaxa(tse_taxa, method="median", top=10) #Pick top 10 taxa
rowData(tse_taxa)[top_features, taxonomyRanks(tse_taxa)]     #Get row data for top 10 taxa

#Library size distribution
tse_taxa <- addPerCellQC(tse_taxa) #Using scater package to calculate QC metrics

#Histogram 
p1 <- ggplot(as.data.frame(colData(tse_taxa))) +
  geom_histogram(aes(x = sum), color = "black", fill = "gray", bins = 20) +
  labs(x = "Library size", y = "Frequency (n)") + 
  # scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), 
  # labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), # Removes the grid
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) # Adds y-axis

#Scatter
df <- as.data.frame(colData(tse_taxa)) %>%
  arrange(sum) %>%
  mutate(index = 1:n())
p2 <- ggplot(df, aes(y = index, x = sum/1e6)) +
  geom_point() +  
  labs(x = "Library size (million reads)", y = "Sample index") +  
  theme_bw() +
  theme(panel.grid.major = element_blank(), # Removes the grid
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) # Adds y-axis

p1 + p2 #Patchwork to plot together


#Plot lib size against variables in the metadata
plotColData(tse_taxa,"sum","Location", colour_by = "Sampler") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(y ="Number of species level assigned reads")

#Taxonomy---
#Observed features as a measure of richness
tse_taxa <- mia::estimateRichness(tse_taxa, 
                                  abund_values  = "counts", 
                                  index = "observed", 
                                  name="observed")

plotColData(tse_taxa, 
            "observed", 
            "Sample_ID", 
            colour_by = "Location",
            shape_by = "Sampler",
            point_size = 4) +
  theme(axis.text.x = element_text(angle=45,hjust=1)) + 
  ylab(expression(Richness[Observed]))
# Alpha diversity --
# Using the counts/rel abund to estimate shannon index 
tse_taxa <- mia::estimateDiversity(tse_taxa, 
                                   abund_values= "counts",
                                   index = "shannon", 
                                   name = "shannon")

# Alpha Diversity ---------------------------------------------------------


tse_taxa <- mia::estimateDiversity(tse_taxa, 
                                   abund_values = "relabundance",
                                   index = "shannon", 
                                   name = "shannon_ra")

df <- as.data.frame(colData(tse_taxa)) # Making a dataframe with col data

#Plotting diversity comparison of Location
df$Location <- factor(df$Location)

# For significance testing, all different combinations are determined
comb_loc <- split(t(combn(levels(df$Location), 2)), 
                  seq(nrow(t(combn(levels(df$Location), 2)))))

ggplot(df, aes(x = Location, y = shannon)) +
  geom_boxplot(outlier.shape = NA) +   # Outliers are removed
  geom_jitter(width = 0.2) + 
  geom_signif(comparisons = comb_loc, map_signif_level = FALSE,
              textsize = 3, vjust = -0.5, step_increase = 0.1) +
  theme(text = element_text(size = 10))


tse_taxa <- estimateEvenness(tse_taxa, 
                             assay.type =  "counts", 
                             index="simpson")

tse_taxa <- estimateDominance(tse_taxa, 
                              assay.type = "counts", 
                              index="relative")

tse_taxa <- mia::estimateDiversity(tse_taxa, 
                                   assay.type = "counts",
                                   index = "log_modulo_skewness")

tse_taxa <- mia::estimateDivergence(tse_taxa,
                                    assay.type = "counts",
                                    reference = "median",
                                    FUN = vegan::vegdist)

plots <- lapply(c("observed", "shannon", "simpson", "relative", "log_modulo_skewness"), #Plot all diversity measures
                plotColData,
                object = tse_taxa,
                x = "Sampler",
                colour_by = "Sampler")

plots <- lapply(plots, "+", 
                theme(axis.text.x = element_blank(),
                      axis.title.x = element_blank(),
                      axis.ticks.x = element_blank()))

((plots[[1]] | plots[[2]] | plots[[3]]) / 
    (plots[[4]] |  plots[[5]])) +
  plot_layout(guides = "collect")



#Beta Diversity -----------------------------------------------------------------

tse_taxa <- runMDS(tse_taxa,    
                   FUN = vegan::vegdist, #using vegdist function of vegan package
                   method = "bray",
                   name = "PCoA_BC",
                   exprs_values = "relabundance")

#sample dissimilarity visualised on 2D scale using plotReducedDim function 
p <- plotReducedDim(tse_taxa, "PCoA_BC",
                    colour_by = "Location",
                    shape_by = "Sampler",
                    point_size = 5)

# Calculate explained variance
e <- attr(reducedDim(tse_taxa, "PCoA_BC"), "eig")
rel_eig <- e / sum(e[e > 0])

# Add explained variance for each axis
p <- p + labs(x = paste("PCoA 1 (", round(100 * rel_eig[[1]], 1), "%", ")", sep = ""),
              y = paste("PCoA 2 (", round(100 * rel_eig[[2]], 1), "%", ")", sep = "")) +
  ggtitle("PCoA MARTi data (relative abundance, BC)")
p

#Comparing different beta diversity measures ------
tse_taxa <- runMDS(tse_taxa,
                   FUN = vegan::vegdist,
                   name = "MDS_euclidean",
                   method = "euclidean",
                   exprs_values = "counts")

tse_taxa <- runNMDS(tse_taxa,
                    FUN = vegan::vegdist,
                    name = "NMDS_BC")
tse_taxa <- runNMDS(tse_taxa,
                    FUN = vegan::vegdist,
                    name = "NMDS_euclidean",
                    method = "euclidean")

plots <- lapply(c("PCoA_BC", "MDS_euclidean", "NMDS_BC", "NMDS_euclidean"),
                plotReducedDim,
                object = tse_taxa,
                colour_by = "Location",
                shape_by = "Sampler",
                point_size = 2,
                point_alpha = 1)

((plots[[1]] | plots[[2]]) / (plots[[3]] | plots[[4]])) +
  plot_layout(guides = "collect")
