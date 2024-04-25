library(mia)
library(miaViz)
library(ggplot2)
library(BiocManager)
library(Biostrings)
library(dplyr)
library(knitr)
library(scater)
library(patchwork)
library(ggsignif)
library(vegan)
library(stringr)

#Following this tutorial https://microbiome.github.io/OMA/ 
#Using some hacks from last time I set up the Phyloseq object in Phyloseq_analysis.R script


#Loading data ----

#Abundance table, 1st column is feature IDs, first row is sample IDs, other vakues are numeric
# Abundance table (e.g. ASV data; to assay data)
counts  <- read.csv("read_data/MARTI_samp_comp_read_data/Phyloseq_data/OTU_table.csv") 
rownames(counts) <- counts[,1]
counts <- counts[-1]

#Row data first column provides feature IDs, first row provides column headers, taxonomic mapping
# Taxonomy table (to rowData)
tax    <- read.csv("read_data/MARTI_samp_comp_read_data/Phyloseq_data/taxa_table.csv") 
rownames(tax) <- rownames(counts)

#Column data, first col sample IDsm first row cikumn headers, contains metadata
# Sample data (to colData)
samples <- read.csv("read_data/MARTI_samp_comp_read_data/Phyloseq_data/Sample_table.csv")
rownames(samples) <- samples[,1]
#weird empty row on the bottom
samples <- samples %>%  na.omit

# Let us ensure that the data is in correct (numeric matrix) format:
counts <- as.matrix(counts)

#Other checks
# coldata rownames match assay colnames
all(rownames(samples) == colnames(counts)) # our dataset should be true
class(samples) # should be data.frame or DataFrame
# rowdata rownames match assay rownames
all(rownames(tax) == rownames(counts)) # our dataset should be true
class(tax) # should be data.frame or DataFrame
# Counts 
class(counts) # should be a numeric matrix

#Constructing TreeSummarizedExperiment--------

# Create a TreeSE
tse_taxa <- TreeSummarizedExperiment(assays =  SimpleList(counts = counts),
                                     colData = DataFrame(samples),
                                     rowData = DataFrame(tax))

tse_taxa

#Subsetting 4.1.2----------
# Show dimensions (features x samples)
dim(tse_taxa)

# Inspect possible values for flow rate (can use this technique to look for anything in the metadata)
unique(tse_taxa$Flow_rate)
unique(tse_taxa$Sampler)

# Show the frequency of each value
tse_taxa$Flow_rate %>% table()


#Say we are only interested in samples from NHM we can subset them (could subset by anything in the metadata)
# Subset by sample
tse_subset_NHM <- tse_taxa[ , tse_taxa$Location %in% "NHM"]
tse_subset_CF <- tse_taxa[ , tse_taxa$Location %in% "CFarm"]

# Show dimensions - now only 20 samples
dim(tse_subset_NHM)

#Can also subset by feature 
# Inspect possible values for phylum
unique(rowData(tse_taxa)$Phylum)

#Can get Phylum frequency 
#Would be helpful to get this from various subsets to compare
rowData(tse_taxa)$Kingdom %>% table()

# Subset by feature
tse_subset_Fungi<- tse_taxa[rowData(tse_taxa)$Kingdom %in% "Fungi" & !is.na(rowData(tse_taxa)$Kingdom), ]

# Show dimensions
dim(tse_subset_Fungi)

#Subset by sample & feature
# Subset by sample and feature and remove NAs
tse_subset_by_Viridiplantae_CFarm <- tse_taxa[rowData(tse_taxa)$Kingdom %in% "Viridiplantae" 
                                         & !is.na(rowData(tse_taxa)$Kingdom),
                                         tse_taxa$Location %in% "Cfarm"]
dim(tse_subset_by_Viridiplantae_CFarm)

#Chp 5 - Exploration & QC----------

#Abundance
#add relative abundances - MARGIN = samples by default

tse_taxa <- transformCounts(tse_taxa, method = "relabundance")

#plot these relative abundances
#n = 40 - looking at top 40
plotAbundanceDensity(tse_taxa, layout = "jitter", abund_values = "relabundance",
                     n = 40, point_size=1, point_shape=19, point_alpha=0.1) + 
  scale_x_log10(label=scales::percent)


#Plotting relative abundance of top 5 taxa as a density plot, with sampler indicated by colour
plotAbundanceDensity(tse_taxa, layout = "density", abund_values = "relabundance",
                     n = 5, colour_by="Sampler", point_alpha=1/10) +
  scale_x_log10()

#Same but coloured by location and showing top 10
plotAbundanceDensity(tse_taxa, layout = "density", abund_values = "relabundance",
                     n = 10, colour_by="Location", point_alpha=1/10) +
  scale_x_log10()

#5.2 Prevalence-----

#freq of samples where certain microbes were detected, given as sample size or percentage 
#investigating prevalence allows focus on changes which affect most samples or identify rare microbes 

#Looking at frequency at 1% relative abundance 
head(getPrevalence(tse_taxa, detection = 1/100, sort = TRUE, as_relative = TRUE))

#Prevalence analysis - investigate microbiome prevalence at a selected taxonomic level 

# Agglomerate taxa abundances to Phylum level, and add the new table, to the altExp slot
altExp(tse_taxa,"Phylum") <- agglomerateByRank(tse_taxa, "Phylum")
# Check prevalence for the Phylum abundance table from the altExp slot
head(getPrevalence(altExp(tse_taxa,"Phylum"), detection = 1/100, sort = TRUE,
                   abund_values = "counts", as_relative = TRUE))

#Getting rare taxa - ither functions to do this 
getRareTaxa(tse_taxa, detection = 0, prevalence = 50/100,
                         rank = "Phylum", sort = TRUE)

#Plotting prevalence

#analysing the Phylum level abundances
rowData(altExp(tse_taxa,"Phylum"))$prevalence <- 
  getPrevalence(altExp(tse_taxa,"Phylum"), detection = 1/100, sort = FALSE,
                abund_values = "counts", as_relative = TRUE)
#5.3 QC---------
# Pick the top taxa
top_features <- getTopTaxa(tse_taxa, method="median", top=10)

# Check the information for these
rowData(tse_taxa)[top_features, taxonomyRanks(tse_taxa)]

#total counts calculated using function from scater package - only including classified reads not total reads
perCellQCMetrics(tse_taxa)

#Can add this to the tse object and keep it in sample data 
tse_taxa <- addPerCellQC(tse_taxa)
colData(tse_taxa)

#distribution of calculated library sizes can be visualised as a histogram 
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

#Or as a scatter plot
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

library(patchwork)
p1 + p2

#Can also plot library size against other variables in the metadata e.g sampler
plotColData(tse_taxa,"sum","Sampler", colour_by = "Air_volume") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(y ="Number of species level assigned reads")

plotColData(tse_taxa,"sum","Location", colour_by = "Sampler") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(y ="Number of species level assigned reads")

#Remove insect samples first 
#Didn't seem to affect the graph, not sure why
tse_subset_not_insect<- tse_taxa[rowData(tse_taxa)$Class != "Insecta" , ]
#removes about 100 species
dim(tse_subset_not_insect)
dim(tse_taxa)
#Calculating the sum values
tse_subset_not_insect <- addPerCellQC(tse_subset_not_insect)
colData(tse_subset_not_insect)
#Now plotting this data
plotColData(tse_subset_not_insect,"sum","Sampler", colour_by = "Sampler") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(y ="Number of species level assigned reads")


## Chapter 6 - Taxonomic information ## --------
# checking the taxa data
checkTaxonomy(tse_taxa)
taxonomyRanks(tse_taxa)
#find proportion of taxa we have the info for - in this case none are empty 
table(taxonomyRankEmpty(tse_taxa, rank = "Species"))

#6.3 Data agglomeration -----
#agglomerate count data on taxonomic levels and track influence of changing conditions through these level s
#Use the agglomerateByRank function and store data as alternative experiment 

tse_taxa <- transformCounts(tse_taxa, abund_values = "counts", method = "relabundance")
altExp(tse_taxa, "Family") <- agglomerateByRank(tse_taxa, rank = "Family",
                                           agglomerateTree = TRUE)
altExp(tse_taxa, "Family")

#6.5 Data transformation --
#Abundances of all taxa in one sample 
taxa.abund.NHM_Micro_2 <- getAbundanceSample(tse_taxa, 
                                     sample_id = "NHM_Micro_2",
                                     abund_values = "counts")
taxa.abund.NHM_Micro_2[1:10]

#6.5 Abundance of specific taxa in all samples 
taxa.abundances <- getAbundanceFeature(tse_taxa, 
                                       feature_id = "Aspergillus wentii",
                                       abund_values = "counts")
taxa.abundances[1:10]


#Chapter 7 - Community Diversity ------------

#7.1.1 Richness ------
#Calculating the number of features present within a community 
#Observed features as a measure of richness
#Results are added to colData of the tse object so we can plot this output
tse_taxa <- mia::estimateRichness(tse_taxa, 
                             abund_values  = "counts", 
                             index = "observed", 
                             name="observed")

head(colData(tse_taxa)$observed)

#Plotting observed richness
plotColData(tse_taxa, 
            "observed", 
            "Sample_ID", 
            colour_by = "Sampler",
            shape_by = "Sample_length") +
  theme(axis.text.x = element_text(angle=45,hjust=1)) + 
  ylab(expression(Richness[Observed]))


#7.1.2 - Diversity ------

#Using the counts data to estimate shannon index 
tse_taxa <- mia::estimateDiversity(tse_taxa, 
                              abund_values = "counts",
                              index = "shannon", 
                              name = "shannon")
#same but using relative abundance values
tse_taxa <- mia::estimateDiversity(tse_taxa, 
                                   abund_values = "relabundance",
                                   index = "shannon", 
                                   name = "shannon_ra")
#relative abundance & counts produce same shannon diversity 
tail(colData(tse_taxa)$shannon)
tail(colData(tse_taxa)$shannon_ra)

#Plotting Diversity --------
#Comparing samplers 
#Making a dataframe with the col data
df <- as.data.frame(colData(tse_taxa))

df$Sampler <- factor(df$Sampler)

# For significance testing, all different combinations are determined
comb <- split(t(combn(levels(df$Sampler), 2)), 
              seq(nrow(t(combn(levels(df$Sampler), 2)))))

#Too many comparisons for this to look good as a graph better for comparing just 2 things 
ggplot(df, aes(x = Sampler, y = shannon)) +
  # Outliers are removed, because otherwise each data point would be plotted twice; 
  # as an outlier of boxplot and as a point of dotplot.
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2) + 
  geom_signif(comparisons = comb, map_signif_level = FALSE,
              textsize = 3, vjust = -0.5, step_increase = 0.1) +
  theme(text = element_text(size = 10))

#Same again but comparing location 
df$Location <- factor(df$Location)

# For significance testing, all different combinations are determined
comb_loc <- split(t(combn(levels(df$Location), 2)), 
              seq(nrow(t(combn(levels(df$Location), 2)))))

#Too many comparisons for this to look good as a graph better for comparing just 2 things 
ggplot(df, aes(x = Location, y = shannon)) +
  # Outliers are removed, because otherwise each data point would be plotted twice; 
  # as an outlier of boxplot and as a point of dotplot.
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2) + 
  geom_signif(comparisons = comb_loc, map_signif_level = FALSE,
              textsize = 3, vjust = -0.5, step_increase = 0.1) +
  theme(text = element_text(size = 10))

#Same again but comparing sample length
df$Sample_length <- factor(df$Sample_length)

# For significance testing, all different combinations are determined
comb_len <- split(t(combn(levels(df$Sample_length), 2)), 
                  seq(nrow(t(combn(levels(df$Sample_length), 2)))))

#Too many comparisons for this to look good as a graph better for comparing just 2 things 
ggplot(df, aes(x = Sample_length, y = shannon)) +
  # Outliers are removed, because otherwise each data point would be plotted twice; 
  # as an outlier of boxplot and as a point of dotplot.
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2) + 
  geom_signif(comparisons = comb_len, map_signif_level = FALSE,
              textsize = 3, vjust = -0.5, step_increase = 0.1) +
  theme(text = element_text(size = 10))

#7.1.4 Evenness - Simpsons diversity ---------

tse_taxa <- estimateEvenness(tse_taxa, 
                        abund_values = "counts", 
                        index="simpson")
head(colData(tse_taxa)$simpson)

#Dominance 
tse_taxa <- estimateDominance(tse_taxa, 
                             abund_values = "counts", 
                             index="relative")
head(colData(tse_taxa)$relative)

#Rarity
tse_taxa <- mia::estimateDiversity(tse_taxa, 
                              abund_values = "counts",
                              index = "log_modulo_skewness")
#Divergence
tse_taxa <- mia::estimateDivergence(tse_taxa,
                               abund_values = "counts",
                               reference = "median",
                               FUN = vegan::vegdist)
#Plot all diversity measures -------
plots <- lapply(c("observed", "shannon", "simpson", "relative", "log_modulo_skewness"),
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

#Chapter 8 - Beta Diversity ---------
#Unsupervised ordination 
#analyse beta diversity in the data, observing vatiation between CFarm samples and NHM
#calculate rel abundance 
tse_taxa <- transformCounts(tse_taxa,
                       method = "relabundance")

#Comparing communities by beta analysis 
#visual representation of the groups by 2D ordination , then estimate rel abundances and MDS ordination from Bray-Curtis 
#dissimilarity is calculated with function supplied to FUN argument 
#This whole thing isn't working, ned to troubleshoot when I'm more alert
# Perform PCoA
tse_taxa <- runMDS(tse_taxa,    
              FUN = vegan::vegdist, #using vegdist function of vegan package
              #using bray-curtis index to calculate
              method = "bray",
              name = "PCoA_BC",
              exprs_values = "relabundance")


#sample dissimilarity visualised on 2D scale using plotReducedDim function 
#this function allows us to colour by/ shape by things
# Create ggplot object
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
  ggtitle("PCoA Species level MARTi data (relative abundance, BC)")
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
                 point_alpha = 1) #stops the points being transparent
 
 ((plots[[1]] | plots[[2]]) / (plots[[3]] | plots[[4]])) +
   plot_layout(guides = "collect")
 
 ##Straight PCA ordination - Similar to Euclidean distance 
 tse_taxa <- runPCA(tse_taxa,
               name = "PCA",
               exprs_values = "counts",
               ncomponents = 10)
 plotReducedDim(tse_taxa, "PCA",
                colour_by = "Location",
                shape_by = "Sampler",
                point_size = 5)

 #8.2 Supervised ordination----------- 
 # Covariates that are being analyzed - currently the way this script is set up it only works with 3 variables - not sure why yet
variable_names <- c("Sampler", "Location", "Sample_length")
 
#Apply relative transform
tse_taxa <- transformCounts(tse_taxa, method = "relabundance")
 
# Create a formula - all variables will be included as predictors in the formula
formula <- as.formula(paste0("assay ~ ", str_c(variable_names, collapse = " + ")) )
 
 # # Perform RDA- this doesn't work

rda <- calculateRDA(tse_taxa,
                     abund_values = "relabundance",
                    formula = formula,
                    distance = "bray",
                    na.action = na.exclude)

 
#Get the rda object
rda <- attr(rda, "rda")
# Calculate p-value and variance for whole model
set.seed(436) #setting a random number that means you will get reproducible results each time you run the model
permanova <- anova.cca(rda, permutations = 999)

# Create a data.frame for results
rda_info <- as.data.frame(permanova)["Model", ]

# Calculate p-value and variance for each variable
# by = "margin" --> the order or variables does not matter
set.seed(4585)
permanova <- anova.cca(rda, by = "margin",  permutations = 999)
# Add results to data.frame
rda_info <- rbind(rda_info, permanova)

# Add info about total variance
rda_info[ , "Total variance"] <- rda_info["Model", 2] +
  rda_info["Residual", 2]

# Add info about explained variance
rda_info[ , "Explained variance"] <- rda_info[ , 2] / 
  rda_info[ , "Total variance"]

# Loop through variables, calculate homogeneity
homogeneity <- list()
# Get colDtaa
coldata <- colData(tse_taxa)
# Get assay
assay <- t(assay(tse_taxa, "relabundance"))


for( variable_name in rownames(rda_info) ){ #looping through variable names in rda_info object
  # If data is continuous or discrete
  if( variable_name %in% c("Model", "Residual") || #checking if variable name data is continuous or discrete
      length(unique(coldata[[variable_name]])) / #ratio of uniq values to total values in colData
      length(coldata[[variable_name]]) > 0.2 ){ #is greater than 0.2
    # Do not calculate homogeneity for continuous data - setting at NA if these things are true
    temp <- NA
  } else{ #if these things aren't true and the data is discrete then...
    # Calculate homogeneity for discrete data
    # Calculate homogeneity
    set.seed(413)
    temp <- anova( #applying anova to the results of BC on groups
      betadisper(  #this function computes multivariate homogenity based on BC
        vegdist(assay, method = "bray"),
        group = coldata[[variable_name]] ), #defining the groups for analysis
      permutations = permutations )["Groups", "Pr(>F)"] #specifying no of permutations and indexing the p-value
  }
  # Add info to the list -storing homogeneity info
  homogeneity[[variable_name]] <- temp
}
# Add homogeneity to information
rda_info[["Homogeneity p-value (NULL hyp: distinct/homogeneous --> permanova suitable)"]] <-
  homogeneity

knitr::kable(rda_info)

# Load ggord for plotting
library(ggord)

# Since na.exclude was used, if there were rows missing information, they were 
# dropped off. Subset coldata so that it matches with rda.
coldata <- coldata[ rownames(rda$CCA$wa), ]

# Adjust names
# Get labels of vectors
vec_lab_old <- rownames(rda$CCA$biplot)

# Loop through vector labels
vec_lab <- sapply(vec_lab_old, FUN = function(name){
  # Get the variable name
  variable_name <- variable_names[ str_detect(name, variable_names) ]
  # If the vector label includes also group name
  if( !any(name %in% variable_names) ){
    # Get the group names
    group_name <- unique( coldata[[variable_name]] )[ 
      which( paste0(variable_name, unique( coldata[[variable_name]] )) == name ) ]
    # Modify vector so that group is separated from variable name
    new_name <- paste0(variable_name, " \U2012 ", group_name)
  } else{
    new_name <- name
  }
  # Add percentage how much this variable explains, and p-value
  new_name <- expr(paste(!!new_name, " (", 
                         !!format(round( rda_info[variable_name, "Explained variance"]*100, 1), nsmall = 1), 
                         "%, ",italic("P"), " = ", 
                         !!gsub("0\\.","\\.", format(round( rda_info[variable_name, "Pr(>F)"], 3), 
                                                     nsmall = 3)), ")"))
  
  return(new_name)
})
# Add names
names(vec_lab) <- vec_lab_old

# Create labels for axis
xlab <- paste0("RDA1 (", format(round( rda$CCA$eig[[1]]/rda$CCA$tot.chi*100, 1), nsmall = 1 ), "%)")
ylab <- paste0("RDA2 (", format(round( rda$CCA$eig[[2]]/rda$CCA$tot.chi*100, 1), nsmall = 1 ), "%)")

# Create a plot        
plot_location <- ggord(rda, grp_in = coldata[["Location"]], vec_lab = vec_lab,
              alpha = 0.5,
              size = 4, addsize = -4,
              #ext= 0.7, 
              txt = 3.5, repel = TRUE, 
              #coord_fix = FALSE
) + 
  # Adjust titles and labels
  guides(colour = guide_legend("Location"),
         fill = guide_legend("Location"),
         group = guide_legend("Location"),
         shape = guide_legend("Location"),
         x = guide_axis(xlab),
         y = guide_axis(ylab)) +
  theme( axis.title = element_text(size = 10) )
plot_location

plot_sampler <- ggord(rda, grp_in = coldata[["Sampler"]], vec_lab = vec_lab,
                       alpha = 0.5,
                       size = 4, addsize = -4,
                       #ext= 0.7, 
                       txt = 3.5, repel = TRUE, 
                       #coord_fix = FALSE
) + 
  # Adjust titles and labels
  guides(colour = guide_legend("Sampler"),
         fill = guide_legend("Sampler"),
         group = guide_legend("Sampler"),
         shape = guide_legend("Sampler"),
         x = guide_axis(xlab),
         y = guide_axis(ylab)) +
  theme( axis.title = element_text(size = 10) )
plot_sampler

#8.3.1 Visualising most dominant genus in PCoA ---------
#For each sample plots which genus is most abundant
# Agglomerate to genus level
tse_genus <- agglomerateByRank(tse_taxa,
                               rank = "Genus")

# Convert to relative abundances
tse_genus <- transformCounts(tse_taxa,
                             method = "relabundance",
                             abund_values = "counts")

# Add info on dominant genus per sample
tse_genus <- addPerSampleDominantTaxa(tse_genus,
                                      abund_values = "relabundance",
                                      name = "dominant_taxa")
#perform PCoA with BC
tse_genus <- runMDS(tse_genus,
                    FUN = vegan::vegdist,
                    name = "PCoA_BC_genus",
                    exprs_values = "relabundance")

# Getting the top taxa
top_taxa <- getTopTaxa(tse_genus,
                       top = 6,
                       abund_values = "relabundance")

# Naming all the rest of non top-taxa as "Other"
most_abundant <- lapply(colData(tse_genus)$dominant_taxa,
                        function(x) {if (x %in% top_taxa) {x} else {"Other"}})

# Storing the previous results as a new column within colData
colData(tse_genus)$most_abundant <- as.character(most_abundant)

# Calculating percentage of the most abundant
most_abundant_freq <- table(as.character(most_abundant))
most_abundant_percent <- round(most_abundant_freq / sum(most_abundant_freq) * 100, 1)

# Retrieving the explained variance
e <- attr(reducedDim(tse_genus, "PCoA_BC_genus"), "eig")
var_explained <- e / sum(e[e > 0]) * 100

# Define colors for visualization
my_colors <- c("black", "blue", "lightblue", "darkgray", "magenta", "darkgreen", "red")

# Visualization
plot <-plotReducedDim(tse_genus, "PCoA_BC_genus",
                      colour_by = "most_abundant", point_size = 4) +
  scale_colour_manual(values = my_colors,
                      labels = paste0(names(most_abundant_percent), "(", most_abundant_percent, "%)")) +
  labs(x = paste("PC 1 (", round(var_explained[1], 1), "%)"),
       y = paste("PC 2 (", round(var_explained[2], 1), "%)"),
       color = "")

plot

##-------
## Chapter 9 - Community Composition ## ---------
##-------

#Composition barchart--------

# Computing relative abundance
tse <- transformCounts(tse_taxa, abund_values = "counts", method = "relabundance")

# Getting top taxa on a Phylum level - just top 5
tse_phylum <- agglomerateByRank(tse, rank ="Phylum", onRankOnly=TRUE)
top_taxa <- getTopTaxa(tse_phylum,top = 5, assay.type = "relabundance")

# Renaming the "Phylum" rank to keep only top taxa and the rest to "Other"
phylum_renamed <- lapply(rowData(tse)$Phylum,
                         function(x){if (x %in% top_taxa) {x} else {"Other"}})
rowData(tse)$Phylum <- as.character(phylum_renamed)

# Visualizing the composition barplot, with samples order by "Bacteroidetes"
#need to add x axis labels
plotAbundance(tse, assay_name="relabundance", rank = "Phylum",
              order_rank_by="abund", 
              order_sample_by = "Location", features = "Sampler", add_legend = FALSE,
              add_x_text = TRUE)[[1]] +
  theme(axis.text.x = element_text(angle = 90))

#Composition heatmap-----
# Add clr-transformation on samples
tse_phylum <- transformCounts(tse_phylum, abund_values = "counts",
                              method = "relabundance", pseudocount = 1)

tse_phylum <- transformCounts(tse_phylum, assay.type = "relabundance",
                              method = "clr", pseudocount = 1)

# Add z-transformation on features (taxa)
tse_phylum <- transformCounts(tse_phylum, assay.type = "clr", 
                              MARGIN = "features",
                             method = "z", name = "clr_z")

#Plotting heatmap
# Melt the assay for plotting purposes
df <- meltAssay(tse_phylum, assay.type = "clr_z")

# Determines the scaling of colours
maxval <- round(max(abs(df$clr_z)))
limits <- c(-maxval, maxval)
breaks <- seq(from = min(limits), to = max(limits), by = 0.5)
colours <- c("darkblue", "blue", "white", "red", "darkred")

# Creates a ggplot object
ggplot(df, aes(x = SampleID, y = FeatureID, fill = clr_z)) +
  geom_tile() +
  scale_fill_gradientn(name = "CLR + Z transform", 
                       breaks = breaks, limits = limits, colours = colours) + 
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=45, hjust=1),
        legend.key.size = unit(1, "cm")) +
  labs(x = "Samples", y = "Taxa")


#Chapter 11 Differential abundance analysis-----------
#For this example I am going compare the two locations
#In the future could expand to compare the samplers

#Packages
library(ALDEx2)
library(Maaslin2)
library(MicrobiomeStat)
library(ANCOMBC)
library(GUniFrac)
library(tibble)

# set random seed because some tools can randomly vary and then produce different results:
set.seed(13253)

# Make location a factor
tse_taxa$Location <- factor(tse_taxa$Location)

# how many observations do we have per group?
as.data.frame(colData(tse_taxa)) %>% 
  count(Location) %>%
  kable()

#11.1.2 Prevalence filtering-----
#Data manipulation pre analysis, can alter these thresholds in the future 
#Applying them now means they will be used in all data

tse_genus_10 <- agglomerateByRank(tse_taxa, rank = "Genus") %>% #agglomerating to Genus
  transformCounts(abund_values = "counts",
                  method = "relabundance",
                  MARGIN = "samples") %>%
  # subset based on the relative abundance assay              
  subsetByPrevalentTaxa(detection = 0,
                        prevalence = 10/100, #only including taxa at 10% or greater
                        abund_values = "relabundance")

# Add also clr abundances
tse_genus_10 <- transformCounts(tse_taxa, method="clr", pseudocount=1) 

#ALDEx2------

# Generate Monte Carlo samples of the Dirichlet distribution for each sample.
# Convert each instance using the centered log-ratio transform.
# This is the input for all further analyses.
set.seed(254)
x <- aldex.clr(assay(tse_genus_10), tse_genus_10$Location)     

#The t-test
# calculates expected values of the Welch's t-test and Wilcoxon rank
# test on the data returned by aldex.clr
x_tt <- aldex.ttest(x, paired.test = FALSE, verbose = FALSE)

#Effect sizes
# Determines the median clr abundance of the feature in all samples and in groups,
# the median difference between the two groups,
# the median variation within each group and the effect size,
# which is the median of the ratio of the between group difference and the larger of the variance within groups

x_effect <- aldex.effect(x, CI = TRUE, verbose = FALSE)

# combine all outputs 
aldex_out <- data.frame(x_tt, x_effect)

#Plotting
#Create a Bland-Altman plot which shows association between relative abundance and magnitude of diff per sample (left)
#Also plot dispersion on x axis against log-ratio abundance (right)
#red dots represent genera that are differently abundant between the 2 groups
#black points are rare taxa
#grey points are abundant taxa 

par(mfrow = c(1, 2))
aldex.plot(
  aldex_out, 
  type = "MA", 
  test = "welch", 
  xlab = "Log-ratio abundance",
  ylab = "Difference",
  cutoff = 0.05
)
aldex.plot(
  aldex_out, 
  type = "MW", 
  test = "welch",
  xlab = "Dispersion",
  ylab = "Difference",
  cutoff = 0.05
)

#Evaluation of differential abundant comes from the p-value 
#developers say statistically different features are those wherer 95% of the CI of the effect sie does not cross 0
#This is indicated in the overlap col, which indicates proportion of overlap

rownames_to_column(aldex_out, "Genus") %>%
  filter(wi.eBH <= 0.05)  %>% # here we chose the wilcoxon output rather than tt
  dplyr::select(Genus, we.eBH, wi.eBH, effect, overlap) %>%
  kable()

#ANOM-BC--------

# Agglomerate data to genus level and add this new abundance table to the altExp slot
altExp(tse_taxa, "Genus") <- agglomerateByRank(tse_taxa, "Genus")

# Identify prevalent genera, setting at 20% for this analysis
prevalent.genera <- getPrevalentFeatures(altExp(tse_taxa, "Genus"), detection = 0, prevalence = 20/100)

# Run ANCOM-BC at the genus level and only including the prevalent genera
#getting an error I don't understand - skipping for now
out <- ancombc2(
  data = altExp(tse_taxa, "Genus")[prevalent.genera, ],
  assay_name = "counts", 
  fix_formula = "Location", 
  p_adj_method = "fdr", 
  prv_cut = 0,
  group = "Locationn", 
  struc_zero = TRUE, 
  neg_lb = TRUE,
  global = TRUE # multi group comparison will be deactivated automatically 
)


#MaAsLin2
# maaslin expects features as columns and samples as rows 
# for both the abundance table as well as metadata 

# We can specify different GLMs/normalizations/transforms.
# Let us use similar settings as in Nearing et al. (2021):
maaslin2_out <- Maaslin2(
  t(assay(tse_taxa)),
  data.frame(colData(tse_taxa)),
  output = "DAA example",
  transform = "AST",
  fixed_effects = "Location",
  # random_effects = c(...), # you can also fit MLM by specifying random effects
  # specifying a ref is especially important if you have more than 2 levels
  reference = "Location,NHM",  
  normalization = "TSS",
  standardize = FALSE,
  min_prevalence = 0 # prev filterin already done
)

#Creates a table for every taxon, head shows top 5
kable(head(filter(maaslin2_out$results, qval <= 0.05)))


#LinDA
meta <- as.data.frame(colData(tse_taxa)) %>% dplyr::select(Location)
linda.res <- linda(
  as.data.frame(assay(tse_taxa)), 
  meta, 
  formula = '~Location', 
  alpha = 0.05, 
  prev.filter = 0, 
  mean.abund.filter = 0)

linda_out <- linda.res$output$LocationNHM

# to scan the table for genera where H0 could be rejected:
kable(head(filter(as.data.frame(linda_out), reject)))

#ZicoSeq--------
set.seed(123)
meta <- as.data.frame(colData(tse_taxa))
zicoseq.obj <- GUniFrac::ZicoSeq(meta.dat = meta, 
                                 feature.dat = as.matrix(assay(tse_taxa)),
                                 grp.name = 'Location',
                                 adj.name = NULL, 
                                 feature.dat.type = 'count',
                                 prev.filter = 0,
                                 perm.no = 999,
                                 mean.abund.filter = 0,
                                 max.abund.filter = 0,
                                 return.feature.dat = T)

zicoseq_out <- cbind.data.frame(p.raw=zicoseq.obj$p.raw, p.adj.fdr=zicoseq.obj$p.adj.fdr) 
kable(head(filter(zicoseq_out, p.adj.fdr<0.05)))

#Confounding variables ------
#LinDA Confounding --------
linda_cov <- linda(
  as.data.frame(assay(tse_taxa, "counts")), 
  as.data.frame(colData(tse_taxa)), 
  formula = '~ Location + Sampler + Air_volume + Sample_length' , 
  prev.filter = 0, 
  mean.abund.filter = 0)

linda.res <- linda_cov$output$LocationNHM
kable(head(filter(linda.res, reject==T)))

#ZicoSeq confounding-----
set.seed(123)
zicoseq.obj <- GUniFrac::ZicoSeq(meta.dat = as.data.frame(colData(tse_taxa)) , 
                                 feature.dat = as.matrix(assay(tse_taxa)),
                                 grp.name = 'Location',
                                 adj.name = 'Air_volume', 
                                 feature.dat.type = 'count',
                                 prev.filter = 0,
                                 perm.no = 999,
                                 mean.abund.filter = 0,
                                 max.abund.filter = 0,
                                 return.feature.dat = T)


zicoseq_out <- cbind.data.frame(p.raw=zicoseq.obj$p.raw,
                                p.adj.fdr=zicoseq.obj$p.adj.fdr) 
kable(head(filter(zicoseq_out, p.adj.fdr<0.05)))

#Chapter 14 - Visualisation ---------
#Plot QC data -----------

# obtain QC data
tse_taxa <- addPerCellQC(tse_taxa)
tse_taxa <- addPerFeatureQC(tse_taxa)


# plot QC Mean against Species
plotRowData(tse_taxa, "mean", "Species") +
  theme(axis.text.x = element_blank()) +
  labs(x = "Species", y = "QC Mean")

# plot QC Sum against Sample ID, colour-labeled by Sample Type
plotColData(tse_taxa, "sum", "Sample_ID", colour_by = "Sampler") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Sample ID", y = "QC Sum")

#Instead convert to a df to plot with ggplot
# store colData into a data frame
coldata <- as.data.frame(colData(tse_taxa))
# plot Number of Samples against Sampler used- equal in my case so not needed
ggplot(coldata, aes(x = Sampler)) +
  geom_bar(width = 0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Air Sampler",
       y = "Number of Samples")

#Diversity estimation------

#Alpha diversity 
#estimate shannon diversity index
tse_taxa <- mia::estimateDiversity(tse_taxa, 
                              abund_values = "counts",
                              index = "shannon", 
                              name = "shannon")
# plot shannon diversity index, colour-labeled by Sampler
plotColData(tse_taxa, "shannon", colour_by = "Sampler")



