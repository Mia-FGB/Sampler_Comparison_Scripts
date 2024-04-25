
#loading libraries
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(stringi)
library(tidyr)
#Reading in the data file
untouched_data <- read.csv("sampler_comparison_PHIbase_taxaID_readcount.csv", header=TRUE)

#--Normalising the data---------------

#Need to normalise the reads for each location, e.g. hits per million reads
#Therefore we need the total number of reads for each sample, not just the classified reads
data <- 
  untouched_data %>% 
  #only include species woth more than 5 reads
  filter(read_count >=5) %>% 
  mutate(hits_per_100.000 = (read_count * 100000)/(total_read)) %>% 
  #log scale-Mutate the data to add a new row, making R see the barcodes as factors not numbers 
  mutate(log_hits = log10(hits_per_100.000), barcode = factor(barcode))

#Graph
sampler.comparison.heatmap.phi <- ggplot(data = data, mapping = aes(x = barcode,
                                                       y = species,
                                                       fill = hits_per_100.000)) +
  geom_tile() +
  xlab(label = "Barcode") + ylab(label = "Species") +
  scale_fill_distiller(name = "Hits per 100,000", palette = "RdBu") +
  theme_bw() +
  ggtitle(label = "Sampler comparison Species Abundance (hits per 100,000 >5)") +
  #Editing axis title size
  theme(axis.text.y = element_text(vjust = 0.5, hjust=1, size = 15,)) +
  theme(axis.title.y = element_text(size = 20, vjust = 0.5)) +
  theme(axis.title.x = element_text(size = 20, hjust = 0.5)) +
  theme(axis.text.x = element_text(vjust=0.5, size = 15,)) +
  # Change legend key size and key width
  theme(legend.key.size = unit(6, "cm"),legend.key.width = unit(1.5,"cm"))+
  #Change legend font size
  theme(legend.title = element_text(size = 15)) +
  theme(legend.text = element_text(size =15)) +
  #Changing the title size
  theme(plot.title = element_text(size = 20, hjust = 0.5))

sampler.comparison.heatmap.phi
#Saving the plot to right dimensions, much taller y axis
ggsave(filename= "sampler_comparison_heatmap_phibase.png", plot = sampler.comparison.heatmap.phi, width=20, height=20)


