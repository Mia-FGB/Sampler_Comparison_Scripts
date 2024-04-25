library(ggplot2)
library(tidyr)
library(dplyr)
library(reshape2)

data <- read.csv("read_data/Marti_megan_comparison_barcode-14/marti_megan_read_counts_barcode14.csv")


#need to reshape the data
data_long <- data %>% 
  pivot_longer(cols = -Species, names_to = "Software", values_to = "Read_Count") %>% 
  filter(Species != "Unclassified")

#Calculate the average read counts for each species
species_averages <- data_long %>%
  group_by(Species) %>%
  summarise(Average_Read_Count = mean(Read_Count))

# Order the species based on the highest average read count
#don't want this for heatmap
ordered_species <- species_averages %>%
  arrange(desc(Average_Read_Count)) %>%
  pull(Species)

# Reorder the species factor based on the highest average read count
data_long$Species <- factor(data_long$Species, levels = ordered_species)

# Define custom color palette
my_colors <- c("#449C75", "#5AE8A9", "#25509C", "#437EE8", "#B5713E", "#ECBC97", "#8B57A1", "#DBAFED"
               , "#ED9151", "#F7AF54")


#barplot for each species, ordered by read count
ggplot(data_long, aes(x = Species, y = Read_Count, fill = Software)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = my_colors, name = "Software")
  

#Trying other graphs to better show the parameters 
#logging the data for better spread
data_long$log_Read_Count <- log(data_long$Read_Count + 1)

data_long$Software <- factor(data_long$Software ,
                             levels = c("MA_LCA_20_MI_70_As", "MA_LCA_20_MI_70_Su","MA_LCA_20_MI_60_As", "MA_LCA_20_MI_60_Su", 
                                        "MA_LCA_100_MI_70_As", "MA_LCA_100_MI_70_Su", "MA_LCA_100_MI_60_As", "MA_LCA_100_MI_60_Su", 
                                        "ME_As", "ME_Su"))

#Heatmap
#think about how to group the data along the y axis
ggplot(data_long, aes(x = Software, y = Species, fill = log_Read_Count)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#449C75") +
  labs(title = "Metagenomic Taxonomic Classification",
       x = "Software and Thresholds",
       y = "Species") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





