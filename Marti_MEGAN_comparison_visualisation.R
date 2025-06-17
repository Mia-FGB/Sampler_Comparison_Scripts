library(ggplot2)
library(tidyr)
library(dplyr)
library(reshape2)
library(viridis)

#Script for MARTi and MEGAN comparison
# notes here - https://docs.google.com/document/d/1KHRCVzLtTlWvmLMS4uyc59MQRDNI5T8OwjvozabSrfc/edit?tab=t.0
#  Updated this June 2025

custom_theme <- theme_minimal(base_size = 12) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.3),
    panel.grid.minor = element_blank(),
  )

data <- read.csv("Marti_megan_comparison_barcode-14/marti_megan_read_counts_barcode14.csv")


# Filter the data to remove assigned read counts and just use summarised
data_filtered <- data[, !grepl("As", colnames(data))]

#Rename the data
colnames(data_filtered) <- gsub("_Su$", "", colnames(data_filtered))  # remove _Su at the end
colnames(data_filtered) <- gsub("^ME$", "MEGAN", colnames(data_filtered))  # rename MEGAN

#need to reshape the data
data_long <- data_filtered %>% 
  pivot_longer(cols = -Species, names_to = "RunName", values_to = "Read_Count") %>% 
  filter(Species != "Unclassified")

#Custim order
custom_order <- c("MEGAN", "MA_LCA_100_MI_60", "MA_LCA_100_MI_70", "MA_LCA_20_MI_60", "MA_LCA_20_MI_70")
data_long$RunName <- factor(data_long$RunName, levels = custom_order)


#Calculate the average read counts for each species
species_averages <- data_long %>%
  group_by(Species) %>%
  summarise(Average_Read_Count = mean(Read_Count))

# Order the species based on the highest average read count
# don't want this for heatmap
ordered_species <- species_averages %>%
  arrange(desc(Average_Read_Count)) %>%
  pull(Species)

# Reorder the species factor based on the highest average read count
data_long$Species <- factor(data_long$Species, levels = ordered_species)

# Define custom color palette
my_colors <- c("#449C75", "#B5713E", "#ECBC97", "#25509C", "#437EE8") 

#barplot for each species, ordered by read count
ggplot(data_long, aes(x = Species, y = Read_Count, fill = RunName)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = my_colors, name = "RunName") +
  custom_theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) 

# Plotting just top 15 species
top_species <- data_long %>%
  group_by(Species) %>%
  summarise(total_reads = sum(Read_Count)) %>%
  top_n(15, total_reads) %>%
  pull(Species)

top_15 <- data_long %>% filter(Species %in% top_species)

top_15_plot <- ggplot(top_15, aes(x = Species, y = Read_Count, fill = RunName)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = my_colors, name = "RunName") +
  custom_theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1))  +
  labs(
    y = "Read Count",
  )
  
ggsave("images/marti_megan_comparison/top15_readcount.pdf", top_15_plot, width=8, height=6, dpi = 300 )

#Heatmap
heatmap <- ggplot(data_long, aes(x = RunName, y = Species, fill = Read_Count)) +
  geom_tile() +
  scale_fill_viridis_c(
    trans = "log1p",  # log1p handles 0 safely (log1p(x) = log(1 + x))
    name = "Read Count (log)",
    breaks = c(1, 10, 100, 1000, 10000),
    labels = scales::comma_format()(c(1, 10, 100, 1000, 10000)),
    option = "D"  # Must be a string, e.g. "D" not just D
  ) +
  custom_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("images/marti_megan_comparison/heatmap_0625.pdf", heatmap, width=8, height=6, dpi = 300 )


