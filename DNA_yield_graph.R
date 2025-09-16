# R Script to plot DNA Yield data from the sampler comparison experiment 

# packages 
library(dplyr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(tidyr)
library(patchwork)
library(showtext)

# Use a font that supports Unicode
font_add_google("Roboto")
showtext_auto()

#read in data
data <- read.csv("metadata/sampler_comparison_DNA_yield.csv")
#rename sampler
data$Air_Sampler <- gsub("Micro", "μ", data$Air_Sampler)

# volume of air sampled
data_air_vol <- data %>% 
  mutate(air_volume = Collection_time * Collection_rate) %>% 
  # dna_yield numeric
  mutate(DNA_yield = as.numeric(DNA_yield), 
         genera_per_1000_reads = as.numeric(genera_per_1000_reads)) %>% 
  # Replace NA values with 0.1 so they log okay
  replace_na(list(
    DNA_yield = 0.1,
    genera_per_1000_reads = 0.1
  )) %>%
  # log the dna yield and air volume
  mutate(log_DNA_yield = log10(DNA_yield), log_air_volume = log10(air_volume),
         log_genera_per_1000_reads = log10(genera_per_1000_reads)) 
  

# Create Insect label col
data_air_vol$Insect_Label <- ifelse(data_air_vol$Insect == "Yes", "*", NA)

#Group data 7 remove those containing Insect 
insect_remove_normalised <- subset(data_air_vol, !grepl("Yes", Insect)) %>% 
  mutate(normalised_yield = DNA_yield/air_volume)


# Theme for plots ------
custom_theme <- theme_minimal(base_size = 12) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.3),
    panel.grid.minor = element_blank(),
  )

# air_sampler_palette <- c(
#   "InnovaPrep Bobcat" = "#FF7F50", 
#   "Coriolis Compact" = "#008080", 
#   "Coriolis μ" = "#FFD700", 
#   "InnovaPrep Cub" = "#8A2BE2", 
#   "SASS 3100" = "#00CED1"
# )

air_sampler_palette <- setNames(brewer.pal(5, "Set2"), 
                                c("InnovaPrep Bobcat", "Coriolis Compact", 
                                  "Coriolis μ", "InnovaPrep Cub", 
                                  "SASS 4100"))

# Get some summary stats to report
summary <- data_air_vol %>% 
  group_by(Air_Sampler, Collection_time) %>% 
  summarise(
    n = n(),
    mean_yield = mean(DNA_yield, na.rm = TRUE),
    sd_yield   = sd(DNA_yield, na.rm = TRUE),
    se_yield   = sd(DNA_yield, na.rm = TRUE) / sqrt(n())
  )

# Insect removed
summary_no_insect <- insect_remove_normalised %>% 
  group_by(Air_Sampler, Collection_time) %>% 
  summarise(
    n = n(),
    mean_yield = mean(DNA_yield, na.rm = TRUE),
    se_yield   = sd(DNA_yield, na.rm = TRUE) / sqrt(n()),
    mean_norm_yield = mean(normalised_yield, na.rm = TRUE),
    se_norm_yield   = sd(normalised_yield, na.rm = TRUE) / sqrt(n())
  )

# Insect samples 
insect_samples <- subset(data_air_vol, grepl("Yes", Insect))

# Plots ---------
# Scatter graph of DNA yield against volume of air sampled 
ggplot(data_air_vol, aes(x = air_volume, y= DNA_yield, colour = Air_Sampler)) +
  geom_point() + custom_theme

#DNA Yield against Air samples -----------------

# # Samplers are different shapes
#  ggplot(data_air_vol, aes(x = air_volume, y= DNA_yield)) +
#   geom_point(aes(shape = Air_Sampler, colour = Location), size = 5) +
#   #so colours match my other graphs
#   scale_color_manual(values = c("#FF7F50", "#008080")) + 
#   #change the point shapes to match other graphs
#   scale_shape_manual(values = c(19, 17, 15, 3, 4)) +
#   theme_bw() +
#   #Also find a way to highlight which samples contained insects
#   geom_text(aes(label = Insect_Label, hjust = -0.25, size = 5)) +
#   #axis labels
#   labs(x = "Log Volume of air sampled (L)", y = "Log DNA yield (ng)", title = "DNA yield (ng)") +
#   custom_theme +
#   scale_y_log10() + scale_x_log10(labels = scales::comma)


# Samplers diff colours ----
dna_yield <- ggplot(data_air_vol, aes(x = air_volume, y = DNA_yield)) +
  
  # Translucent points with coloured border
  geom_point(aes(fill = Air_Sampler, shape = Location, colour = Air_Sampler),
             size = 5, alpha = 0.7, stroke = 1) +
  
  # Manual fill and border colour palette
  scale_fill_manual(name = "Air Sampler", values = air_sampler_palette) +
  
  # Add visible asterisk where insect was present (no legend here)
  geom_text(data = subset(data_air_vol, Insect == "Yes"),
            aes(label = Insect_Label),
            size = 8, hjust = -0.4, show.legend = FALSE) +
  
  # Colour and shape scales
  scale_color_manual(name = "Air Sampler", values = air_sampler_palette) +
  scale_shape_manual(name = "Location", values = c(
    "Church Farm" = 21, "NHM" = 22, "Yes" = 8),
    labels = c("Church Farm", "NHM")) +
  
  # Labels and theme
  labs(x = "Volume of air sampled (L), log scale", 
       y = "DNA yield (ng)",
       colour = "Air Sampler") +
  custom_theme +
  scale_x_log10(labels = scales::comma)

no_insect_dna_yield <- ggplot(insect_remove_normalised, aes(x = air_volume, y = DNA_yield)) +
  
  # Translucent points with coloured border
  geom_point(aes(fill = Air_Sampler, shape = Location, colour = Air_Sampler),
             size = 5, alpha = 0.7, stroke = 1) +
  
  # Manual fill and border colour palette
  scale_fill_manual(name = "Air Sampler", values = air_sampler_palette) +
  
  # Colour and shape scales
  scale_color_manual(name = "Air Sampler", values = air_sampler_palette) +
  scale_shape_manual(name = "Location", values = c(
    "Church Farm" = 21, "NHM" = 22, "Yes" = 8),
    labels = c("Church Farm", "NHM")) +
  
  # Labels and theme
  labs(x = "Volume of air sampled (L), log scale", 
       y = "DNA yield (ng)",
       colour = "Air Sampler") +
  custom_theme +
  scale_x_log10(labels = scales::comma)

no_insect_dna_yield

#Combine plots
com_dna_yield <-  dna_yield + no_insect_dna_yield +
  plot_annotation(tag_levels = "A") & 
  plot_layout(guides = "collect") &
  theme(legend.position = "right")


# Save plots
ggsave("images/DNA_yield/no_inse_yield_air_vol_log.pdf", no_insect_dna_yield, width=10, height=6, dpi = 300 )
ggsave("images/DNA_yield/yield_air_vol_log.pdf", dna_yield, width=10, height=6, dpi = 300 )
ggsave("images/DNA_yield/Comb_DNA_yield_scatter_log.pdf", com_dna_yield, width=12, height=6, dpi = 300 )



#Plotting genera count per 1000 reads----------------------

#Simple graph
ggplot(data_air_vol, aes(x = air_volume, y= genera_per_1000_reads, colour = Air_Sampler)) +
  geom_point(aes(shape = Air_Sampler, colour = Location), size = 5) + theme_bw() 

#Log genera count  graph 
# sampler = shape, location = colour
ggplot(data_air_vol, aes(x = air_volume, y= genera_per_1000_reads)) +
  geom_point(aes(shape = Air_Sampler, colour = Location), size = 5) + 
  #so colours match my other graphs
  scale_color_manual(values = c("#FF7F50", "#008080")) + 
  #change the point shapes to match other graphs
  scale_shape_manual(values = c(19, 17, 15, 3, 4)) +
  theme_bw() +
  #Also find a way to highlight which samples contained insects
  #geom_text(aes(label = Insect, hjust = -0.25)) +
  #axis labels
  labs(x = "Volume of air sampled (L)", y = "Genera per 1000 reads", title = "Diversity accumulation") +
  theme(plot.title = element_text(hjust = 0.5))



# Thesis ready graph - sampler = colour, location = shape
genera_plot <- ggplot(data_air_vol, aes(x = air_volume, y = genera_per_1000_reads)) +
  
  # Translucent points with coloured border
  geom_point(aes(fill = Air_Sampler, shape = Location, colour = Air_Sampler),
             size = 5, alpha = 0.7, stroke = 1) +
  
  # Manual fill and border colour palette
  scale_fill_manual(name = "Air Sampler", values = air_sampler_palette) +
  
  # Add visible asterisk where insect was present (no legend here)
  geom_text(data = subset(data_air_vol, Insect == "Yes"),
            aes(label = Insect_Label),
            size = 8, hjust = -0.4, show.legend = FALSE) +
  
  # Colour and shape scales
  scale_color_manual(name = "Air Sampler", values = air_sampler_palette) +
  scale_shape_manual(name = "Location", values = c(
    "Church Farm" = 21, "NHM" = 22, "Yes" = 8),
    labels = c("Church Farm", "NHM")) +
  
  # Labels and theme
  labs(x = "Volume of air sampled (L)", 
       y = "Genera per 1000 reads",
       colour = "Air Sampler") +
  
  custom_theme
  # + scale_x_log10(labels = scales::comma)

ggsave("images/DNA_yield/genera_per_1k_air.pdf", genera_plot, width=10, height=6, dpi = 300 )


#Boxplots ----------------------

#Boxplot- not considering location 
yield_boxplot <- ggplot(insect_remove_normalised, aes(x=Air_Sampler, y = DNA_yield)) +
  geom_boxplot(aes(fill = Air_Sampler), colour = "black" , alpha = 0.7) + 
  scale_fill_manual(name = "Air Sampler", values = air_sampler_palette) +
  custom_theme +
  labs(y = "DNA yield (ng)", x = "Air Sampler") +
  scale_x_discrete(labels = c(
    "Coriolis Compact" = "Coriolis\nCompact",
    "Coriolis μ" = "Coriolis μ",
    "InnovaPrep Bobcat" = "InnovaPrep\nBobcat",
    "InnovaPrep Cub" = "InnovaPrep\nCub",
    "SASS 4100" = "SASS\n4100"
  ))

# DNA yield normalised by air volume
norm_yield_boxplot <- ggplot(insect_remove_normalised, aes(x=Air_Sampler, y = normalised_yield)) +
  geom_boxplot(aes(fill = Air_Sampler), colour = "black" , alpha = 0.7) + 
  scale_fill_manual(name = "Air Sampler", values = air_sampler_palette) +
  custom_theme +
  labs(y = " Normalised DNA yield (ng/L)", x = "Air Sampler") +
  scale_x_discrete(labels = c(
    "Coriolis Compact" = "Coriolis\nCompact",
    "Coriolis μ" = "Coriolis μ",
    "InnovaPrep Bobcat" = "InnovaPrep\nBobcat",
    "InnovaPrep Cub" = "InnovaPrep\nCub",
    "SASS 4100" = "SASS\n4100"
  ))

comb_boxplot <- yield_boxplot + norm_yield_boxplot +
  plot_layout(guides = "collect") &
  plot_annotation(tag_levels = "A") & 
  theme(legend.position = "right")

comb_boxplot

ggsave("images/DNA_yield/DNA_yield_boxplot.pdf", yield_boxplot,width=6, height=6, dpi = 300 )
ggsave("images/DNA_yield/Norm_DNA_yield_boxplot.pdf", norm_yield_boxplot,width=6, height=6, dpi = 300 )
ggsave("images/DNA_yield/Comb_DNA_yield_boxplot.pdf", comb_boxplot,width=14, height=6, dpi = 300 )



# Yield vs Time -----------------------------------------------------------
data <- data %>% 
  mutate(DNA_yield = as.numeric(DNA_yield))

#mean, min & max data for each Airsampler DNA yield 
data_mean_sd<- data %>% 
  group_by(Air_Sampler, Collection_time) %>% 
  mutate(mean = mean(DNA_yield), sd = sd(DNA_yield))

limits <- c(-100, 300)
breaks <- seq(-100, 300, by=100)

#Standard data
standard <- ggplot(data_mean_sd,
       (aes(x = Collection_time, y = mean, col = Air_Sampler))) +
  geom_line(position=position_dodge(width=2)) +
  geom_point(position=position_dodge(width=2)) + 
  geom_errorbar(aes(ymin = (mean - sd), ymax = (mean + sd)),
                width = 5, alpha = .2,
                position=position_dodge(width=2)) +
  theme_bw() +
  labs(y = " Mean DNA yield (ng)", , x = "") +
  scale_x_continuous( breaks = c(25,50)) +
  scale_y_continuous(limits=limits, breaks=breaks)
  
#Insect removed
insect_remove_mean <- insect_remove_normalised %>% 
  group_by(Air_Sampler, Collection_time) %>% 
  mutate(mean = mean(DNA_yield), sd = sd(DNA_yield))
  
no_insect <- ggplot(insect_remove_mean,
      (aes(x = Collection_time, y = mean, col = Air_Sampler))) +
  geom_line(position=position_dodge(width=2)) +
  geom_point(position=position_dodge(width=2)) + 
  geom_errorbar(aes(ymin = (mean - sd), ymax = (mean + sd)),
                width=5, alpha = 0.2,
                position=position_dodge(width=2)) +
  theme_bw() +
  labs(y = " Mean DNA yield (ng)", x = "Collection time (mins)") +
  scale_x_continuous( breaks = c(25,50)) +
  scale_y_continuous(limits=limits, breaks=breaks)
  
#Yield normalised by air volume & no insect
insect_remove_norm_mean <- insect_remove_normalised %>% 
  group_by(Air_Sampler, Collection_time) %>% 
  mutate(mean = mean(normalised_yield), sd = sd(normalised_yield))

norm <- ggplot(insect_remove_norm_mean,
                    (aes(x = Collection_time, y = mean, col = Air_Sampler))) +
  geom_line(position=position_dodge(width=2)) +
  geom_point(position=position_dodge(width=2)) + 
  geom_errorbar(aes(ymin = (mean - sd), ymax = (mean + sd)),
                width=5, alpha = 0.2,
                position=position_dodge(width=2)) +
  theme_bw() +
  labs(y = " Mean normalised DNA yield (ng)", x = "") +
  scale_x_continuous( breaks = c(25,50))

#Using patchwork to combine 
#Currently all have diff y axis scales - not nice 
Comb_DY_line <- standard + no_insect + norm + plot_layout(guides = "collect")
ggsave("images/DNA_Yield/Comb_DNA_yield_Line.pdf", Comb_DY_line,width=15, height=6, dpi = 300 )

