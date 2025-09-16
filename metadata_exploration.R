# R script to look at metadata

#Packages
library(dplyr)
library(ggplot2)
library(tidyr)
library(reshape2)
library(scales)
library(RColorBrewer)

# Use a font that supports Unicode
font_add_google("Roboto")
showtext_auto()

# Plot aesthetics ---------
sampler_colours <- setNames(brewer.pal(5, "Set2"), 
                            c("Bobcat", "Compact", 
                              "μ", "Cub", 
                              "Sass"))
# Theme for plots ------
custom_theme <- theme_minimal(base_size = 12) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.3),
    panel.grid.minor = element_blank(),
  )

location_labels <- c(
  "NHM" = "Natural History Museum",
  "Cfarm" = "Church Farm"
)


duration_labels <- c(
  "25" = "25 minutes",
  "50" = "50 minutes"
)

sampler_levels <- c(
  "Compact",
  "μ",
  "Bobcat",
  "Cub",
  "Sass"
)


# Load in data ---------
# Metadata
meta <- read.csv("metadata/sample_table.csv")

# Remove commas and convert to numeric
meta$NumReads <- as.numeric(gsub(",", "", meta$NumReads))
meta$Sampler <- gsub("Micro", "μ", meta$Sampler)

# Order 
meta$Sampler <- factor(meta$Sampler, levels = sampler_levels)


# Look at read count per sample -------
meta_summary_samp <- meta %>%
  group_by(Sampler) %>%
  summarise(
    mean_reads = mean(NumReads, na.rm = TRUE),
    se_reads   = sd(NumReads, na.rm = TRUE) / sqrt(n()),
    n_samples  = n(),
    .groups = "drop"
  )

meta_summary_len <- meta %>%
  group_by(Sample_length) %>%
  summarise(
    mean_reads = mean(NumReads, na.rm = TRUE),
    se_reads   = sd(NumReads, na.rm = TRUE) / sqrt(n()),
    n_samples  = n(),
    .groups = "drop"
  )

# Summarise mean and SE of reads by sampler and duration
meta_summary <- meta %>%
  group_by(Sampler, Sample_length) %>%
  summarise(
    mean_reads = mean(NumReads, na.rm = TRUE),
    se_reads   = sd(NumReads, na.rm = TRUE) / sqrt(n()),
    n_samples  = n(),
    .groups = "drop"
  )

# Plot
reads_bar <- ggplot(meta_summary,
                    aes(x = Sampler, y = mean_reads, fill = Sampler)) +
  geom_bar(stat = "identity", width = 0.7, position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = mean_reads - se_reads,
                    ymax = mean_reads + se_reads),
                width = 0.2, colour = "black",
                position = position_dodge(width = 0.8)) +
  facet_grid(cols = vars(Sample_length),
             labeller = labeller(Sample_length = duration_labels)) +
  scale_fill_manual(values = sampler_colours) +
  labs(
    x = "Sampler",
    y = "Mean reads per sample ± SE"
  ) +
  scale_x_discrete(labels = c(
    "Compact" = "Coriolis\nCompact",
    "μ"       = "Coriolis μ",
    "Bobcat"  = "InnovaPrep\nBobcat",
    "Cub"     = "InnovaPrep\nCub",
    "Sass"    = "SASS\n4100"
  )) +
  scale_y_continuous(labels = comma) + 
  custom_theme +
  theme(legend.position = "none") 

reads_bar 

ggsave("Images/meta/Reads_Per_Sampler.pdf",
       plot = reads_bar, width = 10, height = 7)


# DNA yield vs read count ----------------------
dna_reads_plot <- ggplot(meta, aes(x = DNA_yield, y = NumReads, colour = Sampler)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, colour = "black", linetype = "dashed") +
  scale_colour_manual(values = sampler_colours) +
  labs(
    x = "DNA yield (ng)",
    y = "Read count",
    colour = "Sampler"
  ) +
  custom_theme

dna_reads_plot

ggsave("Images/meta/DNA_vs_Reads.pdf",
       plot = dna_reads_plot, width = 10, height = 7)


# Air volume vs read count ----------------------
# Not a strong relationship
air_reads_plot <- ggplot(meta, aes(x = Air_volume, y = NumReads, colour = Sampler)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, colour = "black", linetype = "dashed") +
  scale_colour_manual(values = sampler_colours) +
  labs(
    x = "Air volume sampled (L)",
    y = "Read count",
    colour = "Sampler"
  ) +
  custom_theme

air_reads_plot

ggsave("Images/meta/AirSampled_vs_Reads.pdf",
       plot = air_reads_plot, width = 10, height = 7)
