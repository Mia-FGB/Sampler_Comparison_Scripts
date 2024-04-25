library(dplyr)
library(ggplot2)
library(patchwork)

#read in data
data <- read.csv("metadata/sampler_comparison_DNA_yield.csv")

#adding a column for volume of air sampled
data_air_vol <- data %>% 
  mutate(air_volume = Collection_time * Collection_rate) %>% 
  #makind dna_yield numeric
  mutate(DNA_yield = as.numeric(DNA_yield), genera_per_1000_reads = as.numeric(genera_per_1000_reads)) %>% 
  #log the dna yield and air volume
  mutate(log_DNA_yield = log10(DNA_yield), log_air_volume = log10(air_volume),
         log_genera_per_1000_reads = log10(genera_per_1000_reads)) 


#plot a scatter graph of DNA yield against volume of air sampled - simple
ggplot(data_air_vol, aes(x = air_volume, y= DNA_yield, colour = Air_Sampler)) +
  geom_point() +theme_bw()

#Normal graph -----------------
#taking the normal data and then logging after so the axis remain the same
ggplot(data_air_vol, aes(x = air_volume, y= DNA_yield)) +
  geom_point(aes(shape = Air_Sampler, colour = Location), size = 5) +
  #so colours match my other graphs
  scale_color_manual(values = c("#FF7F50", "#008080")) + 
  #change the point shapes to match other graphs
  scale_shape_manual(values = c(19, 17, 15, 3, 4)) +
  theme_bw() +
  #Also find a way to highlight which samples contained insects
  geom_text(aes(label = Insect, hjust = -0.25)) +
  #axis labels
  labs(x = "Log Volume of air sampled (L)", y = "Log DNA yield (ng)", title = "DNA yield (ng)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_log10() + scale_x_log10(labels = scales::comma)

#Log data graph ---------------
ggplot(data_air_vol, aes(x = log_air_volume, y= log_DNA_yield)) +
  geom_point(aes(shape = Air_Sampler, colour = Location), size = 5) + theme_bw() +
  #Also find a way to highlight which samples contained insects
  geom_text(aes(label = Insect, hjust = -0.25)) +
  #axis labels
  labs(x = "Log Volume of air sampled (L)", y = "Log DNA yield (ng)", title = "DNA yield (ng)") +
  theme(plot.title = element_text(hjust = 0.5))

#--------------------------------------------------------------
#Plotting genera count per 1000 reads----------------------
#not sure if since it's for 1000 reads it's considered normalised?
#--------------------------------------------------------------
#plot a scatter graph of genera count against volume of air sampled - simple
ggplot(data_air_vol, aes(x = air_volume, y= genera_per_1000_reads, colour = Air_Sampler)) +
  geom_point(aes(shape = Air_Sampler, colour = Location), size = 5) + theme_bw() 

#need to log the axis for better seperation 
#Log genera count  graph ---------------
ggplot(data_air_vol, aes(x = log_air_volume, y= genera_per_1000_reads)) +
  geom_point(aes(shape = Air_Sampler, colour = Location), size = 5) + 
  #so colours match my other graphs
  scale_color_manual(values = c("#FF7F50", "#008080")) + 
  #change the point shapes to match other graphs
  scale_shape_manual(values = c(19, 17, 15, 3, 4)) +
  theme_bw() +
  #Also find a way to highlight which samples contained insects
  #geom_text(aes(label = Insect, hjust = -0.25)) +
  #axis labels
  labs(x = "Log Volume of air sampled (L)", y = "Genera per 1000 reads", title = "Diversity accumulation") +
  theme(plot.title = element_text(hjust = 0.5))

#--------------------------------------------------------------
#Trying to plot a graph of averaged data----------------------
#--------------------------------------------------------------

#Remove nd's and read the DNA_yield variable as numeric
clean_data <- data_air_vol %>% 
  #need to remove any nd's
  filter(DNA_yield != "nd") %>% 
  #makind dna_yield numeric
  mutate(DNA_yield = as.numeric(DNA_yield))

#check the strucutre of the dataset
str(clean_data)


#grouping the data - removing those with insect 
insect_remove_normalised <- subset(data_air_vol, !grepl("Insect", Insect)) %>% 
  filter(DNA_yield != "nd") %>% 
  mutate(normalised_yield = DNA_yield/air_volume)


#Boxplot
yield_boxplot <- ggplot(insect_remove_normalised, aes(x=Air_Sampler, y = DNA_yield)) +
  geom_boxplot() + theme_bw() +
  labs(y = "DNA yield (ng)", x = "Air Sampler") 

norm_yield_boxplot <- ggplot(insect_remove_normalised, aes(x=Air_Sampler, y = normalised_yield)) +
  geom_boxplot() + theme_bw() +
  labs(y = " Normalised DNA yield (ng)", x = "Air Sampler") 

comb_boxplot <- yield_boxplot + norm_yield_boxplot

ggsave("images/DNA_yield_boxplot.pdf", yield_boxplot,width=6, height=6, dpi = 300 )
ggsave("images/Norm_DNA_yield_boxplot.pdf", norm_yield_boxplot,width=6, height=6, dpi = 300 )
ggsave("images/Comb_DNA_yield_boxplot.pdf", comb_boxplot,width=12, height=6, dpi = 300 )



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

