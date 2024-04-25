setwd("~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Sampler comparison")

#loading libraries
library(dplyr)
library(ggplot2)
#Reading in the data file
untouched_data <- read.csv("sampler_comparison_PHIbase_taxaID_readcount.csv", header=TRUE)