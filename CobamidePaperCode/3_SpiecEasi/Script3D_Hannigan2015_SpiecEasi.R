
# Script 3D ----------------------------------------------------------------
# Analysis for Figure 4 

# Cobamide Skin Microbiome Manuscript
# Spiec-Easi - Hannigan 2015 samples
# Written by Mary Hannah Swaney
# University of Wisconsin - Madison

# Load required packages and data -----------------------------------------

library(tidyr)
library(SpiecEasi)

# set working directory to the code directory
setwd("/Users/mhswaney/GitHub/scripts/mhswaney/cobamide_paper/CobamidePaperCode/")

# load in the species to be used for analysis
source("3_SpiecEasi/Script3A_SpiecEasiTaxaFilter.R")

# extract sample IDs for each microenvironment
sebaceous <- filter(metadata, Study == "Hannigan" & microenvironment == "sebaceous")$SampleID
moist <- filter(metadata, Study == "Hannigan" & microenvironment == "moist")$SampleID
dry <- filter(metadata, Study == "Hannigan" & microenvironment == "dry")$SampleID
foot <- filter(metadata, Study == "Hannigan" & microenvironment == "foot")$SampleID

sites <- list(sebaceous,moist,dry,foot)

# Run SpiecEasi on each microenvironment 
Hannigan_SE <- lapply(sites,FUN=function(x) {
  
  # filter the abundance table to include the identified species for analysis 
  species_count <- Hannigan[which(Hannigan$Species %in% taxa),c("Species",x)]
  
  # separate the taxonomy
  species_count <- separate(species_count, "Species", into= c("D","P","C","O","F","G","S"),sep=";")
  
  # create a matrix from the count data
  species_matrix <- t(as.matrix(species_count[,8:ncol(species_count)]))
  
  # name the columns as the species names
  colnames(species_matrix) <- species_count$S
  
  ##### mb #######
  
  set.seed(100)
  
  # run Spiec Easi
  Hannigan.se.mb <- spiec.easi(species_matrix, method='mb', lambda.min.ratio=1e-2,
                               nlambda=20, pulsar.params=list(rep.num=50, seed = 100))
  return(Hannigan.se.mb)
})

names(Hannigan_SE) <- c("sebaceous","moist","dry","foot")

#save(Hannigan_SE, file="3_SpiecEasi/Hannigan.se.mb.RData")

