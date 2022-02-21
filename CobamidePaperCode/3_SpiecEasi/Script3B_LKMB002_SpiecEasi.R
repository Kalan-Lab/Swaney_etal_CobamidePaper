
# Script 3B ----------------------------------------------------------------
# Analysis for Figure 4 

# Cobamide Skin Microbiome Manuscript
# Spiec-Easi - LKMB002 samples
# Written by Mary Hannah Swaney
# University of Wisconsin - Madison

# Load required packages and data -----------------------------------------

library(tidyverse)
library(SpiecEasi)

# set working directory to the code directory
setwd("/Users/mhswaney/GitHub/scripts/mhswaney/cobamide_paper/CobamidePaperCode/")

# load in the species to be used for analysis
source("3_SpiecEasi/Script3A_SpiecEasiTaxaFilter.R")

# extract sample IDs for each microenvironment
sebaceous <- filter(metadata, Study == "LKMB002" & microenvironment == "sebaceous")$SampleID
moist <- filter(metadata, Study == "LKMB002" & microenvironment == "moist")$SampleID
dry <- filter(metadata, Study == "LKMB002" & microenvironment == "dry")$SampleID
foot <- filter(metadata, Study == "LKMB002" & microenvironment == "foot")$SampleID

sites <- list(sebaceous,moist,dry,foot)

# Run SpiecEasi on each microenvironment 
LKMB002_SE <- lapply(sites,FUN=function(x) {
  
  # filter the abundance table to include the identified species for analysis 
  species_count <- LKMB002[which(LKMB002$Species %in% taxa),c("Species",x)]
  
  # separate the taxonomy
  species_count <- separate(species_count, "Species", into= c("D","P","C","O","F","G","S"),sep=";")
  
  # create a matrix from the count data
  species_matrix <- t(as.matrix(species_count[,8:ncol(species_count)]))
  
  # name the columns as the species names
  colnames(species_matrix) <- species_count$S
  
  ##### mb #######
  
  set.seed(100)
  
  # run Spiec Easi
  LKMB002.se.mb <- spiec.easi(species_matrix, method='mb', lambda.min.ratio=1e-2,
                            nlambda=20, pulsar.params=list(rep.num=50, seed = 100))
  return(LKMB002.se.mb)
})

names(LKMB002_SE) <- c("sebaceous","moist","dry","foot")

#save(LKMB002_SE, file="3_SpiecEasi/LKMB002.se.mb.RData")
