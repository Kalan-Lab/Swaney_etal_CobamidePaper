
# Script 5B ----------------------------------------------------------------
# Figure 6A

# Cobamide Skin Microbiome Manuscript
# Corynebacterium phylogeny 
# Written by Mary Hannah Swaney
# University of Wisconsin - Madison

# Load required packages and data -----------------------------------------

library(tidyverse)
library(ggtree)
library(readxl)
library(ggtree)
library(ape)
library(ggnewscale)

# set working directory to the code directory
setwd("/Users/mhswaney/GitHub/scripts/mhswaney/cobamide_paper/CobamidePaperCode/")

KO_db <- read_excel("data/SupplementalMaterial.xlsx",sheet = "S6",skip = 1,)[,1:3]

tree <- read.tree("5_Corynebacterium_comparative_genomics/data/Corynebacterium_Phylogenetic_Tree.txt")

pangenome <- read.delim("5_Corynebacterium_comparative_genomics/data/Corynebacterium_Pangenome_Summary.txt")

KOFamScan_results <- read_csv("5_Corynebacterium_comparative_genomics/data/Corynebacterium_B12_KOFamScan_Results.csv")

for_heatmap <- KOFamScan_results[,c("species","Gene ID")] # subset only the species name and gene name (cobA, etc.)
for_heatmap$present <- 1 # add a value of 1 for each item in the data frame to represent gene presence
for_heatmap$`Gene ID` <- factor(for_heatmap$`Gene ID`, levels= KO_db$`Gene ID`)

# spread out the values into a matrix-like format - species vs genes. Fill any blanks with 0 (this will represent gene absence)
matrix <- spread(for_heatmap,key="Gene ID",value="present",fill=0)

# rename some of the species to match up with phylogenetic tree naming - fixes some underscores/dashes
matrix[32,1] <- "Corynebacterium_lactis_RW2_5"
matrix[67,1] <- "Corynebacterium_ureicelerivoransstrain_IMMIB_RIV_2301"
matrix[13,1] <- "Corynebacterium_deserti_GIMN1010"
matrix[16,1] <- "Corynebacterium_efficiens_YS_314"
matrix[3,1] <- "Corynebacterium_aquilae_DSM_44791strain_S_613"
matrix[63,1] <- "Corynebacterium_terpenotabidum_Y_11"

# combine similar genes and gene fusions for pathway simplification
matrix$hemD <- ifelse(matrix$`cobA-hemD` == 1, 1, matrix$hemD)
matrix$cobA <- ifelse(matrix$`cobA-hemD` == 1, 1, matrix$cobA)
matrix$`cobA` <- matrix$cobA + matrix$cysG
matrix$`cobIJ` <- matrix$cobIJ + matrix$`cobI/cbiL`

# exclude the extra columns that were combined
matrix <- matrix[,c(1:6,8,10,12:28)]

# because I only want to show presence/absence in my heatmap, I'm changing any values above 1 back to 1 
matrix$hemD <- ifelse(matrix$hemD > 1, 1, matrix$hemD)
matrix$cobA <- ifelse(matrix$cobA > 1, 1, matrix$cobA)
matrix$cobIJ <- ifelse(matrix$cobIJ > 1, 1, matrix$cobIJ)

# Phylogenetic tree plotting ----------------------------------------------

# root phylogenetic tree with T. paurometabola
rooted <- root(tree, which(tree$tip.label == "Tsukamurella_paurometabola_DSM_20162"))

# create tree object - color text by ecosystem association
# if there is the error: "Error in DataMask$new(.data, caller_env) : argument "caller_env" is missing, with no default"
# need to install dplyr v1.0.5

p <- ggtree(rooted) %<+% pangenome + geom_tiplab(size=1.5, align=TRUE) + 
  scale_color_manual(values=c("#d1453b","#008f94")) + xlim(0,2)

# flip nodes
p <- flip(p, 67, 131)
p <- flip(p, 88, 12)
p <- flip(p, 89,13)
p <- flip(p, 90,14)

# read in ecosystem association data for each gennome
pangenome <- pangenome %>% rename("species" = "Strain")

# join together ordered tree species names with ecosystem association
order <- rooted$tip.label           
order <- data.frame(order)
colnames(order) = c("species")
ecosystem <- dplyr:: left_join(order, pangenome)
ecosystem2 <- as.matrix(ecosystem[,"Ecosystem"])
rownames(ecosystem2) <- ecosystem$species
colnames(ecosystem2) <- c("Ecosystem")
ecosystem2 <- as.data.frame(ecosystem2)

# plot tree with ecosystem association
p1 <- gheatmap(p, ecosystem2[, "Ecosystem", drop = FALSE], width = 0.05, offset=0.5, color = "black",
               colnames_angle = 90,
               colnames_offset_y = 76,
               font.size = 2) +
  scale_fill_manual(breaks = c("Environment","Host"), values = c("#F8766D", "#00BFC4")) +
  labs(fill="Ecosystem")

# join together ordered tree species names with b12 pathway gene presence/absence
pathway <- dplyr::left_join(order, matrix)
order <- rooted$tip.label           
order <- data.frame(order)
colnames(order) = c("species")
pathway <- dplyr:: left_join(order, matrix)
pathway2 <- as.matrix(pathway[,2:25])
rownames(pathway2) <- pathway$species
colnames(pathway2) <- colnames(pathway[2:25])
pathway2 <- as.data.frame(pathway2)

p2 <- p1 + new_scale_fill()

# plot tree with B12 pathway heatmap
p3 <- gheatmap(p2, pathway2, width=1,offset=0.62, color="white", 
               colnames_position="top", low="#DDC7D6", high="#A11856",
               colnames_angle = 90,
               colnames_offset_y = 2,
               font.size = 2) +
  theme(legend.position = "none")

# join together ordered tree species names with genome length
genomeLength <- as.matrix(pangenome[,"total_length"])
rownames(genomeLength) <- pangenome$species
colnames(genomeLength) <- c("Genome Length")
genomeLength <- as.data.frame(genomeLength)

p4 <- p3 + new_scale_fill()

# plot final tree - Figure 6A
final <- gheatmap(p4, genomeLength, 
                  width = 0.05, offset=0.55, color = "black", 
                  high= "#203f87",low="white",
                  colnames_angle = 90,
                  colnames_offset_y = 77,
                  font.size = 2) +
  theme(legend.position = "left") +
  labs(fill="Genome length")

final




