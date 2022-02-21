
# Script 3A ----------------------------------------------------------------
# Analysis for Figure 4

# Cobamide Skin Microbiome Manuscript
# Spiec-Easi Taxa Filter
# Written by Mary Hannah Swaney
# University of Wisconsin - Madison

# Load required packages -----------------------------------------

library(dplyr)
library(readr)
library(readxl)

# set working directory to the code directory
setwd("/Users/mhswaney/GitHub/scripts/mhswaney/cobamide_paper/CobamidePaperCode/")

# Read in abundance data ----------------------------------------------------------

# LKMB002
LKMB002 <- read_csv("data/LKMB002_AbundanceTable.csv")
LKMB002_samples <- read_delim("data/SamplesForAnalysisLKMB002.txt", delim = "\t", col_names = TRUE) # list of samples for analysis (>1500000 reads)
LKMB002 <- LKMB002[,c("Species",LKMB002_samples$K_for_analysis.SampleID)] # samples above 1500000 reads

# Oh
Oh <- read_csv("data/Oh2016_AbundanceTable.csv")
Oh_samples <- read_delim("data/SamplesForAnalysisOh2016.txt",delim = "\t", col_names = TRUE) # list of samples for analysis (>1500000 reads)
Oh <- Oh[,c("Species",Oh_samples$colnames.Oh_rarefy_counts.)] # samples above 1500000 reads

# Hannigan
Hannigan <- read_csv("data/Hannigan2015_AbundanceTable.csv")
Hannigan_samples <- read_delim("data/SamplesForAnalysisHannigan2015.txt",delim = "\t", col_names = TRUE) # list of samples for analysis (>125000 reads)
Hannigan <- Hannigan[,c("Species",Hannigan_samples$keep.Run)] # samples above 125000 reads


# Get metadata for included samples ---------------------------------------

FormatMetadata <- function(){
  
  # read in metadata for LKMB02 metagenomes
  meta_LKMB002 <- read_excel("data/SupplementalMaterial.xlsx", sheet = "S1", skip = 1)
  
  # add full skin site names and microenvironment to metadata
  SkinSite <- c("Af","Al","Ax","Ba","Ch","Ea","Fh","Gb","Hp","Ic","Id","Mb","Na","Neg","Oc","Pa","Pc","Ph","Ra", "Sc","Tn","Tw","Um","Vf","Mock")
  full <- c("antecubital fossa","alar crease","axilla","back","cheek","external_auditory_canal","forehead",
            "glabella","hypothenar palm","inguinal_crease","interdigital_web_space","manubrium","nare",
            "negative_control","occiput","hypothenar palm","popliteal_fossa","plantar_heel",
            "retroaricular_crease","scalp","toenail","toe_web_space","umbilicus","volar_forearm","Mock")
  microenvironment <- c("moist","sebaceous","moist","sebaceous","sebaceous","sebaceous","sebaceous",
                        "sebaceous","dry","moist","moist","sebaceous","moist","negative",
                        "sebaceous","dry","moist","foot","sebaceous","sebaceous","foot","foot","moist","dry","Mock")
  sites_map <- data.frame(SkinSite,full,microenvironment)
  meta_LKMB002 <- merge(meta_LKMB002,sites_map,by="SkinSite")
  meta_LKMB002 <- meta_LKMB002[,c(1,2,19,20)]
  meta_LKMB002$Study <- "LKMB002"
  rownames(meta_LKMB002) <- meta_LKMB002$SampleID
  meta_LKMB002 <- meta_LKMB002[LKMB002_samples$K_for_analysis.SampleID,]
  
  meta_Oh <- read_excel("data/SupplementalMaterial.xlsx", sheet = "S2", skip = 1) %>% filter(Study == "Oh2016")
  meta_Oh <- meta_Oh[!duplicated(meta_Oh$SampleID), ]
  meta_Oh <- meta_Oh[,c("SampleID","SkinSite")]
  meta_Oh <- merge(meta_Oh,sites_map,by="SkinSite")
  meta_Oh$Study <- "Oh"
  rownames(meta_Oh) <- meta_Oh$SampleID
  meta_Oh <- meta_Oh[Oh_samples$colnames.Oh_rarefy_counts.,]
  
  meta_Hannigan <- read_excel("data/SupplementalMaterial.xlsx", sheet = "S2", skip = 1) %>% filter(Study == "Hannigan2015")
  meta_Hannigan <- meta_Hannigan[,c("RunID","SkinSite")]
  colnames(meta_Hannigan) <- c("SampleID","SkinSite")
  meta_Hannigan <- merge(meta_Hannigan,sites_map,by="SkinSite")
  meta_Hannigan$Study <- "Hannigan"
  rownames(meta_Hannigan) <- meta_Hannigan$SampleID
  meta_Hannigan <- meta_Hannigan[Hannigan_samples$keep.Run,]
  
  metadata <- rbind(meta_LKMB002, meta_Oh)
  metadata <- rbind(metadata, meta_Hannigan)
  return(metadata)}

metadata <- FormatMetadata()
metadata$SampleID

# Filter taxa -------------------------------------------------------------

# identify common species between all three datasets
common_species <- intersect(intersect(Oh$Species, Hannigan$Species),LKMB002$Species)

# join abundance data
combined <- full_join(LKMB002,Oh)
combined <- full_join(combined,Hannigan)
# remove NAs from full join
combined[is.na(combined)] <- 0
rownames(combined) <- combined$Species

# get the mean relative abundance for each species
abundance <- combined
abundance[,-1] = apply(abundance[,-1],2,function(x){x/sum(x)*100})
mean <- data.frame(Species = abundance$Species, Mean =rowMeans(abundance[,2:ncol(abundance)]))
#colSums(abundance[,2:ncol(abundance)])

# filter to only include species above 0.015% relative abundance
mean <- filter(mean,Mean > 0.015)
combined <- combined %>% filter(Species %in% mean$Species)

# include species found in all three datasets
combined <- combined[common_species,]
combined <- combined[!is.na(combined$Species),]

# identify taxa present in greater than 55% of all samples
combined_taxfilt <- combined[rowSums(combined[2:ncol(combined)] != 0) > ncol(combined[2:ncol(combined)])*0.55,]
taxa <- combined_taxfilt$Species 


# Get mean relative abundances for selected taxa -----------------------------------

# group samples by microenvironment
seb <- metadata %>% filter(microenvironment == "sebaceous")
moist <- metadata %>% filter(microenvironment == "moist")
dry <- metadata %>% filter(microenvironment == "dry")
foot <- metadata %>% filter(microenvironment == "foot")

# get abundances for selected taxa by microenvironment
seb_avg <- abundance[abundance$Species %in% taxa,c("Species",seb$SampleID)]
moist_avg <- abundance[abundance$Species %in% taxa,c("Species",moist$SampleID)]
dry_avg <- abundance[abundance$Species %in% taxa,c("Species",dry$SampleID)]
foot_avg <- abundance[abundance$Species %in% taxa,c("Species",foot$SampleID)]

# calculate mean relative abundance of each species within each microenvironment
seb_avg <- data.frame(Species = seb_avg$Species, Sebaceous_Mean =rowMeans(seb_avg[,2:ncol(seb_avg)]))
moist_avg <- data.frame(Species = moist_avg$Species, Moist_Mean =rowMeans(moist_avg[,2:ncol(moist_avg)]))
dry_avg <- data.frame(Species = dry_avg$Species, Dry_Mean =rowMeans(dry_avg[,2:ncol(dry_avg)]))
foot_avg <- data.frame(Species = foot_avg$Species, Foot_Mean =rowMeans(foot_avg[,2:ncol(foot_avg)]))

avgs <- left_join(seb_avg,moist_avg)
avgs <- left_join(avgs,dry_avg)
avgs <- left_join(avgs,foot_avg) # 

#write_csv(avgs, "SpiecEasi/data/taxa_mean_abundances_by_microenvironment_SpiecEasi.csv")
#write_csv(mean, "SpiecEasi/data/taxa_mean_abundances_SpiecEasi.csv")
