
# Script 4 ----------------------------------------------------------------
# Figure 5A, Supplemental Figures 7-9, Supplemental Tables 3-4

# Cobamide Skin Microbiome Manuscript
# Healthy skin microbiome abundance analysis
# Written by Mary Hannah Swaney
# University of Wisconsin - Madison

# Load required packages and data -----------------------------------------

library(readr)
library(tidyverse)
library(phyloseq)
library(readxl)
library(MMUPHin)
library(taxa)
library(vegan)
library(cowplot)
library(ggpubr)

# set working directory to the code directory
setwd("/Users/mhswaney/GitHub/scripts/mhswaney/cobamide_paper/CobamidePaperCode/")

# read in the final count table - LKMB002
LKMB002 <- read_csv("data/LKMB002_AbundanceTable.csv")

# read in the final count table - Oh
Oh <- read_csv("data/Oh2016_AbundanceTable.csv")
Oh_samples <- read_excel("data/SupplementalMaterial.xlsx", sheet="S2", skip=1) %>% filter(Study == "Oh2016")
Oh_unique_METIDs <- unique(Oh_samples$SampleID)
Oh <- Oh[,c("Species",Oh_unique_METIDs)] # only include 594 original samples in Oh 2016 paper

# Rarefy samples to 1,500,000 reads ---------------------------------------

## Create phyloseq object for LKMB002 metagenomes

# make a taxonomy table by splitting the relative abundance table Taxa column (Domain;phylum;class;etc)
K_taxonomy_table <- as.matrix(separate(LKMB002, Species , into = c("Domain","Phylum","Class","Order","Family","Genus","Species"), sep = ";", remove = FALSE,
                                       convert = FALSE, extra = "warn", fill = "warn"))
K_taxonomy_table <- K_taxonomy_table[,2:8] # taxonomy
rownames = LKMB002$Species
K_phyloseq_abundance <- data.matrix(LKMB002[,2:ncol(LKMB002)]) # abundances
rownames(K_phyloseq_abundance) = rownames
rownames(K_taxonomy_table) = rownames
K_ABUND = otu_table(K_phyloseq_abundance, taxa_are_rows = TRUE) # abundance table
K_TAX = tax_table(K_taxonomy_table) # taxonomy table
K_physeq = phyloseq(K_ABUND, K_TAX) # phyloseq object

# rarefy LKMB002 samples - 1,500,000 reads, without replacement
#LKMB002_rarefy <- rarefy_even_depth(physeq = K_physeq,sample.size = 1500000,rngseed = 100, replace = FALSE,verbose = TRUE)
#save(LKMB002_rarefy, file = "Corynebacterium_Microbiome_Abundance_Analysis/data/LKMB002_rarefy.RData")

load(file = "4_Corynebacterium_Microbiome_Abundance_Analysis/data/LKMB002_rarefy.RData")

## Create phyloseq object for Oh metagenomes

# make a taxonomy table by splitting the relative abundance table Taxa column (Domain;phylum;class;etc)
O_taxonomy_table <- as.matrix(separate(Oh, Species , into = c("Domain","Phylum","Class","Order","Family","Genus","Species"), sep = ";", remove = FALSE,
                                       convert = FALSE, extra = "warn", fill = "warn"))
O_taxonomy_table <- O_taxonomy_table[,2:8] # taxonomy
rownames = Oh$Species
O_phyloseq_abundance <- data.matrix(Oh[,2:ncol(Oh)])
rownames(O_phyloseq_abundance) = rownames
rownames(O_taxonomy_table) = rownames
O_ABUND = otu_table(O_phyloseq_abundance, taxa_are_rows = TRUE)
O_TAX = tax_table(O_taxonomy_table) # taxonomy table
O_physeq = phyloseq(O_ABUND, O_TAX) # phyloseq object

# rarefy Oh 2016 samples - 1,500,000 reads, without replacement
#Oh_rarefy <- rarefy_even_depth(physeq = O_physeq,sample.size = 1500000,rngseed = 100, replace = FALSE,verbose = TRUE)
#save(Oh_rarefy, file = "Corynebacterium_Microbiome_Abundance_Analysis/data/Oh_rarefy.RData")

load(file = "4_Corynebacterium_Microbiome_Abundance_Analysis/data/Oh_rarefy.RData")

# Create relative abundance tables ----------------------------------------

# extract LKMB002 counts table from phyloseq object
LKMB002_rarefy_counts <- as.data.frame(otu_table(LKMB002_rarefy))
LKMB002_rarefy_counts$Species <- rownames(LKMB002_rarefy_counts)
LKMB002_rarefy_counts <- LKMB002_rarefy_counts[,c(220,1:219)]

# extract Oh2016 counts table from phyloseq object
Oh_rarefy_counts <- as.data.frame(otu_table(Oh_rarefy))
#write_delim(data.frame(colnames(Oh_rarefy_counts)), "data/SamplesForAnalysisOh2016.txt") # Oh samples included in this analysis and SPIEC-EASI analysis
Oh_rarefy_counts$Species <- rownames(Oh_rarefy_counts)
Oh_rarefy_counts <- Oh_rarefy_counts[,c(493,1:492)] # re-order columns

combined_counts <- full_join(Oh_rarefy_counts,LKMB002_rarefy_counts) # join count tables for the two data sets
combined_counts[is.na(combined_counts)] <- 0

combined_counts_matrix <- data.matrix(combined_counts[,2:ncol(combined_counts)]) # convert to matrix
rownames(combined_counts_matrix) <- combined_counts$Species

# Format metadata ---------------------------------------------------------

# read in metadata for LKMB002 metagenomes
meta_LKMB002 <- read_excel("data/SupplementalMaterial.xlsx", sheet="S1", skip=1)

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

# get all non-negative and non-mock samples for analysis that are >1500000 read count
K_for_analysis <- meta_LKMB002 %>% filter(SampleID %in% colnames(LKMB002_rarefy_counts))
K_for_analysis <- K_for_analysis %>% filter(microenvironment != "Mock" & microenvironment != "negative")
#write_delim(data.frame(K_for_analysis$SampleID), "data/SamplesForAnalysisLKMB002.txt") # LKMB002 samples included in this analysis and SPIEC-EASI analysis

# read in metadata for Oh2016 metagenomes

meta_Oh <- read_excel("data/SupplementalMaterial.xlsx", sheet="S2", skip=1) %>% filter(Study == "Oh2016")
meta_Oh <- meta_Oh %>% filter(Study == "Oh2016")
meta_Oh <- meta_Oh[!duplicated(meta_Oh$SampleID),]
meta_Oh <- left_join(meta_Oh, sites_map)
meta_Oh <- meta_Oh[,c("SkinSite","SampleID","full","microenvironment","Study")]

meta_all <- rbind(meta_LKMB002, meta_Oh) # join metadata from the two studies

# Study effect adjustment using MMUPHin ------------------------------------------

# extract unique samples in the rarefied counts table
samples <- colnames(combined_counts_matrix)

# attach metadata (skin site characteristic and skin site) for each sample
metadata <- meta_all %>% filter(SampleID %in% samples & microenvironment != "Mock" & microenvironment != "negative") # samples for analysis
colnames(combined_counts_matrix) <- as.character(colnames(combined_counts_matrix))
rownames(metadata)<- as.character(unique(metadata$SampleID))
order <- rownames(metadata)
combined_counts_matrix <- combined_counts_matrix[,order] # order matrix columns by metadata order

# correct for study effect using MMUPHin
fit_adjust_batch2 <- adjust_batch(feature_abd = combined_counts_matrix,
                                  batch = "Study",
                                  data = metadata,
                                  control = list(verbose = TRUE))

batch_adjusted <- data.frame(fit_adjust_batch2$feature_abd_adj) # extract adjusted counts table
batch_adjusted$Species <- rownames(batch_adjusted)
colnames(batch_adjusted) <- gsub("\\.", "-", colnames(batch_adjusted))
batch_adjusted <- batch_adjusted[,c(705,1:704)] # re-order columns

# convert to relative abundances
batch_adjusted_relabund <- batch_adjusted
batch_adjusted_relabund[,-1] = apply(batch_adjusted_relabund[,-1],2,function(x){x/sum(x)*1})
#colSums(batch_adjusted_relabund[,2:ncol(batch_adjusted_relabund)]) # should add up to 1

# gather data into long format
gathered <- batch_adjusted_relabund %>% gather(SampleID, Abundance, 2:705)

# Calculate alpha diversity metrics ---------------------------------------

## Create phyloseq object (counts) from batch-adjusted table
taxonomy_table <- as.matrix(separate(batch_adjusted, Species , into = c("Domain","Phylum","Class","Order","Family","Genus","Species"), sep = ";", remove = FALSE,
                                     convert = FALSE, extra = "warn", fill = "warn"))
taxonomy_table <- taxonomy_table[,2:8]
rownames = batch_adjusted$Species
phyloseq_abundance <- data.matrix(batch_adjusted[,2:ncol(batch_adjusted)]) # using counts
rownames(phyloseq_abundance) = rownames
rownames(taxonomy_table) = rownames
ABUND = otu_table(phyloseq_abundance, taxa_are_rows = TRUE)
TAX = tax_table(taxonomy_table)
physeq = phyloseq(ABUND, TAX)

diversity <- estimate_richness(physeq = physeq) # calculate alpha diversity metrics
diversity$SampleID <- rownames(diversity)
diversity$SampleID <- gsub("\\.", "-", diversity$SampleID)

# Function to plot Bray Curtis dissimilarity - by microenvironment -----------------

plot_NMDS <- function(site, meta, abundance_long, abundance_table, diversity, taxa){
  
  # filter metadata for samples from specified site
  meta <- meta %>% filter(microenvironment == site)
  # filter abundance data for samples from specified site
  SampleIDs <- meta$SampleID
  abundance_long <- abundance_long %>% filter(SampleID %in% SampleIDs)
  
  if (taxa == "Corynebacterium"){
    # Identify which sample IDs have Corynebacterium B12 producers and add up their relative abundances in each sample
    B12 <- abundance_long %>% filter(Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium amycolatum" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium kroppenstedtii" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium ulcerans" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium pseudotuberculosis" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium diphtheriae" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium argentoratense" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium choanis" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium epidermidicanis" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium geronticis" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium glucuronolyticum" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium kutscheri" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium lactis" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium matruchotii" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium mustelae" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium pelargi" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium pseudopelargi" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium rouxii" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. 1959" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. ATCC 6931" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sphenisci" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium vitaeruminis" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. sy039" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium aquilae" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. NML93-0612" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium xerosis"
    )} 
  
  else if (taxa == "Cutibacterium"){
    B12 <- abundance_long %>% filter(Species == "Bacteria;Actinobacteria;Actinomycetia;Propionibacteriales;Propionibacteriaceae;Cutibacterium;Cutibacterium acnes" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Propionibacteriales;Propionibacteriaceae;Cutibacterium;Cutibacterium avidum" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Propionibacteriales;Propionibacteriaceae;Cutibacterium;Cutibacterium granulosum")}
  
  # add up the abundance of each B12 producer in each sample
  B12 <- aggregate(B12$Abundance, by=list(SampleID=B12$SampleID), FUN=sum)
  colnames(B12) <- c("SampleID","Abundance")
  
  # merge metadata with B12 producer abundance data
  B12 <- merge(meta, B12[c("SampleID","Abundance")], by="SampleID")
  B12$Abundance <- B12$Abundance * 100 # convert to relative abundance out of 100%
  B12$log10abundance <- log10(B12$Abundance) # take log of coryne producer abundance
  
  abundance <- abundance_table[,c("Species",SampleIDs)] # relative abundance table
  
  print(site)
  
  B12 <- left_join(B12, diversity) # add in alpha diversity
  
  abundances_sqrt <- sqrt(abundance[,2:ncol(abundance)]) # square root transformation
  abundances_sqrt$Species <- abundance$Species
  abundances_sqrt[is.na(abundances_sqrt)] <- 0
  
  # prepare data for beta diversity calculations
  obj <- parse_tax_data(abundances_sqrt,
                        class_cols = "Species",
                        class_sep = ";",
                        class_regex = "^([a-z]{0,1})_{0,2}(.*)$",
                        class_key = c("tax_rank" = "taxon_rank", "name" = "taxon_name"))
  obj$data$class_data <- NULL
  names(obj$data) <- "OTU_abundance"
  
  # bray curtis dissimilarity calculation
  beta_dist <- vegdist(t(obj$data$OTU_abundance[,SampleIDs]),
                       method = "bray")
  set.seed(100)
  # NMDS using Vegan package
  mds <- metaMDS(beta_dist, k=2, trymax=100)
  mds_data <- as.data.frame(mds$points)
  mds_data$SampleID <- rownames(mds_data)
  mds_data <- dplyr::left_join(B12,mds_data)
  
  # plot NMDS
  plot <- ggplot(mds_data, aes(x=MDS1, y=MDS2, color=log10abundance,size=Shannon)) + geom_point() +
    labs(color = paste0("Log10 ", taxa,"\ncobamide producer abundance")) +
    theme_classic() +
    ggtitle(site) +
    scale_color_gradient(low="blue",high="orange") +
    guides(color = "colorbar", size= "legend") +
    scale_size(range=c(1,7))
  
  print(cor.test(mds_data$log10abundance,mds_data$Shannon,method="spearman", exact=FALSE))
  
  plot
  
}

# Generate NMDS plots for each microenvironment

# gathered = batch-adjusted relative abundances in long format
# batch_adjusted_relabund = batch-adjusted relative abundance table
# diversity = diversity metrics computed on batch-adjusted counts table

# Corynebacterium abundances
# Data for Supplemental Table 3
sebaceous_plot_Co <- plot_NMDS("sebaceous",metadata, gathered, batch_adjusted_relabund, diversity, taxa = "Corynebacterium")
moist_plot_Co <- plot_NMDS("moist",metadata, gathered, batch_adjusted_relabund, diversity, taxa = "Corynebacterium")
dry_plot_Co <- plot_NMDS("dry",metadata, gathered, batch_adjusted_relabund, diversity, taxa = "Corynebacterium")
foot_plot_Co <- plot_NMDS("foot",metadata, gathered, batch_adjusted_relabund, diversity, taxa = "Corynebacterium")

# final plot - Figure 5A
plot_grid(sebaceous_plot_Co, moist_plot_Co, dry_plot_Co, foot_plot_Co, nrow=2)

# Cutibacterium abundances
sebaceous_plot_Cu <- plot_NMDS("sebaceous",metadata, gathered, batch_adjusted_relabund, diversity, taxa = "Cutibacterium")
moist_plot_Cu <- plot_NMDS("moist",metadata, gathered, batch_adjusted_relabund, diversity, taxa = "Cutibacterium")
dry_plot_Cu <- plot_NMDS("dry",metadata, gathered, batch_adjusted_relabund, diversity, taxa = "Cutibacterium")
foot_plot_Cu <- plot_NMDS("foot",metadata, gathered, batch_adjusted_relabund, diversity, taxa = "Cutibacterium")

# final plot - Cutibacterium - Supplemental Figure 7
plot_grid(sebaceous_plot_Cu, moist_plot_Cu, dry_plot_Cu, foot_plot_Cu, nrow=2)

# High vs low CPC abundance analysis --------------------------------------

## Identify which SampleIDs have B12 producers and add up their relative abundances in each sample
B12 <- gathered %>% filter(Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium amycolatum" |
                             Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium kroppenstedtii" |
                             Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium ulcerans" |
                             Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium pseudotuberculosis" |
                             Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium diphtheriae" |
                             Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium argentoratense" |
                             Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium choanis" |
                             Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium epidermidicanis" |
                             Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium geronticis" |
                             Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium glucuronolyticum" |
                             Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium kutscheri" |
                             Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium lactis" |
                             Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium matruchotii" |
                             Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium mustelae" |
                             Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium pelargi" |
                             Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium pseudopelargi" |
                             Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium rouxii" |
                             Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. 1959" |
                             Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. ATCC 6931" |
                             Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sphenisci" |
                             Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium vitaeruminis" |
                             Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. sy039" |
                             Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium aquilae" |
                             Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. NML93-0612" |
                             Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium xerosis"
)

# add up the abundance of each B12 producer in each sample
B12 <- aggregate(B12$Abundance, by=list(SampleID=B12$SampleID), FUN=sum)
colnames(B12) <- c("SampleID","Abundance")

# merge metadata with B12 producer abundance data and diversity data
B12 <- merge(metadata, B12[c("SampleID","Abundance")], by="SampleID")
B12 <- left_join(B12, diversity)
B12$Abundance <- B12$Abundance * 100

B12$log10abundance <- log10(B12$Abundance) # take log of coryne producer abundance
B12 <- B12 %>% filter(!is.infinite(log10abundance)) # exclude a few samples with no coryne B12 producers 

summary(B12$Abundance)

# first and third quartiles for CPC abundance
low <- B12 %>% filter(Abundance < 0.05)
low_samples <- low$SampleID
high <- B12 %>% filter(Abundance > 0.75)
high_samples <- high$SampleID

# extract samples with "low" and "high" CPC abundance
low_data <- gathered %>% filter(SampleID %in% low_samples)
low_data$Group <- "Low"
high_data <- gathered %>% filter(SampleID %in% high_samples)
high_data$Group <- "High"

data <- rbind(low_data,high_data)
data <- data %>% separate(Species, into=c("Domain","Phylum","Class","Order","Family","Genus","Species"), sep=";")

# assign taxa under 10% abundance as "Other"
data$name_other <- ifelse(data$Abundance <= 0.1, "Other", as.character(data$Species))

# add all "Other" abundances together
dat_plot <- data %>% group_by(SampleID, name_other, Group) %>% summarise(Abundance = sum(Abundance))

mhs45rainbow = c("#771155", "#AA4488", "#CC99BB", "#e6c5e0","#fae6f8",
                 "#114477", "#4477AA", "#77AADD", "#b9cfeb","#e6ebfa",
                 "#117777", "#44AAAA", "#77CCCC", "#bae6e0","#e1f4f5",
                 "#117744", "#44AA77", "#88CCAA", "#bce3d0","#dff2e5",
                 "#777711", "#AAAA44", "#DDDD77", "#e7eba4","#f3f5cb",
                 "#774411", "#AA7744", "#DDAA77", "#edc99f","#f7e1c8",
                 "#771122", "#AA4455", "#DD7788", "#e8aeae","#f7dada",
                 "#4c285c","#9167a3","#d0a3e3","#dac2ed","#ebdff5",
                 "#636363","#8f8f8f","#c9c9c9","#e3e3e3","#f0f0f0")

dat_plot$Group <- factor(dat_plot$Group, levels=c("Low","High"))

# plot - Supplemental Figure 9
ggplot(dat_plot, aes(x=SampleID, y=Abundance*100,fill=name_other)) + geom_col() +
  theme_classic() + 
  facet_wrap(~Group, scales = "free_x") + 
  theme(axis.text.x = element_blank()) + ylab("Relative abundance") + 
  scale_fill_manual(values = mhs45rainbow) + xlab("Sample") + 
  labs(fill="Species")

# Function to correlate non-B12 producing Corynebacteria vs Shannon Diversity -----------------

corr <- function(site, meta, abundance_long, diversity, taxa){
  
  # filter metadata for samples from specified site
  meta <- meta %>% filter(microenvironment == site)
  # filter abundance data for samples from specified site
  SampleIDs <- meta$SampleID
  abundance_long <- abundance_long %>% filter(SampleID %in% SampleIDs)
  
  if (taxa == "Corynebacterium_non"){
    # Identify which sample IDs have Corynebacterium non-B12 producers and add up their relative abundances in each sample
    B12 <- abundance_long %>% filter(Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium ammoniagenes" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium anserum" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium atypicum" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium aurimucosum" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium bovis" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium callunae" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium camporealensis" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium casei" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium crudilactis" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium cystitidis" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium doosanense" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium efficiens" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium endometrii" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium flavescens" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium frankenforstense" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium genitalium" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium deserti" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium glaucum" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium glutamicum" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium glyciniphilium" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium halotolerans" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium humireducens" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium imitans" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium jeikeium" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium kefirresidentii" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium macginleyii" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium marinum" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium maris" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium minutissimum" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium mycetoides" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium nuruki" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium phocae" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium propinquum" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium provencense" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium renale" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium resistens" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium riegelii" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sanguinis" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium segmentosum" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium simulans" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium singulare" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. 2019" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. 2039" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. 2183" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. 2184" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. 4H37-19" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. FDAARGOS 1242" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. LMM-1652" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. Marseille-Q3630" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. MC1420" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. NML98-0116" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. zg-320" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. zg-917" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. ZJ-599" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. stationis" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. striatum" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. suranareeae" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. terpenotabidum" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. testudinoris" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. timonense" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. tuberculostearicum" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. urealyticum" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. ureicelerivorans" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. uterequi" |
                                       Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. variabile" 
    )} 
  
  
  # add up the abundance of each B12 producer in each sample
  B12 <- aggregate(B12$Abundance, by=list(SampleID=B12$SampleID), FUN=sum)
  colnames(B12) <- c("SampleID","Abundance")
  
  # merge metadata with B12 producer abundance data
  B12 <- merge(meta, B12[c("SampleID","Abundance")], by="SampleID")
  B12$Abundance <- B12$Abundance * 100 # convert to relative abundance out of 100%
  B12$log10abundance <- log10(B12$Abundance) # take log of coryne producer abundance
  
  print(site)
  
  B12 <- left_join(B12, diversity) # add in alpha diversity
  
  print(cor.test(B12$log10abundance,B12$Shannon,method="spearman", exact=FALSE))
  
}

# Calculate correlation between non-B12 producing Corynebacteria and Shannon Diversity
# Data for Supplemental Table 4
sebaceous_corr <- corr("sebaceous",metadata, gathered, diversity, taxa = "Corynebacterium_non")
moist_corr <- corr("moist",metadata, gathered, diversity, taxa = "Corynebacterium_non")
dry_corr <- corr("dry",metadata, gathered, diversity, taxa = "Corynebacterium_non")
foot_corr <- corr("foot",metadata, gathered, diversity, taxa = "Corynebacterium_non")


# Producer vs Non-Producer Correlation Analysis  --------------------------------------

nonB12 <- gathered %>% filter(Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium ammoniagenes" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium anserum" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium atypicum" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium aurimucosum" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium bovis" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium callunae" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium camporealensis" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium casei" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium crudilactis" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium cystitidis" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium deserti" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium doosanense" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium efficiens" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium endometrii" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium flavescens" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium frankenforstense" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium genitalium" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium glaucum" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium glutamicum" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium glyciniphilium" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium halotolerans" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium humireducens" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium imitans" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium jeikeium" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium kefirresidentii" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium macginleyii" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium marinum" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium maris" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium minutissimum" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium mycetoides" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium nuruki" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium phocae" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium propinquum" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium provencense" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium renale" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium resistens" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium rieglii" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sanguinis" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium segmentosum" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium simulans" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium singulare" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. 2019" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. 2039" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. 2183" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. 2184" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. 4H37-19" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. FDAARGOS 1242" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. LMM-1652" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. Marseille-Q3630" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. MC1420" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. NML98-0116" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. zg-320" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. zg-917" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. ZJ-599" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. stationis" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. striatum" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. suranareeae" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. terpenotabidum" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. testudinoris" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. timonense" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. tuberculostearicum" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. urealyticum" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. ureicelerivorans" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. uterequi" |
                                Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. variabile" )

# add up the abundance of each nonB12 producer in each sample
nonB12 <- aggregate(nonB12$Abundance, by=list(SampleID=nonB12$SampleID), FUN=sum)
colnames(nonB12) <- c("SampleID","Non-Abundance")

# merge metadata with nonB12 producer abundance data and diversity data
nonB12 <- merge(metadata, nonB12[c("SampleID","Non-Abundance")], by="SampleID")
nonB12 <- left_join(nonB12, diversity)
nonB12$`Non-Abundance` <- nonB12$`Non-Abundance` * 100

nonB12$log10nonabundance <- log10(nonB12$`Non-Abundance`) # take log of coryne non-producer abundance
nonB12 <- nonB12 %>% filter(!is.infinite(log10nonabundance)) # exclude a few samples with no non-coryne B12 producers 

all <- left_join(B12,nonB12[,c(1,16)])

# Supplemental Figure 8
ggplot(all,aes(x=log10abundance, y=log10nonabundance)) + geom_point() + facet_wrap(~microenvironment,scales="free") + theme_classic() +
  stat_cor(method = "spearman", label.x = -0.5, label.y = -1.5) + ylab("Log10 Corynebacterium non-producer abundance") + 
  xlab("Log10 Corynebacterium producer abundance")

