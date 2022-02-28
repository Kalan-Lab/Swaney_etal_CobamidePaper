
# Script 4B  ----------------------------------------------------------------
# Figure 6B-C, Supplemental Figure 5

# Cobamide Skin Microbiome Manuscript
# Atopic dermatitis skin microbiome abundance analysis
# Written by Mary Hannah Swaney
# University of Wisconsin - Madison

# Load required packages and data -----------------------------------------

library(readr)
library(tidyverse)
library(phyloseq)
library(readxl)
library(taxa)
library(vegan)
library(cowplot)
library(FSA)

setwd("/Users/mhswaney/GitHub/scripts/mhswaney/cobamide_paper/CobamidePaperCode/")

# read metadata for samples
meta <- read_excel("data/SupplementalMaterial.xlsx", skip = 1, sheet = "S7")
meta <- meta %>% separate(SkinSite, into=c("SkinSite","Side"), sep="-")
meta <- meta %>%
  separate(Subject_ID, into = c("Group", "Subject_Num"), sep = 2,remove = FALSE)

# read in final count table - Byrd 2017
Byrd <- read_csv("4_Corynebacterium_Microbiome_Abundance_Analysis/data/Byrd2017_AbundanceTable.csv")

# create relative abundance table
relativeAbundances <- Byrd
relativeAbundances[,-1] = apply(relativeAbundances[,-1],2,function(x){x/sum(x)*100})
#colSums(relativeAbundances[,2:ncol(relativeAbundances)]) # should add up to 1

# gather data into long format
gather <- relativeAbundances %>% gather(SampleID, Abundance, 2:ncol(relativeAbundances))

# Format metadata ---------------------------------------------------------

# add full skin site names and microenvironment to metadata
SkinSite <- c("Ac","Al","Ax","Ba","Ch","Ea","Fh","Gb","Hp","Ic","Id","Mb","Na","Neg","Oc","Pa","Pc","Ph","Ra", "Sc","Tn","Tw","Um","Vf","Mock")
full <- c("antecubital fossa","alar crease","axilla","back","cheek","external_auditory_canal","forehead",
          "glabella","hypothenar palm","inguinal_crease","interdigital_web_space","manubrium","nare",
          "negative_control","occiput","hypothenar palm","popliteal_fossa","plantar_heel",
          "retroaricular_crease","scalp","toenail","toe_web_space","umbilicus","volar_forearm","Mock")
microenvironment <- c("moist","sebaceous","moist","sebaceous","sebaceous","sebaceous","sebaceous",
                      "sebaceous","dry","moist","moist","sebaceous","moist","negative",
                      "sebaceous","dry","moist","foot","sebaceous","sebaceous","foot","foot","moist","dry","Mock")
sites_map <- data.frame(SkinSite,full,microenvironment)
meta <- merge(meta,sites_map,by="SkinSite")

# Calculate alpha diversity metrics ---------------------------------------

## Create phyloseq object for Byrd 2017 metagenomes

# make a taxonomy table by splitting the relative abundance table Taxa column (Domain;phylum;class;etc)
taxonomy_table <- as.matrix(separate(Byrd, Species , into = c("Domain","Phylum","Class","Order","Family","Genus","Species"), sep = ";", remove = FALSE,
                                     convert = FALSE, extra = "warn", fill = "warn"))
taxonomy_table <- taxonomy_table[,2:8] # taxonomy
rownames = Byrd$Species
phyloseq_abundance <- data.matrix(Byrd[,2:ncol(Byrd)]) # abundances
rownames(phyloseq_abundance) = rownames
rownames(taxonomy_table) = rownames
ABUND = otu_table(phyloseq_abundance, taxa_are_rows = TRUE) # abundance table
TAX = tax_table(taxonomy_table) # taxonomy table
physeq = phyloseq(ABUND, TAX) # phyloseq object

diversity <- estimate_richness(physeq = physeq) # calculate alpha diversity metrics
diversity$SampleID <- rownames(diversity)

# Function to plot Bray Curtis dissimilarity - by microenvironment -----------------

plot_NMDS <- function(site, meta, abundance_long, abundance_table, diversity){
  
  # filter metadata for samples from specified site
  meta <- meta %>% filter(microenvironment == site)
  # filter abundance data for samples from specified site
  SampleIDs <- unique(meta$SampleID)
  abundance_long <- abundance_long %>% filter(SampleID %in% SampleIDs)
  
  ## Identify which SampleIDs have Corynebacterium B12 producers and add up their relative abundances in each sample
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
                                     Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. sy09" |
                                     Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium aquilae" |
                                     Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. NML93-0612" |
                                     Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium xerosis"
  )
  
  
  # add up the abundance of each B12 producer in each sample
  B12 <- aggregate(B12$Abundance, by=list(SampleID=B12$SampleID), FUN=sum)
  colnames(B12) <- c("SampleID","Abundance")
  
  # merge metadata with B12 producer abundance data
  B12 <- merge(meta, B12[c("SampleID","Abundance")], by="SampleID")
  
  abundance <- abundance_table[,c("Species",SampleIDs)] # relative abundance table
  
  print(site)
  
  B12 <- left_join(B12, diversity) # add in alpha diversity
  B12$log10abundance <- log10(B12$Abundance) # take log of coryne producer abundance
  B12 <- B12 %>% filter(!is.na(log10abundance)) # exclude a few samples with no coryne B12 producers 
  
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
  #NMDS using Vegan package
  mds <- metaMDS(beta_dist, k=2, trymax=100)
  mds_data <- as.data.frame(mds$points)
  mds_data$SampleID <- rownames(mds_data)
  mds_data <- dplyr::left_join(B12,mds_data)
  mds_data$SCORAD <- as.numeric(mds_data$SCORAD)
  
  # plot NMDS
  plot <- ggplot(mds_data, aes(x=MDS1, y=MDS2, color=log10abundance,size=Shannon,shape=Timepoint)) + geom_point() +
    labs(shape = "Timepoint", color = "Log10 Corynebacterium\ncobamide producer abundance") +
    theme_classic() +
    ggtitle(site) +
    scale_color_gradient(low="blue",high="orange") +
    guides(color = "colorbar", size= "legend") +
    scale_size(range=c(1,7))

  plot

}

# Generate NMDS plots for each microenvironment
# gathered = relative abundances in long format
# relativeAbundances = relative abundance table
# diversity = diversity metrics computed on counts table
sebaceous_plot <- plot_NMDS("sebaceous",meta = meta, abundance_long = gather, abundance_table = relativeAbundances, diversity = diversity)
moist_plot <- plot_NMDS("moist", meta = meta, abundance_long = gather, abundance_table = relativeAbundances, diversity = diversity)
dry_plot <- plot_NMDS("dry", meta = meta, abundance_long = gather, abundance_table = relativeAbundances, diversity = diversity)

# final plot - Supplemental Figure 5
plot_grid(sebaceous_plot, moist_plot, dry_plot)


# Coryne cobamide producer abundance - disease state analysis ----------------------

## Identify which SampleIDs have B12 producers and add up their relative abundances in each sample
B12 <- gather %>% filter(Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium amycolatum" |
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
                             Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. sy09" |
                             Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium aquilae" |
                             Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium sp. NML93-0612" |
                             Species == "Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium xerosis"
)

# add up the abundance of each B12 producer in each sample
B12 <- aggregate(B12$Abundance, by=list(SampleID=B12$SampleID), FUN=sum)
colnames(B12) <- c("SampleID","Abundance")

# merge metadata with B12 producer abundance data and diversity data
B12 <- merge(meta, B12[c("SampleID","Abundance")], by="SampleID")
B12 <- left_join(B12, diversity)
B12 <- B12 %>% distinct(SampleID, .keep_all = TRUE)

B12$log10abundance <- log10(B12$Abundance) # take log of coryne producer abundance
B12 <- B12 %>% filter(!is.infinite(log10abundance)) # exclude a few samples with no coryne B12 producers 
avg <- B12 %>% group_by(Subject_Num,SkinSite,Timepoint) %>% summarise(mean(log10abundance)) # average log abundance
avg <- left_join(avg,B12)

# plot cobamide-producing coryne (CPC) abundance by disease state - Figure 6B
ggplot(B12, aes(x=Timepoint,y=log10abundance, fill=Group)) + geom_boxplot() + geom_point() +
  theme_classic() + ylab("Log10 Corynebacterium cobamide producer abundance") + 
  scale_fill_manual(values=c("#DEEFB7","#5FB49C"))

B12 %>% count(Timepoint)

# statistics
kruskal.test(log10abundance ~ Timepoint, data = B12) # p = 0.000204
dunnTest(log10abundance ~ Timepoint, data = B12, method="bonferroni")

avg$SkinSite <- factor(avg$SkinSite, levels = c("Pc","Ac","Ic","Ra","Oc","Gb","Fh","Vf"))

# plot cobamide-producing coryne (CPC) abundance by skin site and disease state - Figure 6C
ggplot(avg, aes(x=Timepoint,y=avg$`mean(log10abundance)`, color=Group, group=Subject_ID)) + geom_line(color="black") +
  geom_point(size=2) +   geom_point(shape = 1,size = 2,colour = "black") + 
  facet_wrap(~SkinSite) + theme_classic() +  ylab("Log10 Corynebacterium cobamide producer abundance") + 
  scale_color_manual(values=c("#DEEFB7","#5FB49C", "red")) 
