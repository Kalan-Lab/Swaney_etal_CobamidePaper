
# Script 7 ----------------------------------------------------------------
# Supplemental Figure 2

# Cobamide Skin Microbiome Manuscript
# Mock community analysis
# Written by Mary Hannah Swaney
# University of Wisconsin - Madison

# Load required packages and data -----------------------------------------

library(readr)
library(dplyr)
library(tidyr)
library(readxl)
library(data.table)
library(cowplot)
library(pheatmap)

# set working directory to the code directory
setwd("/Users/mhswaney/GitHub/scripts/mhswaney/cobamide_paper/CobamidePaperCode/")

mock <- read_csv("7_MockAnalysis/data/data_for_mock.csv")

# cobamide-dependent enzymes  20-member community --------------------------

# filter only cobamide-dependent or single copy gene results
cob_dep <- mock %>% filter( Protein.Description == "Cobamide dependent enzymes" | 
                              Protein.Description == "Single copy genes")

cob_dep <- cob_dep %>% filter(new_est_reads > 1)

# add up the hits for each taxa, for each gene
Mock20 <- cob_dep %>% filter(SampleID == "LKMN002-3-94_S293") %>% group_by(genus, Gene.Ortholog.Group, SampleID) %>% summarise(sum(NormHits))

colnames(Mock20) <- c("Genus","Gene","SampleID","NormHits")

# filter out single-copy core gene hits and non-specific HMMs
Mock20 <- Mock20 %>% filter( Gene != "nusA" &
                               Gene != "smpB" & 
                               Gene != "rptR")

# re-factor the levels for visual purposes
Mock20$Gene <- factor(Mock20$Gene, levels = c("rpoB","D-ornithine aminomutase","mutA","nrdJ_Z","metH","queG","eutB","pduC","D-lysine aminomutase","glmE",
                                              "B12-dep radical SAM","btuB"))

# color palette for plot
mhs45rainbow = c("#771155", "#AA4488", "#CC99BB", "#e6c5e0","#fae6f8",
                 "#114477", "#4477AA", "#77AADD", "#b9cfeb","#e6ebfa",
                 "#117777", "#44AAAA", "#77CCCC", "#bae6e0","#e1f4f5",
                 "#117744", "#44AA77", "#88CCAA", "#bce3d0","#dff2e5",
                 "#777711", "#AAAA44", "#DDDD77", "#e7eba4","#f3f5cb",
                 "#774411", "#AA7744", "#DDAA77", "#edc99f","#f7e1c8",
                 "#771122", "#AA4455", "#DD7788", "#e8aeae","#f7dada",
                 "#4c285c","#9167a3","#d0a3e3","#dac2ed","#ebdff5",
                 "#636363","#8f8f8f","#c9c9c9","#e3e3e3","#f0f0f0")

# composition of 20-member staggered mock community
mock_real <- data.frame(Genus = c("Acinetobacter","Bacillus","Phocaeicola",
                                  "Bifidobacterium","Clostridium",
                                  "Cutibacterium", "Deinococcus","Enterococcus",
                                  "Escherichia","Helicobacter","Lactobacillus",
                                  "Neisseria","Porphyromonas","Pseudomonas",
                                  "Luteovulum","Schaalia","Staphylococcus",
                                  "Staphylococcus","Streptococcus","Streptococcus"),
                        `NormHits` = c(0.18,1.80,0.02,0.02,1.80,0.18,0.02,0.02,18,0.18,.18,0.18,
                                       18,1.8,18,0.02,1.8,18,1.8,18))
mock_real$SampleID <- "20 strain staggered mix"
mock_real$Gene <- "True abundance"

# composition of 6-member even mock community
mock_real_skin <- data.frame(Genus=c("Acinetobacter johnsonii","Corynebacterium striatum","Micrococcus luteus",
                                     "Cutibacterium acnes","Staphylococcus epidermidis","Streptococcus mitis"),
                             `NormHits` = c(16.7,16.7,16.7,16.7,16.7,16.7))
mock_real_skin$SampleID <- "Skin microbiome genomic mix"
mock_real_skin$Gene <- "True abundance"


mock_communities_real <- rbind(mock_real,mock_real_skin)
colnames(mock_communities_real) <- c("Genus","NormHits","SampleID", "Gene")

plot <- rbind(mock_communities_real,Mock20)

Mock20 <- plot %>% filter(SampleID == "20 strain staggered mix" | SampleID == "LKMN002-3-94_S293")

Mock20 <-Mock20 %>% group_by(Genus,Gene, SampleID) %>% summarise(sum(NormHits))
colnames(Mock20) <- c("Genus","Gene","SampleID","NormHits")

# re-factor the levels for visual purposes
Mock20$Gene <- factor(Mock20$Gene, levels = c("rpoB","D-ornithine aminomutase","mutA","nrdJ_Z","metH","queG","eutB","pduC","D-lysine aminomutase","glmE",
                                             "B12-dep radical SAM","btuB"))


p1 <- ggplot(Mock20, aes(x=Gene, y=NormHits, fill=Genus)) +
  geom_bar(stat="identity", color = "black") +
  xlab("") +  ylab("") + 
  theme_classic() +
  labs(fill = "Genus") +
  scale_fill_manual(values=mhs45rainbow) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust=1),
    panel.background = element_rect(colour = "black", size=1),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    plot.title = element_text(hjust = 0.5)) + ggforce::facet_row(~SampleID, scales = 'free', space = 'free')

p2 <- ggplot(Mock20, aes(x=Gene, y=NormHits, fill=Genus)) +
  geom_bar(stat="identity", position="fill", color="black") +
  xlab("") +  ylab("") + 
  theme_classic() +
  labs(fill = "Genus") +
  scale_fill_manual(values = mhs45rainbow) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust=1),
    panel.background = element_rect(colour = "black", size=1),
    strip.background = element_blank(),
    strip.text.x = element_blank()) + ggforce::facet_row(~SampleID, scales = 'free', space = 'free')

# Supplemental Figure 2C
plot_grid(p2,p1, nrow=2)

# cobamide-dependent enzymes  6-member community --------------------------

# add up the hits for each taxa, for each gene
Mock6 <- cob_dep %>% filter(SampleID == "LKMN002-3-95_S294") %>% group_by(species, Gene.Ortholog.Group, SampleID) %>% summarise(sum(NormHits))

colnames(Mock6) <- c("Species","Gene","SampleID","NormHits")

# filter out single-copy core gene hits
Mock6 <- Mock6 %>% filter( Gene != "nusA" &
                             Gene != "smpB" & 
                             Gene != "rptR")

# re-factor the levels for visual purposes
Mock6$Gene <- factor(Mock6$Gene, levels = c("rpoB","D-ornithine aminomutase","mutA","nrdJ_Z","metH","queG","eutB","pduC","D-lysine aminomutase","glmE",
                                            "B12-dep radical SAM","btuB"))

mock_real_skin <- data.frame(Species=c("Acinetobacter johnsonii","Corynebacterium striatum","Micrococcus luteus",
                                       "Cutibacterium acnes","Staphylococcus epidermidis","Streptococcus mitis"),
                             `NormHits` = c(16.7,16.7,16.7,16.7,16.7,16.7))
mock_real_skin$SampleID <- "Skin microbiome genomic mix"
mock_real_skin$Gene <- "True abundance"

colnames(mock_real_skin) <- c("Species","NormHits","SampleID", "Gene")

plot <- rbind(mock_real_skin,Mock6)

Mock6 <- plot %>% filter( SampleID == "LKMN002-3-95_S294" | SampleID == "Skin microbiome genomic mix")

Mock6 <-Mock6 %>% group_by(Species,Gene, SampleID) %>% summarise(sum(NormHits))
colnames(Mock6) <- c("Species","Gene","SampleID","NormHits")

# re-factor the levels for visual purposes
Mock6$Gene <- factor(Mock6$Gene, levels = c("True abundance", "rpoB","D-ornithine aminomutase","mutA","nrdJ_Z","metH","queG","eutB","pduC","D-lysine aminomutase","glmE",
                                            "B12-dep radical SAM","btuB"))

p3 <- ggplot(Mock6, aes(x=Gene, y=NormHits, fill=Species)) +
  geom_bar(stat="identity", color = "black") +
  xlab("") +  ylab("") + 
  theme_classic() +
  labs(fill = "Species") +
  scale_fill_manual(values=mhs45rainbow) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust=1),
    panel.background = element_rect(colour = "black", size=1),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    plot.title = element_text(hjust = 0.5)) + ggforce::facet_row(~SampleID, scales = 'free', space = 'free')

p4 <- ggplot(Mock6, aes(x=Gene, y=NormHits, fill=Species)) +
  geom_bar(stat="identity", position="fill", color="black") +
  xlab("") +  ylab("") + 
  theme_classic() +
  labs(fill = "Family") +
  scale_fill_manual(values = mhs45rainbow) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust=1),
    panel.background = element_rect(colour = "black", size=1),
    strip.background = element_blank(),
    strip.text.x = element_blank()) + ggforce::facet_row(~SampleID, scales = 'free', space = 'free')

# Supplemental Figure 2A
plot_grid(p4,p3, nrow=2)

# Heatmaps ----------------------------------------------------------------

## 20-member staggered 

mock_real_species <- c("Acinetobacter baumannii","Bacillus cereus","Phocaeicola vulgatus",
                       "Bifidobacterium adolescentis","Clostridium beijerinckii",
                       "Cutibacterium acnes", "Deinococcus radiodurans","Enterococcus faecalis",
                       "Escherichia coli","Helicobacter pylori","Lactobacillus gasseri",
                       "Neisseria meningitidis","Porphyromonas gingivalis","Pseudomonas aeruginosa",
                       "Luteovulum sphaeroides","Schaalia odontolytica","Staphylococcus aureus",
                       "Staphylococcus spidermidis","Streptococcus mutans","Streptococcus agalactiae")

filtered_staggered <- mock %>% filter(name %in% mock_real_species & SampleID == "LKMN002-3-94_S293")
filtered_staggered <- filtered_staggered[,c(3,8,9)]

filtered_staggered <-filtered_staggered %>% group_by(name,Gene.Ortholog.Group) %>% summarise(new_est_reads=sum(new_est_reads))

matrix <- filtered_staggered %>% spread(Gene.Ortholog.Group, new_est_reads)
matrix[is.na(matrix)] <- 0 
colnames(matrix)
Ba <- data.frame(name="Bifidobacterium adolescentis") # no reads identified for Bifidobacterium adolescentis, create empty row
matrix <- rbind(matrix,Ba)
Sa <- data.frame(name="Streptococcus agalactiae") # no reads identified for Streptococcus agalactiae, create empty row
matrix <- rbind(matrix,Sa)
matrix[is.na(matrix)] <- 0 

matrix_order <- c("name","rpoB","cobI_cbiL","cobM_cbiF","cobJ_cbiH","cobL","cobH_cbiC","cobB_cbiA","cobQ_cbiP","cobD_cbiB","cobP_cobU","cobS","bluB","cbiZ",
                  "D-ornithine aminomutase","mutA","nrdJ_Z","metH","queG","eutB","D-lysine aminomutase","rptR","btuB")

matrix <- matrix[,matrix_order]
species <- matrix$name
matrix <- as.matrix(matrix[,2:ncol(matrix)])
rownames(matrix) <- species
matrix[matrix > 0] <- 1 # any species with more than 1 read, set as 1

# Supplemental Figure 2D
pheatmap::pheatmap(mat=matrix,cluster_rows = FALSE, cluster_cols = FALSE, col=c("grey80","black"))

## 6-member even 

filtered_even <- mock %>% filter(name %in% mock_real_skin$Species & SampleID == "LKMN002-3-95_S294")

filtered_even <- filtered_even[,c(3,8,9)]
filtered_even <-filtered_even %>% group_by(name,Gene.Ortholog.Group) %>% summarise(new_est_reads=sum(new_est_reads))

matrix <- filtered_even %>% spread(Gene.Ortholog.Group, new_est_reads)
matrix[is.na(matrix)] <- 0 

matrix_order <- c("name","rpoB","cobI_cbiL","cobM_cbiF","cobJ_cbiH","cobL","cobH_cbiC","cobB_cbiA","cobQ_cbiP","cobD_cbiB","cobP_cobU","cobS","bluB",
                  "D-ornithine aminomutase","mutA","nrdJ_Z","metH","queG","eutB")

matrix <- matrix[,matrix_order]
species <- matrix$name
matrix <- as.matrix(matrix[,2:ncol(matrix)])
rownames(matrix) <- species
matrix[matrix > 0] <- 1 # any species with more than 1 read, set as 1

# Supplemental Figure 2B
pheatmap::pheatmap(mat=matrix,cluster_rows = FALSE, cluster_cols = FALSE, col=c("grey80","black"))


# False positive rate -----------------------------------------------------

# taxonomic classification of all sequencing reads from mock community samples
mock_reads_unfiltered <- read_csv("7_MockAnalysis/data/MockCommunityCounts_Unfiltered.csv")

# 20 member community
all_mock <- mock %>% filter(SampleID == "LKMN002-3-94_S293") %>% group_by(name) %>% summarise(all_reads = sum(new_est_reads))

Mock20_ForComparison <- mock_reads_unfiltered[,c("Species","LKMN002-3-94_S293")]
Mock20_ForComparison <- Mock20_ForComparison %>% separate(Species, into=c("Domain","Phylum","Class","Order","Family","Genus","name"), sep=";")


Mock20Compare <-  left_join(Mock20_ForComparison,all_mock)
Mock20Compare[is.na(Mock20Compare)] <- 0

# find which species are likely contaminants - present in both samples
# any species that is present from the hmm hits and not in the actual sample can assumed to be a false positive hit
Mock20Compare$contaminant <- ifelse(!Mock20Compare$name %in% mock_real_species & !Mock20Compare$`LKMN002-3-94_S293` > 0, TRUE,FALSE)

false_pos_20 <- Mock20Compare %>% filter(contaminant == TRUE)
false_rate_20 <- (sum(false_pos_20$all_reads) / sum(Mock20Compare$all_reads)) * 100

# 6 member community
all_mock <- mock %>% filter(SampleID == "LKMN002-3-95_S294") %>% group_by(name) %>% summarise(all_reads = sum(new_est_reads))

Mock6_ForComparison <- mock_reads_unfiltered[,c("Species","LKMN002-3-95_S294")]
Mock6_ForComparison <- Mock6_ForComparison %>% separate(Species, into=c("Domain","Phylum","Class","Order","Family","Genus","name"), sep=";")

Mock6Compare <-  left_join(all_mock, Mock6_ForComparison)
Mock6Compare[is.na(Mock6Compare)] <- 0

# find which species are likely contaminants - present in both samples
# any species that is present from the hmm hits and not in the actual sample can assumed to be a false positive hit
Mock6Compare$contaminant <- ifelse(!Mock6Compare$name %in% mock_real_skin$Species & !Mock6Compare$`LKMN002-3-95_S294` > 0, TRUE,FALSE)

false_pos_6 <- Mock6Compare %>% filter(contaminant == TRUE)
false_rate_6 <- (sum(false_pos_6$all_reads) / sum(Mock6Compare$all_reads)) * 100

# false positive rate
mean(c(false_rate_20, false_rate_6))




