
# Script 1 ----------------------------------------------------------------
# Figures 1-3A, Supplemental Figures 1,3-4

# Cobamide Skin Microbiome Manuscript
# Cobamide Metagenome Survey
# Written by Mary Hannah Swaney
# University of Wisconsin - Madison

# Load required packages and data -----------------------------------------

library(readr)
library(tidyverse)
library(cowplot)
library(pheatmap)

# set working directory to the code directory
setwd("/Users/mhswaney/GitHub/scripts/mhswaney/cobamide_paper/CobamidePaperCode/")

# HMM results from LKMB002, Oh 2016, and Hannigan 2015 metagenomes
hmm <- read_csv("1_CobamideMetagenomeSurvey/data/HMM_Results.csv")
length(unique(hmm$SampleID)) # total metagenomes analyzed

infernal <- read_csv("1_CobamideMetagenomeSurvey/data/INFERNAL_Results.csv")

# COBAMIDE DEPENDENCE #####################################################

# filter only cobamide-dependent or single copy gene results
cob_dep <- hmm %>% filter( Protein.Description == "Cobamide dependent enzymes" | 
                             Protein.Description == "Single copy genes")

# cobamide-dependent enzymes by microenvironment --------------------------

# add up the hits for each taxa, for each gene, within each microenvironment
summed <- cob_dep %>% group_by(family, Gene.Ortholog.Group, microenvironment) %>% summarise(sum(NormHits))

# call taxa with fewer than a certain # of hits "Other"- cutoff depends on how many taxa can be shown in the legend
summed$Name_filtered <- ifelse(summed$`sum(NormHits)` <= 0.00025, "Other", as.character(summed$family))
summed <- summed %>% group_by(family) %>% mutate(Plot_name = case_when(
  any(Name_filtered != "Other") ~ family, 
  any(Name_filtered == "Other") ~ Name_filtered))

# add up the hits for the "Other" taxa
summed <- summed %>% group_by(Plot_name,Gene.Ortholog.Group, microenvironment) %>% summarise(sum(`sum(NormHits)`))
colnames(summed) <- c("Plot_name","Gene","microenvironment","NormHits")

# filter out single-copy core gene hits
summed <- summed %>% filter( Gene != "nusA" &
                               Gene != "smpB" &
                               Gene!= "rpoB")

## add in rpoB - for plotting purposes
# all taxa for cob-dep hits
taxa <- unique(summed$Plot_name)

# add up rpoB hits for each taxa, for each gene, in each microenvironment
rpoB_summed <- cob_dep %>% group_by(family,Gene.Ortholog.Group,microenvironment) %>% summarise(sum(NormHits))
colnames(rpoB_summed) <- c("Name","Gene","microenvironment","NormHits")
rpoB_summed <- rpoB_summed %>% filter(Gene == "rpoB")

# rename other taxa that aren't in cob-dep list to "Other"
rpoB_summed$Plot_name <- ifelse(rpoB_summed$Name %in% taxa, as.character(rpoB_summed$Name),"Other")

# normalized rpoB hits in each microenvironment
rpoB_summed <- rpoB_summed %>% group_by(Plot_name,Gene, microenvironment) %>% summarise(sum(NormHits))
colnames(rpoB_summed) <- c("Plot_name","Gene","microenvironment","NormHits")

# order rpoB taxa from most to least # of hits
rpoB_summed <- rpoB_summed[order(-rpoB_summed$NormHits),]
unique(rpoB_summed$Plot_name)

# bind together the cob-dep and rpoB dataframes
summed <- bind_rows(summed,rpoB_summed)
unique(summed$Gene)

# re-factor the levels for visual purposes
summed$Gene <- factor(summed$Gene, levels = c("rpoB","D-ornithine aminomutase","mutA","nrdJ_Z","metH","queG","eutB","pduC","D-lysine aminomutase","glmE",
                                              "B12-dep radical SAM","btuB"))

factors <- unique(rpoB_summed$Plot_name) # taxa for plotting
factors <- factors[c(1:3,5:38)] # taxa names minus other
factors <- c("Other",factors) # re-order taxa names

# specify the factor levels to go in the order of hits from highest to lowest 
summed$Plot_name <- factor(summed$Plot_name, levels = factors)

# color palette for plot
mhs40rainbow = c("#C6C6C6","#771155", "#AA4488", "#CC99BB", "#e6c5e0","#fae6f8",
                 "#114477", "#4477AA", "#77AADD", "#b9cfeb","#e6ebfa",
                 "#117777", "#44AAAA", "#77CCCC", "#bae6e0","#e1f4f5",
                 "#117744", "#44AA77", "#88CCAA", "#bce3d0","#dff2e5",
                 "#777711", "#AAAA44", "#DDDD77", "#e7eba4","#f3f5cb",
                 "#774411", "#AA7744", "#DDAA77", "#edc99f","#f7e1c8",
                 "#771122", "#AA4455", "#DD7788", "#e8aeae","#f7dada",
                 "#4c285c","#9167a3","#d0a3e3","#dac2ed","#ebdff5")

total_sebaceous <- summed %>% filter(microenvironment == "sebaceous")
total_moist <- summed %>% filter(microenvironment == "moist")
total_dry <- summed %>% filter(microenvironment == "dry")
total_foot <- summed %>% filter(microenvironment == "foot")

# plot sebaceous total cobamide-dependent hits
p1 <- ggplot(total_sebaceous, aes(x=Gene, y=NormHits, fill=Plot_name)) +
  geom_bar(stat="identity") +
  xlab("") +  ylab("") + 
  facet_wrap(~microenvironment) +
  theme_classic() +
  labs(fill = "Family") +
  scale_fill_manual(values=mhs40rainbow) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust=1),
    panel.background = element_rect(colour = "black", size=1),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.position="none",
    plot.title = element_text(hjust = 0.5)) 

# plot moist total cobamide-dependent hits
p2 <- ggplot(total_moist, aes(x=Gene, y=NormHits, fill=Plot_name)) +
  geom_bar(stat="identity") +
  xlab("") + ylab("") + 
  facet_wrap(~microenvironment) +
  theme_classic() +
  labs(fill = "Family") +
  scale_fill_manual(values=mhs40rainbow) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust=1),
    panel.background = element_rect(colour = "black", size=1),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.position="none",
    plot.title = element_text(hjust = 0.5))

# plot dry total cobamide-dependent hits
p3 <- ggplot(total_dry, aes(x=Gene, y=NormHits, fill=Plot_name)) +
  geom_bar(stat="identity") +
  xlab("") + ylab("") + 
  facet_wrap(~microenvironment) +
  theme_classic() +
  labs(fill = "Family") +
  scale_fill_manual(values=mhs40rainbow) +
  theme(
    axis.text.x = element_text(angle = 45, hjust=1),
    panel.background = element_rect(colour = "black", size=1),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.position="none",
    plot.title = element_text(hjust = 0.5)) 

# plot foot total cobamide-dependent hits
p4 <- ggplot(total_foot, aes(x=Gene, y=NormHits, fill=Plot_name)) +
  geom_bar(stat="identity") +
  xlab("") + ylab("") + 
  facet_wrap(~microenvironment) +
  theme_classic() +
  labs(fill = "Family") +
  scale_fill_manual(values=mhs40rainbow) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust=1),
    panel.background = element_rect(colour = "black", size=1),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.position="none",
    plot.title = element_text(hjust = 0.5))

# plot sebaceous relative proportion cobamide-dependent hits
p5 <- ggplot(total_sebaceous, aes(x=Gene, y=NormHits, fill=Plot_name)) +
  geom_bar(stat="identity", position="fill") +
  xlab("") +  ylab("") + 
  facet_wrap(~microenvironment) +
  theme_classic() +
  labs(fill = "Family") +
  scale_fill_manual(values=mhs40rainbow) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust=1),
    panel.background = element_rect(colour = "black", size=1),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.position="none") 

# plot moist relative proportion cobamide-dependent hits
p6 <- ggplot(total_moist, aes(x=Gene, y=NormHits, fill=Plot_name)) +
  geom_bar(stat="identity", position="fill") +
  xlab("") + ylab("") + 
  facet_wrap(~microenvironment) +
  theme_classic() +
  labs(fill = "Family") +
  scale_fill_manual(values=mhs40rainbow) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust=1),
    panel.background = element_rect(colour = "black", size=1),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.position="none")

# plot dry relative proportion cobamide-dependent hits
p7 <- ggplot(total_dry, aes(x=Gene, y=NormHits, fill=Plot_name)) +
  geom_bar(stat="identity", position="fill") +
  xlab("") + ylab("") + 
  facet_wrap(~microenvironment) +
  theme_classic() +
  labs(fill = "Family") +
  scale_fill_manual(values=mhs40rainbow) +
  theme(
    axis.text.x = element_text(angle = 45, hjust=1),
    panel.background = element_rect(colour = "black", size=1),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.position="none") 

# plot foot relative proportion cobamide-dependent hits
p8 <- ggplot(total_foot, aes(x=Gene, y=NormHits, fill=Plot_name)) +
  geom_bar(stat="identity", position="fill") +
  xlab("") + ylab("") + 
  facet_wrap(~microenvironment) +
  theme_classic() +
  labs(fill = "Family") +
  scale_fill_manual(values=mhs40rainbow) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust=1),
    panel.background = element_rect(colour = "black", size=1),
    strip.background = element_blank(),
    strip.text.x = element_blank())

# get legend from p8 - legend is the same for p1-p8
legend<- get_legend(
  p8 + 
    guides(color = guide_legend(nrow = 1))+
    theme(legend.box.margin = margin(0, 0, 0, 12))
)

top <- plot_grid(p5, p6, p7, p8 + theme(legend.position = "none"), ncol=4, align="vh")
bottom <- plot_grid(p1, p2, p3, p4, ncol=4, align="vh")

## plot total normalized hits and relative proportions for cobamide-dependent enzymes - Figure 2
plot_grid(top,bottom, rel_widths=c(4,1),nrow=2, align="vh") 

## legend for plot
plot(legend)

# Unique species for cob-dependent gene hits ------------------------------

# group by cobamide dependent enzymes or rpoB
cob_dep <- hmm %>% filter((Protein.Description == "Cobamide dependent enzymes" | Gene.Ortholog.Group == "rpoB"))

unique_taxa <- cob_dep %>% group_by(Gene.Ortholog.Group,Protein.Description) %>% summarise(unique_taxa = n_distinct(name))

taxa_order <- unique_taxa[order(-unique_taxa$unique_taxa),]$Gene.Ortholog.Group
unique_taxa$Gene.Ortholog.Group <- factor(unique_taxa$Gene.Ortholog.Group, levels=taxa_order)

unique_taxa$Protein.Description <- factor(unique_taxa$Protein.Description, levels = c("Single copy genes","Cobamide dependent enzymes"))

## plot unique species per cobamide-dependent gene - Supplemental Figure 4
ggplot(unique_taxa, aes(x=Gene.Ortholog.Group, y=unique_taxa, fill=Protein.Description)) + 
  geom_col() + 
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  ylab("Unique species") +
  geom_text(aes(label=unique_taxa), vjust=-0.25) + 
  scale_fill_manual(values=c("#3034ab", "#8d94e3")) + 
  theme(legend.position = "bottom",
        axis.title.x = element_blank(),
        legend.title = element_blank())

# number of unique species with cobamide-dependent enzymes
cob_dep_num_species <- hmm %>% filter(Protein.Description == "Cobamide dependent enzymes") %>% summarise(unique_taxa = n_distinct(name))

# number of unique species encoding rpoB
rpoB_num_species <- hmm %>% filter(Gene.Ortholog.Group == "rpoB") %>% summarise(unique_taxa = n_distinct(name))


# COBAMIDE BIOSYNTHESIS ###################################################

# identify top 40 taxa in overall dataset - by # of rpoB hits
cob_biosynth_rpoB <- hmm %>% filter((Gene.Ortholog.Group == "rpoB"))
summed <- cob_biosynth_rpoB %>% group_by(species) %>% summarise(sum(NormHits))
top <- summed[order(-summed$`sum(NormHits)`),]
top_40 <- top[1:40,1]  

cob_biosynth <- hmm %>% filter((Protein.Description == "Cobamide biosynthesis" | Gene.Ortholog.Group == "rpoB"))

# identify de novo producers ----------------------------------------------

# count # of reads for each cobamide biosynthesis gene-species pairing
all_taxa_reads <- cob_biosynth %>% filter(Gene.Ortholog.Group != "cobA_cysG" & Gene.Ortholog.Group != "rpoB" & Gene.Ortholog.Group != "bluB" & Gene.Ortholog.Group != "cbiZ") %>%
  group_by(name,Gene.Ortholog.Group) %>% summarise(total_reads = sum(new_est_reads))

# count how many biosynthesis genes are present for each species
taxa_gene_count <- all_taxa_reads%>%
  group_by(name) %>%
  summarise(distinct_genes = n_distinct(Gene.Ortholog.Group))

# filter out species that have at least 7 identified cob biosynthesis genes
greaterthan5 <- taxa_gene_count %>% filter(distinct_genes >= 5)


# microenvironment heatmap - species-level --------------------------------

# group all metagenomic samples together by microenvironment and sum up hits per gene
summed_hm <- cob_biosynth %>% group_by(species,Gene.Ortholog.Group,microenvironment) %>% summarise(sum(NormHits))
# only include hits to the top 40 taxa
summed_filter <- summed_hm %>% filter(species %in% top_40$species)

# re-name species not in the top 40 taxa as "Other"
other_genes <- summed_hm %>% filter(!(species %in% top_40$species)) %>% group_by(Gene.Ortholog.Group,microenvironment) %>% summarise(sum(`sum(NormHits)`))
other_genes$name <- "Other"
colnames(other_genes) <- c("Gene.Ortholog.Group","microenvironment","sum(NormHits)","species")
summed_filter <- rbind(summed_filter, other_genes)

# get hits for rpoB in top 40 species
rpoB <- hmm %>% filter((Gene.Ortholog.Group == "rpoB")) %>% group_by(species,Gene.Ortholog.Group,microenvironment) %>% summarise(sum(NormHits))
rpoB_filter <- rpoB %>% filter(species %in% top_40$species)

# rename species not in top 40 list as "Other"
other <- rpoB %>% filter(!(species %in% top_40$species)) %>% group_by(Gene.Ortholog.Group,microenvironment) %>% summarise(sum(`sum(NormHits)`))
other$name <- "Other"
colnames(other) <- c("Gene.Ortholog.Group","microenvironment","sum(NormHits)","species")
rpoB_filter <- rbind(rpoB_filter,other)

# join rpoB hits and cob biosynthesis hits
summed_filter <- rbind(summed_filter,rpoB_filter)

# order of microenvironments for heatmap
sites <- unique(cob_biosynth$microenvironment)

# function to plot heatmap for each micr oenvironment
microenvironment_heatmap_df <- lapply(sites, FUN = function(site){
  
  site_df <- summed_filter %>% filter(microenvironment == site)
  for_hm <- site_df[,c(1,2,4)]
  for_hm <- for_hm %>% group_by(Gene.Ortholog.Group,species) %>% summarise(sum(`sum(NormHits)`)) # sum all hits per gene for each species
  colnames(for_hm) <- c("Gene","Name","Abundance")
  # determine if any top 40 species have no cob biosynthesis hits and if so, add 0 to total hits
  missing <- lapply(top_40$species, FUN = function(taxa){
    ifelse((!(taxa %in% for_hm$Name)),taxa,NA)})
  missing_taxa <- as.data.frame(missing[!is.na(missing)])
  if (nrow(missing_taxa)>0){
    df <- data.frame(Gene = unique(cob_biosynth$Gene.Ortholog.Group), Name = missing_taxa[1,1], Abundance = as.double(0))
    for_hm <- bind_rows(for_hm[,1:3],df[,1:3])}
  colnames(for_hm) <- c("Gene","Name","Abundance")
  matrix <- spread(for_hm,Gene, Abundance) # create wide data 
  matrix_new <- matrix[,c("Name","rpoB","cobI_cbiL","cobJ_cbiH","cobM_cbiF","cobL","cobH_cbiC","cobB_cbiA","cobQ_cbiP","cobD_cbiB","cobP_cobU","cobS","bluB")] #order of genes
  matrix_new[,"cbiZ"] <- ifelse("cbiZ" %in% colnames(matrix), matrix[,"cbiZ"], NA) # code for possible error with cbiZ if it is absent from a given microenvironment
  matrix_new <- matrix_new[,c("Name","rpoB","cobI_cbiL","cobJ_cbiH","cobM_cbiF","cobL","cobH_cbiC","cobB_cbiA","cobQ_cbiP","cobD_cbiB","cobP_cobU","cobS","bluB","cbiZ")]
  matrix_new[is.na(matrix_new)] <- as.numeric(0)
  taxa <- as.list(matrix_new[,1]) # get taxa order
  matrix_new <- as.matrix(matrix_new[,2:ncol(matrix_new)]) # convert to matrix type
  relabund <- t(t(matrix_new)/rowSums(t(matrix_new))) # get relative proportions
  relabund[is.nan(relabund)] <- as.numeric(0) 
  rownames(relabund) <- taxa[[1]]
  return(relabund * 100) # multiply by 100 to get relative abundances
})

# create large matrix containing each microenvironment matrix
large_matrix <- cbind(microenvironment_heatmap_df[[1]],microenvironment_heatmap_df[[2]])
large_matrix <- cbind(large_matrix,microenvironment_heatmap_df[[3]])
large_matrix <- cbind(large_matrix,microenvironment_heatmap_df[[4]])
large_matrix <- large_matrix[,c(14:26, 1:13, 27:52)]

# get family classifications for the species
family_classification <- data.frame("species" = rownames(large_matrix))
family_classification <- left_join(family_classification,hmm[,c("species","family")])
family_classification <- unique(family_classification)
family_classification[is.na(family_classification$family),"family"] <- "NA"
family_classification_2 <- data.frame(family_classification[,2])
rownames(family_classification_2) <- family_classification$species

family_classification_2$family_classification...2. <- factor(family_classification_2$family_classification...2., levels = unique(family_classification_2$family_classification...2.))

# plotsettings
my.breaks <- c(seq(0,0.01, by=0.001),seq(0.010001, 100, by=0.3)) 
my.colors <- c(colorRampPalette(colors=c("black","black"))(length(my.breaks)/50), colorRampPalette(colors = c("#290b91","#ff4475","#fec021","#ffe9b0"))(length(my.breaks)/1.01))
ann_colors = list(family_classification...2. = c("NA"="black", Carnobacteriaceae="#AA4488", Corynebacteriaceae="#CC99BB", Lactobacillaceae="#77AADD", Lawsonellaceae="#117777", Micrococcaceae="#77CCCC", Moraxellaceae="#AAAA44", Pasteurellaceae="#DDDD77", Peptoniphilaceae= "#774411", Propionibacteriaceae="#AA7744", Staphylococcaceae="#E8CCB5", Streptococcaceae="#AA4455", Veillonellaceae="#DD7788"))

# order of heatmaps
sites

## plot species-level heatmap - Figure 1C
pheatmap(large_matrix, cluster_rows = TRUE, cluster_cols = FALSE, color = my.colors,
         breaks = my.breaks,gaps_col=c(13,26,39), 
         annotation_row = family_classification_2,
         annotation_colors = ann_colors,
         border_color = "black", )


# microenvironment heatmap - family-level ---------------------------------


# identify top 20 families in overall dataset - by # of rpoB hits
cob_biosynth_rpoB <- hmm %>% filter((Gene.Ortholog.Group == "rpoB"))
summed <- cob_biosynth_rpoB %>% group_by(family) %>% summarise(sum(NormHits))
top <- summed[order(-summed$`sum(NormHits)`),]
top_20 <- top[1:20,1]  

cob_biosynth <- hmm %>% filter((Protein.Description == "Cobamide biosynthesis" | Gene.Ortholog.Group == "rpoB"))

# group all metagenomic samples together by microenvironment and sum up hits per gene
summed_hm <- cob_biosynth %>% group_by(family,Gene.Ortholog.Group,microenvironment) %>% summarise(sum(NormHits))
# only include hits to the top 20 taxa
summed_filter <- summed_hm %>% filter(family %in% top_20$family)

other_genes <- summed_hm %>% filter(!(family %in% top_20$family)) %>% group_by(Gene.Ortholog.Group,microenvironment) %>% summarise(sum(`sum(NormHits)`))
other_genes$name <- "Other"
colnames(other_genes) <- c("Gene.Ortholog.Group","microenvironment","sum(NormHits)","family")
summed_filter <- rbind(summed_filter, other_genes)

# get hits for rpoB in top 20 family

rpoB <- hmm %>% filter((Gene.Ortholog.Group == "rpoB")) %>% group_by(family,Gene.Ortholog.Group,microenvironment) %>% summarise(sum(NormHits))
rpoB_filter <- rpoB %>% filter(family %in% top_20$family)

# rename family not in top 20 list as "Other"
other <- rpoB %>% filter(!(family %in% top_20$family)) %>% group_by(Gene.Ortholog.Group,microenvironment) %>% summarise(sum(`sum(NormHits)`))
other$name <- "Other"
colnames(other) <- c("Gene.Ortholog.Group","microenvironment","sum(NormHits)","family")
rpoB_filter <- rbind(rpoB_filter,other)

# join rpoB hits and cob biosynthesis hits
summed_filter <- rbind(summed_filter,rpoB_filter)

sites <- unique(cob_biosynth$microenvironment)
sites

# function to plot heatmap for each microenvironment
microenvironment_heatmap_df <- lapply(sites, FUN = function(site){
  
  site_df <- summed_filter %>% filter(microenvironment == site)
  for_hm <- site_df[,c(1,2,4)]
  for_hm <- for_hm %>% group_by(Gene.Ortholog.Group,family) %>% summarise(sum(`sum(NormHits)`)) # sum all hits per gene for each family
  colnames(for_hm) <- c("Gene","Name","Abundance")
  # determine if any top 20 family have no cob biosynthesis hits and if so, add 0 to total hits
  missing <- lapply(top_20$family, FUN = function(taxa){
    ifelse((!(taxa %in% for_hm$Name)),taxa,NA)})
  missing_taxa <- as.data.frame(missing[!is.na(missing)])
  if (nrow(missing_taxa)>0){
    df <- data.frame(Gene = unique(cob_biosynth$Gene.Ortholog.Group), Name = missing_taxa[1,1], Abundance = as.double(0))
    for_hm <- bind_rows(for_hm[,1:3],df[,1:3])}
  colnames(for_hm) <- c("Gene","Name","Abundance")
  matrix <- spread(for_hm,Gene, Abundance) # create wide data 
  matrix_new <- matrix[,c("Name","rpoB","cobI_cbiL","cobJ_cbiH","cobM_cbiF","cobL","cobH_cbiC","cobB_cbiA","cobQ_cbiP","cobD_cbiB","cobP_cobU","cobS","bluB")] #order of genes
  matrix_new[,"cbiZ"] <- ifelse("cbiZ" %in% colnames(matrix), matrix[,"cbiZ"], NA) # code for possible error with cbiZ if it is absent from a given microenvironment
  matrix_new <- matrix_new[,c("Name","rpoB","cobI_cbiL","cobJ_cbiH","cobM_cbiF","cobL","cobH_cbiC","cobB_cbiA","cobQ_cbiP","cobD_cbiB","cobP_cobU","cobS","bluB","cbiZ")]
  matrix_new[is.na(matrix_new)] <- as.numeric(0)
  taxa <- as.list(matrix_new[,1]) # get taxa order
  matrix_new <- as.matrix(matrix_new[,2:ncol(matrix_new)]) # convert to matrix type
  relabund <- t(t(matrix_new)/rowSums(t(matrix_new))) # get relative proportions
  relabund[is.nan(relabund)] <- as.numeric(0) 
  rownames(relabund) <- taxa[[1]]
  return(relabund * 100) # multiply by 100 to get relative abundances
})

# create large matrix containing each microenvironment matrix
large_matrix <- cbind(microenvironment_heatmap_df[[1]],microenvironment_heatmap_df[[2]])
large_matrix <- cbind(large_matrix,microenvironment_heatmap_df[[3]])
large_matrix <- cbind(large_matrix,microenvironment_heatmap_df[[4]])
large_matrix <- large_matrix[,c(14:26, 1:13, 27:52)]

# plot settings
my.breaks <- c(seq(0,0.01, by=0.001),seq(0.010001, 100, by=0.3)) 
my.colors <- c(colorRampPalette(colors=c("black","black"))(length(my.breaks)/50), colorRampPalette(colors = c("#290b91","#ff4475","#fec021","#ffe9b0"))(length(my.breaks)/1.01))
ann_colors = list(family_classification...2. = c("NA"="black", Bifidobacteriaceae="#771155", Carnobacteriaceae="#AA4488", Corynebacteriaceae="#CC99BB", Dermacoccaceae= "#114477",Enterobacteriaceae= "#4477AA", Lactobacillaceae="#77AADD", Lawsonellaceae="#117777", Micrococcaceae="#44AAAA", Moraxellaceae="#77CCCC", Pasteurellaceae="#AAAA44", Peptoniphilaceae= "#DDDD77", Propionibacteriaceae="#774411", Staphylococcaceae="#AA7744", Streptococcaceae="#DDAA77", Streptomycetaceae="#771122", Veillonellaceae="#AA4455", Xanthomonadaceae="#DD7788"))
# order of heatmaps
sites

## plot family-level heatmap - Supplemental Figure 3
pheatmap(large_matrix, cluster_rows = TRUE, cluster_cols = FALSE, color = my.colors,
         breaks = my.breaks,gaps_col=c(13,26,39), 
         annotation_colors = ann_colors,
         border_color=NA)


# per taxon cobamide biosynthesis gene distribution -----------------------


# SampleIDs and their skin site information
samples <- unique(hmm[,c("SampleID","SkinSite","microenvironment")])

# filter hits for only cobamide biosynthesis genes
cob_biosynth <- hmm %>% filter((Protein.Description == "Cobamide biosynthesis"))
cob_biosynth <- cob_biosynth[!is.na(cob_biosynth$NormHits),]

# summarise total biosynthesis hits in each microenvironment
cob_biosynth %>% group_by(microenvironment) %>% summarise(sum(new_est_reads))

# summarise total biosynthesis hits in each family
tax <- cob_biosynth %>% group_by(family) %>% summarise(sum(new_est_reads))

# add up all normalized hits for the cob biosynth genes within each sample (remove SCG and bza genes)
total_hits <- cob_biosynth %>% filter(Gene.Ortholog.Group != "smpB" &
                                        Gene.Ortholog.Group != "nusA") %>% group_by(SampleID) %>% summarise(sum(NormHits))
colnames(total_hits) <- c("SampleID","TotalCobHits")

# add up hits for each taxa within each sample
per_sample_taxon <- cob_biosynth %>% filter(Gene.Ortholog.Group != "smpB" &
                                              Gene.Ortholog.Group != "nusA") %>%
  group_by(SampleID,family) %>% summarise(sum(NormHits))
colnames(per_sample_taxon) <- c("SampleID","Taxon","TaxonHits")

# calculate the contribution of each taxon to biosynthesis hits - divide taxon hits by total biosynthesis hits within sample
per_sample_taxon <- left_join(per_sample_taxon,total_hits)
per_sample_taxon$TaxonContribution <- per_sample_taxon$TaxonHits / per_sample_taxon$TotalCobHits
per_sample_taxon <- left_join(per_sample_taxon,samples)

# get the total biosynthesis hits for each taxon within the dataset and order from most to least - for plotting purposes
sum <- cob_biosynth %>% filter(Gene.Ortholog.Group != "smpB" &
                                 Gene.Ortholog.Group != "nusA") %>%
  group_by(family) %>% summarise(sum(NormHits))
top_sum <- sum[order(-sum$`sum(NormHits)`),]

# top 6 taxa to plot - most biosynthesis hits
taxa_plot <- unique(top_sum[1:6,]$family)

# filter to only include top 6 taxa
filtered <- per_sample_taxon %>% filter(Taxon %in% taxa_plot)
# re-factor levels skin site - plotting order
filtered$SkinSite <- factor(filtered$SkinSite, levels = c("alar crease","back","cheek","external_auditory_canal",
                                                          "forehead","glabella","manubrium","occiput","retroaricular_crease","scalp","antecubital fossa","axilla","inguinal_crease",
                                                          "interdigital_web_space","nare","popliteal_fossa",
                                                          "umbilicus","hypothenar palm","volar_forearm","plantar_heel","toe_web_space","toenail"))

## plot taxon contribution - Figure 1B
ggplot(filtered, aes(x=SkinSite, y=TaxonContribution,fill=microenvironment)) + 
  geom_boxplot(outlier.size = 0.5) + facet_grid(~Taxon) + 
  theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() + geom_point(size=0.5) + scale_fill_manual(values = c("#f3eac2","#9ad3bc","#f5b461","#ec524b"))

# COBALAMIN RIBOSWITCHES ##############################################

# cobalamin riboswitches by microenvironment ------------------------------

# add up the riboswitch hits for each taxa within each microenvironment
summed <- infernal %>% group_by(family,microenvironment) %>% summarise(sum(NormHits))
colnames(summed) <- c("Plot_name","microenvironment","NormHits")

summed <- summed[order(-summed$NormHits),] # order from most to least hits
factors <- unique(summed$Plot_name)

# specify the factor levels to go in the order of hits from highest to lowest 
summed$Plot_name <- factor(summed$Plot_name, levels = factors)

# color palette for plot
mhs40rainbow = c("#771155", "#AA4488", "#CC99BB", "#e6c5e0","#fae6f8",
                 "#114477", "#4477AA", "#77AADD", "#b9cfeb","#e6ebfa",
                 "#117777", "#44AAAA", "#77CCCC", "#bae6e0","#e1f4f5",
                 "#117744", "#44AA77", "#88CCAA", "#bce3d0","#dff2e5",
                 "#777711", "#AAAA44", "#DDDD77", "#e7eba4","#f3f5cb",
                 "#774411", "#AA7744", "#DDAA77", "#edc99f","#f7e1c8",
                 "#771122", "#AA4455", "#DD7788", "#e8aeae","#f7dada",
                 "#4c285c","#9167a3","#d0a3e3","#dac2ed","#ebdff5")

## total riboswitch count by microenvironment
rib_counts <- infernal %>% group_by(microenvironment) %>% summarise(sum(new_est_reads))

total_sebaceous <- summed %>% filter(microenvironment == "sebaceous")
total_moist <- summed %>% filter(microenvironment == "moist")
total_dry <- summed %>% filter(microenvironment == "dry")
total_foot <- summed %>% filter(microenvironment == "foot")

# plot sebaceous relative proportion riboswitch hits
rp1 <- ggplot(total_sebaceous, aes(x=factor(1), y=NormHits, fill=Plot_name)) +
  geom_bar(stat="identity",position = "fill") +
  xlab(paste0("n=",rib_counts[which(rib_counts$microenvironment == "sebaceous"),"sum(new_est_reads)"])) +
  ylab("Relative proportion") + 
  facet_wrap(~microenvironment) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        legend.position="none",  
        strip.background = element_blank(),
        strip.text.x = element_blank())+
  labs(fill = "Family") +
  scale_fill_manual(values=mhs40rainbow, drop=FALSE)

# plot moist relative proportion riboswitch hits
rp2 <- ggplot(total_moist, aes(x=factor(1), y=NormHits, fill=Plot_name)) +
  geom_bar(stat="identity",position = "fill") +
  xlab(paste0("n=",rib_counts[which(rib_counts$microenvironment == "moist"),"sum(new_est_reads)"])) +
  ylab("Relative proportion") + 
  facet_wrap(~microenvironment) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        legend.position="none", strip.background = element_blank(),
        strip.text.x = element_blank())+
  labs(fill = "Family") +
  scale_fill_manual(values=mhs40rainbow, drop=FALSE)

# plot dry relative proportion riboswitch hits
rp3 <- ggplot(total_dry, aes(x=factor(1), y=NormHits, fill=Plot_name)) +
  geom_bar(stat="identity",position = "fill") +
  xlab(paste0("n=",rib_counts[which(rib_counts$microenvironment == "dry"),"sum(new_est_reads)"])) +
  ylab("Relative proportion") + 
  facet_wrap(~microenvironment) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        legend.position="none", 
        strip.background = element_blank(),
        strip.text.x = element_blank())+
  labs(fill = "Family") +
  scale_fill_manual(values=mhs40rainbow, drop=FALSE)

# plot foot relative proportion riboswitch hits
rp4 <- ggplot(total_foot, aes(x=factor(1), y=NormHits, fill=Plot_name)) +
  geom_bar(stat="identity",position = "fill") +
  xlab(paste0("n=",rib_counts[which(rib_counts$microenvironment == "foot"),"sum(new_est_reads)"])) +
  ylab("Relative proportion") + 
  facet_wrap(~microenvironment) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        legend.position="none", 
        strip.background = element_blank(),
        strip.text.x = element_blank())+
  labs(fill = "Family") +
  scale_fill_manual(values=mhs40rainbow, drop=FALSE)

# plot sebaceous zoomed in riboswitch hits
rp5 <- ggplot(total_sebaceous, aes(x=factor(1), y=NormHits, fill=Plot_name)) +
  geom_bar(stat="identity",position = "fill") +
  facet_wrap(~microenvironment) +
  theme_void()+
  theme(legend.position="none")+
  labs(fill = "Family") +
  scale_fill_manual(values=mhs40rainbow, drop=FALSE) + 
  coord_cartesian(y=c(0,0.006)) + 
  ggtitle("0.006") + 
  xlab(paste0("n=",rib_counts[which(rib_counts$microenvironment == "sebaceous"),"sum(new_est_reads)"]))

# plot moist zoomed in riboswitch hits
rp6 <- ggplot(total_moist, aes(x=factor(1), y=NormHits, fill=Plot_name)) +
  geom_bar(stat="identity",position = "fill") +
  facet_wrap(~microenvironment) +
  theme_void()+
  theme(legend.position="none")+
  labs(fill = "Family") +
  scale_fill_manual(values=mhs40rainbow, drop=FALSE) + 
  coord_cartesian(y=c(0,0.04)) +
  ggtitle("0.04") +
  xlab(paste0("n=",rib_counts[which(rib_counts$microenvironment == "moist"),"sum(new_est_reads)"]))

# plot dry zoomed in riboswitch hits
rp7 <- ggplot(total_dry, aes(x=factor(1), y=NormHits, fill=Plot_name)) +
  geom_bar(stat="identity",position = "fill") +
  facet_wrap(~microenvironment) +
  theme_void()+
  theme(legend.position="none")+
  labs(fill = "Family") +
  scale_fill_manual(values=mhs40rainbow, drop=FALSE) + 
  coord_cartesian(y=c(0,0.03)) + 
  ggtitle("0.03") +
  xlab(paste0("n=",rib_counts[which(rib_counts$microenvironment == "dry"),"sum(new_est_reads)"]))

# plot foot zoomed in riboswitch hits
rp8 <- ggplot(total_foot, aes(x=factor(1), y=NormHits, fill=Plot_name)) +
  geom_bar(stat="identity",position = "fill") +
  facet_wrap(~microenvironment) +
  theme_void()+
  labs(fill = "Family") +
  scale_fill_manual(values=mhs40rainbow, drop=FALSE) + 
  coord_cartesian(y=c(0,0.2)) + 
  ggtitle("0.2") 

legend<- get_legend(
  rp8 + 
    guides(color = guide_legend(nrow = 1))
)

## plot relative proportions for riboswitch hits - Figure 3A 
plot_grid(rp1,rp5,rp2,rp6,rp3,rp7,rp4,rp8 + theme(legend.position= "none"),nrow=1, rel_widths = c(1,1.1,1,1.1,1,1.1,1,1.1)) 

## legend for Figure 3A
plot(legend)

# total riboswitch count by species
riboswitch_species <- infernal %>% group_by(name) %>% summarise(sum(new_est_reads))


# Read counts per sample --------------------------------------------------


# cobamide biosynthesis
cob_biosynth <- hmm %>% filter(Protein.Description == "Cobamide biosynthesis")
cob_biosynth_per_sample <- cob_biosynth %>% group_by(SampleID, microenvironment) %>% summarise(sum(new_est_reads))
cob_biosynth_per_sample$microenvironment <- factor(cob_biosynth_per_sample$microenvironment, levels=c("sebaceous","moist","dry","foot"))

cob_biosynth_per_sample %>% group_by(microenvironment) %>% summarise(median_reads = median(`sum(new_est_reads)`))

#  plot per sample read count - cobamide biosynthesis
biosynth_plot <- ggplot(cob_biosynth_per_sample, aes(x=microenvironment, y=log10(`sum(new_est_reads)`), fill=microenvironment)) + geom_boxplot(outlier.size = 0.5,) +
  geom_jitter(width=0.05,size=0.5) + theme_classic() + ylab("Log 10 cobamide biosynthesis read count") + 
  xlab("Microenvironment") + 
  theme(legend.position = "none")

# cobamide dependence
cob_dep <- hmm %>% filter(Protein.Description == "Cobamide dependent enzymes")
cob_dep_per_sample <- cob_dep %>% group_by(SampleID, microenvironment) %>% summarise(sum(new_est_reads))
cob_dep_per_sample$microenvironment <- factor(cob_dep_per_sample$microenvironment, levels=c("sebaceous","moist","dry","foot"))

cob_dep_per_sample %>% group_by(microenvironment) %>% summarise(median_reads = median(`sum(new_est_reads)`))

#  plot per sample read count - cobamide dependence
dep_plot <- ggplot(cob_dep_per_sample, aes(x=microenvironment, y=log10(`sum(new_est_reads)`), fill=microenvironment)) + geom_boxplot(outlier.size = 0.5) +
  geom_jitter(width=0.05,size=0.5) + theme_classic() + ylab("Log 10 cobamide dependent enzyme read count") + 
  xlab("Microenvironment") + 
  theme(legend.position = "none")

# rpoB
rpob <- hmm %>% filter(Gene.Ortholog.Group == "rpoB")
rpoB_per_sample <- rpob %>% group_by(SampleID, microenvironment) %>% summarise(sum(new_est_reads))
rpoB_per_sample$microenvironment <- factor(rpoB_per_sample$microenvironment, levels=c("sebaceous","moist","dry","foot"))

rpoB_per_sample %>% group_by(microenvironment) %>% summarise(median_reads = median(`sum(new_est_reads)`))

#  plot per sample read count - rpoB
rpoB_plot <- ggplot(rpoB_per_sample, aes(x=microenvironment, y=log10(rpoB_per_sample$`sum(new_est_reads)`), fill=microenvironment)) + geom_boxplot(outlier.size = 0.5) +
  geom_jitter(width=0.05,size=0.5) + theme_classic() + ylab("Log 10 rpoB read count") + 
  xlab("Microenvironment") + 
  theme(legend.position = "none")

## final plot - read counts per sample - Supplemental Figure 1
plot_grid(biosynth_plot, dep_plot, rpoB_plot, nrow=1)

## Archaea 

archaea <- hmm %>% filter(superkingdom == "Archaea") %>% summarise(sum(new_est_reads))
