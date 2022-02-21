
# Script 5A ----------------------------------------------------------------
# Figure 6B-E and Supplemental Figure 12

# Cobamide Skin Microbiome Manuscript
# Corynebacterium pangenome analysis
# Written by Mary Hannah Swaney
# University of Wisconsin - Madison


# Load required packages and data -----------------------------------------

library(readxl)
library(tidyverse)
library(cowplot)

# set working directory to the code directory
setwd("/Users/mhswaney/GitHub/scripts/mhswaney/cobamide_paper/CobamidePaperCode/")

#read in the data
pangenome <- read.delim("5_Corynebacterium_comparative_genomics/data/Corynebacterium_Pangenome_Summary.txt")
enrich <- read.delim("5_Corynebacterium_comparative_genomics/data/Corynebacterium_Enriched_Functions.txt")

# Figure 6B-E plotting -------------------------------------------------------

### Box plot of the total number of gene glusters in each genome
p2 <- ggplot(pangenome, aes(x=Ecosystem, y=num_gene_clusters, fill = Ecosystem)) + 
  geom_boxplot() + 
  geom_point() + 
  theme_classic(base_size = 12) + 
  ylab("Number of gene clusters per genome") + 
  theme(legend.position = "none")

# stats between host vs environment total num gene clusters
wilcox.test(num_gene_clusters ~ Ecosystem, data = pangenome) # p = 1.767e-07

# average num gene clusters by Ecosystem
pangenome %>% group_by(Ecosystem) %>% summarise(median(num_gene_clusters))

### Box plot of total genome length

pangenome$Mbp <- pangenome$total_length/1000000

p1 <- ggplot(pangenome, aes(x=Ecosystem, y=Mbp, fill = Ecosystem)) + 
  geom_boxplot() + 
  geom_point() + 
  theme_classic(base_size = 12) + 
  ylab("Genome length (Mbp)") + 
  theme(legend.position = "none")

# stats between host vs environment genome length
wilcox.test(Mbp ~ Ecosystem, data = pangenome) # p = 3.885e-06

# average genome length by ecosystem
pangenome %>% group_by(Ecosystem) %>% summarise(median(total_length/1E6))

### Enriched functions by ecosystem association

host <- enrich %>% filter(associated_groups=="Host" & adjusted_q_value<=0.05) # get host-associated significant COGs

host$percent_Host <- host$N_Host*host$p_Host/(host$N_Host+host$N_Environment) # percent of all genomes
host$percent_Environment <- host$N_Environment*host$p_Environment/(host$N_Host+host$N_Environment) # percent of all genomes

cog <- host[,c("COG_FUNCTION","adjusted_q_value","percent_Host","percent_Environment","enrichment_score")]
colnames(cog) <- c("COG Function","Adjusted q value","Host","Environment","Enrichment score")

cog <- cog %>% gather("Ecosystem","Percentage",3:4)

cog$`Adjusted q value` <- format(cog$`Adjusted q value`, scientific = FALSE, digits = 1)
cog$`Enrichment score` <- format(cog$`Enrichment score`, scientific = FALSE, digits = 3)
cog$`COG Function`<- factor(cog$`COG Function`, levels=rev(cog$`COG Function`[1:19])) # order by significance 

### Bar plot of host-associated significantly enriched COGs

p3 <- ggplot(cog,aes(x=`COG Function`,y=Percentage,fill=Ecosystem)) + 
  geom_col() +
  theme_classic() + 
  ylab("Proportion of genomes") + 
  xlab("Significantly enriched COG Functions in host-associated genomes (q<0.05)") +
  ylim(0,0.75) + 
  theme(axis.text.x = element_text(angle=40, hjust=1),
      text = element_text(size=8),
      legend.position="none") 

environment <- enrich %>% filter(associated_groups=="Environment" & adjusted_q_value<=0.05) # get environment-associated significant COGs

environment$percent_Host <- environment$N_Host*environment$p_Host/(environment$N_Host+environment$N_Environment) # percent of all genomes
environment$percent_Environment <- environment$N_Environment*environment$p_Environment/(environment$N_Host+environment$N_Environment) # percent of all genomes

cog <- environment[,c("COG_FUNCTION","adjusted_q_value","percent_Host","percent_Environment","enrichment_score")]
colnames(cog) <- c("COG Function","Adjusted q value","Host","Environment","Enrichment score")

cog <- cog[order(cog$`Adjusted q value`),] # order by significance
cog <- cog[1:20,] # get top 20

cog <- cog %>% gather("Ecosystem","Percentage",3:4)

cog$`Adjusted q value` <- format(cog$`Adjusted q value`, scientific = FALSE, digits = 1)
cog$`Enrichment score` <- format(cog$`Enrichment score`, scientific = FALSE, digits = 3)
cog$`COG Function`<- factor(cog$`COG Function`, levels=rev(cog$`COG Function`[1:20])) # order by significance 

p4 <- ggplot(cog,aes(x=`COG Function`,y=Percentage,fill=Ecosystem)) + 
  geom_col() +
  theme_classic() + 
  ylab("Proportion of genomes") + 
  xlab("Significantly enriched COG Functions in environment-associated genomes (q<0.05)") +
  ylim(0,0.75) + 
  theme(axis.text.x = element_text(angle=40, hjust=1),
      text = element_text(size=8),
      legend.position="none")  

top <- plot_grid(p1,p2, rel_widths = c(1,1.02))

plot_grid(top,p3, p4,nrow=3, rel_heights = c(0.7,1,1.2))


# Supplemental Figure 12 plotting ------------------------------------------

Strain <- lapply(pangenome$Strain,FUN=function(X){str_replace(X,"strain","")}) # remove 'strain' string from strain names
Strain <- t(data.frame(Strain))
pangenome$Strain <- as.vector(Strain)

# relabel C. matruchotii with correct strain ID
pangenome[pangenome$Strain == "Corynebacterium_matruchotii",][,"Strain"] <- "Corynebacterium_matruchotii_NCTC10206"

# plot dot plot - number of singleton gene clusters per genome - Supplemental Figure 12A
ggplot(pangenome, aes(x=reorder(Strain, singleton_gene_clusters),y=singleton_gene_clusters, color=Ecosystem)) + 
  geom_point(size=3) + 
  geom_segment(aes(x=Strain, 
                   xend=Strain,
                   y=min(singleton_gene_clusters),
                   yend= max(singleton_gene_clusters)),
                linetype="longdash",
               color="grey",
               size=0.25) + 
  coord_flip() + 
  theme(text=element_text(size=6)) +
  xlab("Strain") + 
  ylab("Number of singleton gene clusters") +
  theme_classic(base_size = 9) + 
  labs(color="Ecosystem")

# plot box plot - number of singleton gene clusters per genome - Supplemental Figure 12B
ggplot(pangenome, aes(x=Ecosystem, y=singleton_gene_clusters, fill=Ecosystem)) + 
  geom_boxplot()+
  theme_classic() + 
  xlab("Ecosystem") + 
  ylab("Number of singleton gene clusters") + 
  geom_point()

# statistical significance
wilcox.test(singleton_gene_clusters ~ Ecosystem, data = pangenome) # p = 0.9147


