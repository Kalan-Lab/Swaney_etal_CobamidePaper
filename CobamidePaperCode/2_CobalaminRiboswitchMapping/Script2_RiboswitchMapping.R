
# Script 2 ----------------------------------------------------------------
# Figure 3B-G

# Cobamide Skin Microbiome Manuscript
# Riboswitch Mapping
# Written by Mary Hannah Swaney
# University of Wisconsin - Madison

# Load required packages and data -----------------------------------------

library(readxl)
library(tidyverse)
library(gggenes)
library(rtracklayer)
library(GenomicAlignments)
library(ggbio)

# set working directory to the code directory
setwd("/Users/mhswaney/GitHub/scripts/mhswaney/cobamide_paper/CobamidePaperCode/")

# genes surrounding riboswitches in the six analyzed genomes
genes_near_riboswitches <- read_excel("data/SupplementalMaterial.xlsx",sheet="S4",skip=1)


# Cutibacterium acnes -----------------------------------------------------

Cacnes <- genes_near_riboswitches %>% filter(Species == "Cutibacterium acnes")

# plot genes - Figure 3B
ggplot(Cacnes, aes(xmin=start, xmax=end, y=region, forward=strand,fill=Category)) +
  theme_genes() + 
  geom_gene_arrow() + 
  facet_wrap(~ region, scales="free",ncol=1) + 
  scale_fill_manual(values=c("#3B76E3","#70D5E5","#6cc0a6","#FAA088","#EB4962","#BF19B7","darkgray","white","black"))

ca <- import.gff3("2_CobalaminRiboswitchMapping/data/genomes/Cutibacterium_acnes_KPA171202.gff")

# riboswitch sequences aligned to C. acnes genome with bowtie2 
bam <- GenomicAlignments::readGAlignments("2_CobalaminRiboswitchMapping/data/riboswitch_alignments/riboswitch_C_acnes.sorted.bam")

alignments <- autoplot(bam, color="black")
genome <- autoplot(ca[1], color = "#CECECE", fill = "#CECECE")

# plot aligned riboswitches along genome - Figure 3B
tracks(alignments,genome, label.bg.color = "white", track.bg.color = "white",heights=c(1, 0.25)) +
  theme_classic() + theme(axis.line = element_blank())


# Veillonella parvula -----------------------------------------------------

Vpar <- genes_near_riboswitches %>% filter(Species == "Veillonella_parvula")

# plot genes - Figure 3C
ggplot(Vpar, aes(xmin=start, xmax=end, y=region, forward=strand, fill=Category)) +
  theme_genes() + 
  geom_gene_arrow() + 
  facet_wrap(~ region, scales="free",ncol=1) + 
  scale_fill_manual(values=c("#3B76E3", "#FAA088","darkgray","white"))

vp <- import.gff3("2_CobalaminRiboswitchMapping/data/genomes/Veillonella_parvula_DSM2008.gff")

# riboswitch sequences aligned to V. parvula genome with bowtie2 
bam <- GenomicAlignments::readGAlignments("2_CobalaminRiboswitchMapping/data/riboswitch_alignments/riboswitch_V_parvula.sorted.bam")

alignments <- autoplot(bam, color="black")
genome <- autoplot(vp[1], color = "#CECECE", fill = "#CECECE")

# plot aligned riboswitches along genome - Figure 3C
tracks(alignments,genome, label.bg.color = "white", track.bg.color = "white",heights=c(1, 0.25)) +
  theme_classic() + theme(axis.line = element_blank())


# Pseudomonas putida ------------------------------------------------------

Pput <- genes_near_riboswitches %>% filter(Species == "Pseudomonas_putida")

# plot genes - Figure 3D
ggplot(Pput, aes(xmin=start, xmax=end, y=region, forward=strand,fill=Category)) +
  theme_genes() + 
  geom_gene_arrow() + 
  facet_wrap(~ region, scales="free",ncol=1) + 
  scale_fill_manual(values=c("#3B76E3","#70D5E5","#6CC0A6","white"))

ppu <- import.gff3("2_CobalaminRiboswitchMapping/data/genomes/Pseudomonas_putida_KT2440.gff")

# riboswitch sequences aligned to P. putida genome with bowtie2 
bam <- GenomicAlignments::readGAlignments("2_CobalaminRiboswitchMapping/data/riboswitch_alignments/riboswitch_P_putida.bam")

alignments <- autoplot(bam, color="black")
genome <- autoplot(ppu[1], color = "#CECECE", fill = "#CECECE")

# plot aligned riboswitches along genome - Figure 3D
tracks(alignments,genome, label.bg.color = "white", track.bg.color = "white",heights=c(1, 0.25)) +
  theme_classic() + theme(axis.line = element_blank())


# Corynebacterium kroppenstedtii ------------------------------------------

Ckrop <- genes_near_riboswitches %>% filter(Species == "Corynebacterium_kroppenstedtii")

# plot genes - Figure 3E
ggplot(Ckrop, aes(xmin=start, xmax=end, y=region, forward=strand,fill=Category)) +
  theme_genes() + 
  geom_gene_arrow() + 
  facet_wrap(~ region, scales="free",ncol=1) + 
  scale_fill_manual(values=c("#3B76E3","#70D5E5","#F4E3BC","white"))

ck <- import.gff3("2_CobalaminRiboswitchMapping/data/genomes/Corynebacterium_kroppenstedtii_DSM44385.gff")

# riboswitch sequences aligned to C. kroppenstedtii genome with bowtie2 
bam <- GenomicAlignments::readGAlignments("2_CobalaminRiboswitchMapping/data/riboswitch_alignments/riboswitch_C_kroppenstedtii.sorted.bam")

alignments <- autoplot(bam, color="black")
genome <- autoplot(ck[1], color = "#CECECE", fill = "#CECECE")

# plot aligned riboswitches along genome - Figure 3E
tracks(alignments,genome, label.bg.color = "white", track.bg.color = "white",heights=c(1, 0.25)) +
  theme_classic() + theme(axis.line = element_blank())


# Corynebacterium amycolatum ----------------------------------------------

Camyc <- genes_near_riboswitches %>% filter(Species == "Corynebacterium_amycolatum")

# plot genes - Figure 3F
ggplot(Camyc, aes(xmin=start, xmax=end, y=region, forward=strand,fill=Category)) +
  theme_genes() + 
  geom_gene_arrow() + 
  facet_wrap(~ region, scales="free",ncol=1) + 
  scale_fill_manual(values=c("#3b76e3","#70D5E5","#6CC0A6","#F4E3BC", "white"))

cam <- import.gff3("2_CobalaminRiboswitchMapping/data/genomes/Corynebacterium_amycolatum_FDAARGOS1107.gff")

# riboswitch sequences aligned to C. amycolatum genome with bowtie2 
bam <- GenomicAlignments::readGAlignments("2_CobalaminRiboswitchMapping/data/riboswitch_alignments/riboswitch_C_amycolatum.sorted.bam")

alignments <- autoplot(bam, color="black")
genome <- autoplot(cam[1], color = "#CECECE", fill = "#CECECE")

# plot aligned riboswitches along genome - Figure 3F
tracks(alignments,genome, label.bg.color = "white", track.bg.color = "white",heights=c(1, 0.25)) +
  theme_classic() + theme(axis.line = element_blank())


# Streptococcus sanguinis -------------------------------------------------

Ssang <- genes_near_riboswitches %>% filter(Species == "Streptococcus_sanguinis")

# plot genes - Figure 3G
ggplot(Ssang, aes(xmin=start, xmax=end, y=region, forward=strand,fill=Category)) +
  theme_genes() + 
  geom_gene_arrow() + 
  facet_wrap(~ region, scales="free",ncol=1) + 
  scale_fill_manual(values=c("#6CC0A6","#C6DDA8","darkgray","white"))

ss <- import.gff3("2_CobalaminRiboswitchMapping/data/genomes/Streptococcus_sanguinis_SK36.gff")

# riboswitch sequences aligned to S. sanguinis genome with bowtie2 
bam <- GenomicAlignments::readGAlignments("2_CobalaminRiboswitchMapping/data/riboswitch_alignments/riboswitch_S_sanguinis.sorted.bam")

alignments <- autoplot(bam, color="black")
genome <- autoplot(ss[1], color = "#CECECE", fill = "#CECECE")

# plot aligned riboswitches along genome - Figure 3G
tracks(alignments,genome, label.bg.color = "white", track.bg.color = "white",heights=c(1, 0.25)) +
  theme_classic() + theme(axis.line = element_blank())
