
# Script 3E ----------------------------------------------------------------
# Figure 4, Supplemental Table 2, and Supplemental Figure 6

# Cobamide Skin Microbiome Manuscript
# Spiec-Easi - Consensus Network Analysis
# Written by Mary Hannah Swaney
# University of Wisconsin - Madison

# Load required packages and data -----------------------------------------

library(readr)
library(tidyverse)
library(igraph)
library(SpiecEasi)
library(Matrix)
library(data.table)
library(ggplot2)
library(cowplot)

# set working directory to the code directory
setwd("/Users/mhswaney/GitHub/scripts/mhswaney/cobamide_paper/CobamidePaperCode/")

# load the 12 networks - 4 microenvironments, 3 studies
load("3_SpiecEasi/data/LKMB002.se.mb.RData")
load("3_SpiecEasi/data/Oh.se.mb.RData")
load("3_SpiecEasi/data/Hannigan.se.mb.RData")

# load cobamide-related species info
cobamide_info <- read_csv("3_SpiecEasi/data/CobamideInfoForSpiecEasi.csv")
rownames(cobamide_info) <- cobamide_info$Species


# Create new triangle shape for igraph ------------------------------------

# triangle vertex shape
mytriangle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  
  symbols(x=coords[,1], y=coords[,2], bg=vertex.color,
          stars=cbind(vertex.size, vertex.size, vertex.size),
          add=TRUE, inches=FALSE)
}
# clips as a circle
add_shape("triangle", clip=shapes("circle")$clip,
          plot=mytriangle)



# Create consensus network for each microenvironment ----------------------

order <- c("LKMB002","Oh","Hannigan")

# join together the spiec easi objects from the different studies for each microenvironment
sebaceous.mb <- list(LKMB002_SE$sebaceous, Oh_SE$sebaceous, Hannigan_SE$sebaceous)
names(sebaceous.mb) <- order
moist.mb <- list(LKMB002_SE$moist, Oh_SE$moist, Hannigan_SE$moist)
names(moist.mb) <- order
dry.mb <- list(LKMB002_SE$dry, Oh_SE$dry, Hannigan_SE$dry)
names(dry.mb) <- order
foot.mb <- list(LKMB002_SE$foot, Oh_SE$foot, Hannigan_SE$foot)
names(foot.mb) <- order
# create one large list that contains all of the networks
networks <- list(sebaceous.mb, moist.mb, dry.mb, foot.mb)
names(networks) <- c("sebaceous","moist","dry","foot")

# function to create consensus networks
# adapted from: https://github.com/zdk123/SpiecEasiSLR_manuscript

consensusSE <- function(network,microenvironment_for_mean){
  
  # get optimal beta matrix 
  se.mb <- lapply(network,FUN=function(x){
    sebeta <- symBeta(getOptBeta(x), mode='maxabs')
    diag(sebeta) <- 0
    colnames(sebeta) <- rownames(sebeta) <- cobamide_info$Species
    sebeta})
  nets_li <- list("se.mb" =se.mb)
  
  # get edges and convert to direction of the sign (pos or neg)
  gr2edges <- function(x) {
    edges <- Matrix::summary(x)
    edges[,1] <- cobamide_info$Species[edges[,1]]
    edges[,2] <- cobamide_info$Species[edges[,2]]
    edges[,1:2] <- t(apply(edges[,1:2], 1, sort))
    edges <- as.data.frame(edges)
    edges$edge <- paste(edges[,1], edges[,2], sep=".")
    edges$sign <- ifelse(edges$x > 0, "pos","neg")
    tmp <- split(edges, edges$sign)
    names(tmp) <- ifelse(names(tmp)=="pos", "pos", "neg")
    tmp
  }
  
  edges_li <- lapply(nets_li, function(nets) purrr::transpose(lapply(nets, gr2edges)))
  
  
  pos_edges <- table(unlist(lapply(edges_li$se.mb$pos, function(x) x$edge)))
  neg_edges <- table(unlist(lapply(edges_li$se.mb$neg, function(x) x$edge)))
  
  ## build consensus weighted adjancency matrix
  taxa <- sort(unique(unname(
    c(unlist(lapply(edges_li$se.mb$pos, function(x) x[,1:2])),
      unlist(lapply(edges_li$se.mb$neg, function(x) x[,1:2])))
  )))
  
  rebuild_adj <- function(x) {
    if (is.null(x)) {
      return(Matrix::Matrix(data=0, ncol=length(taxa), nrow=length(taxa)))
    }
    
    sign(sparseMatrix(i=match(x$i, taxa),
                      j=match(x$j, taxa),
                      x=x$x, dims=c(length(taxa), length(taxa))))
  }
  
  edges_full <- list(pos=lapply(edges_li$se.mb$pos, rebuild_adj),
                     neg=lapply(edges_li$se.mb$neg, rebuild_adj))
  
  pos_ <- Reduce(`+`, x=edges_full$pos)/
    length(edges_full$pos)
  neg_ <- Reduce(`+`, x=edges_full$neg)/
    length(edges_full$neg)
  
  adj <- forceSymmetric(pos_ + neg_)
  rownames(adj) <- colnames(adj) <- taxa
  
  # exclude edges that are not present in at least 2 of 3 networks
  adj[abs(adj) < 2/3] <- 0
  # remove if degree is 0 
  deg <- Matrix::rowSums(adj)
  adj_fli <- adj[deg!=0, deg!=0]
  
  tmpadj <- forceSymmetric(adj_fli)
  
  # get taxa names/order
  subtax <- rownames(adj_fli)
  
  ### Add vertex metadata for igraph ###

  # mean category = group species by their mean relative abundance in each microenvironment
  cobamide_info$Mean_category <- ""
  
  # add the mean relative abundance for each species, depending on the microenvironment, to cobamide_info
  if(microenvironment_for_mean == "sebaceous"){
    cobamide_info$Mean_category <- ifelse(cobamide_info$Sebaceous_Mean >= 10, ">10%", cobamide_info$Mean_category)
    cobamide_info$Mean_category <- ifelse(cobamide_info$Sebaceous_Mean < 10 & cobamide_info$Sebaceous_Mean >=1, "1-10%", cobamide_info$Mean_category)
    cobamide_info$Mean_category <- ifelse(cobamide_info$Sebaceous_Mean < 1 & cobamide_info$Sebaceous_Mean >=0.1, "0.1-1%", cobamide_info$Mean_category)
    cobamide_info$Mean_category <- ifelse(cobamide_info$Sebaceous_Mean < 0.1, "<0.1%", cobamide_info$Mean_category)
  }
  if(microenvironment_for_mean == "moist"){
    cobamide_info$Mean_category <- ifelse(cobamide_info$Moist_Mean >= 10, ">10%", cobamide_info$Mean_category)
    cobamide_info$Mean_category <- ifelse(cobamide_info$Moist_Mean < 10 & cobamide_info$Moist_Mean >=1, "1-10%", cobamide_info$Mean_category)
    cobamide_info$Mean_category <- ifelse(cobamide_info$Moist_Mean < 1 & cobamide_info$Moist_Mean >=0.1, "0.1-1%", cobamide_info$Mean_category)
    cobamide_info$Mean_category <- ifelse(cobamide_info$Moist_Mean < 0.1, "<0.1%", cobamide_info$Mean_category)
  }
  if(microenvironment_for_mean == "dry"){
    cobamide_info$Mean_category <- ifelse(cobamide_info$Dry_Mean >= 10, ">10%", cobamide_info$Mean_category)
    cobamide_info$Mean_category <- ifelse(cobamide_info$Dry_Mean < 10 & cobamide_info$Dry_Mean >=1, "1-10%", cobamide_info$Mean_category)
    cobamide_info$Mean_category <- ifelse(cobamide_info$Dry_Mean < 1 & cobamide_info$Dry_Mean >=0.1, "0.1-1%", cobamide_info$Mean_category)
    cobamide_info$Mean_category <- ifelse(cobamide_info$Dry_Mean < 0.1, "<0.1%", cobamide_info$Mean_category)
  }
  if(microenvironment_for_mean == "foot"){
    cobamide_info$Mean_category <- ifelse(cobamide_info$Foot_Mean >= 10, ">10%", cobamide_info$Mean_category)
    cobamide_info$Mean_category <- ifelse(cobamide_info$Foot_Mean < 10 & cobamide_info$Foot_Mean >=1, "1-10%", cobamide_info$Mean_category)
    cobamide_info$Mean_category <- ifelse(cobamide_info$Foot_Mean < 1 & cobamide_info$Foot_Mean >=0.1, "0.1-1%", cobamide_info$Mean_category)
    cobamide_info$Mean_category <- ifelse(cobamide_info$Foot_Mean < 0.1, "<0.1%", cobamide_info$Mean_category)
  }
  
  merged <- left_join(data.frame("Species" = subtax),cobamide_info) # attach cobamide info to the species list
  
  # categories
  prod <- merged$BiosynthesisCategory # cobamide biosynthesis category
  B12user <- merged$B12_dependence # true/false b12 dependence
  depgenes <- merged$B12_dep_gene_count # number of B12-dependent genes
  phylum <- merged$P # phylum classification
  Mean_category <- merged$Mean_category # discrete mean relative abundances
  
  # generate igraph object
  spiec.graph <- SpiecEasi::adj2igraph(abs(tmpadj),
                                       rmEmptyNodes=TRUE,
                                       vertex.attr = list(name=subtax,B12producer=prod, B12user = B12user, depgenes = depgenes, phylum = phylum,
                                                          Mean_category = Mean_category)
  )
  
  ### Assign igraph attributes ###
  
  # shape - specifies biosynthesis category
  V(spiec.graph)$shape = "square"                             
  V(spiec.graph)[which(B12producer =="Non-producer" )]$shape <- "square"  
  V(spiec.graph)[which(B12producer =="Precursor salvager") ]$shape <- "triangle"  
  V(spiec.graph)[which(B12producer =="Producer") ]$shape <- "circle"  
  
  # color - specifies phylum
  V(spiec.graph)$color = "gray"    
  V(spiec.graph)[which(phylum == "Actinobacteria") ]$color <- "#d17584"  
  V(spiec.graph)[which(phylum =="Bacteroidetes") ]$color <- "#AA4499"  
  V(spiec.graph)[which(phylum =="Firmicutes" )]$color <- "#999933"
  V(spiec.graph)[which(phylum =="Proteobacteria") ]$color <- "#66ccbb"  
  V(spiec.graph)[which(phylum =="Basidiomycota" )]$color <- "#DDcc77"  
  
  # arbitrary numbers assigned to phyla - for statistics
  V(spiec.graph)$color_num = 0   
  V(spiec.graph)[which(phylum == "Actinobacteria") ]$color_num <- 1 
  V(spiec.graph)[which(phylum =="Bacteroidetes") ]$color_num <- 2 
  V(spiec.graph)[which(phylum =="Firmicutes" )]$color_num <- 3
  V(spiec.graph)[which(phylum =="Proteobacteria") ]$color_num <- 4 
  V(spiec.graph)[which(phylum =="Basidiomycota" )]$color_num <- 5  
  V(spiec.graph)[which(phylum =="Fusobacteria" )]$color_num <- 6
  
  # frame color - specifies cobamide dependence
  V(spiec.graph)$frame.color = NA                            
  V(spiec.graph)[which(B12user == FALSE )]$frame.color <- NA  
  V(spiec.graph)[which(B12user == TRUE )]$frame.color <- "black"  
  
  # size - specifies mean relative abundance category
  V(spiec.graph)$size = 0                            
  V(spiec.graph)[which(Mean_category == ">10%")]$size <-16
  V(spiec.graph)[which(Mean_category == "1-10%")]$size <- 12
  V(spiec.graph)[which(Mean_category == "0.1-1%")]$size <- 8
  V(spiec.graph)[which(Mean_category == "<0.1%")]$size <- 4
  
  # edges color - positive = green, negative = red
  # adapted from: http://psbweb05.psb.ugent.be/conet/microbialnetworks/solutions.php#spieceasiedge
  edges=E(spiec.graph)
  edge.colors=c()
  for(e.index in 1:length(edges)){
    adj.nodes=ends(spiec.graph,edges[e.index])
    xindex=which(taxa==adj.nodes[1])
    yindex=which(taxa==adj.nodes[2])
    beta=adj[xindex,yindex]
    if(beta>0){
      edge.colors=append(edge.colors,"forestgreen")
    }else if(beta<0){
      edge.colors=append(edge.colors,"#FF298F")
    }
  }
  E(spiec.graph)$color=edge.colors
  
  return(spiec.graph)
}


# Plot the consensus networks with igraph - Figure 4A -------------------------------------------

sebaceous <- consensusSE(networks$sebaceous,"sebaceous")
set.seed(99)
sb.coord <- layout.fruchterman.reingold(sebaceous)
#pdf("SpiecEasi/consensus_sebaceous.pdf", width=6, height=6)
plot(x=sebaceous, vertex.label=NA,edge.width=0.75, main = "Sebaceous")
dev.off()

# moist consensus network
moist <- consensusSE(networks$moist,"moist")
set.seed(99)
sb.coord <- layout.fruchterman.reingold(moist)
#pdf("SpiecEasi/consensus_moist.pdf", width=6, height=6)
plot(x=moist, vertex.label = NA, edge.width=0.75, main = "Moist")
dev.off()

# dry consensus network
dry <- consensusSE(networks$dry,"dry")
set.seed(99)
sb.coord <- layout.fruchterman.reingold(dry)
#pdf("SpiecEasi/consensus_dry.pdf", width=6, height=6)
plot(x=dry, vertex.label=NA, edge.width=0.75, main = "Dry")
dev.off()

# foot consensus network
foot <- consensusSE(networks$foot,"foot")
set.seed(99)
sb.coord <- layout.fruchterman.reingold(foot)
#pdf("SpiecEasi/consensus_foot.pdf", width=6, height=6)
plot(x=foot, vertex.label=NA, edge.width=0.75, main = "Foot")
dev.off()


# Compute network summary statistics ---------------------------------------------

# compute some standard network statistics
# https://github.com/ramellose/NetworkUtils/blob/master/R/computeStat.R
computeStat = function(graph, mode){
  score = assortativity.degree(graph)
  if (mode == "transitivity"){
    score = transitivity(graph, type=c("global")) # transitivity - probability that two adjacent nodes are connected - can reveal existence of tightly connected communities
    # also known as clustering coefficient
  }
  else if (mode == "degree"){
    score = igraph::degree(graph) # number of degrees for each node
  }
  else if (mode == "assortativity_phylum"){
    score = assortativity_nominal(graph,types=V(graph)$color_num) # assortativity by phylum - preference for nodes to connect to other nodes of the same phylum
  }
  else if (mode == "assortativity_degree"){
    score = assortativity_degree(graph) # assortativity by degree - preference for nodes to connect to other nodes with similar degrees
  }
  else if (mode == "modularity"){
    groups = igraph::cluster_walktrap(graph) # tries to find densely connected communities
    score = igraph::modularity(graph, membership(groups)) #  measures the strength of division of a network into modules 
  }
  
  return(score)
}

NetworkStatistics <- function(microenvironment){
  stat <- list(computeStat(microenvironment,"transitivity"),
               computeStat(microenvironment,"assortativity_phylum"),
               computeStat(microenvironment,"assortativity_degree"),
               computeStat(microenvironment,"modularity"))
  names(stat) <- c("Transitivity","Assortativity_phylum","Assortativity_degree","Modularity")
  return(stat)
}

# network statistics for each microenvironment network
sebaceous_stat <- NetworkStatistics(sebaceous)
moist_stat <- NetworkStatistics(moist)
dry_stat <- NetworkStatistics(dry)
foot_stat <- NetworkStatistics(foot)

Transitivity <- data.frame(sebaceous = sebaceous_stat$Transitivity, moist = moist_stat$Transitivity, dry = dry_stat$Transitivity, foot = foot_stat$Transitivity) %>%
  gather(microenvironment,Transitivity)
Assortativity_phylum <- data.frame(sebaceous = sebaceous_stat$Assortativity_phylum, moist = moist_stat$Assortativity_phylum, dry = dry_stat$Assortativity_phylum, foot = foot_stat$Assortativity_phylum) %>%
  gather(microenvironment,Assortativity_phylum)
Assortativity_degree <- data.frame(sebaceous = sebaceous_stat$Assortativity_degree, moist = moist_stat$Assortativity_degree, dry = dry_stat$Assortativity_degree, foot = foot_stat$Assortativity_degree) %>% 
  gather(microenvironment, Assortativity_degree)
Modularity <- data.frame(sebaceous = sebaceous_stat$Modularity, moist = moist_stat$Modularity, dry = dry_stat$Modularity, foot = foot_stat$Modularity) %>%
  gather(microenvironment, Modularity)

# join each microenvironment together
summary_stat <- left_join(Transitivity, Assortativity_phylum)
summary_stat <- left_join(summary_stat, Assortativity_degree)
summary_stat <- left_join(summary_stat, Modularity)

# for Supplemental Table 1
summary_stat <- summary_stat %>% gather(statistic, value, 2:5)

# count number of edges and nodes for each network and calculate density
edge_vertex <- list("sebaceous" = data.frame(edges = gsize(sebaceous), vertex=gorder(sebaceous), density = edge_density(sebaceous), microenvironment="sebaceous"),
                    "moist" = data.frame(edges = gsize(moist), vertex=gorder(moist), density = edge_density(moist), microenvironment="moist"),
                    "dry" = data.frame(edges = gsize(dry), vertex=gorder(dry), density = edge_density(dry), microenvironment="dry"),
                    "foot" = data.frame(edges = gsize(foot), vertex=gorder(foot), density = edge_density(foot), microenvironment="foot"))

# for Supplemental Table 1
edge_vertex <- rbindlist(edge_vertex)

# Integration of cobamide biosynthesis and use ----------------------------

# count the edges between species of each cobamide biosynthesis category
edgesQuantify <- function(igraph){
  
  edges_quantify <- as.data.frame(get.edgelist(igraph))
  
  
  V1 <- data.frame(Species = edges_quantify$V1)
  V2 <- data.frame(Species = edges_quantify$V2)
  
  V1 <- left_join(V1,cobamide_info[,c("Species","BiosynthesisCategory")])
  V2 <- left_join(V2,cobamide_info[,c("Species","BiosynthesisCategory")])
  
  V <- data.frame(V1$BiosynthesisCategory,V2$BiosynthesisCategory)
  
  V$combined <- with(V, paste0(V1.BiosynthesisCategory, "--",V2.BiosynthesisCategory))
  
  a <- V %>% group_by(combined) %>% summarise(n())
  
  is.integer0 <- function(x)
  {
    is.integer(x) && length(x) == 0L
  }
  count <- data.frame(0)
  count$salv_prod <- ifelse("Precursor salvager--Producer" %in% a$combined & "Producer--Precursor salvager" %in% a$combined, a[which(a$combined == "Precursor salvager--Producer"),2][[1]] + a[which(a$combined == "Producer--Precursor salvager"),2][[1]],a[which(a$combined == "Producer--Precursor salvager"),2][[1]]) 
  count$prod_prod <- a[which(a$combined == "Producer--Producer"),2][[1]]
  count$non_prod <- a[which(a$combined == "Non-producer--Producer"),2][[1]] + a[which(a$combined == "Producer--Non-producer"),2][[1]]
  count$non_non <- a[which(a$combined == "Non-producer--Non-producer"),2][[1]]
  count$salv_non <- a[which(a$combined == "Non-producer--Precursor salvager"),2][[1]] + a[which(a$combined == "Precursor salvager--Non-producer"),2][[1]]
  count <- count[,2:ncol(count)]
  return(as.data.frame(count))
  
}

igraphs <- list(sebaceous,moist,dry,foot)
counts <- lapply(igraphs, FUN=edgesQuantify) # count edges function
names(counts) <- c("sebaceous","moist","dry","foot")
final <- data.frame(rbindlist(counts))
final$microenvironment <- c("sebaceous","moist","dry","foot")
colnames(final) <- c("Producer-Salvager","Producer-Producer", "Producer-Non-producer","Non-producer-Non-producer","Salvager-Non-producer","microenvironment" )
gathered <- final %>% gather(A,B, 1:5)

gathered$microenvironment <- factor(gathered$microenvironment, levels = c("sebaceous","moist","dry","foot"))

edge_categories <- left_join(gathered, edge_vertex)
edge_categories$normalized <- edge_categories$B / edge_categories$edges
edge_categories$microenvironment <- factor(edge_categories$microenvironment, levels = c("sebaceous","moist","dry","foot"))

# plot edge categories
p1 <- ggplot(edge_categories, aes(x=A, y=normalized * 100, fill=microenvironment)) + geom_col(position="dodge2",width = 0.9) + 
  theme_classic() + xlab("Edge category") + theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  ylab("Percentage of total edge count (%)") + scale_fill_manual(values=c("#882255","#CC6677","#DDCC77","#44AA99"))

average_edge <- gathered %>% group_by(A) %>% summarise(mean = mean(B), sd = sd(B)) 

# count the edges between de novo producers and non-producers/salvagers with cobamide dependence
edgesQuantifyDependence <- function(igraph){
  
  edges_quantify <- as.data.frame(get.edgelist(igraph))
  
  V1 <- data.frame(Species = edges_quantify$V1)
  V2 <- data.frame(Species = edges_quantify$V2)
  
  V1 <- left_join(V1,cobamide_info[,c("Species","BiosynthesisCategory","B12_dependence")])
  V2 <- left_join(V2,cobamide_info[,c("Species","BiosynthesisCategory","B12_dependence")])
  
  V <- data.frame(V1$Species, V1$B12_dependence,V1$BiosynthesisCategory, V2$Species, V2$B12_dependence, V2$BiosynthesisCategory)
  
  V <- V %>% filter((V1.B12_dependence == TRUE & 
                       (V1.BiosynthesisCategory == "Non-producer" | V1.BiosynthesisCategory == "Precursor salvager") & V2.BiosynthesisCategory == "Producer") |
                      (V2.B12_dependence == TRUE & (V2.BiosynthesisCategory == "Non-producer" | V1.BiosynthesisCategory == "Precursor salvager") & V1.BiosynthesisCategory == "Producer"))
  
  
  producer_to_user <- nrow(V)
  
  return(producer_to_user)
}

igraphs <- list(sebaceous,moist,dry,foot)
counts <- lapply(igraphs, FUN=edgesQuantifyDependence) # count edges between producer and non-producer/dependent function
names(counts) <- c("sebaceous","moist","dry","foot")
final <- data.frame(unlist(counts))
final$microenvironment <- rownames(final)
colnames(final) <- c("Producer to User", "microenvironment")

prod_to_user <- left_join(final, edge_vertex)
prod_to_user$normalized <- prod_to_user$`Producer to User` / prod_to_user$edges
prod_to_user$Category <- "Producer-to-User"
prod_to_user$microenvironment <- factor(prod_to_user$microenvironment, levels = c("sebaceous","moist","dry","foot"))

# plot number of producer to non-producer/dependent edges
p2 <- ggplot(prod_to_user, aes(x=Category, y=normalized*100, fill=microenvironment)) + geom_col(position="dodge2", width=0.375) + 
  theme_classic() + xlab("Edge category") + theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  ylab("Percentage of total edge count (%)") + scale_fill_manual(values=c("#882255","#CC6677","#DDCC77","#44AA99"))

# get species in each network and attach cobamide dependence/biosynthesis info
UniqueTaxa <- function(igraph){
  
  edges_quantify <- as.data.frame(get.edgelist(igraph))
  
  V1 <- data.frame(Species = edges_quantify$V1)
  V2 <- data.frame(Species = edges_quantify$V2)
  both <- rbind(V1,V2)
  unique(both$Species)
}

igraphs <- list(sebaceous,moist,dry,foot)
counts <- lapply(igraphs, FUN=UniqueTaxa) # species present in each microenvironment network
names(counts) <- c("sebaceous","moist","dry","foot")
counts <- lapply(counts, FUN=function(x)
  left_join(data.frame(Species = x),cobamide_info)
)
counts[[1]]$microenvironment <- "sebaceous"
counts[[2]]$microenvironment <- "moist"
counts[[3]]$microenvironment <- "dry"
counts[[4]]$microenvironment <- "foot"
species_plot <- rbindlist(counts)

species_plot$B12_dependence <- factor(species_plot$B12_dependence, levels =c(TRUE,FALSE))
species_plot$BiosynthesisCategory <- factor(species_plot$BiosynthesisCategory, levels =c("Producer","Non-producer","Precursor salvager"))
species_plot$microenvironment <- factor(species_plot$microenvironment, levels = c("sebaceous","moist","dry","foot"))

# plot number of species withcobamide dependence
p3 <- ggplot(species_plot, aes(x=B12_dependence, fill=microenvironment)) + geom_bar(position="dodge2",width=0.4) +theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
  xlab("Cobamide dependence") + ylab("Number of species") +
  scale_fill_manual(values=c("#882255","#CC6677","#DDCC77","#44AA99"))

# plot number of species within each biosynthesis category 
p4 <- ggplot(species_plot, aes(x=BiosynthesisCategory, fill=microenvironment)) + geom_bar(position="dodge2",width=0.6) +theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  xlab("Cobamide biosynthesis category") + ylab("Number of species") +
  scale_fill_manual(values=c("#882255","#CC6677","#DDCC77","#44AA99")) + 
  labs(fill= "Microenvironment")

legend <- get_legend(
  # create some space to the left of the legend
  p4 + guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "right")
)


p <- plot_grid(
  p4 + theme(legend.position="none"),
  p3 + theme(legend.position="none"),
  p1 + theme(legend.position="none"),
  p2,align = "hv"
)


# Figure 4B-E
plot_grid(p, ncol = 1)

# average number of edge categories
average_category <- species_plot %>% group_by(microenvironment,BiosynthesisCategory) %>% summarise(count = n()) %>%
  group_by(BiosynthesisCategory) %>% summarise(mean=mean(count), sd = sd(count))

# average number of dependent/not dependent species
average_dependence <- species_plot %>% group_by(microenvironment,B12_dependence) %>% summarise(count = n()) %>%
  group_by(B12_dependence) %>% summarise(mean=mean(count), sd = sd(count))


# Degree distribution -----------------------------------------------------

# plot degree distribution - Supplemental Figure 6
dev.off()
plot(0:(length(degree.distribution(sebaceous))-1), degree.distribution(sebaceous), ylim=c(0,.4), xlim=c(0,11), type='b',
     ylab="Frequency", xlab="Degree", col = "#882255",lwd=2)
lines(0:(length(degree.distribution(moist))-1), degree.distribution(moist),col="#CC6677",type='b', lwd=2)
lines(0:(length(degree.distribution(dry))-1), degree.distribution(dry),col="#DDCC77",type='b', lwd=2)
lines(0:(length(degree.distribution(foot))-1), degree.distribution(foot),col="#44AA99",type='b', lwd=2)
legend("topright", legend=c("Sebaceous","Moist","Dry","Foot"),
       col=c("#882255","#CC6677","#DDCC77","#44AA99"),
       lty=1, cex=1,lwd = 2)

# count number of degrees for each node (species) in each network
sebaceous_degree <- data.frame(degree_count = computeStat(sebaceous,"degree"), microenvironment = "sebaceous")
moist_degree <- data.frame(degree_count = computeStat(moist,"degree"), microenvironment = "moist")
dry_degree <- data.frame(degree_count = computeStat(dry,"degree"), microenvironment = "dry")
foot_degree <- data.frame(degree_count = computeStat(foot, "degree"), microenvironment = "foot")

degree <- rbind(sebaceous_degree,moist_degree)
degree <- rbind(degree, dry_degree)
degree <- rbind(degree, foot_degree)

degree %>% group_by(microenvironment) %>% summarise(mean = mean(degree_count), sd = sd(degree_count))

