# Title     : Main file for data processing
# Objective : Create a workflow for the data processing
# Created by: metz.seba91@gmail.com
# Created on: 18/03/2021


# Libraries ----

library(stringi)
library(stringr)
library(pr2database)
library(ggplot2)
library(treemapify)
library(Biostrings)
library(dplyr)
library(RColorBrewer)

# Functions ----
source("scripts/Functions/ReadStampa.R")
source("scripts/Functions/GetNetwork.R")
source("scripts/Functions/GetSeqs.R")
source("scripts/Functions/StampaPlots.R")
source("scripts/Functions/ColorPallets.R")

#path to folder to save results
path <- file.path("02-results/")

# Load inputs data ----
 # - ASVs table (rownames = ASVs, colnames = Samples)
 # - Taxonomy - STAMPA results: taxonomy(ASV  abundance identity  Taxonomy[separated by "|"]  IDs)
 # - CIM and CEM results - STAMPA results
 # - All to all - VSEARCH results = query+target+id

ASVTABLE <- readRDS("01-data/table_formated.rds")

TAX <- readStampa(readfile="01-data/stampa_ASVs/stampa.results", 
                  taxnames = c("S", "D", "C", "O", "F", "G", "Sp"), 
                  taxcorrect = T)

CIM <- readStampa(readfile="01-data/stampa_CI/stampa.results", 
                  taxnames = c("S", "D", "C", "O", "F", "G", "Sp"), 
                  taxcorrect = T)

CEM <- readStampa(readfile="01-data/stampa_ENV/stampa.results", 
                 taxnames = c("S", "D", "C", "O", "F", "G", "Sp"), 
                 taxcorrect = T)

ALLTOALL <-read.table("01-data/allpairs_vsearch.cl", sep = "\t", 
                      stringsAsFactors = F, 
                      col.names = c("ASV1", "ASV2", "identity"))

ALLTOALL$ASV1 <- stri_replace_all_regex(ALLTOALL$ASV1, "(ASV_[0-9]+)(_[0-9]+)", "$1")
ALLTOALL$ASV2 <- stri_replace_all_regex(ALLTOALL$ASV2, "(ASV_[0-9]+)(_[0-9]+)", "$1")

FASTA <- readBStringSet("01-data/ASVs.fasta", format = "fasta") # Read fasta for phylogeny
names(FASTA) <- unlist(str_extract_all(names(FASTA), "ASV_[0-9]+"))

################################################################################
##### Table curation ----

TAX <- TAX %>% dplyr::filter(!D %in% c("Metazoa", "Streptophyta"), 
                             !ASV %in% c("ASV_3130", "ASV_3465", "ASV_439",
                                         "ASV_579", "ASV_751", "ASV_722", "ASV_804",
                                         "ASV_840", "ASV_1252", "ASV_1288", "ASV_1493",
                                         "ASV_1562", "ASV_1644", "ASV_1711","ASV_1798",
                                         "ASV_1815", "ASV_1942", "ASV_2207", "ASV_2227",
                                         "ASV_2314", "ASV_2349", "ASV_2350", "ASV_2418",
                                         "ASV_2420", "ASV_2484", "ASV_2549", "ASV_2582",
                                         "ASV_2637", "ASV_2647", "ASV_2739", "ASV_2778",
                                         "ASV_2804", "ASV_2832", "ASV_2840", "ASV_2854",
                                         "ASV_2893", "ASV_3024", "ASV_3097", "ASV_3123",
                                         "ASV_3127", "ASV_3146", "ASV_3178",
                                         "ASV_3200", "ASV_3214", "ASV_3376", "ASV_3386",
                                         "ASV_3465", "ASV_3467", "ASV_3471", "ASV_3472",
                                         "ASV_3492", "ASV_3565", "ASV_3308") ) # ASV_840, ASV_1711, ASV_1942 XCELLIDAE

#Manually checked
TAX$C[grep(T, TAX$ASV %in% c("ASV_1191", "ASV_1821", "ASV_2530"))] <- "Ulvophyceae"
TAX$C[grep(T, TAX$ASV == "ASV_1351")] <- "Trebouxiophyceae"
TAX$C[grep(T, TAX$ASV == "ASV_1757")] <- "Perkinsida"
TAX$D[grep(T, TAX$ASV == "ASV_1757")] <- "Perkinsea"
TAX$C[grep(T, TAX$ASV %in% c("ASV_2004", "ASV_3136"))] <- "Pedinophyceae"
TAX$C[grep(T, TAX$ASV %in% c("ASV_3324", "ASV_3339", "ASV_3424"))] <- "Chlorophyceae"
TAX$C[grep(T, TAX$ASV == "ASV_3364")] <- "Dinophyceae"
TAX$D[grep(T, TAX$ASV == "ASV_3364")] <- "Dinoflagellata"

###
# Filtering STAMPA DATA
STAMPA <- read.table("01-data/stampa_ASVs/stampa.results", sep = "\t")
STAMPA <- STAMPA %>% dplyr::filter(V1 %in% TAX$ASV)
write.table(STAMPA, "02-results/STAMPA.FILTER.tsv", sep = "\t", row.names = F, quote = F)

################################################################################
###### Exploration analysis------
    # N of reads, ASVs
    # ASVs Contribution
    # Stampa plot for more abundant groups

aux <- ASVTABLE[grep(T, rownames(ASVTABLE) %in% TAX$ASV), ]
sum(aux) # Reads 1227460
colSums(aux) # Reads per samples

# P1_2014  P1_2016  P2_2014  P2_2016  P3_2015  P3_2016  P5_2014  P5_2015 P12_2015 C16_2015 
# 121534   238547    42032    53894    59465   212011    13930   173981   267674    44392
# 3167 ASVs

# TreeMAP --- FUNCTION
pdf(file.path(path, "treemap.pdf"), width = 8.5, height = 7)
  teemapPlot(Stampa = TAX, level1 = "S", level2 = "C")
dev.off()

# Stampa plot
tax.aux <- TAX  %>% dplyr::filter(S %in% c("Alveolata", "Archaeplastida", "Stramenopiles", 
                                          "Opisthokonta", "Cryptista", "Rhizaria"))
tax.aux <- data.frame(taxonomy =  tax.aux$S, pident = tax.aux$pident)

pdf(file.path(path, "densityPlot.pdf"), width = 6.5, height = 5)
  print(densityplot_stampa(tax.aux))
dev.off()

###### Novelty analysis ----
  # Data preparation
  # Global network analysis
  # Novelty barplot
  # 2D Novelty network plot for each Supergroup

### Data preparation

# merge tables
aux <- merge(CIM[,c(1,3)], CEM[,c(1,3)], by = "ASV")
colnames(aux) <- c("ASV", "CIM", "CEM")

# create nodes
nodes <- merge(TAX, aux, all.x = T)
colnames(nodes)[1]<-"name"

# create relationships
ALLTOALL <- ALLTOALL %>% dplyr::filter(ASV1 %in% TAX$ASV &  ASV2 %in% TAX$ASV)
relations <- data.frame(from = ALLTOALL$ASV1, to = ALLTOALL$ASV2, 
                        identity = ALLTOALL$identity)

### Global network analysis
#Complete network for Supergroup novelty comparison
network <- get_network(edges = 95, # identity value
                       NODES = nodes, RELATIONS = relations)

#plot network at specific taxonomic level (color level)
mcim <- mean(network$NODES$CIM)
mcem <- mean(network$NODES$CEM)

aux <- prepare_for_plot(network, color_level = "S")

pdf(file.path(path, "Supergroup_novelty.pdf"), width = 5, height = 5)
  plot_network(aux, plot_lines = T, plot_mem = F, mcim =  mean(network$NODES$CIM), 
               mcem = mean(network$NODES$CEM))
dev.off()

### Novelty barplot
#plot barplots with % of ASVs in each quadrant accordint to the mean of the total datasets
aux2 <- network$NODE
# labels
aux2$quadrant <- unlist(lapply(1:nrow(aux2), function(x)
  ifelse(aux2$CIM[x] <= mean(aux2$CIM) & aux2$CEM[x] <= mean(aux2$CEM), "3th",
         ifelse(aux2$CIM[x] <= mean(aux2$CIM) & aux2$CEM[x] > mean(aux2$CEM), "4th",
                ifelse(aux2$CIM[x] > mean(aux2$CIM) & aux2$CEM[x] > mean(aux2$CEM), "1st", "2nd")))))


#plot supergroups
table_n <- aux2 %>% dplyr::group_by(S, quadrant) %>% summarise(N = n()) %>%
  mutate(percent = (N / sum(N)), cumsum = cumsum(percent), Size=ifelse(
    S %in% c("Amoebozoa", "Apusozoa") & quadrant == "3th", paste0("N=",sum(N)),
    ifelse(quadrant == "4th", paste0("N=",sum(N)), "")), total = sum(N))
colnames(table_n) <- c("Tax", "Quadrant", "N", "per", "cumssum", "Size", "Total")


#plot
g <- ggplot(table_n, aes(fill=Quadrant, y=per, x= reorder(Tax, -Total))) + 
  scale_y_continuous(labels = scales::percent) +
  labs ( y = "Percentage of ASVs", x = "Supergroup") +
  geom_bar(stat="identity", color="black") +
  geom_text(aes(y=cumssum, label=Size), vjust=-0.8) +
  guides(fill=guide_legend(reverse=T))+
  scale_fill_manual(values = c("#88e0e3", "#81e278", "#ff6363", "#b1aaaa"), name = "Quadrant")+
  theme_minimal()+
  theme(
    text = element_text(family="serif", colour = "black", size = 12),
    axis.text.x = element_text(angle = 45, colour = "black"),
    axis.text.y = element_text(colour = "black"))

pdf(file.path(path,"barplot_supergroups.pdf"), width = 8, height = 6)
  print(g)
dev.off()


### 2D Novelty network plot for each Supergroup and prepare for phylogeny
#filter the network to get the novel ASVs or environmental ASVs

for (i in unique(TAX$S)) {
  network_group <- get_network(level = "S", group = i, edges = 95, 
                               NODES = nodes, RELATIONS = relations)
  aux <- prepare_for_plot(network_group, color_level = "C")
  pdf(file.path(path, paste0("Novelty/",i,"_novelty.pdf")), width = 3.5, height = 3.5)
    print(plot_network(aux, plot_lines = T, plot_mem = F, mcimg = mcim, 
                       mcemg = mcem, mcem = mean(network_group$NODES$CEM), mcim = mean(network_group$NODES$CIM), Title = i))
  dev.off()
  
  stat <- get_stat(network_group$NODES, "C")
  write.table(stat, file.path(path,paste0("Novelty/", i,"_novelty_stat.tsv")), sep = "\t", quote = F, row.names = F)
  
  network_novel <- get_novels(network = network_group,
                              table = ASVTABLE, stat = T,
                              maxCIM = mean(network_group$NODES$CIM),
                              maxCEM = mean(network_group$NODES$CEM))
  aux <- prepare_for_plot(network = network_novel, color_level = "C")
  pdf(file.path(path, paste0("Novelty/",i,"_novelty_novels.pdf")), width = 3.5, height = 3.5)
    print(plot_network(aux, plot_lines = T, plot_mem = T, mcimg = mcim, 
                       mcemg = mcem, mcem = mean(network_group$NODES$CEM), mcim = mean(network_group$NODES$CIM)))
  dev.off()
}


#ASVs PARTICULAR GROUPS

#Perkinsea
network_group <- get_network(level = "C", group = "Alveolata", edges = 95, 
                             NODES = nodes, RELATIONS = relations)
aux <- prepare_for_plot(network_group)
print(plot_network(aux, plot_lines = T, plot_mem = F, mcimg = mcim, 
                   mcemg = mcem, mcem = mean(network_group$NODES$CEM), mcim = mean(network_group$NODES$CIM)))

network_novel <- get_novels(level = "D", group = "Perkinsea", network = network_group,
                            table = ASVTABLE, stat = T,
                            maxCIM = mean(network_group$NODES$CIM),
                            maxCEM = mean(network_group$NODES$CEM))

new_fasta <- FASTA[grep(T, names(FASTA) %in% rownames(network_novel$NOVEL))]
aux_seq <- data.frame(ASV = names(new_fasta), seq = as.character(new_fasta))
aux_seq <- merge(aux_seq, network_novel$NODES, all.x = T, by.x = "ASV", by.y ="name")
aux_seq$ASV <- paste0(aux_seq$ASV, "_", str_replace_all(aux_seq$C, "\\W+", "_"), "_cluster_", aux_seq$Cluster)
format_seq <- Biostrings::DNAStringSet(aux_seq$seq)
names(format_seq) <- aux_seq$ASV

Biostrings::writeXStringSet(format_seq, file.path(path, paste0("Phylogeny/Perkinsea/Perkinsea_ASVs.fasta")))


#Bicosoecida
network_group <- get_network(level = "S", group = "Stramenopiles", edges = 95, 
                             NODES = nodes, RELATIONS = relations)
aux <- prepare_for_plot(network_group)
print(plot_network(aux, plot_lines = T, plot_mem = F, mcimg = mcim, 
                   mcemg = mcem, mcem = mean(network_group$NODES$CEM), mcim = mean(network_group$NODES$CIM)))

network_novel <- get_novels(level = "C", group = "Bicoecea", network = network_group,
                            table = ASVTABLE, stat = T,
                            maxCIM = mean(network_group$NODES$CIM),
                            maxCEM = mean(network_group$NODES$CEM))

new_fasta <- FASTA[grep(T, names(FASTA) %in% rownames(network_novel$NOVEL))]
aux_seq <- data.frame(ASV = names(new_fasta), seq = as.character(new_fasta))
aux_seq <- merge(aux_seq, network_novel$NODES, all.x = T, by.x = "ASV", by.y ="name")
aux_seq$ASV <- paste0(aux_seq$ASV, "_", str_replace_all(aux_seq$C, "\\W+", "_"), "_cluster_", aux_seq$Cluster)
format_seq <- Biostrings::DNAStringSet(aux_seq$seq)
names(format_seq) <- aux_seq$ASV

Biostrings::writeXStringSet(format_seq, file.path(path, paste0("Phylogeny/Bicosoecida/Bicosoecida_ASVs.fasta")))


#Opisthokonta

network_group <- get_network(level = "S", group = "Opisthokonta", edges = 95, 
                             NODES = nodes, RELATIONS = relations)

aux <- prepare_for_plot(network_group)
print(plot_network(aux, plot_lines = T, plot_mem = F, mcimg = mcim, 
                   mcemg = mcem, mcem = mean(network_group$NODES$CEM), mcim = mean(network_group$NODES$CIM)))

network_novel <- get_novels(level = "G", group = "Choanoflagellata", network = network_group,
                            table = ASVTABLE, stat = T,
                            maxCIM = mean(network_group$NODES$CIM),
                            maxCEM = mean(network_group$NODES$CEM))

new_fasta <- FASTA[grep(T, names(FASTA) %in% rownames(network_novel$NOVEL))]
aux_seq <- data.frame(ASV = names(new_fasta), seq = as.character(new_fasta))
aux_seq <- merge(aux_seq, network_novel$NODES, all.x = T, by.x = "ASV", by.y ="name")
aux_seq$ASV <- paste0(aux_seq$ASV, "_", str_replace_all(aux_seq$C, "\\W+", "_"), "_cluster_", aux_seq$Cluster)
format_seq <- Biostrings::DNAStringSet(aux_seq$seq)
names(format_seq) <- aux_seq$ASV

Biostrings::writeXStringSet(format_seq, file.path(path, paste0("Phylogeny/Opisthokonta/Choanoflagellata_ASVs.fasta")))


#Archaeplastida <-  Mamiellophyceae + Pedinophyceae

network_group <- get_network(level = "S", group = "Archaeplastida", edges = 95, 
                             NODES = nodes, RELATIONS = relations)
aux <- prepare_for_plot(network_group)
print(plot_network(aux, plot_lines = T, plot_mem = F, mcimg = mcim, 
                   mcemg = mcem, mcem = mean(network_group$NODES$CEM), mcim = mean(network_group$NODES$CIM)))

#Mamiellophyceae
network_novel <- get_novels(level = "C", group = "Mamiellophyceae", network = network_group,
                            table = ASVTABLE, stat = T,
                            maxCIM = mean(network_group$NODES$CIM),
                            maxCEM = mean(network_group$NODES$CEM))

new_fasta <- FASTA[grep(T, names(FASTA) %in% rownames(network_novel$NOVEL))]
aux_seq <- data.frame(ASV = names(new_fasta), seq = as.character(new_fasta))
aux_seq <- merge(aux_seq, network_novel$NODES, all.x = T, by.x = "ASV", by.y ="name")
aux_seq$ASV <- paste0(aux_seq$ASV, "_", str_replace_all(aux_seq$C, "\\W+", "_"), "_cluster_", aux_seq$Cluster)
format_seq <- Biostrings::DNAStringSet(aux_seq$seq)
names(format_seq) <- aux_seq$ASV

Biostrings::writeXStringSet(format_seq, file.path(path, paste0("Phylogeny/Mamiellophyceae/Mamiellophyceae_ASVs.fasta")))


#Pedinophyceae
network_novel <- get_novels(level = "C", group = "Pedinophyceae", network = network_group,
                            table = ASVTABLE, stat = T,
                            maxCIM = mean(network_group$NODES$CIM),
                            maxCEM = mean(network_group$NODES$CEM))

new_fasta <- FASTA[grep(T, names(FASTA) %in% rownames(network_novel$NOVEL))]
aux_seq <- data.frame(ASV = names(new_fasta), seq = as.character(new_fasta))
aux_seq <- merge(aux_seq, network_novel$NODES, all.x = T, by.x = "ASV", by.y ="name")
aux_seq$ASV <- paste0(aux_seq$ASV, "_", str_replace_all(aux_seq$C, "\\W+", "_"), "_cluster_", aux_seq$Cluster)
format_seq <- Biostrings::DNAStringSet(aux_seq$seq)
names(format_seq) <- aux_seq$ASV

Biostrings::writeXStringSet(format_seq, file.path(path, paste0("Phylogeny/Pedinophyceae/Pedinophyceae_ASVs.fasta")))
