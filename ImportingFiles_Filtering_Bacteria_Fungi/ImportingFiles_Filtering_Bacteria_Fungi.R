####################################################################################
####                                                                            ####
#### Host Species and Environment Shape the Skin Microbiota of Mexican Axolotls ####
####                         Soto-Cortés et al., 2024                           ####
####                           esotocor@gmail.com                               ####
####                                                                            ####
####################################################################################

#Loading libraries

library(phyloseq) 
library(vegan) 
library(qiime2R)
library(microbiome)

#Bacterial data

phy.bac <- qza_to_phyloseq("table_16S.qza", "rooted-tree_16S.qza", "taxa_16S.qza","metadata_16S_rename.txt")

#Fungal data 
phy.its <- qza_to_phyloseq("table3.qza", "rooted-tree.qza", "taxa.qza","metadata_its_rare.txt")

#Removing contaminants in bacterial dataset

filterPhyla2 <- c("Chloroplast", "Mitochondria", "Eukaryota", "d__Eukaryota", "Vertebrata", "Archaea", "Protista", "d__Archaea")
phy2 <- subset_taxa(phy.bac, !Kingdom %in% filterPhyla2)
phy2 <- subset_taxa(phy2, !Phylum %in% filterPhyla2)
phy2 <- subset_taxa(phy2, !Class %in% filterPhyla2)
phy2 <- subset_taxa(phy2, !Order %in% filterPhyla2)
phy2 <- subset_taxa(phy2, !Family %in% filterPhyla2)
phy2 <- subset_taxa(phy2, !Genus %in% filterPhyla2)
phy2 <- subset_taxa(phy2, Phylum != "unidentified")
phy2 <- subset_taxa(phy2, Phylum != "Unknown")
phy2.bac  <- subset_taxa(phy2, !Species %in% filterPhyla2)


#Removing contaminants in fungal dataset

filterPhyla2 <- c("Chloroplast", "Mitochondria", "Vertebrata", "Archaea", "Protista", "d__Archaea")
filterPhyla3 <- c(" ")
phy2 <- subset_taxa(phy.its, !Kingdom %in% filterPhyla2)
phy2 <- subset_taxa(phy2, !Phylum %in% filterPhyla2)
phy2 <- subset_taxa(phy2, !Class %in% filterPhyla2)
phy2 <- subset_taxa(phy2, !Order %in% filterPhyla2)
phy2 <- subset_taxa(phy2, !Family %in% filterPhyla2)
phy2 <- subset_taxa(phy2, !Genus %in% filterPhyla2)
phy2 <- subset_taxa(phy2, !Species %in% filterPhyla2)
phy2 <- subset_taxa(phy2, !Phylum %in% filterPhyla3)
phy2 <- subset_taxa(phy2, Phylum != "unidentified")
phy2 <- subset_taxa(phy2, Phylum != "Unknown")
phy2.its <- subset_taxa(phy2, Phylum != "Basidiobolomycota")


### Rarefaction ### 

#Bacteria
otu_tab.bac <- t(abundances(phy2.bac))
phy.bac.rar <- rarefy_even_depth(phy2.bac, rngseed=1, sample.size = 10000, replace = F) #NORMALISED SAMPLES
phy.bac.rar

saveRDS(phy.bac.rar, "Archivos_Core_Final/phy.bac.rar.rds")
phy.bac.rar<- readRDS("Archivos_Core_Final/phy.bac.rar.rds")

#Fungi
otu_tab.its <- t(abundances(phy2.its))
phy.its.rar <- rarefy_even_depth(phy2.its, rngseed=1, sample.size = 10000, replace = F) #NORMALISED SAMPLES
phy.its.rar 

saveRDS(phy.its.rar, "Archivos_Core_Final/phy.its.rar.rds")
phy.its.rar<- readRDS("Archivos_Core_Final/phy.its.rar.rds")
phy.its.rar
