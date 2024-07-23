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
library(microbiome)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(FSA)
library(agricolae)

### SKIN VS WATER SAMPLES COMPARISONS ###

##Importing phyloseq objects that contain all samples 

#Bacteria dataset. In the paper is called it "all.data.bac"
phy.bac.rar<- readRDS("Archivos_Core_Final/phy.bac.rar.rds")

#Fungal dataset. In the paper is called it ""all.data.its"
phy.its.rar<- readRDS("Archivos_Core_Final/phy.its.rar.rds")

#Subsetting "all.data.bac"

phy.bac.axo.rar <- subset_samples(phy.bac.rar, SampleType == "Axolotl")
phy.rar.mex.bac <- subset_samples(phy.bac.axo.rar, Species == "A. mexicanum")
phy.bac.axo.water <- subset_samples(phy.bac.rar, SampleType == "Water")

mex.bac <- subset_samples(phy.bac.rar, Species == "A. mexicanum")
mex.bac.water <- subset_samples(phy.bac.axo.water, Species == "A. mexicanum")
mex.axo <- subset_samples(phy.bac.axo.rar, Species == "A. mexicanum")
and.bac <- subset_samples(phy.bac.rar, Species == "A. andersoni")
and.axo <- subset_samples(phy.bac.axo.rar, Species == "A. andersoni")
and.bac.water <- subset_samples(phy.bac.axo.water, Species == "A. andersoni")
dum.bac <- subset_samples(phy.bac.rar, Species == "A. dumerillii")
dum.bac.water <- subset_samples(phy.bac.axo.water, Species == "A. dumerillii")
dum.axo <- subset_samples(phy.bac.axo.rar, Species == "A. dumerillii")
tay.bac <- subset_samples(phy.bac.rar, Species == "A. taylori")
tay.bac.water <- subset_samples(phy.bac.axo.water, Species == "A. taylori")
tay.axo <- subset_samples(phy.bac.axo.rar, Species == "A. taylori")

#Subsetting "all.data.bac"
phy.its.axo.rar <- subset_samples(phy.its.rar, SampleType == "Axolotl")
phy.rar.mex.its <- subset_samples(phy.its.axo.rar, Species == "A. mexicanum")
phy.its.axo.water <- subset_samples(phy.its.rar, SampleType == "Water")

mex.its <- subset_samples(phy.its.rar, Species == "A. mexicanum")
mex.its.water <- subset_samples(phy.its.axo.water, Species == "A. mexicanum")
mex.axo.its <- subset_samples(phy.its.axo.rar, Species == "A. mexicanum")
and.its <- subset_samples(phy.its.rar, Species == "A. andersoni")
and.its.water <- subset_samples(phy.its.axo.water, Species == "A. andersoni")
and.axo.its <- subset_samples(phy.its.axo.rar, Species == "A. andersoni")
dum.its <- subset_samples(phy.its.rar, Species == "A. dumerilii")
dum.its.water <- subset_samples(phy.its.axo.water, Species == "A. dumerilii")
dum.axo.its <- subset_samples(phy.its.axo.rar, Species == "A. dumerilii")
tay.its <- subset_samples(phy.its.rar, Species == "A. taylori")
tay.its.water <- subset_samples(phy.its.axo.water, Species == "A. taylori")
tay.axo.its <- subset_samples(phy.its.axo.rar, Species == "A. taylori")


## ALPHA DIVERSITY##

#Calculating metrics
phy.bac.div <- microbiome::alpha(phy.bac.rar, index = "all")
phy.bac.div$Kingdom <- c("Bacteria")
phy.its.div <- microbiome::alpha(phy.its.rar, index = "all")
phy.its.div$Kingdom <- c("Fungi")
phy.alpha  <- rbind(phy.bac.div, phy.its.div)
phy.bac.axo.rar

# get the metadata out as seprate object
phy.bac.meta <- meta(phy.bac.rar) #modifying data order to merge the 2 dataframes
phy.its.meta <- meta(phy.its.rar)
phy.alpha.meta  <- rbind(phy.bac.meta, phy.its.meta)

phy.alpha$sam_name <- rownames(phy.alpha)
phy.alpha.meta$sam_name <- rownames(phy.alpha.meta)

# merge these two data frames into one
alpha.metrics <- merge(phy.alpha,phy.alpha.meta, by = "sam_name")

#Saving dataframe as .csv
write.csv(alpha.metrics, "AlphaDiversity.metrics.csv", row.names = FALSE)


#Subsetting data for each analysis
#Bacteria
alpha.bact <- alpha.metrics[alpha.metrics$Kingdom=="Bacteria" , ]
alpha.bac.and <- alpha.bact[alpha.bact$Species=="A. andersoni" , ]
alpha.bac.dum <- alpha.bact[alpha.bact$Species=="A. dumerillii" , ]
alpha.bac.mex <- alpha.bact[alpha.bact$Species=="A. mexicanum" , ]
alpha.bac.tay <- alpha.bact[alpha.bact$Species=="A. taylori" , ]


alpha.bact.axo <- alpha.bact[alpha.bact$SampleType=="Axolotl" , ]
write.csv(alpha.bact.axo, "alpha.axo.bac.all.csv", row.names = FALSE)
alpha.bact.water <- alpha.bact[alpha.bact$SampleType=="Water" , ]

#Fungi
alpha.its <- alpha.metrics[alpha.metrics$Kingdom=="Fungi" , ]
alpha.its.and <- alpha.its[alpha.its$Species=="A. andersoni" , ]
alpha.its.dum <- alpha.its[alpha.its$Species=="A. dumerillii" , ]
alpha.its.mex <- alpha.its[alpha.its$Species=="A. mexicanum" , ]
alpha.its.tay <- alpha.its[alpha.its$Species=="A. taylori" , ]
alpha.its.axo <- alpha.its[alpha.its$SampleType=="Axolotl" , ]
write.csv(alpha.its.axo, "alpha.axo.its.all.csv", row.names = FALSE)
alpha.its.water <- alpha.its[alpha.its$SampleType=="Water" , ]


##Evaluating data distribution of variable "observed ASVs"

hist(alpha.bact$observed)
summary(alpha.bact$observed)
shapiro.test(alpha.bact$observed) # p = 0.0012, < 0.05, non-normal data

hist(alpha.its$observed)
summary(alpha.its$observed)
shapiro.test(alpha.its$observed) # p=1.504e-07, <0.05, non-normal

#Evaluating statistical differences between groups (sample type)

wcx.st.bac <- wilcox.test(observed ~ SampleType, data = alpha.bact)
wcx.st.bac 
wcx.st.its <- wilcox.test(observed ~ SampleType, data = alpha.its)
wcx.st.its 

## Boxplots ##
#Bacteria
bx.bac <- ggplot(alpha.bact,aes(x=SampleType,y=observed, fill=SampleType)) + 
  geom_violin(lwd=1) +scale_fill_manual(values = c("#a97f2fff","#2e77abff")) +
  theme_bw()+geom_jitter(position = position_jitter(seed = 1, width = 0.2))+
  theme(legend.position="bottom", legend.text = element_text(size = 16, face = "bold"),
        legend.title = element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_text(colour = "black", size = 17),
        axis.title.y=element_text(face="bold",colour = "black", size = 17),
        axis.title.x=element_blank())+
  ylab("Observed ASVs  (Bacteria)") + facet_grid(~Species)
bx.bac <- bx.bac + theme(strip.text = element_text(face = "bold.italic",size = 16))
bx.bac
bx.bac.p <- bx.bac + stat_compare_means(label =  "p.signif", label.x = 1.5)+
  theme(legend.position = "none")
bx.bac.p

#Fungi
bx.its <- ggplot(alpha.its,aes(x=SampleType,y=observed, fill=SampleType)) + 
  geom_violin(lwd=1) +scale_fill_manual(values = c("#a97f2fff","#2e77abff")) +
  theme_bw()+ geom_jitter(position = position_jitter(seed = 1, width = 0.2))+
  theme(legend.position="bottom", legend.text = element_text(size = 16, face = "bold"),
        legend.title = element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_text(colour = "black", size = 17),
        axis.title.y=element_text(face="bold",colour = "black", size = 17),
        axis.title.x=element_blank())+
  ylab("Observed ASVs (Fungi)") + facet_grid(~Species)
bx.its <- bx.its + theme(strip.text = element_text(face = "bold.italic",size = 16))
bx.its
bx.its.p <- bx.its + stat_compare_means(label =  "p.signif", label.x = 1.5)
bx.its.p


bx.obs.st <- ggarrange(bx.bac.p, bx.its.p, ncol=1, nrow = 2,
                       font.label = list(size = 13, face = "bold"),
                       common.legend = T,legend = "bottom",
                       labels = c("A","B"))
bx.obs.st
ggsave("Diversity.obs.st.viol.svg", bx.obs.st, width = 8, height = 10,
       limitsize = F, units = "in", path = "./Plots/Articulo/")


### BETA DIVERSITY

#Beta diversity per each host species comparing between sample type
nmds.bray.bac.and <- ordinate(and.bac, "PCoA", "bray")
plot.bray.bac.and <- plot_ordination(and.bac, nmds.bray.bac.and) +
  geom_point(size=8, pch=21, colour="black", aes(fill=SampleType),stroke=2)+ 
  scale_fill_manual(values = c("#A6696F", "#7A85B7")) + theme_bw()+
  theme(legend.text = element_text(size = 18, face = "bold"),
        legend.title = element_text(size = 18,face = "bold"),
        axis.text.y=element_text(colour = "black", size = 18),
        axis.text.x=element_text(colour = "black", size = 19),
        axis.title.y=element_text(face = "bold",colour = "black", size = 19),
        axis.title.x=element_text(face = "bold",colour = "black", size = 19))

nmds.bray.bac.dum <- ordinate(dum.bac, "PCoA", "bray")
plot.bray.bac.dum <- plot_ordination(dum.bac, nmds.bray.bac.dum, color= "SampleType") +
  geom_point(size=8, pch=21, colour="black", aes(fill=SampleType),stroke=2)+ 
  scale_fill_manual(values = c("#A6696F", "#7A85B7")) + theme_bw()+
  theme(legend.text = element_text(size = 18, face = "bold"),
        legend.title = element_text(size = 18,face = "bold"),
        axis.text.y=element_text(colour = "black", size = 18),
        axis.text.x=element_text(colour = "black", size = 19),
        axis.title.y=element_text(face = "bold",colour = "black", size = 19),
        axis.title.x=element_text(face = "bold",colour = "black", size = 19))

nmds.bray.bac.mex <- ordinate(mex.bac, "PCoA", "bray")
plot.bray.bac.mex <- plot_ordination(mex.bac, nmds.bray.bac.mex, color= "SampleType") +
  geom_point(size=8, pch=21, colour="black", aes(fill=SampleType),stroke=2)+ 
  scale_fill_manual(values = c("#A6696F", "#7A85B7")) + theme_bw()+
  theme(legend.text = element_text(size = 18, face = "bold"),
        legend.title = element_text(size = 18,face = "bold"),
        axis.text.y=element_text(colour = "black", size = 18),
        axis.text.x=element_text(colour = "black", size = 19),
        axis.title.y=element_text(face = "bold",colour = "black", size = 19),
        axis.title.x=element_text(face = "bold",colour = "black", size = 19))

nmds.bray.bac.tay <- ordinate(tay.bac, "PCoA", "bray")
plot.bray.bac.tay <- plot_ordination(tay.bac, nmds.bray.bac.tay, color= "SampleType") +
  geom_point(size=8, pch=21, colour="black", aes(fill=SampleType),stroke=2)+
  scale_fill_manual(values = c("#A6696F", "#7A85B7")) + theme_bw()+
  theme(legend.text = element_text(size = 18, face = "bold"),
        legend.title = element_text(size = 18,face = "bold"),
        axis.text.y=element_text(colour = "black", size = 18),
        axis.text.x=element_text(colour = "black", size = 19),
        axis.title.y=element_text(face = "bold",colour = "black", size = 19),
        axis.title.x=element_text(face = "bold",colour = "black", size = 19))


ord.bac.plot <- ggarrange(plot.bray.bac.and,plot.bray.bac.dum,plot.bray.bac.mex, plot.bray.bac.tay, ncol=4, nrow = 1,
                          labels = c("A","B","C","D"), font.label = list(size = 13, face = "bold"),
                          common.legend = F,legend = "none")
ord.bac.plot


ggsave("bac.nmds.png", ord.bac.plot, width = 16, height = 6,
       limitsize = F, units = "in", path = "./Plots/Articulo/")

#Fungi
nmds.bray.its.and <- ordinate(and.its, "PCoA", "bray")
plot.bray.its.and <- plot_ordination(and.its, nmds.bray.its.and, color= "SampleType") +
  geom_point(size=8, pch=21, colour="black", aes(fill=SampleType),stroke=2)+
  scale_fill_manual(values = c("#A6696F", "#7A85B7")) + theme_bw()+
  theme(legend.text = element_text(size = 18, face = "bold"),
        legend.title = element_text(size = 18,face = "bold"),
        axis.text.y=element_text(colour = "black", size = 18),
        axis.text.x=element_text(colour = "black", size = 19),
        axis.title.y=element_text(face = "bold",colour = "black", size = 19),
        axis.title.x=element_text(face = "bold",colour = "black", size = 19))

nmds.bray.its.dum <- ordinate(dum.its, "PCoA", "bray")
plot.bray.its.dum <- plot_ordination(dum.its, nmds.bray.its.dum, color= "SampleType") +
  geom_point(size=8, pch=21, colour="black", aes(fill=SampleType),stroke=2)+ 
  scale_fill_manual(values = c("#A6696F", "#7A85B7")) + theme_bw()+
  theme(legend.text = element_text(size = 18, face = "bold"),
        legend.title = element_text(size = 18,face = "bold"),
        axis.text.y=element_text(colour = "black", size = 18),
        axis.text.x=element_text(colour = "black", size = 17),
        axis.title.y=element_text(face = "bold",colour = "black", size = 19),
        axis.title.x=element_text(face = "bold",colour = "black", size = 19)) 

nmds.bray.its.mex <- ordinate(mex.its, "PCoA", "bray")
plot.bray.its.mex <- plot_ordination(mex.its, nmds.bray.its.mex, color= "SampleType") +
  geom_point(size=8, pch=21, colour="black", aes(fill=SampleType),stroke=2)+
  scale_fill_manual(values = c("#A6696F", "#7A85B7")) + theme_bw()+
  theme(legend.text = element_text(size = 18, face = "bold"),
        legend.title = element_text(size = 18,face = "bold"),
        axis.text.y=element_text(colour = "black", size = 18),
        axis.text.x=element_text(colour = "black", size = 19),
        axis.title.y=element_text(face = "bold",colour = "black", size = 19),
        axis.title.x=element_text(face = "bold",colour = "black", size = 19))

nmds.bray.its.tay <- ordinate(tay.its, "PCoA", "bray")
plot.bray.its.tay <- plot_ordination(tay.its, nmds.bray.its.tay, color= "SampleType") +
  geom_point(size=8, pch=21, colour="black", aes(fill=SampleType),stroke=2)+
  scale_fill_manual(values = c("#A6696F", "#7A85B7")) + theme_bw()+
  theme(legend.text = element_text(size = 18, face = "bold"),
        legend.title = element_text(size = 18,face = "bold"),
        axis.text.y=element_text(colour = "black", size = 18),
        axis.text.x=element_text(colour = "black", size = 19),
        axis.title.y=element_text(face = "bold",colour = "black", size = 19),
        axis.title.x=element_text(face = "bold",colour = "black", size = 19))

ord.its.plot <- ggarrange(plot.bray.its.and,plot.bray.its.dum,plot.bray.its.mex, plot.bray.its.tay, ncol=4, nrow = 1,
                          labels = c("E","F","G","H"), font.label = list(size = 13, face = "bold"),
                          common.legend = T,legend = "bottom")

ord.div.plot <- ggarrange(ord.bac.plot, ord.its.plot, ncol=1, nrow = 2,
                          font.label = list(size = 13, face = "bold"),
                          common.legend = T,legend = "bottom")
ord.div.plot


ggsave("beta.st.final.svg", ord.div.plot, width = 19, height = 13,
       limitsize = F, units = "in", path = "./Plots/Articulo/")


#Beta diversity all samples
#Bacteria
#Bray-Curtis
nmds.bray.bac.st <- ordinate(phy.bac.rar, "PCoA", "bray")
plot.bray.bac.st <- plot_ordination(phy.bac.rar, nmds.bray.bac.st, color = "SampleType",shape="Species") +
  geom_point(size=5)+ scale_color_manual(values=c("#A6696F", "#7A85B7"))+ theme_bw()+scale_shape_manual(values=c(15, 16, 17,18))+
  theme(legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        axis.text.y=element_text(colour = "black", size = 18),
        axis.text.x=element_text(colour = "black", size = 19),
        axis.title.y=element_text(colour = "black", size = 19),
        axis.title.x=element_text(colour = "black", size = 19))
plot.bray.bac.st 
ggsave("Bray.pcoa.bac.st.png", plot.bray.bac.st, width = 10, height = 8,
       limitsize = F, units = "in", path = "./Plots/Articulo/")

#Jaccard
nmds.jac.bac.st <- ordinate(phy.bac.rar, "PCoA", "jaccard",binary=T)
plot.jac.bac.st <- plot_ordination(phy.bac.rar, nmds.jac.bac.st, color = "SampleType",shape = "Species") +
  geom_point(size=5)+ scale_color_manual(values=c("#A6696F", "#7A85B7"))+ theme_bw()+scale_shape_manual(values=c(15, 16, 17,18))+
  theme(legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        axis.text.y=element_text(colour = "black", size = 18),
        axis.text.x=element_text(colour = "black", size = 19),
        axis.title.y=element_text(colour = "black", size = 19),
        axis.title.x=element_text(colour = "black", size = 19))
plot.jac.bac.st
ggsave("Jac.pcoa.bac.st.png", plot.jac.bac.st, width = 10, height = 8,
       limitsize = F, units = "in", path = "./Plots/Articulo/")

#Fungi
nmds.bray.its.st <- ordinate(phy.its.rar, "PCoA", "bray")
plot.bray.its.st <- plot_ordination(phy.its.rar, nmds.bray.its.st, color = "SampleType",shape="Species") +
  geom_point(size=5)+ scale_color_manual(values=c("#A6696F", "#7A85B7"))+ theme_bw()+scale_shape_manual(values=c(15, 16, 17,18))+
  theme(legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        axis.text.y=element_text(colour = "black", size = 18),
        axis.text.x=element_text(colour = "black", size = 19),
        axis.title.y=element_text(colour = "black", size = 19),
        axis.title.x=element_text(colour = "black", size = 19))
plot.bray.its.st
ggsave("Bray.pcoa.its.st.png", plot.bray.its.st, width = 10, height = 8,
       limitsize = F, units = "in", path = "./Plots/Articulo/")

nmds.jac.its.st <- ordinate(phy.its.rar, "PCoA", "jaccard",binary=T)
plot.jac.its.st <- plot_ordination(phy.its.rar, nmds.jac.its.st, color = "SampleType",shape="Species") +
  geom_point(size=5)+ scale_color_manual(values=c("#A6696F", "#7A85B7"))+ theme_bw()+scale_shape_manual(values=c(15, 16, 17,18))+
  theme(legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        axis.text.y=element_text(colour = "black", size = 18),
        axis.text.x=element_text(colour = "black", size = 19),
        axis.title.y=element_text(colour = "black", size = 19),
        axis.title.x=element_text(colour = "black", size = 19))
plot.jac.its.st
ggsave("Jac.pcoa.its.st.png", plot.jac.its.st, width = 10, height = 8,
       limitsize = F, units = "in", path = "./Plots/Articulo/")

fig.ord.all.st <- ggarrange(plot.bray.bac.st, plot.bray.its.st,
                            plot.jac.bac.st, plot.jac.its.st, 
                            ncol=2, nrow = 2,
                            labels = c("A","B","C","D"), font.label = list(size = 15, face = "bold"),
                            common.legend = T,legend = "right")

fig.ord.all.st

### HOST SPECIES COMPARISONS ###

##Importing phyloseq objects that contain skin associated microbiota  

#Bacteria dataset. In the paper is called it "core.data.bac"
phy.axo.core.bac<- readRDS("Archivos_Core_Final/phy.axo.core.bac.last")

#Fungal dataset. In the paper is called it "core.data.its"
phy.axo.core.its<- readRDS("Archivos_Core_Final/phy.axo.core.its.last")

#Calculating alpha diversity

phy.bac.div.core <- microbiome::alpha(phy.axo.core.bac, index = "all")
phy.bac.div.core$Kingdom <- c("Bacteria")
phy.its.div.core <- microbiome::alpha(phy.axo.core.its, index = "all")
phy.its.div.core$Kingdom <- c("Fungi")
phy.alpha.core  <- rbind(phy.bac.div.core, phy.its.div.core)

# get the metadata out as separate object
phy.bac.meta.core <- meta(phy.axo.core.bac) #modifying data order to merge the 2 dataframes
phy.its.meta.core <- meta(phy.axo.core.its)
phy.alpha.meta.core  <- rbind(phy.bac.meta.core, phy.its.meta.core)

phy.alpha.core$sam_name <- rownames(phy.alpha.core)
phy.alpha.meta.core$sam_name <- rownames(phy.alpha.meta.core)

# merge these two data frames into one
alpha.metrics.core <- merge(phy.alpha.core,phy.alpha.meta.core, by = "sam_name")
alpha.metrics.core <- read.csv("AlphaDiversityCore.metrics.csv", header = T)

alpha.metrics.core.bac <- subset(alpha.metrics.core, Kingdom == "Bacteria")
alpha.metrics.core.its <- subset(alpha.metrics.core, Kingdom == "Fungi")



#Alpha diversity core
#Alpha diversity boxplots core

bac.summarized.core = alpha.metrics.core.bac %>% group_by(Species) %>% summarize(observed=max(observed))
alpha.metrics.core.bac <- alpha.metrics.core.bac[order(alpha.metrics.core.bac$Species),]
hsd.bac=HSD.test(aov(observed~Species,data=alpha.metrics.core.bac), "Species", group=T)
hsd.bac
hsd.df.bac <- as.data.frame(hsd.bac$groups)
hsd.df.bac$sample <- row.names(hsd.df.bac)
hsd.df.bac <- hsd.df.bac[order(hsd.df.bac$sample),]

axo.core.bac <- ggplot(alpha.metrics.core.bac, aes(x=Species, y=observed, fill=Species)) + 
  geom_violin(lwd=1) + geom_jitter(position = position_jitter(seed = 1, width = 0.2))+
  scale_fill_locuszoom() + theme_classic2() +
  geom_text(data=bac.summarized.core,aes(x=Species,y=0.2+observed,label=hsd.df.bac$groups),vjust=0,fontface="bold", size = 5)+
  theme(legend.position="none",
        legend.text = element_text(face = "bold.italic",colour = "black", size = 11),
        legend.title = element_text(face = "bold",colour = "black", size = 11),
        axis.text.x=element_blank(),
        axis.text.y=element_text(colour = "black", size = 12),
        axis.title.y=element_text(face = "bold",colour = "black", size = 12),
        axis.title.x=element_blank())+
  ylab("Observed ASVs (Bacteria)")
axo.core.bac

its.summarized.core = alpha.metrics.core.its %>% group_by(Species) %>% summarize(observed=max(observed))
alpha.metrics.core.its <- alpha.metrics.core.its[order(alpha.metrics.core.its$Species),]
hsd=HSD.test(aov(observed~Species,data=alpha.metrics.core.its), "Species", group=T)
hsd
hsd.df <- as.data.frame(hsd$groups)
hsd.df$sample <- row.names(hsd.df)
hsd.df <- hsd.df[order(hsd.df$sample),]

axo.core.its <- ggplot(alpha.metrics.core.its, aes(x=Species, y=observed, fill=Species)) + 
  geom_violin(lwd=1) + geom_jitter(position = position_jitter(seed = 1, width = 0.2))+
  scale_fill_locuszoom() + theme_classic2() +
  geom_text(data=its.summarized.core,aes(x=Species,y=0.2+observed,label=hsd.df$groups),vjust=0,fontface="bold", size = 5)+
  theme(legend.position="bottom",
        legend.text = element_text(face = "bold.italic",colour = "black", size = 11),
        legend.title = element_text(face = "bold",colour = "black", size = 11),
        axis.text.x=element_blank(),
        axis.text.y=element_text(colour = "black", size = 12),
        axis.title.y=element_text(face = "bold",colour = "black", size = 12),
        axis.title.x=element_blank())+
  ylab("Observed ASVs (Fungi)")
axo.core.its

fig2.1 <- ggarrange(axo.core.bac, axo.core.its, ncol=2, nrow = 1,
                    labels = c("A","B"), font.label = list(size = 15, face = "bold"),
                    common.legend = F)
fig2.1


## Beta diversity of skin associated microbiota, comparing between host species

#Bacteria

nmds.bray.bac.core <- ordinate(phy.axo.core.bac, "PCoA", "bray")
plot.bray.bac.core <- plot_ordination(phy.axo.core.bac, nmds.bray.bac.core, color= "Species") +
  geom_point(size=5, pch=21, colour="black", aes(fill=Species),stroke=1.5)+
  scale_fill_locuszoom() + theme_classic2() +
  theme(legend.text = element_text(size = 13, face = "bold.italic"),
        legend.title = element_text(size = 13,face = "bold"),
        axis.text.y=element_text(colour = "black", size = 12),
        axis.text.x=element_text(colour = "black", size = 13),
        axis.title.y=element_text(face = "bold",colour = "black", size = 13),
        axis.title.x=element_text(face = "bold",colour = "black", size = 13))
plot.bray.bac.core

#Fungi

nmds.bray.its.core <- ordinate(phy.axo.core.its, "PCoA", "bray")
plot.bray.its.core <- plot_ordination(phy.axo.core.its, nmds.bray.its.core, color = "Species") +
  geom_point(size=5, pch=21, colour="black", aes(fill=Species),stroke=1.5)+
  scale_fill_locuszoom() + theme_classic2() + 
  theme(legend.text = element_text(size = 13, face = "bold.italic"),
        legend.title = element_text(size = 13,face = "bold"),
        axis.text.y=element_text(colour = "black", size = 12),
        axis.text.x=element_text(colour = "black", size = 13),
        axis.title.y=element_text(face = "bold", colour = "black", size = 13),
        axis.title.x=element_text(face = "bold", colour = "black", size = 13))
plot.bray.its.core




fig.sp <- ggarrange(axo.core.bac,axo.core.its,
                    plot.bray.bac.core, plot.bray.its.core, ncol=2, nrow = 2,
                    labels = c("A","B","C","D"), font.label = list(size = 15, face = "bold"),
                    common.legend = T, legend = "bottom")
fig.sp


#STATISTICAL ANALYSIS

#Kruskal Wallis

kw.sp.bac <- kruskal.test(observed ~ Species, data = alpha.metrics.core.bac)
kw.sp.bac 

#Dunn test
dunn.sp.bac <- dunnTest(observed ~ Species, data = alpha.metrics.core.bac, method = "bonferroni")
dunn.sp.bac

kw.disp.bac <- kruskal.test(dispersion ~ group, data = dc_spp)
kw.disp.bac 

kw.sp.its <- kruskal.test(observed ~ Species, data = alpha.metrics.core.its)
kw.sp.its 
dunn.sp.its <- dunnTest(observed ~ Species, data = alpha.metrics.core.its, method = "bonferroni")
dunn.sp.its 

### PERMANOVA AND PERMUTEST ANALYSYS TO EVALUATE MICROBIAL COMPOSITION BETWEEN SAMPLE TYPE ##

## BACTERIA

## PERMANOVA AND PERMUTEST CONSIDERING ALL SAMPLES 

#Bray-Curtis

metadf.bac <- data.frame(sample_data(phy.bac.rar))
dist.bac.bray <- distance(phy.bac.rar, method = "bray")
bac.st.betadis <- betadisper(d = dist.bac.bray, group = metadf.bac$SampleType, type = "centroid")
plot(bac.st.betadis)
permutest(bac.st.betadis, pairwise = TRUE) 

perma.bac <- adonis(t(phy.bac.rar@otu_table) ~ SampleType, 
                    data = metadf.bac, permutations=999, method = "bray")
perma.bac 

ps.disper.sp <- betadisper(dist.sp, meta.bac.axo$Species)
permutest(ps.disper.sp, pairwise = TRUE)

#Jaccard
dist.bac.jac <- distance(phy.bac.rar, method = "jaccard",binary=T)
bac.st.betadis.jac <- betadisper(d = dist.bac.jac, group = metadf.bac$SampleType, type = "centroid")
plot(bac.st.betadis.jac)
permutest(bac.st.betadis.jac, pairwise = TRUE) 

metadf.bac <- data.frame(sample_data(phy.bac.rar))
perma.bac.jac <- adonis(t(phy.bac.rar@otu_table) ~ SampleType, 
                        data = metadf.bac, permutations=999, method = "jaccard",binary=T)
perma.bac.jac 

## PERMANOVA AND PERMUTEST  BETWEEN SAMPLE TYPE FOR A. ANDERSONI
#Bray-Curtis
dist.bac.bray.and <- distance(and.bac, method = "bray")
bac.st.betadis.and <- betadisper(d = dist.bac.bray.and, group = metadf.bac.and$SampleType, type = "centroid")
plot(bac.st.betadis.and)
permutest(bac.st.betadis.and, pairwise = TRUE) 

metadf.bac.and <- data.frame(sample_data(and.bac))
perma.bac.and <- adonis(t(and.bac@otu_table) ~ SampleType, 
                        data = metadf.bac.and, permutations=999, method = "bray")
perma.bac.and #

#Jaccard
dist.bac.jac.and <- distance(and.bac, method = "jaccard",binary=T)
bac.st.betadis.and.jac <- betadisper(d = dist.bac.jac.and, group = metadf.bac.and$SampleType, type = "centroid")
plot(bac.st.betadis.and.jac)
permutest(bac.st.betadis.and.jac, pairwise = TRUE) 

metadf.bac.and <- data.frame(sample_data(and.bac))
perma.bac.and.jac <- adonis(t(and.bac@otu_table) ~ SampleType, 
                            data = metadf.bac.and, permutations=999, method = "jaccard",binary=T)
perma.bac.and.jac 

## PERMANOVA AND PERMUTEST  BETWEEN SAMPLE TYPE FOR A. dumerilii
#Bray-Curtis
dist.bac.bray.dum <- distance(dum.bac, method = "bray")
bac.st.betadis.dum <- betadisper(d = dist.bac.bray.dum, group = metadf.bac.dum$SampleType, type = "centroid")
plot(bac.st.betadis.dum)
permutest(bac.st.betadis.dum, pairwise = TRUE)

metadf.bac.dum <- data.frame(sample_data(dum.bac))
perma.bac.dum <- adonis(t(dum.bac@otu_table) ~ SampleType, 
                        data = metadf.bac.dum, permutations=999, method = "bray")
perma.bac.dum 

#Jaccard
dist.bac.jac.dum <- distance(dum.bac, method = "jaccard",binary=T)
bac.st.betadis.dum.jac <- betadisper(d = dist.bac.jac.dum, group = metadf.bac.dum$SampleType, type = "centroid")
plot(bac.st.betadis.dum.jac)
permutest(bac.st.betadis.dum.jac, pairwise = TRUE) 

metadf.bac.dum <- data.frame(sample_data(dum.bac))
perma.bac.dum.jac <- adonis(t(dum.bac@otu_table) ~ SampleType, 
                            data = metadf.bac.dum, permutations=999, method = "jaccard",binary=T)
perma.bac.dum.jac 

## PERMANOVA AND PERMUTEST  BETWEEN SAMPLE TYPE FOR A. mexicanum
#Bray-Curtis
dist.bac.bray.mex <- distance(mex.bac, method = "bray")
bac.st.betadis.mex <- betadisper(d = dist.bac.bray.mex, group = metadf.bac.mex$SampleType, type = "centroid")
plot(bac.st.betadis.mex)
permutest(bac.st.betadis.mex, pairwise = TRUE)

metadf.bac.mex <- data.frame(sample_data(mex.bac))
perma.bac.mex <- adonis(t(mex.bac@otu_table) ~ SampleType, 
                        data = metadf.bac.mex, permutations=999, method = "bray")
perma.bac.mex 

#Jaccard
dist.bac.jac.mex <- distance(mex.bac, method = "jaccard",binary=T)
bac.st.betadis.mex.jac <- betadisper(d = dist.bac.jac.mex, group = metadf.bac.mex$SampleType, type = "centroid")
plot(bac.st.betadis.mex.jac)
permutest(bac.st.betadis.mex.jac, pairwise = TRUE)

metadf.bac.mex <- data.frame(sample_data(mex.bac))
perma.bac.mex.jac <- adonis(t(mex.bac@otu_table) ~ SampleType, 
                            data = metadf.bac.mex, permutations=999, method = "jaccard",binary=T)
perma.bac.mex.jac 

## PERMANOVA AND PERMUTEST  BETWEEN SAMPLE TYPE FOR A. taylori
#Bray-Curtis
dist.bac.bray.tay <- distance(tay.bac, method = "bray")
bac.st.betadis.tay <- betadisper(d = dist.bac.bray.tay, group = metadf.bac.tay$SampleType, type = "centroid")
plot(bac.st.betadis.tay)
permutest(bac.st.betadis.tay, pairwise = TRUE) 

metadf.bac.tay <- data.frame(sample_data(tay.bac))
perma.bac.tay <- adonis(t(tay.bac@otu_table) ~ SampleType, 
                        data = metadf.bac.tay, permutations=999, method = "bray")
perma.bac.tay 

#Jaccard
dist.bac.jac.tay <- distance(tay.bac, method = "jaccard",binary=T)
bac.st.betadis.mex.jac <- betadisper(d = dist.bac.jac.tay, group = metadf.bac.tay$SampleType, type = "centroid")
plot(bac.st.betadis.mex.jac)
permutest(bac.st.betadis.mex.jac, pairwise = TRUE)

metadf.bac.tay <- data.frame(sample_data(tay.bac))
perma.bac.tay.jac <- adonis(t(tay.bac@otu_table) ~ SampleType, 
                            data = metadf.bac.tay, permutations=999, method = "jaccard",binary=T)
perma.bac.tay.jac

## FUNGI

## PERMANOVA AND PERMUTEST CONSIDERING ALL SAMPLES 

#Bray-Curtis
dist.its.bray <- distance(phy.its.rar, method = "bray")
its.st.betadis <- betadisper(d = dist.its.bray, group = metadf.its$SampleType, type = "centroid")
plot(its.st.betadis)
permutest(its.st.betadis, pairwise = TRUE) 

metadf.its <- data.frame(sample_data(phy.its.rar))
perma.its <- adonis(t(phy.its.rar@otu_table) ~ SampleType, 
                    data = metadf.its, permutations=999, method = "bray")
perma.its 

#Jaccard
dist.its.jac <- distance(phy.its.rar, method = "jaccard",binary=T)
its.st.betadis.jac <- betadisper(d = dist.its.jac, group = metadf.its$SampleType, type = "centroid")
plot(its.st.betadis.jac)
permutest(its.st.betadis.jac, pairwise = TRUE)

perma.its.jac <- adonis(t(phy.its.rar@otu_table) ~ SampleType, 
                        data = metadf.its, permutations=999, method = "jaccard",binary=T)
perma.its.jac 

## PERMANOVA AND PERMUTEST  BETWEEN SAMPLE TYPE FOR A. andersoni
#Bray-Curtis
dist.its.bray.and <- distance(and.its, method = "bray")
its.st.betadis.and <- betadisper(d = dist.its.bray.and, group = metadf.its.and$SampleType, type = "centroid")
plot(its.st.betadis.and)
permutest(its.st.betadis.and, pairwise = TRUE) 

metadf.its.and <- data.frame(sample_data(and.its))
perma.its.and <- adonis(t(and.its@otu_table) ~ SampleType, 
                        data = metadf.its.and, permutations=999, method = "bray")
perma.its.and

#Jaccard
dist.its.jac.and <- distance(and.its, method = "jaccard",binary=T)
its.st.betadis.and.jac <- betadisper(d = dist.its.jac.and, group = metadf.its.and$SampleType, type = "centroid")
plot(its.st.betadis.and.jac)
permutest(its.st.betadis.and.jac, pairwise = TRUE) 

perma.its.and.jac <- adonis(t(and.its@otu_table) ~ SampleType, 
                            data = metadf.its.and, permutations=999, method = "jaccard",binary=T)
perma.its.and.jac 

## PERMANOVA AND PERMUTEST  BETWEEN SAMPLE TYPE FOR A. dumerilii

#Bray-Curtis
dist.its.bray.dum <- distance(dum.its, method = "bray")
its.st.betadis.dum <- betadisper(d = dist.its.bray.dum, group = metadf.its.dum$SampleType, type = "centroid")
plot(its.st.betadis.dum)
permutest(its.st.betadis.dum, pairwise = TRUE) 

metadf.its.dum <- data.frame(sample_data(dum.its))
perma.its.dum <- adonis(t(dum.its@otu_table) ~ SampleType, 
                        data = metadf.its.dum, permutations=999, method = "bray")
perma.its.dum 

#Jaccard
dist.its.jac.dum <- distance(dum.its, method = "jaccard",binary=T)
its.st.betadis.dum.jac <- betadisper(d = dist.its.jac.dum, group = metadf.its.dum$SampleType, type = "centroid")
plot(its.st.betadis.dum.jac)
permutest(its.st.betadis.dum.jac, pairwise = TRUE) 

perma.its.dum.jac <- adonis(t(dum.its@otu_table) ~ SampleType, 
                            data = metadf.its.dum, permutations=999, method = "jaccard",binary=T)
perma.its.dum.jac

## PERMANOVA AND PERMUTEST  BETWEEN SAMPLE TYPE FOR A. mexicanum

#Bray-Curtis
dist.its.bray.mex <- distance(mex.its, method = "bray")
its.st.betadis.mex <- betadisper(d = dist.its.bray.mex, group = metadf.its.mex$SampleType, type = "centroid")
plot(its.st.betadis.mex)
permutest(its.st.betadis.mex, pairwise = TRUE) 

metadf.its.mex <- data.frame(sample_data(mex.its))
perma.its.mex <- adonis(t(mex.its@otu_table) ~ SampleType, 
                        data = metadf.its.mex, permutations=999, method = "bray")
perma.its.mex 

#Jaccard

dist.its.jac.mex <- distance(mex.its, method = "jaccard",binary=T)
its.st.betadis.mex.jac <- betadisper(d = dist.its.jac.mex, group = metadf.its.mex$SampleType, type = "centroid")
plot(its.st.betadis.mex.jac)
permutest(its.st.betadis.mex.jac, pairwise = TRUE)

perma.its.mex.jac <- adonis(t(mex.its@otu_table) ~ SampleType, 
                            data = metadf.its.mex, permutations=999, method = "jaccard",binary=T)
perma.its.mex.jac 

## PERMANOVA AND PERMUTEST  BETWEEN SAMPLE TYPE FOR A. taylori
#Bray-Curtis
dist.its.bray.tay <- distance(tay.its, method = "bray")
its.st.betadis.tay <- betadisper(d = dist.its.bray.tay, group = metadf.its.tay$SampleType, type = "centroid")
plot(its.st.betadis.tay)
permutest(its.st.betadis.tay, pairwise = TRUE) 

metadf.its.tay <- data.frame(sample_data(tay.its))
perma.its.tay <- adonis(t(tay.its@otu_table) ~ SampleType, 
                        data = metadf.its.tay, permutations=999, method = "bray")
perma.its.tay 

#Jaccard
dist.its.jac.tay <- distance(tay.its, method = "jaccard",binary=T)
its.st.betadis.tay.jac <- betadisper(d = dist.its.jac.tay, group = metadf.its.tay$SampleType, type = "centroid")
plot(its.st.betadis.tay.jac)
permutest(its.st.betadis.tay.jac, pairwise = TRUE) 

perma.its.tay.jac <- adonis(t(tay.its@otu_table) ~ SampleType, 
                            data = metadf.its.tay, permutations=999, method = "jaccard",binary=T)
perma.its.tay.jac 


### PERMANOVA AND PERMUTEST ANALYSYS TO EVALUATE DIFFERENCES IN SKIN MICROBIOTA COMPOSITION BETWEEN HOST SPECIES ##

## Bacteria 

#Bray-Curtis
meta.bac.axo <- data.frame(sample_data(phy.axo.core.bac))
dist.sp <- distance(phy.axo.core.bac, method = "bray")
permanova.sp <- adonis(t(phy.axo.core.bac@otu_table) ~ Species, 
                       data = meta.bac.axo, permutations=999, method = "bray")
permanova.sp #
print(as.data.frame(permanova_loc$aov.tab)["SampleType", "Pr(>F)"])
ps.disper.sp <- betadisper(dist.sp, meta.bac.axo$Species,type = "centroid")
permutest(ps.disper.sp, pairwise = TRUE) 


#Jaccard
dist.sp.jac <- distance(phy.axo.core.bac, method = "jaccard",binary=T)
permanova.sp.jac <- adonis(t(phy.axo.core.bac@otu_table) ~ Species, 
                           data = meta.bac.axo, permutations=999, method = "jaccard",binary=T)
permanova.sp.jac 
ps.disper.sp.jac <- betadisper(dist.sp.jac, meta.bac.axo$Species,type = "centroid")
permutest(ps.disper.sp.jac, pairwise = TRUE) 


## Fungi

#Bray-Curtis
meta.its.axo <- data.frame(sample_data(phy.axo.core.its))
dist.sp.its <- distance(phy.axo.core.its, method = "bray")
permanova.sp.its <- adonis(t(phy.axo.core.its@otu_table) ~ Species, 
                           data = meta.its.axo, permutations=999, method = "bray")
permanova.sp.its
print(as.data.frame(permanova_loc$aov.tab)["SampleType", "Pr(>F)"])

ps.disper.sp.its <- betadisper(dist.sp.its, meta.its.axo$Species,type = "centroid")
permutest(ps.disper.sp.its, pairwise = TRUE) 

dc_spp.its <-as.data.frame(ps.disper.sp.its$distances)
dc_spp.its$group <- ps.disper.sp.its$group
names(dc_spp.its)

#Jaccard

dist.sp.jac.its <- distance(phy.axo.core.its, method = "jaccard",binary=T)
permanova.sp.jac.its <- adonis(t(phy.axo.core.its@otu_table) ~ Species, 
                               data = meta.its.axo, permutations=999, method = "jaccard",binary=T)
permanova.sp.jac.its 
ps.disper.sp.jac.its <- betadisper(dist.sp.jac.its, meta.its.axo$Species,type = "centroid")
permutest(ps.disper.sp.jac.its, pairwise = TRUE) 
