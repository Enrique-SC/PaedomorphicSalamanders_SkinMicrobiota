####################################################################################
####                                                                            ####
#### Host Species and Environment Shape the Skin Microbiota of Mexican Axolotls ####
####                         Soto-Cortés et al., 2024                           ####
####                           esotocor@gmail.com                               ####
####                                                                            ####
####################################################################################

#Loading libraries
library(MuMIn)
library(ggfortify)
library("FactoMineR")
library(corrplot)
library("Hmisc")
library(tidyverse)
library(caret)
library(leaps) 
library(MASS)
library(car)
library(lme4)
library(nlme)
library(vegan)
library(ggplot2)
library(ggsci)
library(jtools)
library(sjPlot)
library(mia)
library(QsRutils)
library(factoextra)

#Loading skin-asccociated microbiota dataset
#Bacteria
phy.axo.core.bac <- readRDS("Archivos_Core_Final/phy.axo.core.bac.last")
#Fungi
phy.axo.core.its<- readRDS("Archivos_Core_Final/phy.axo.core.its.last")

phy.bac.div.core <- microbiome::alpha(phy.axo.core.bac, index = "all")
phy.bac.div.core$Kingdom <- c("Bacteria")
phy.its.div.core <- microbiome::alpha(phy.axo.core.its, index = "all")
phy.its.div.core$Kingdom <- c("Fungi")
phy.alpha.core  <- rbind(phy.bac.div.core, phy.its.div.core)

# get the metadata out as seprate object
phy.bac.meta.core <- meta(phy.axo.core.bac) #modifying data order to merge the 2 dataframes
phy.its.meta.core <- meta(phy.axo.core.its)
phy.alpha.meta.core  <- rbind(phy.bac.meta.core, phy.its.meta.core)

phy.alpha.core$sam_name <- rownames(phy.alpha.core)
phy.alpha.meta.core$sam_name <- rownames(phy.alpha.meta.core)

# merge these two data frames into one
alpha.metrics.core <- merge(phy.alpha.core,phy.alpha.meta.core, by = "sam_name")

alpha.metrics.core.bac <- subset(alpha.metrics.core, Kingdom == "Bacteria")
alpha.metrics.core.its <- subset(alpha.metrics.core, Kingdom == "Fungi")

# Generating PCA to evaluate envrionmental differences between localities
#Loading matrix data 
matrix.bac <- read.csv("metadata_modelos_bac_final.csv", header = T)
matrix.its <- read.csv("metadata_modelos_its_final.csv", header = T)
matrix.pca <- matrix.bac[7:32]

pca.bac <- prcomp(matrix.pca, scale = T)
pca.bac.scores <- as.data.frame(pca.bac$x)
pca.bac.rotations <- as.data.frame(pca.bac$rotation)
pca.bac.scores$Locality <- matrix.bac$Locality
pca.bac.scores$Species <- matrix.bac$Species

#Variance explained calculation for PC1 and PC2

k <-  100*pca.bac$sdev^2/sum(pca.bac$sdev^2)
var_pc1 <- as.character(round(k[1], digits = 2))
var_pc2 <- as.character(round(k[2], digits = 2)) 
var_pc1 <- paste0("(", var_pc1, "%)")
var_pc2 <- paste0("(", var_pc2, "%)")

pc1 <- "PC1"
pc2 <- "PC2"
pca1_var <- paste(pc1, var_pc1)
pca2_var <- paste(pc2, var_pc2)

#Ploting PCA

pca.loadings <- data.frame(variables = rownames(pca.bac$rotation), pca.bac$rotation)
pca.bac.scores$Locality <- factor(pca.bac.scores$Locality, levels = c("Zacapu", "Patzcuaro","CDMX", "Alchichica"))

plotpca.bac <- ggplot(pca.bac.scores, aes(x=PC1, y=PC2, colour=Locality)) +
  geom_point(size=8)+ scale_color_locuszoom() + theme_classic()+
  theme(legend.text = element_text(size = 13,face = "bold"),
        legend.title = element_text(size = 15),
        axis.text.y=element_text(colour = "black", size = 14),
        axis.text.x=element_text(colour = "black", size = 15),
        axis.title.y=element_text(face = "bold",colour = "black", size = 15),
        axis.title.x=element_text(face = "bold",colour = "black", size = 15))+
  geom_segment(data = pca.loadings, aes(x = 0, y = 0, xend = (PC1*15),
                                        yend = (PC2*15)), arrow = arrow(length = unit(1/2, "picas")),
               color = "black") +  annotate("text", x = (pca.loadings$PC1*15), y = (pca.loadings$PC2*15),
                                            label = pca.loadings$variables)
plotpca.bac
plotpca.bac <- plotpca.bac + xlab(pca1_var) + ylab(pca2_var)
plotpca.bac

#PERMANOVA to evaluate environmental differences between localities

perma.env <- adonis((matrix.pca) ~ Locality, 
                    data = matrix.bac, permutations=999, method = "euclidean")
perma.env

#Identifying variables that contribute to PC1 and PC2
#Checking  the amount of influence that each predictor variable has on each principal component
#For PC 1
loading_Scores_PC_1 <- pca.bac$rotation[,1]
fac_scores_PC_1 <- abs(loading_Scores_PC_1)
fac_scores_PC_1_ranked <- names(sort(fac_scores_PC_1,decreasing = T))
fac_scores_PC_1_ranked
#For PC 2
loading_Scores_PC_2 <- pca.bac$rotation[,2]
fac_scores_PC_2 <- abs(loading_Scores_PC_2)
fac_scores_PC_2_ranked <- names(sort(fac_scores_PC_2,decreasing = T))
fac_scores_PC_2_ranked

### MODELS
#Two-step approach selecting the variables that remained as predictors of microbial alpha diversity
#1) Pairwise Pearson correlations among selected variables to identify and discard those with a pairwise correlation higher than r>0.7
matrix.bac.models <- matrix.bac[5:32]
matrix.bac.models <- decostand(matrix.bac.models, method='standardize')

cor_matrix_bac <- cor(matrix.bac.models) # Correlation matrix
cor_matrix_bac_rm <- cor_matrix_bac                  # Modify correlation matrix
cor_matrix_bac_rm[upper.tri(cor_matrix_bac_rm)] <- 0
diag(cor_matrix_bac_rm) <- 0
matrix.bac.scores.new <- matrix.bac.models[ , !apply(cor_matrix_bac_rm,    # Remove highly correlated variables
                                                     2,
                                                     function(x) any(x > 0.7))]


#Fungi
matrix.its.models <- matrix.its[5:32]
matrix.its.models <- decostand(matrix.its.models, method='standardize')

cor_matrix_its <- cor(matrix.its.models) # Correlation matrix
cor_matrix_its_rm <- cor_matrix_its                 # Modify correlation matrix
cor_matrix_its_rm[upper.tri(cor_matrix_its_rm)] <- 0
diag(cor_matrix_its_rm) <- 0

matrix.its.scores.new <- matrix.its.models[ , !apply(cor_matrix_its_rm,    # Remove highly correlated variables
                                                     2,
                                                     function(x) any(x > 0.7))]

#2) The least-correlated variables, together with host species, were included in a stepwise forward and backward regression model to select variables with significant effects on observed ASVs

##Bacteria
matrix.bac.scores.new.mod <-matrix.bac.scores.new
matrix.bac.scores.new.mod$Species <-div.bac$Species
matrix.bac.scores.new.mod$observed <-div.bac$observed

#Full model
full.bac <- lm(observed ~., data = matrix.bac.scores.new.mod)
summary(full.bac)

# Stepwise regression model
step.model.bac <- stepAIC(full.bac, direction = "both")
summary(step.model.bac) 

## Fungi
matrix.its.scores.new.mod <-matrix.its.scores.new
matrix.its.scores.new.mod$Species <-div.its$Species
matrix.its.scores.new.mod$observed <-div.its$observed

#Full model
full.its <- lm(observed ~., data = matrix.its.scores.new.mod)
summary(full.its)

# Stepwise regression model
step.model.its <- stepAIC(full.its, direction = "both")
summary(step.model.its)

#Bacterial model including the variables that remained after stepwise selection
matrix.bac.scores.new.mod$Species<- div.bac$Species

bac.lm.core <- lm(observed ~ Bio_7+Pp+Tmin, data = matrix.bac.scores.new.mod)
summary(bac.lm.core)

##Fungal model including the variables that remained after stepwise selection

its.lm.core <- lm(observed ~ Temperature+Bio_6+Bio_7+Pp, data = matrix.its.scores.new.mod)
summary(its.lm.core)

#Distance-based redundancy analysis (dbRDA)
# Extract the OTU table
bac.veg <- as.data.frame(veganotu(phy.axo.core.bac))
its.veg <- as.data.frame(veganotu(phy.axo.core.its))

#Selecting the least-correlated variables previously chosen. Then, including Host species variable in the dataframe

matrix.rda.bac <- matrix.bac.scores.new
matrix.rda.bac$Species <- matrix.bac$Species

matrix.rda.its <- matrix.its.scores.new
matrix.rda.its$Species <- matrix.its$Species

# dbRDA models
#Bacteria
res.bac.core <- capscale(bac.veg ~ ., data=matrix.rda.bac,dist="bray")

#Fungi
res.its.core <- capscale(its.veg ~ ., data=matrix.rda.its,dist="bray")

#Stepwise selection
# Set up full and null models for 'ordistep'
# Full model

rda1.bac <- capscale(bac.veg ~ ., data=matrix.rda.bac)
rda1.its <- capscale(its.veg ~ ., data=matrix.rda.its)

# Intercept-only (null) model

rda0.bac <- capscale(bac.veg ~ 1, data=matrix.rda.bac)
rda0.its <- capscale(its.veg ~ 1, data=matrix.rda.its)
head(matrix.its.scores.new)

# Perform forward and backward selection of explanatory variables
#Bacteria
step.env.bac <- ordistep(rda0.bac, scope=formula(rda1.bac), direction='both')
#Fungi
step.env.its <- ordistep(rda0.its, scope=formula(rda1.its), direction='both')

# code to get variable names from 'ordistep' and 'envfit' results

vars.ordistep.bac <- gsub('^. ', '', rownames(step.env.bac$anova))
vars.envfit.bac <- names(which(vif.cca(res.bac.core) <= 10))
vars.bac <- unique(c(vars.ordistep.bac, vars.envfit.bac))

vars.ordistep.its <- gsub('^. ', '', rownames(step.env.its$anova))
vars.envfit.its <- names(which(vif.cca(res.its.core) <= 10))
vars.its <- unique(c(vars.ordistep.its, vars.envfit.its))

# select variables to keep from table 'Y'
matrix.bac.kept <- matrix.rda.bac[, vars.bac]
res.bac.kept <- capscale(bac.veg ~ ., data=matrix.bac.kept,dist="bray")

matrix.its.kept <- matrix.rda.its[, vars.its]
res.its.kept <- capscale(its.veg ~ ., data=matrix.its.kept,dist="bray")

# set up dataframes for plotting the results
sit <- as.data.frame(cbind(bac.veg, scores(res.bac.kept, display='sites')))
sit$Species <- phy.axo.core.bac@sam_data$Species
vec <- data.frame(scores(res.bac.kept, display='bp'))

sit.its <- as.data.frame(cbind(its.veg, scores(res.its.kept, display='sites')))
sit.its$Species <- phy.axo.core.its@sam_data$Species
vec.its <- data.frame(scores(res.its.kept, display='bp'))

# use these to adjust length of arrows and position of arrow labels
adj.vec <- 1.5
adj.txt <- 1.5

#  dbrda plot 
dbrda.bac <- ggplot(sit, aes(x=CAP1, y=CAP2, color=Species)) +
  geom_point(size=5, pch=21, colour="black", aes(fill=Species),stroke=1.5) +
  geom_segment(data=vec, inherit.aes=F, 
               mapping=aes(x=0, y=0, xend=adj.vec*CAP1, yend=adj.vec*CAP2), 
               arrow=arrow(length=unit(0.2, 'cm'))) + 
  geom_text(data=vec, inherit.aes=F, size=5, fontface="bold",
            mapping=aes(x=adj.txt*CAP1, y=adj.txt*CAP2, 
                        label=c('A. dumerilii','A. mexicanum','A. taylori','Water pH',"Water temperature", "Pp", "Tmin")))+
  theme_bw() + scale_fill_locuszoom() +
  theme(legend.text = element_text(size = 15, face = "bold.italic"),
        legend.title = element_text(size = 15,face = "bold"),
        axis.text.y=element_text(colour = "black", size = 14),
        axis.text.x=element_text(colour = "black", size = 15),
        axis.title.y=element_text(face = "bold",colour = "black", size = 15),
        axis.title.x=element_text(face = "bold",colour = "black", size = 15))
dbrda.bac$labels$x <- ord_labels(res.bac.kept)[1]
dbrda.bac$labels$y <- ord_labels(res.bac.kept)[2]
dbrda.bac

dbrda.its <- ggplot(sit.its, aes(x=CAP1, y=CAP2, color=Species)) +
  geom_point(size=5, pch=21, colour="black", aes(fill=Species),stroke=1.5) +
  geom_segment(data=vec.its, inherit.aes=F,
               mapping=aes(x=0, y=0, xend=adj.vec*CAP1, yend=adj.vec*CAP2), 
               arrow=arrow(length=unit(0.2, 'cm'))) + 
  geom_text(data=vec.its, inherit.aes=F, size=5,fontface="bold",
            mapping=aes(x=adj.txt*CAP1, y=adj.txt*CAP2, 
                        label=c('A. dumerilii','A. mexicanum', 'A. taylori','Water temperature',"Water pH","Pp")))+
  theme_bw() + scale_fill_locuszoom() +
  theme(legend.text = element_text(size = 15, face = "bold.italic"),
        legend.title = element_text(size = 15,face = "bold"),
        axis.text.y=element_text(colour = "black", size = 14),
        axis.text.x=element_text(colour = "black", size = 15),
        axis.title.y=element_text(face = "bold",colour = "black", size = 15),
        axis.title.x=element_text(face = "bold",colour = "black", size = 15))
dbrda.its$labels$x <- ord_labels(res.its.kept)[1]
dbrda.its$labels$y <- ord_labels(res.its.kept)[2]
dbrda.its

#PERMANOVAs

matrix.bac.perma <- matrix.bac.models
matrix.bac.perma$Species <- matrix.bac.scores.new.mod$Species
matrix.its.perma <- matrix.its.models
matrix.its.perma$Species <- matrix.its.scores.new.mod$Species

#PERMANOVA bacteria
perma.rda.bac<-anova.cca(res.bac.kept, by = "terms", permu=999)
perma.rda.bac
head(summary(res.bac.kept)) #46%

#PERMANOVA fungi
perma.rda.its <- anova.cca(res.its.kept, by = "terms", permu=999)
perma.rda.its
head(summary(res.its.kept)) #20%
