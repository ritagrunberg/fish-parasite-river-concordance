rm(list = ls())  # clear Rs brain 
###########################################################################################
library(readr)   # read csv file
library(tidyr)   # tidying data
library(dplyr)   # manipulating df
library(ggplot2) # graphing
library(ggpubr)  # nice package for figures
library(vegan)   # nmds 
###########################################################################################
# AIM TESTING FOR CONCORDANCE BETWEEN ALL PARASITES PRESENT IN ECOSYSTEM AND HOST + ENVIR VARIABLES 
# import dataset
river <- read_csv("C:/Users/grunberg/Dropbox/DISSERTATION_CHAPTERS/1_parasite_community_spatial/data/final_river_data.csv")
river <- river[,-c(38,57)] # remove notes that were in excel file... trash 
set.seed(1234567890)
river$plot <- as.factor(river$plot)
river$subplot <- as.factor(river$subplot)
################################################################################################
# ORDINATION PARASITE COMMUNITIES 
##########################################################################################
parasite <- river  %>% filter(river == "Passaic") %>% # filter(season=='fall') %>%
  group_by(river, season, plot, subplot) %>%
  summarise_all(funs(if(is.numeric(.)) sum(., na.rm = TRUE) else first(.))) %>%
  group_by(river, plot, season) %>%
  # summarise_if(.predicate = function(x) is.numeric(x), .funs = funs(mean="mean"))
  summarise_all(funs(if(is.numeric(.)) mean(., na.rm = TRUE) else first(.)))%>%
  mutate(dummy_parasite = rep(1)) # used to add observation without parasites 


# create community matrix and envi data matrix 
p_matrix<-parasite[,c(28:65 #,67
                      )] # makes new matrix with parasite data (columns 28-65)
p_envi<-parasite[,1:27] # makes new matrix with envi data (columns 1-27)

# remove columns 
p_matrix<-p_matrix[, colSums(abs(p_matrix)) !=0] # removes columns with all zeros

# remove rows without parasites
p01<-p_matrix[rowSums(abs(p_matrix))!=0,]  # remove rows with no parasites
envi01<-p_envi[rowSums(abs(p_matrix))!=0,]  # remove rows with no parasites 

#standardize values to unit maximum 
pt<-apply(p01,2, function(x) x/max(x)) 

# global nmds
ordO_p<-metaMDS(pt, distance="bray", trymax = 100, autotransform=FALSE) 
#ordO_p<-metaMDS(p01, distance="bray", trymax = 100, autotransform=FALSE) 
ordO_p

# PERMANOVA
p_para<-adonis(formula=pt ~ as.factor(plot), data=envi01,  permutations=999, method="bray") 
#p_para<-adonis(formula=p01 ~ as.factor(plot) * river, data=envi01,  permutations=999, method="bray") 
p_para

hcoordO<-as.data.frame(scores(ordO_p, display="sites"))#extracts coordinates for plot
pcoordO<-scores(ordO_p, display="species")#extracts coordinates for parasite vectors

parasite_species <- data.frame(parasite = c('Nematode sp1', 'Eocolis', 'C. margninatum', 'Eustronglyoides', 'Leptorhynchoides (la.)',
                      'Neoechinorhynchus (la)', 'Proteocephalus', 'Bothriocephalus', 'Plagioporus',
                      'Strigeoidae (mc)', 'Alloglossidium', 'Isoglaridacris', 'Fessesentis', 'C. cooperi',
                      'Caecincola', 'Acanthocephalus', 'N. cristatus', 'Triganodistomum', 'Echinorhynchus',
                      'Camallanus', "L. thecatus", 'Rhabdochona', 'P. minimum', 'Heterophyidae (mc)', 'U. ambloplitis',
                      'Philonema', 'Raphidascaris', 'Phyllodistomum', 'N. cylindratus'))
pcoordO <-cbind(pcoordO, parasite_species)

parasite_all <- bind_cols(envi01, hcoordO)

topp<-max(parasite_all[,28:29]) #determines maximum and minimum values for the plot axes
bott<-min(parasite_all[,28:29]) #I draw from both data sets because I make the graph sqaure and

parasite_all$plot <- as.factor(parasite_all$plot)
###########################################################################################
# make nmds graphic 
nmds_parasite <-ggplot(parasite_all, aes(x= NMDS1, y=NMDS2))+
  geom_hline(yintercept = 0, lty=2, color="grey") +
  geom_vline(xintercept = 0, lty=2, color="grey") +
  geom_point( size =2.75, aes(fill=plot, shape=season))+
  xlab("NMDS 1 (B-C)") +
  ylab("NMDS 2 (B-C)") +
  xlim(c( bott-0.2, topp+0.2)) +
  ylim(c( bott-0.2, topp+0.2)) +
  scale_shape_manual(values=c(21, 22,23, 24))+
  scale_fill_manual(values=c('#f7f7f7','#969696','#252525'),
                    labels = c("Upstream", "Midstream", "Downstream")) + 
  theme_pubr() +
  ggtitle("Parasite abundance")+
  theme(legend.text = element_text(size=8), legend.box = "horizontal",
        legend.title = element_text(size=9, face="bold"),
        legend.position=c(0.13,0.85)) +
  guides(colour = guide_legend(override.aes = list(size=70)))+
  theme(plot.title = element_text(hjust = 0.5, size=15, face="bold")) + # guides(shape=FALSE)+ 
  guides(fill=guide_legend("Site",  keywidth=0.12,
                           keyheight=0.12, override.aes = list(shape=21),
                           default.unit="inch"),
         shape= guide_legend("Season",  keywidth=0.12,
                             keyheight=0.12, 
                             default.unit="inch"))+
 # geom_text(data=pcoordO, aes(x=NMDS1,y=NMDS2, label=parasite), fontface =2, colour = "black" , size=1.5) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) # add black square panel around graphic

###########################################################################################
############### Jaccard ####################################################
###########################################################################################
parasite_pa<-apply(parasite[,c(28:65 #,67
                               )],2, function(x) {ifelse(x==0,0,1)})#makes presence-absence matrix
p_matrix_pa<-parasite_pa[, colSums(abs(parasite_pa)) !=0] # removes columns with all zeros
p01_pa<-p_matrix_pa[rowSums(abs(p_matrix_pa))!=0,] 
envi01_pa<-p_envi[rowSums(abs(p_matrix_pa))!=0,] 

ordO_p_pa<-metaMDS(p01_pa, distance="jaccard", trymax = 50, autotransform=FALSE)
ordO_p_pa

r_para_pa<-adonis(formula=p01_pa ~ as.factor(plot), data=envi01_pa,  permutations=999, method="jaccard")
r_para_pa

hcoordOpa<-as.data.frame(scores(ordO_p_pa, display="sites"))#extracts coordinates for plot
pcoordOpa<-scores(ordO_p_pa, display="species")#extracts coordinates for parasite vectors
pcoordOpa <-cbind(pcoordOpa, parasite_species)

parasite_all_pa <- bind_cols(envi01_pa, hcoordOpa)

topp<-max(parasite_all_pa[,28:29]) #determines maximum and minimum values for the plot axes
bott<-min(parasite_all_pa[,28:29]) #I draw from both data sets because I make the graph sqaure and

parasite_all_pa$plot <- as.factor(parasite_all_pa$plot)

# make nmds graphic 
nmds_parasite_PA <-ggplot(parasite_all_pa, aes(x= NMDS1, y=NMDS2))+
  geom_hline(yintercept = 0, lty=2, color="grey") +
  geom_vline(xintercept = 0, lty=2, color="grey") +
  geom_point( size =2.75, aes(fill=plot, shape=season))+
  xlab("NMDS 1 (JAC)") +
  ylab("NMDS 2 (JAC)") +
  xlim(c( bott-0.4, topp+0.4)) +
  ylim(c( bott-0.4, topp+0.4)) +
  scale_shape_manual(values=c(21, 22,23, 24))+
  scale_fill_manual(values=c('#f7f7f7','#969696','#252525'),
                    labels = c("Upstream", "Midstream", "Downstream")) + 
  theme_pubr() +
  ggtitle("Parasite composition")+
  theme(legend.text = element_text(size=8), legend.box = "horizontal",
        legend.title = element_text(size=9, face="bold"),
        legend.position=c(0.13,0.85)) +
  guides(colour = guide_legend(override.aes = list(size=70)))+
  theme(plot.title = element_text(hjust = 0.5, size=15, face="bold")) + # guides(shape=FALSE)+ 
  guides(fill=guide_legend("Site",  keywidth=0.12,
                           keyheight=0.12, override.aes = list(shape=21),
                           default.unit="inch"),
         shape= guide_legend("Season",  keywidth=0.12,
                             keyheight=0.12, 
                             default.unit="inch"))+
#  geom_text(data=pcoordOpa, aes(x=NMDS1,y=NMDS2, label=parasite), fontface =2, colour = "black" , size=1.5) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) # add black square panel around graphic

############################################################################
############################################################################
# ORDINATION HABITAT VARIABLES  
##########################################################################################
############################################################################
habitat01 <- river  %>% filter(river == "Passaic")%>% 
  group_by (river, season, plot, subplot) %>%
  summarise_all(funs(if(is.numeric(.)) mean(., na.rm = TRUE) else first(.))) %>%
  group_by(river, plot, season) %>%
  summarise_all(funs(if(is.numeric(.)) mean(., na.rm = TRUE) else first(.))) #%>%
# filter(!(river == "Passaic" & season == "summer" & plot ==1))   # remove summer plot 1 observation becaues no parasite data for that 

# return all rows from DF1 where there are matching values in DF2,
# keeping just columns from DF1.
# habitat01 <- semi_join(habitat01, parasite_all, by=c("river", "season", "plot", "subplot"))

#all habitat data
hab<-habitat01[,6:16] 
#just physical data
hab_phys<-habitat01[,c(6,8,15,16)] 
#just chemical data 
hab_chem<-habitat01[,10:14] 

#information for ordinations 
hab_2<-habitat01[,c(1:4)] 

#without predictor variable rda is a PCA 
ordO_hab<-rda(hab, scale=TRUE)# pca in vegan 
ordO_hab_phy<-rda(habitat01[,c(6,8,15,16)] , scale=TRUE)# pca in vegan 
ordO_hab_chem<-rda(habitat01[,10:14] , scale=TRUE)# pca in vegan 

#biplot(ordO_hab)
biplot(ordO_hab_phy)
biplot(ordO_hab_chem)
summary(ordO_hab_chem)
summary(ordO_hab_phy)
###########################################################################################
# create PCA plot for water quality variables 
sitecoordO_habchem<-as.data.frame(scores(ordO_hab_chem, display="sites"))#extracts coordinates for plots
habitat_chem <- bind_cols(hab_2, sitecoordO_habchem)

topp<-max(habitat_chem[,5:6]) 
bott<-min(habitat_chem[,5:6]) 

habitat_chem$plot <- as.factor(habitat_chem$plot)

nmds_habitat <-ggplot(habitat_chem, aes(x= PC1, y=PC2))+
  geom_hline(yintercept = 0, lty=2, color="grey") +
  geom_vline(xintercept = 0, lty=2, color="grey") +
  geom_point( size =2.75, aes(fill=plot, shape=season))+
  xlab("PC 1") +
  ylab("PC 2") +
  xlim(c( bott, topp)) +
  ylim(c( bott, topp)) +
  scale_shape_manual(values=c(21, 22,23, 24))+
  scale_fill_manual(values=c('#f7f7f7','#969696','#252525'),
                    labels = c("Upstream", "Midstream", "Downstream")) + 
  theme_pubr() +
  ggtitle("Water quality")+
  theme(legend.text = element_text(size=8), legend.box = "horizontal",
        legend.title = element_text(size=9, face="bold"),
        legend.position=c(0.13,0.85)) +
  guides(colour = guide_legend(override.aes = list(size=70)))+
  theme(plot.title = element_text(hjust = 0.5, size=15, face="bold")) + # guides(shape=FALSE)+ 
  guides(fill=guide_legend("Site",  keywidth=0.12,
                           keyheight=0.12, override.aes = list(shape=21),
                           default.unit="inch"),
         shape= guide_legend("Season",  keywidth=0.12,
                             keyheight=0.12, 
                             default.unit="inch"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))  # add black square panel around graphic
###########################################################################################
# create PCA plot for physical stream variables 
sitecoordO_habphy<-as.data.frame(scores(ordO_hab_phy, display="sites"))#extracts coordinates for plots
habitat_phy <- bind_cols(hab_2, sitecoordO_habphy)

topp<-max(habitat_phy[,5:6]) 
bott<-min(habitat_phy[,5:6]) 

habitat_phy$plot <- as.factor(habitat_phy$plot)

nmds_phy <-ggplot(habitat_phy, aes(x= PC1, y=PC2))+
  geom_hline(yintercept = 0, lty=2, color="grey") +
  geom_vline(xintercept = 0, lty=2, color="grey") +
  xlab("PC 1") +
  ylab("PC 2") +
  geom_point( size =2.75, aes(fill=plot, shape=season))+
  xlim(c( bott, topp)) +
  ylim(c( bott, topp)) +
  scale_shape_manual(values=c(21, 22,23, 24))+
  scale_fill_manual(values=c('#f7f7f7','#969696','#252525'),
                    labels = c("Upstream", "Midstream", "Downstream")) + 
  theme_pubr() +
  ggtitle("Physical")+
  theme(legend.text = element_text(size=8), legend.box = "horizontal",
        legend.title = element_text(size=9, face="bold"),
        legend.position=c(0.13,0.85)) +
  guides(colour = guide_legend(override.aes = list(size=70)))+
  theme(plot.title = element_text(hjust = 0.5, size=15, face="bold")) + # guides(shape=FALSE)+ 
  guides(fill=guide_legend("Site",  keywidth=0.12,
                           keyheight=0.12, override.aes = list(shape=21),
                           default.unit="inch"),
         shape= guide_legend("Season",  keywidth=0.12,
                             keyheight=0.12, 
                             default.unit="inch"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))  # add black square panel around graphic

################################################################################################
################################################################################################
# h4: host community abundance within plot 
##########################################################################################
################################################################################################
# abundance based FISH MATRIX 
fish_01 <- river %>% filter (river == "Passaic") %>% 
  mutate(abun = rep(1)) %>% group_by(river, plot, season, subplot, host_species) %>%
  summarise_all(funs(if(is.numeric(.)) sum(., na.rm = TRUE) else first(.))) %>% 
  spread(host_species, abun)
fish_01[,27:ncol((fish_01))][is.na(fish_01[,27:ncol((fish_01))])] <- 0
fish <- fish_01 %>%
  group_by(river, plot, subplot, season) %>%
  summarise_all(funs(if(is.numeric(.)) sum(., na.rm = TRUE) else first(.))) %>%
  group_by(river, plot, season) %>%  
  # filter(!(river == "Passaic" & season == "summer" & plot ==1)) %>%  # remove summer plot 1 observation becaues no parasite data for that 
  summarise_all(funs(if(is.numeric(.)) mean(., na.rm = TRUE) else first(.)))
fish <- fish[,-c(27:65,88)] # remove parasite data

# return all rows from DF1 where there are matching values in DF2,
# keeping just columns from DF1.
#fish <- semi_join(fish, parasite_all, by=c("river", "season", "plot", "subplot"))

# abundance based matrix 
fish_matrix<-fish[,27:ncol((fish))] # makes new matrix with fish data (columns 18-41)
fish_envi<-fish[,1:26] # makes new matrix with envi data (columns 1-17)

fish_matrix<-fish_matrix[, colSums(abs(fish_matrix)) !=0] # removes columns with all zeros

fish_01<-fish_matrix[rowSums(abs(fish_matrix))!=0,] 
fish_envi01<-fish_envi[rowSums(abs(fish_matrix))!=0,] 

# transform the abundance data here; to allow for a uniform comparison we standardize abundance values for each parasite species maximum 
fish_t<-apply(fish_01,2, function(x) x/max(x))#standardizes columns to unit maximum
#second argument =1 to apply to rows and =2 to apply to columns

#Bray-Curtis NMDS
ordO<-metaMDS(fish_t, distance="bray", trymax = 100, autotransform=FALSE)#global NMDS
ordO
#ordO<-metaMDS(fish_01, distance="bray", trymax = 100, autotransform=FALSE)#global NMDS

stressplot(ordO)

f_abun<-adonis(formula=fish_t ~ as.factor(plot) , data=fish_envi01,  permutations=999, method="bray")
#f_abun<-adonis(formula=fish_01 ~ as.factor(plot) * river , data=fish_envi01,  permutations=999, method="bray")
f_abun

# pairwise PERMANOVAs
# pairwise.adonis(fish_t, fish_envi01$plot)

sitecoordO<-as.data.frame(scores(ordO, display="sites"))#extracts coordinates for hosts
fishcoordO<-scores(ordO, display="species")#extracts coordinates for parasite vectors

fish_species <- data.frame(fish = c('A. rupestris', 'A. melas', 'C. commersonii', 'C. spiloptera', 'E. niger', 'E. olmstedi',
                                    'F. diaphanus', 'G. holbrooki', 'L. auritus', 'L. gibbosus', 'L. macrochirus',
                                    'L. cornutus', 'M. dolomieu', 'N. bifrenatus', 'N procne', 'N. gyrinus', 'R. atratulus',
                                    'R. cataractae', 'S. trutta', 'S. atromaculatus', 'S. corporalis', 'U. pygmaea'))
fishcoordO <-cbind(fishcoordO, fish_species)

fish_all <- bind_cols(fish_envi01, sitecoordO)

topp<-max(fish_all[,27:28]) #determines maximum and minimum values for the plot axes
bott<-min(fish_all[,27:28]) #I draw from both data sets because I make the graph sqaure and
#the scale of the X and Y axes should be the same on the graph for proper visual interpretation

fish_all$plot <- as.factor(fish_all$plot)
###########################################################################################
nmds_fish_abun <-ggplot(fish_all, aes(x= NMDS1, y=NMDS2))+
  geom_hline(yintercept = 0, lty=2, color="grey") +
  geom_vline(xintercept = 0, lty=2, color="grey") +
  geom_point( size =2.75, aes(fill=plot, shape=season))+
  xlab("NMDS 1 (B-C)") +
  ylab("NMDS 2 (B-C)") +
  xlim(c( bott-0.2, topp+0.2)) +
  ylim(c( bott-0.2, topp+0.2)) +
  scale_fill_manual(values=c('#f7f7f7','#969696','#252525'),
                    labels = c("Upstream", "Midstream", "Downstream")) + 
  theme_pubr() +
  scale_shape_manual(values=c(21, 22,23, 24))+
  ggtitle("Fish abundance")+
  theme(legend.text = element_text(size=8), legend.box = "horizontal",
        legend.title = element_text(size=9, face="bold"),
        legend.position=c(0.13,0.85)) +
  guides(colour = guide_legend(override.aes = list(size=70)))+
  theme(plot.title = element_text(hjust = 0.5, size=15, face="bold")) + # guides(shape=FALSE)+ 
  guides(fill=guide_legend("Site",  keywidth=0.12,
                           keyheight=0.12, override.aes = list(shape=21),
                           default.unit="inch"),
         shape= guide_legend("Season",  keywidth=0.12,
                             keyheight=0.12, 
                             default.unit="inch"))+
 # geom_text(data=fishcoordO, aes(x=NMDS1,y=NMDS2, label=fish), fontface =2, colour = "black" , size=1.5) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))  # add black square panel around graphic

##########################################################################################
############### Jaccard ####################################################
##########################################################################################
fish_pa<-apply(fish[,27:ncol(fish)] ,2, function(x) {ifelse(x==0,0,1)})#makes presence-absence matrix
fish_01pa<-fish_pa[rowSums(abs(fish_pa))!=0,] 
fish_envi01pa<-fish_envi[rowSums(abs(fish_pa))!=0,] 

#Bray-Curtis NMDS
ordOpa<-metaMDS(fish_01pa, distance="jaccard", trymax = 100, autotransform=FALSE)#global NMDS
stressplot(ordOpa)
ordOpa

r_fishpa<-adonis(formula=fish_01pa ~ as.factor(plot) , data=fish_envi01pa,  permutations=999, method="bray")
r_fishpa

sitecoordO<-as.data.frame(scores(ordOpa, display="sites"))#extracts coordinates for hosts
fishcoordOpa<-scores(ordOpa, display="species")#extracts coordinates for parasite vectors
fishcoordOpa <-cbind(fishcoordOpa, fish_species)

fish_all_pa <- bind_cols(fish_envi01pa, sitecoordO)

topp<-max(fish_all_pa[,27:28]) #determines maximum and minimum values for the plot axes
bott<-min(fish_all_pa[,27:28]) #I draw from both data sets because I make the graph sqaure and
#the scale of the X and Y axes should be the same on the graph for proper visual interpretation

fish_all_pa$plot <- as.factor(fish_all_pa$plot)
###########################################################################################
nmds_fish_pa <-ggplot(fish_all_pa, aes(x= NMDS1, y=NMDS2))+
  geom_hline(yintercept = 0, lty=2, color="grey") +
  geom_vline(xintercept = 0, lty=2, color="grey") +
  geom_point( size =2.75, aes(fill=plot, shape=season))+
  xlab("NMDS 1 (JAC)") +
  ylab("NMDS 2 (JAC)") +
  xlim(c( bott-0.4, topp+0.4)) +
  ylim(c( bott-0.4, topp+0.4)) +
  scale_shape_manual(values=c(21, 22,23, 24))+
  scale_fill_manual(values=c('#f7f7f7','#969696','#252525'),
                    labels = c("Upstream", "Midstream", "Downstream")) + 
  theme_pubr() +
  ggtitle("Fish composition")+
  theme(legend.text = element_text(size=8), legend.box = "horizontal",
        legend.title = element_text(size=9, face="bold"),
        legend.position=c(0.13,0.85)) +
  guides(colour = guide_legend(override.aes = list(size=70)))+
  theme(plot.title = element_text(hjust = 0.5, size=15, face="bold")) + # guides(shape=FALSE)+ 
  guides(fill=guide_legend("Site",  keywidth=0.12,
                           keyheight=0.12, override.aes = list(shape=21),
                           default.unit="inch"),
         shape= guide_legend("Season",  keywidth=0.12,
                             keyheight=0.12, 
                             default.unit="inch"))+
 # geom_text(data=fishcoordOpa, aes(x=NMDS1,y=NMDS2, label=fish), fontface =2, colour = "black" , size=1.5) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

################################################################################################
# h5: host community biomass
##########################################################################################
fish_02 <- river %>% filter (river == "Passaic") %>% 
  mutate(abun = rep(1)) %>% 
  group_by(river, plot, season, subplot, host_species) %>%
  summarise_all(funs(if(is.numeric(.)) sum(., na.rm = TRUE) else first(.))) %>% spread(host_species, dry_weight_g)
fish_02[,27:ncol(fish_02)][is.na(fish_02[,27:ncol(fish_02)])] <- 0
fish2 <- fish_02 %>%
  group_by(river,season, plot, subplot) %>%
  summarise_all(funs(if(is.numeric(.)) sum(., na.rm = TRUE) else first(.))) %>%
  group_by(river,plot,season) %>%  
  # filter(!(river == "Passaic" & season == "summer" & plot ==1)) %>%  # remove summer plot 1 observation becaues no parasite data for that 
  summarise_all(funs(if(is.numeric(.)) mean(., na.rm = TRUE) else first(.))) 
fish2 <- fish2[,-c(27:65, 88)]  # remove parasite data

# return all rows from DF1 where there are matching values in DF2,
# keeping just columns from DF1.
#fish_02 <- semi_join(fish2, parasite_all, by=c("river", "season", "plot", "subplot"))

# abundance based matrix 
fish_matrix2<-fish2[,27:ncol(fish2)] # makes new matrix with fish data (columns 18-41)
fish_envi2<-fish2[,1:26] # makes new matrix with envi data (columns 1-17)

fish_matrix2<-fish_matrix2[, colSums(abs(fish_matrix2)) !=0] # removes columns with all zeros

fish_001<-fish_matrix2[rowSums(abs(fish_matrix2))!=0,] 
fish_envi001<-fish_envi2[rowSums(abs(fish_matrix2))!=0,] 

#Bray-Curtis NMDS
ordO_02<-metaMDS(fish_001, distance="bray", trymax = 100, autotransform=FALSE)#global NMDS

stressplot(ordO)
ordO_02

#PERMANOVA
#p_biomass<-adonis(formula=fish_t0 ~ as.factor(plot), data=fish_envi001,  permutations=999, method="bray")
p_biomass<-adonis(formula=fish_001 ~ as.factor(plot), data=fish_envi001,  permutations=999, method="bray")
p_biomass

# pairwise PERMANOVAs
# pairwise.adonis(fish_t0, fish_envi001$plot)

sitecoordO2<-as.data.frame(scores(ordO_02, display="sites"))#extracts coordinates for hosts
fishcoordO2<-scores(ordO_02, display="species")#extracts coordinates for parasite vectors

fish_all_02 <- bind_cols(fish_envi001, sitecoordO2)

topp<-max(fish_all_02[,27:28]) #determines maximum and minimum values for the plot axes
bott<-min(fish_all_02[,27:28]) #I draw from both data sets because I make the graph sqaure and
#the scale of the X and Y axes should be the same on the graph for proper visual interpretation

fish_all_02$plot <- as.factor(fish_all_02$plot)
###########################################################################################
nmds_biomass <-ggplot(fish_all_02, aes(x= NMDS1, y=NMDS2))+
  geom_hline(yintercept = 0, lty=2, color="grey") +
  geom_vline(xintercept = 0, lty=2, color="grey") +
  geom_point( size =2.75, aes(fill=plot, shape=season))+
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  xlim(c( bott, topp)) +
  ylim(c( bott, topp)) +
  scale_shape_manual(values=c(21, 22,23, 24))+
  scale_fill_manual(values=c('#f7f7f7','#969696','#252525'),
                    labels = c("Upstream", "Midstream", "Downstream")) + 
  theme_pubr() +
  ggtitle("Fish Biomass")+
  theme(legend.text = element_text(size=12), # legend.box = "horizontal",
        legend.title = element_text(size=12, face="bold") #,
        #legend.position=c(0.13,0.85)
        ) +
  guides(colour = guide_legend(override.aes = list(size=70)))+
  theme(plot.title = element_text(hjust = 0.5, size=15, face="bold")) + # guides(shape=FALSE)+ 
  guides(fill=guide_legend("Site",  keywidth=0.12,
                           keyheight=0.12, override.aes = list(shape=21, size=5),
                           default.unit="inch"),
         shape= FALSE) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))  # add black square panel around graphic


nmds_biomass2 <-ggplot(fish_all_02, aes(x= NMDS1, y=NMDS2))+
  geom_hline(yintercept = 0, lty=2, color="grey") +
  geom_vline(xintercept = 0, lty=2, color="grey") +
  geom_point( size =2.75, aes(fill=plot, shape=season))+
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  xlim(c( bott, topp)) +
  ylim(c( bott, topp)) +
  scale_shape_manual(values=c(21, 22,23, 24))+
  scale_fill_manual(values=c('#f7f7f7','#969696','#252525'),
                    labels = c("Upstream", "Midstream", "Downstream")) + 
  theme_pubr() +
  ggtitle("Fish Biomass")+
  theme(legend.text = element_text(size=12), # legend.box = "horizontal",
        legend.title = element_text(size=12, face="bold") #,
        #legend.position=c(0.13,0.85)
  ) +
  guides(colour = guide_legend(override.aes = list(size=70)))+
  theme(plot.title = element_text(hjust = 0.5, size=15, face="bold")) + # guides(shape=FALSE)+ 
  guides(fill=FALSE,
         shape= guide_legend("Season",  keywidth=0.12,
                             keyheight=0.12,  override.aes = list(size=5),
                             default.unit="inch")) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))  # add black square panel around graphic

#############################################################################################
# export graphic of all ordinations 
#############################################################################################
get_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
legend.passaic <- get_legend(nmds_biomass)
legend.passaic2 <- get_legend(nmds_biomass2)

library(gridExtra)
library(grid)

jpeg(file="Passaic_River_ordinations.jpeg", width = 180, height = 270, units = 'mm', res = 300)
grid.arrange(arrangeGrob(nmds_habitat +theme(legend.position="none"),  
                          nmds_phy+theme(legend.position="none"),   
                          nmds_fish_abun+theme(legend.position="none"),  
                          nmds_fish_pa +theme(legend.position="none"),  
                          nmds_parasite +theme(legend.position="none"),   
                          nmds_parasite_PA + theme(legend.position="none"), 
                          nrow=3, ncol=2), 
             legend.passaic, legend.passaic2,
                          nrow=3,heights=c(10, 0.4, 0.4))

dev.off()


jpeg(file="Passaic_River_ordinations_biomass.jpeg", width = 120, height = 120, units = 'mm', res = 300)
nmds_biomass
dev.off()
#############################################################################################
# pairwise procrustes and protests
###########################################################################################
# host biomass and abundance
pro_3 <-(procrustes(ordO, ordO_02))
#jpeg(file="Passaic_procrustes_host_abundance_biomass.jpeg", width = 180, height = 180, units = 'mm', res = 300)
#plot(pro_3, kind="1")
#dev.off()
prot_hb_ha <- protest(ordO, ordO_02, permutations = 999)
prot_hb_ha
c<- c(prot_hb_ha$ss, prot_hb_ha$signif, "fish_abundance", "fish_biomass")

# host abundance and paraiste 
pro_4 <- procrustes(ordO, ordO_p)
#plot(pro_1, kind="1")
prot_ha_p <- protest(ordO, ordO_p, permutations = 999)
prot_ha_p
res_host_parasite <-data.frame(residuals(prot_ha_p))
res_host_parasite<-data.frame(res_host_parasite, envi01)# pointwise residuals 

season_order <- data.frame (order_season= c(1,3,4,2), season = c("fall", "spring", "summer", "winter"))
res_host_parasite <- merge(res_host_parasite, season_order, by= "season")

write.csv(res_host_parasite, "residuals_host_parasite_abun_passaic.csv")

res_host_parasite %>% ggplot()+ geom_point(aes(x=season, y = residuals.prot_ha_p., fill=plot), pch=21, size=5)+
  scale_fill_manual(values=c('#f7f7f7','#969696','#252525'),
                    labels = c("Upstream", "Midstream", "Downstream")) +
  theme_pubr() +labs(x="season", y="Procrustes residuals")
#summary(prot_ha_p)
#plot(prot_ha_p)
#plot(prot_ha_p, kind=2)
d<- c(prot_ha_p$ss, prot_ha_p$signif, "parasite", "fish_abundance")

# host biomass and parasite
pro_5 <- (procrustes(ordO_02, ordO_p))
#plot(pro_5, kind="1")
prot_hb_p <- protest(ordO_02, ordO_p, permutations = 999)
prot_hb_p
e<- c(prot_hb_p$ss, prot_hb_p$signif, "parasite", "fish_biomass")

#############################################
# phys and chem ordinations 
#########################################################
# parasite and chemical 
prot_env_pchem <- protest(ordO_hab_chem, ordO_p, permutations = 999)
prot_env_pchem$t
g<- c(prot_env_pchem$ss, prot_env_pchem$signif, "parasite", "chemical")


#parasite and physical 
prot_env_pphy <- protest(ordO_hab_phy, ordO_p, permutations = 999)
prot_env_pphy
h<- c(prot_env_pphy$ss, prot_env_pphy$signif, "parasite", "physical")

# host abundance and chemical
pro_12 <- procrustes(ordO, ordO_hab_chem)
#plot(pro_1, kind="1")
prot_ha_envchem <- protest(ordO, ordO_hab_chem, permutations = 999)
i <- c(prot_ha_envchem$ss, prot_ha_envchem$signif, "fish_abundance", "chemical")

# host abundance and phys
pro_13 <- procrustes(ordO, ordO_hab_phy)
#plot(pro_1, kind="1")
prot_ha_envphy <- protest(ordO, ordO_hab_phy, permutations = 999)
j <- c(prot_ha_envphy$ss, prot_ha_envphy$signif, "fish_abundance", "physical")

# host biomass and chemical
pro_14 <- procrustes(ordO_02, ordO_hab_chem)
#plot(pro_1, kind="1")
prot_bio_envchem <- protest(ordO_02,ordO_hab_chem, permutations = 999)
k <- c(prot_bio_envchem$ss, prot_bio_envchem$signif, "fish_biomass", "chemical")

# host biomass and phys
pro_15 <- procrustes(ordO_02, ordO_hab_phy)
#plot(pro_1, kind="1")
prot_bio_envphy <- protest(ordO_02,ordO_hab_phy, permutations = 999)
l <- c(prot_bio_envphy$ss, prot_bio_envphy$signif,"fish_biomass", "physical")
#########################################################
# adding PA data
#########################################################
# parasite pa and host pa 
pro_16 <- procrustes(ordOpa, ordO_p_pa)
proc_out <- data.frame(NMDS1=pro_16$Yrot[,1], # parasites
                       NMDS2=pro_16$Yrot[,2], # parasites
                       group1 = rep("parasite"),
                       xNMDS1=pro_16$X[,1], # hosts
                       xNMDS2=pro_16$X[,2],
                       group2 = rep("host")) #hosts 
pro_16$rotation
plot(pro_16)

jpeg(file="procrustes_error_host_parasite_pa_PASSAIC.jpeg", width = 90, height = 90, units = 'mm', res = 300)
ggplot(proc_out)+
  geom_hline(yintercept = 0, lty=2, color="grey") +
  geom_vline(xintercept = 0, lty=2, color="grey") +
  geom_point(aes(x=NMDS1, y=NMDS2, fill=group2), # fill="black", 
             pch=21, size=1.5, show.legend  = TRUE)+ # parasite
  geom_point(aes(x=xNMDS1, y=xNMDS2, fill=group1),  #fill="white",
             pch=21, size=1.5, show.legend  = TRUE) + # host 
  geom_segment(aes(x=NMDS1,y=NMDS2,xend=xNMDS1,yend=xNMDS2),
               arrow=arrow(length=unit(0.2,"cm")), color="#1c9099", lwd=0.5)+ theme_pubr()+
  geom_segment(aes(x=pro_16$rotation[1,2],y=-1,
                   xend=-pro_16$rotation[1,2],yend=1), 
               lwd=0.25)+
  geom_segment(aes(x=-1,y=-pro_16$rotation[1,2],
                   xend=1,yend=pro_16$rotation[1,2]), 
               lwd=0.25)+
  labs(x="Axis 1", y = "Axis 2")+
  scale_fill_manual("", values=c('black','white'),
                    labels=c("parasite", "host"), guide='legend')  + 
  #  ylim(c(-0.75,0.75))+ xlim(c(-0.75,0.75))+
  theme(legend.text = element_text(size=8), legend.box = "horizontal",
        legend.title = element_text(size=0, face="bold"),
        legend.position=c(0.15,0.925)) +
  guides(fill=guide_legend(keywidth=0.12,
                           keyheight=0.12, override.aes = list(shape=21),
                           default.unit="inch"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()

plot(pro_16, kind="1")
prot_pp_hopa <- protest(ordOpa, ordO_p_pa, permutations = 999)
m <- c(prot_pp_hopa$ss, prot_pp_hopa$signif, "parasite_pa", "fish_pa")

res_host_parasite_pa <- residuals(prot_pp_hopa)
res_host_parasite_pa<-data.frame(res_host_parasite_pa, envi01)# pointwise residuals 
season_order <- data.frame (order_season= c(1,3,4,2), season = c("fall", "spring", "summer", "winter"))
res_host_parasite_pa <- merge(res_host_parasite_pa, season_order, by= "season")

write.csv(res_host_parasite_pa, "residuals_host_parasite_pa_passaic.csv")

jpeg(file="Passaic_host_parasite_PA_residuals.jpeg", width = 90, height = 90, units = 'mm', res = 300)
res_host_parasite_pa %>% ggplot()+ geom_point(aes(x=order_season, y = res_host_parasite_pa, fill=plot), pch=21, size=4)+
  scale_fill_manual(values=c('#f7f7f7','#969696','#252525'),
                    labels = c("Upstream", "Midstream", "Downstream")) +
  theme_pubr() +labs(x="", y="Procrustes residuals") +
  guides(fill=guide_legend("",  keywidth=0.12,
                           keyheight=0.12, override.aes = list(shape=21),
                           default.unit="inch"))+
  scale_x_discrete(name ="", limits=c("fall", "winter", "spring", "summer"))
dev.off()

# parasite pa and host abun
pro_17 <- procrustes(ordO, ordO_p_pa)
#plot(pro_1, kind="1")
prot_par_host <- protest(ordO_p_pa, ordO, permutations = 999)
n <- c(prot_par_host$ss, prot_par_host$signif, "parasite_pa", "fish_abundance")

#parasite pa and host biomass
pro_18<- procrustes(ordO_02, ordO_p_pa)
#plot(pro_1, kind="1")
prot_par_host_bio <- protest(ordO_p_pa, ordO_02, permutations = 999)
o <- c(prot_par_host_bio$ss, prot_par_host_bio$signif, "parasite_pa", "fish_biomass")

#parasite pa and chemical envi
pro_19<- procrustes(ordO_hab_chem, ordO_p_pa)
#plot(pro_1, kind="1")
prot_par_chem <- protest(ordO_hab_chem, ordO_p_pa, permutations = 999)
p <- c(prot_par_chem$ss, prot_par_chem$signif, "parasite_pa", "chemical")

#parasite pa and phys env
pro_20<- procrustes(ordO_hab_phy, ordO_p_pa)
#plot(pro_1, kind="1")
prot_par_phys <- protest(ordO_hab_phy, ordO_p_pa, permutations = 999)
q <- c(prot_par_phys$ss, prot_par_phys$signif, "parasite_pa", "physical")

#host pa and phys env
pro_21<- procrustes(ordO_hab_phy, ordOpa)
#plot(pro_1, kind="1")
prot_host_phys <- protest(ordO_hab_phy, ordOpa, permutations = 999)
r <- c(prot_host_phys$ss, prot_host_phys$signif, "fish_pa", "physical")

res_host_phys_pa <- residuals(prot_host_phys)
res_host_phys_pa<-data.frame(res_host_phys_pa, envi01)# pointwise residuals 
season_order <- data.frame (order_season= c(1,3,4,2), season = c("fall", "spring", "summer", "winter"))
res_host_phys_pa <- merge(res_host_phys_pa, season_order, by= "season")

write.csv(res_host_phys_pa, "residuals_host_phys_pa_passaic.csv")

#host pa and chem env
pro_22<- procrustes(ordO_hab_chem, ordOpa)
#plot(pro_1, kind="1")
prot_host_chem <- protest(ordO_hab_chem, ordOpa, permutations = 999)
s <- c(prot_host_chem$ss, prot_host_chem$signif, "fish_pa", "chemical")

# host pa and host abun
pro_23 <- procrustes(ordO, ordOpa)
#plot(pro_1, kind="1")
prot_pa_host <- protest(ordO, ordOpa, permutations = 999)
t <- c(prot_pa_host$ss, prot_pa_host$signif, "fish_pa", "fish_abundance")

#host pa and host biomass
pro_24<- procrustes(ordO_02, ordOpa)
#plot(pro_1, kind="1")
prot_pa_bio <- protest(ordO_02, ordOpa, permutations = 999)
u <- c(prot_pa_bio$ss, prot_pa_bio$signif, "fish_pa", "fish_biomass")

#host pa and parasite density 
pro_25<- procrustes(ordOpa, ordO_p)
#plot(pro_1, kind="1")
prot_pa_paras <- protest(ordOpa, ordO_p, permutations = 999)
plot(prot_pa_paras)
v <- c(prot_pa_paras$ss, prot_pa_paras$signif, "parasite", "fish_pa")

res_host_parasite <- residuals(prot_pa_paras)
res_host_parasite<-data.frame(res_host_parasite, envi01)# pointwise residuals 
season_order <- data.frame (order_season= c(1,3,4,2), season = c("fall", "spring", "summer", "winter"))
res_host_parasite <- merge(res_host_parasite, season_order, by= "season")
write.csv(res_host_parasite, "residuals_hostpa_parasite_passaic.csv")

#parasite pa and parasite density 
pro_26<- procrustes(ordO_p_pa, ordO_p)
#plot(pro_1, kind="1")
prot_pa_paras2 <- protest(ordO_p_pa, ordO_p, permutations = 999)
w <- c(prot_pa_paras2$ss, prot_pa_paras2$signif, "parasite", "parasite_pa")

#habitat and chemical
pro_27<- procrustes(ordO_hab_chem, ordO_hab_phy)
#plot(pro_1, kind="1")
prot_chem_phys <- protest(ordO_hab_chem, ordO_hab_phy, permutations = 999)
x <- c(prot_chem_phys$ss, prot_chem_phys$signif, "physical", "chemical")

################################################################################################################
# put PROTEST CONCORDANCE in one df 
################################################################################################################
passaic_river_protest<-rbind(c,d,e,
                           g,h,
                           i,j,k,l, 
                           m, n, o, p, q,
                           r, s, t, u, 
                           v, 
                           w, x)
passaic_river_protest <- as.data.frame(passaic_river_protest)
passaic_river_protest <- passaic_river_protest %>% dplyr::select(1:4) 
meh <-c('ms_square', "p_value", "var1", "var2")
colnames(passaic_river_protest) <- meh
passaic_river_protest <- passaic_river_protest %>% mutate(river =rep("Passaic"))

write.csv(passaic_river_protest, "passaic_river_protest.csv")

###########################################################################################
# put PERMANOVA in one df 
################################################################################################################
parasite_permanova <- c(p_para$aov.tab[1,1], p_para$aov.tab[3,1], p_para$aov.tab[1,6], p_para$aov.tab[1,4], p_para$aov.tab[1,5])
parasite_pa_permanova <- c(r_para_pa$aov.tab[1,1], r_para_pa$aov.tab[3,1], r_para_pa$aov.tab[1,6], r_para_pa$aov.tab[1,4], r_para_pa$aov.tab[1,5])
fish_numerical_permanova <- c(p_abun$aov.tab[1,1], p_abun$aov.tab[3,1], p_abun$aov.tab[1,6], p_abun$aov.tab[1,4], p_abun$aov.tab[1,5])
fish_biomass_permanova <- c(p_biomass$aov.tab[1,1], p_biomass$aov.tab[3,1], p_biomass$aov.tab[1,6], p_biomass$aov.tab[1,4], p_biomass$aov.tab[1,5])
fish_pa_permanova <- c(r_fishpa$aov.tab[1,1], r_fishpa$aov.tab[3,1], r_fishpa$aov.tab[1,6], r_fishpa$aov.tab[1,4], r_fishpa$aov.tab[1,5])

Passaic_permanova <-rbind(parasite_permanova, parasite_pa_permanova, 
                          fish_numerical_permanova, fish_biomass_permanova, fish_pa_permanova)
thing <-c('DF', 'DF2', 'P_value', 'F_model', 'R2')
yeah <-c('parasite_numerical',"parasite_pa", "fish_numerical", "fish_biomass", "fish_pa")
colnames(Passaic_permanova) <- thing
rownames(Passaic_permanova) <- yeah

write.csv(Passaic_permanova, "passaic_river_PERMANOVA.csv")
