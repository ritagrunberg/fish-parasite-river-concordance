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
###########################################################################################
# function for pariwise adonis with bonferroni correction
pairwise.adonis <- function(x,factors, sim.method = 'bray', p.adjust.m
                            ='bonferroni')
{
  library(vegan)
  co = combn(unique(factors),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  for(elem in 1:ncol(co)){
    ad = adonis(x[factors %in% c(co[1,elem],co[2,elem]),] ~ factors[factors
                                                                    %in% c(co[1,elem],co[2,elem])] , method =sim.method);
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}

################################################################################################
# h2: parasite within plots 
##########################################################################################

parasite <- river  %>% filter(river == "Raritan") %>% 
  group_by(river, season, plot, subplot) %>%
  summarise_all(funs(if(is.numeric(.)) sum(., na.rm = TRUE) else first(.))) %>%
  group_by(river, plot, season) %>%
  summarise_all(funs(if(is.numeric(.)) mean(., na.rm = TRUE) else first(.)))

# create community matrix and envi data matrix 
p_matrix<-parasite[,28:65] # makes new matrix with parasite data (columns 28-65)
p_envi<-parasite[,1:27] # makes new matrix with envi data (columns 1-27)

# remove columns 
p_matrix<-p_matrix[, colSums(abs(p_matrix)) !=0] # removes columns with all zeros

# remove rows without parasites
p01<-p_matrix[rowSums(abs(p_matrix))!=0,]  # remove rows with no parasites
envi01<-p_envi[rowSums(abs(p_matrix))!=0,]  # remove rows with no parasites 

#standardize values to unit maximum 
pt<-apply(p01,2, function(x) x/max(x)) 

# global nmds
#ordO_p<-metaMDS(pt, distance="bray", trymax = 100, autotransform=FALSE) 
ordO_p<-metaMDS(p01, distance="bray", trymax = 100, autotransform=FALSE) 

ordO_p
# PERMANOVA
#p_para<-adonis(formula=pt ~ as.factor(plot), data=envi01,  permutations=999, method="bray") 
p_para<-adonis(formula=p01 ~ as.factor(plot), data=envi01,  permutations=999, method="bray") 

p_para

# pairwise PERMANOVAs
# pairwise.adonis(pt, envi01$plot)

hcoordO<-as.data.frame(scores(ordO_p, display="sites"))#extracts coordinates for plot
pcoordO<-scores(ordO_p, display="species")#extracts coordinates for parasite vectors

parasite_all <- bind_cols(envi01, hcoordO)

topp<-max(parasite_all[,28:29]) #determines maximum and minimum values for the plot axes
bott<-min(parasite_all[,28:29]) #I draw from both data sets because I make the graph sqaure and

parasite_all$plot <- as.factor(parasite_all$plot)
###########################################################################################
nmds_parasite_Raritan <-ggplot(parasite_all, aes(x= NMDS1, y=NMDS2))+
  geom_hline(yintercept = 0, lty=2, color="grey") +
  geom_vline(xintercept = 0, lty=2, color="grey") +
  geom_point( size =2.75, aes(fill=plot), shape=21)+
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  xlim(c( bott, topp)) +
  ylim(c( bott, topp)) +
  scale_shape_manual(values=c(21))+
  scale_fill_manual(values=c('#f7f7f7','#969696','#252525')) + 
  theme_pubr() +
  ggtitle("Parasite")+
  theme(legend.position = "top", legend.title = element_blank(), 
        legend.text = element_text(size=15), legend.box = "horizontal") +
  guides(colour = guide_legend(override.aes = list(size=70)))+
  theme(plot.title = element_text(hjust = 0.5, size=15, face="bold")) +
  theme(legend.title = element_text(size=8.5, face="bold"),
        legend.text = element_text(size=8))+ guides(shape=FALSE)+ 
  theme(legend.position=c(0.13,0.85))+
  guides(fill=guide_legend("Plot",  keywidth=0.12,
                           keyheight=0.12, override.aes = list(shape=21),
                           default.unit="inch"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2)) # add black square panel around graphic

############### Jaccard ####################################################
parasite_pa<-apply(parasite[,28:65],2, function(x) {ifelse(x==0,0,1)})#makes presence-absence matrix
p_matrix_pa<-parasite_pa[, colSums(abs(parasite_pa)) !=0] # removes columns with all zeros
p01_pa<-p_matrix_pa[rowSums(abs(p_matrix_pa))!=0,] 
envi01_pa<-p_envi[rowSums(abs(p_matrix_pa))!=0,] 

ordO_p_pa<-metaMDS(p01_pa, distance="jaccard", trymax = 50, autotransform=FALSE)
r_para_pa<-adonis(formula=p01_pa ~ as.factor(plot), data=envi01_pa,  permutations=999, method="bray")
r_para_pa

############################################################################
# h3. habitat variables within subplots  
##########################################################################################
habitat01 <- river  %>% filter(river == "Raritan")%>% group_by (river, season, plot, subplot) %>%
  summarise_all(funs(if(is.numeric(.)) mean(., na.rm = TRUE) else first(.))) %>%
  group_by(river, plot, season) %>%
  summarise_all(funs(if(is.numeric(.)) mean(., na.rm = TRUE) else first(.))) %>%
  filter(!(season == "summer" & plot ==1))   # remove summer plot 1 observation becaues no parasite data for that 
  
# return all rows from DF1 where there are matching values in DF2,
# keeping just columns from DF1.
# habitat01 <- semi_join(habitat01, parasite_all, by=c("river", "season", "plot", "subplot"))

hab<-habitat01[,6:16] 
hab_phys<-habitat01[,c(6,8,15,16)] 
hab_chem<-habitat01[,10:14] 

hab_2<-habitat01[,c(1:4)] 

#ordO_hab<-metaMDS(hab, trymax = 100, autotransform=FALSE)#global NMDS
pca_hab<- prcomp(hab, scale = TRUE)
summary(pca_hab)
pca_axis <-pca_hab$rotation[,1:2]

apply(hab, 2, var)# check variance across variables

ordO_hab<-rda(hab, scale=TRUE)# pca in vegan 
ordO_hab_phy<-rda(hab_phys, scale=TRUE)# pca in vegan 
ordO_hab_chem<-rda(hab_chem, scale=TRUE)# pca in vegan 

summary(ordO_hab_chem)
summary(ordO_hab_phy)
###########################################################################################
# create PCA plot for water quality variables 
sitecoordO_habchem<-as.data.frame(scores(ordO_hab_chem, display="sites"))#extracts coordinates for plots
habitat_chem <- bind_cols(hab_2, sitecoordO_habchem)

topp<-max(habitat_chem[,5:6]) 
bott<-min(habitat_chem[,5:6]) 

habitat_chem$plot <- as.factor(habitat_chem$plot)

nmds_habitat_Raritan <-ggplot(habitat_chem, aes(x= PC1, y=PC2))+
  geom_hline(yintercept = 0, lty=2, color="grey") +
  geom_vline(xintercept = 0, lty=2, color="grey") +
  geom_point( size =2.75, aes(fill=plot), shape=21)+
  xlab("PC 1") +
  ylab("PC 2") +
  xlim(c( bott, topp)) +
  ylim(c( bott, topp)) +
  scale_shape_manual(values=c(21))+
  scale_fill_manual(values=c('#f7f7f7','#969696','#252525')) + 
  theme_pubr() +
  ggtitle("Water quality")+
  theme(legend.position = "top", legend.title = element_blank(), 
        legend.text = element_text(size=15), legend.box = "horizontal") +
  guides(colour = guide_legend(override.aes = list(size=70)))+
  theme(plot.title = element_text(hjust = 0.5, size=15, face="bold")) +
  theme(legend.title = element_text(size=8.5, face="bold"),
        legend.text = element_text(size=8))+ guides(shape=FALSE)+ 
  theme(legend.position=c(0.13,0.17))+
  guides(fill=guide_legend("Plot",  keywidth=0.12,
                           keyheight=0.12, override.aes = list(shape=21),
                           default.unit="inch"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2)) # add black square panel around graphic
###########################################################################################
# create PCA plot for physical stream variables 
sitecoordO_habphy<-as.data.frame(scores(ordO_hab_phy, display="sites"))#extracts coordinates for plots
habitat_phy <- bind_cols(hab_2, sitecoordO_habphy)

topp<-max(habitat_phy[,5:6]) 
bott<-min(habitat_phy[,5:6]) 

habitat_phy$plot <- as.factor(habitat_phy$plot)

nmds_phy_Raritan <-ggplot(habitat_phy, aes(x= PC1, y=PC2))+
  geom_hline(yintercept = 0, lty=2, color="grey") +
  geom_vline(xintercept = 0, lty=2, color="grey") +
  geom_point( size =2.75, aes(fill=plot), shape=21)+
  xlab("PC 1") +
  ylab("PC 2") +
  xlim(c( bott, topp)) +
  ylim(c( bott, topp)) +
  scale_shape_manual(values=c(21))+
  scale_fill_manual(values=c('#f7f7f7','#969696','#252525')) + 
  theme_pubr() +
  ggtitle("Physical habitat")+
  theme(legend.position = "top", legend.title = element_blank(), 
        legend.text = element_text(size=15), legend.box = "horizontal") +
  guides(colour = guide_legend(override.aes = list(size=70)))+
  theme(plot.title = element_text(hjust = 0.5, size=15, face="bold")) +
  theme(legend.title = element_text(size=8.5, face="bold"),
        legend.text = element_text(size=8))+ guides(shape=FALSE)+ 
  theme(legend.position=c(0.13,0.17))+
  guides(fill=guide_legend("Plot",  keywidth=0.12,
                           keyheight=0.12, override.aes = list(shape=21),
                           default.unit="inch"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2)) # add black square panel around graphic

################################################################################################
# h4: host community abundance within subplot 
##########################################################################################
# abundance based FISH MATRIX 
fish_01 <- river %>% filter (river == "Raritan") %>% mutate(abun = rep(1)) %>% group_by(river, plot, season, subplot, host_species) %>%
  summarise_all(funs(if(is.numeric(.)) sum(., na.rm = TRUE) else first(.))) %>% spread(host_species, abun)
fish_01[,20:84][is.na(fish_01[,20:84])] <- 0
fish <- fish_01 %>%
  group_by(river, plot, subplot, season) %>%
  summarise_all(funs(if(is.numeric(.)) sum(., na.rm = TRUE) else first(.))) %>%
  group_by(river, plot, season) %>%  
  filter(!(season == "summer" & plot ==1)) %>%  # remove summer plot 1 observation becaues no parasite data for that 
  summarise_all(funs(if(is.numeric(.)) mean(., na.rm = TRUE) else first(.)))
fish <- fish[,-c(27:65)] # remove parasite data

# return all rows from DF1 where there are matching values in DF2,
# keeping just columns from DF1.
#fish <- semi_join(fish, parasite_all, by=c("river", "season", "plot", "subplot"))

# abundance based matrix 
fish_matrix<-fish[,27:44] # makes new matrix with fish data (columns 18-41)
fish_envi<-fish[,1:26] # makes new matrix with envi data (columns 1-17)

fish_matrix<-fish_matrix[, colSums(abs(fish_matrix)) !=0] # removes columns with all zeros

fish_01<-fish_matrix[rowSums(abs(fish_matrix))!=0,] 
fish_envi01<-fish_envi[rowSums(abs(fish_matrix))!=0,] 

# transform the abundance data here; to allow for a uniform comparison we standardize abundance values for each parasite species maximum 
fish_t<-apply(fish_01,2, function(x) x/max(x))#standardizes columns to unit maximum
#second argument =1 to apply to rows and =2 to apply to columns

#Bray-Curtis NMDS
#ordO<-metaMDS(fish_t, distance="bray", trymax = 100, autotransform=FALSE)#global NMDS
ordO<-metaMDS(fish_01, distance="bray", trymax = 100, autotransform=FALSE)#global NMDS

stressplot(ordO)
ordO

#p_abun<-adonis(formula=fish_t ~ as.factor(plot) , data=fish_envi01,  permutations=999, method="bray")
p_abun<-adonis(formula=fish_01 ~ as.factor(plot) , data=fish_envi01,  permutations=999, method="bray")

p_abun

# pairwise PERMANOVAs
# pairwise.adonis(fish_t, fish_envi01$plot)

sitecoordO<-as.data.frame(scores(ordO, display="sites"))#extracts coordinates for hosts
fishcoordO<-scores(ordO, display="species")#extracts coordinates for parasite vectors

fish_all <- bind_cols(fish_envi01, sitecoordO)

topp<-max(fish_all[,27:28]) #determines maximum and minimum values for the plot axes
bott<-min(fish_all[,27:28]) #I draw from both data sets because I make the graph sqaure and
#the scale of the X and Y axes should be the same on the graph for proper visual interpretation

fish_all$plot <- as.factor(fish_all$plot)
###########################################################################################
nmds_fish_abun_Raritan <-ggplot(fish_all, aes(x= NMDS1, y=NMDS2))+
  geom_hline(yintercept = 0, lty=2, color="grey") +
  geom_vline(xintercept = 0, lty=2, color="grey") +
  geom_point( size =2.75, aes(fill=plot), shape=21)+
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  xlim(c( bott, topp)) +
  ylim(c( bott, topp)) +
  scale_shape_manual(values=c(21))+
  scale_fill_manual(values=c('#f7f7f7','#969696','#252525')) + 
  theme_pubr() +
  ggtitle("Fish numerical density")+
  theme(legend.position = "top", legend.title = element_blank(), 
        legend.text = element_text(size=15), legend.box = "horizontal") +
  guides(colour = guide_legend(override.aes = list(size=70)))+
  theme(plot.title = element_text(hjust = 0.5, size=15, face="bold")) +
  theme(legend.title = element_text(size=8.5, face="bold"),
        legend.text = element_text(size=8))+ guides(shape=FALSE)+ 
  theme(legend.position=c(0.13,0.17))+
  guides(fill=guide_legend("Plot",  keywidth=0.12,
                           keyheight=0.12, override.aes = list(shape=21),
                           default.unit="inch"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2)) # add black square panel around graphic
############### Jaccard ####################################################
fish_pa<-apply(fish[,27:44],2, function(x) {ifelse(x==0,0,1)})#makes presence-absence matrix
fish_01pa<-fish_pa[rowSums(abs(fish_pa))!=0,] 
fish_envi01pa<-fish_envi[rowSums(abs(fish_pa))!=0,] 

#Bray-Curtis NMDS
ordOpa<-metaMDS(fish_01pa, distance="jaccard", trymax = 100, autotransform=FALSE)#global NMDS
stressplot(ordOpa)
ordOpa

r_fishpa<-adonis(formula=fish_01pa ~ as.factor(plot) , data=fish_envi01pa,  permutations=999, method="bray")
r_fishpa

sitecoordO<-as.data.frame(scores(ordOpa, display="sites"))#extracts coordinates for hosts
fishcoordO<-scores(ordOpa, display="species")#extracts coordinates for parasite vectors

fish_all_pa <- bind_cols(fish_envi01pa, sitecoordO)

topp<-max(fish_all_pa[,27:28]) #determines maximum and minimum values for the plot axes
bott<-min(fish_all_pa[,27:28]) #I draw from both data sets because I make the graph sqaure and
#the scale of the X and Y axes should be the same on the graph for proper visual interpretation

fish_all_pa$plot <- as.factor(fish_all_pa$plot)
###########################################################################################
nmds_fish_pa <-ggplot(fish_all_pa, aes(x= NMDS1, y=NMDS2))+
  geom_hline(yintercept = 0, lty=2, color="grey") +
  geom_vline(xintercept = 0, lty=2, color="grey") +
  geom_point( size =2.75, aes(fill=plot), shape=21)+
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  xlim(c( bott, topp)) +
  ylim(c( bott, topp)) +
  scale_shape_manual(values=c(21))+
  scale_fill_manual(values=c('#f7f7f7','#969696','#252525')) + 
  theme_pubr() +
  ggtitle("Fish P/A")+
  theme(legend.position = "top", legend.title = element_blank(), 
        legend.text = element_text(size=15), legend.box = "horizontal") +
  guides(colour = guide_legend(override.aes = list(size=70)))+
  theme(plot.title = element_text(hjust = 0.5, size=15, face="bold")) +
  theme(legend.title = element_text(size=8.5, face="bold"),
        legend.text = element_text(size=8))+ guides(shape=FALSE)+ 
  theme(legend.position=c(0.13,0.17))+
  guides(fill=guide_legend("Plot",  keywidth=0.12,
                           keyheight=0.12, override.aes = list(shape=21),
                           default.unit="inch"))
################################################################################################
# h5: host community biomass
##########################################################################################
fish_02 <- river %>% filter (river == "Raritan") %>% mutate(abun = rep(1)) %>% 
  group_by(river, plot, season, subplot, host_species) %>%
  summarise_all(funs(if(is.numeric(.)) sum(., na.rm = TRUE) else first(.))) %>% spread(host_species, dry_weight_g)
fish_02[,20:84][is.na(fish_02[,20:84])] <- 0
fish2 <- fish_02 %>%
  group_by(river,season, plot, subplot) %>%
  summarise_all(funs(if(is.numeric(.)) sum(., na.rm = TRUE) else first(.))) %>%
  group_by(river,plot,season) %>%  
  filter(!(season == "summer" & plot ==1)) %>%  # remove summer plot 1 observation becaues no parasite data for that 
  summarise_all(funs(if(is.numeric(.)) mean(., na.rm = TRUE) else first(.))) 
fish2 <- fish2[,-c(27:65)]  # remove parasite data

# return all rows from DF1 where there are matching values in DF2,
# keeping just columns from DF1.
#fish_02 <- semi_join(fish2, parasite_all, by=c("river", "season", "plot", "subplot"))

# abundance based matrix 
fish_matrix2<-fish2[,27:44] # makes new matrix with fish data (columns 18-41)
fish_envi2<-fish2[,1:26] # makes new matrix with envi data (columns 1-17)

fish_matrix2<-fish_matrix2[, colSums(abs(fish_matrix2)) !=0] # removes columns with all zeros

fish_001<-fish_matrix2[rowSums(abs(fish_matrix2))!=0,] 
fish_envi001<-fish_envi2[rowSums(abs(fish_matrix2))!=0,] 

# transform the abundance data here; to allow for a uniform comparison we standardize abundance values for each parasite species maximum 
fish_t0<-apply(fish_001,2, function(x) x/max(x))#standardizes columns to unit maximum
#second argument =1 to apply to rows and =2 to apply to columns

#Bray-Curtis NMDS
#ordO_02<-metaMDS(fish_t0, distance="bray", trymax = 100, autotransform=FALSE)#global NMDS
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
nmds_biomass_raritan <-ggplot(fish_all_02, aes(x= NMDS1, y=NMDS2))+
  geom_hline(yintercept = 0, lty=2, color="grey") +
  geom_vline(xintercept = 0, lty=2, color="grey") +
  geom_point( size =2.75, aes(fill=plot), shape=21)+
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  xlim(c( bott, topp)) +
  ylim(c( bott, topp)) +
  scale_shape_manual(values=c(21))+
  scale_fill_manual(values=c('#f7f7f7','#969696','#252525')) + 
  theme_pubr() +
  ggtitle("Fish biomass density")+
  theme(legend.position = "top", legend.title = element_blank(), 
        legend.text = element_text(size=15), legend.box = "horizontal") +
  guides(colour = guide_legend(override.aes = list(size=70)))+
  theme(plot.title = element_text(hjust = 0.5, size=15, face="bold")) +
  theme(legend.title = element_text(size=8.5, face="bold"),
        legend.text = element_text(size=8))+ guides(shape=FALSE)+ 
  theme(legend.position=c(0.13,0.17))+
  guides(fill=guide_legend("Plot",  keywidth=0.12,
                           keyheight=0.12, override.aes = list(shape=21),
                           default.unit="inch"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2)) # add black square panel around graphic

#############################################################################################
# export graphic of all ordinations 
#############################################################################################

jpeg(file="Raritan_river_ordinations.jpeg", width = 180, height = 270, units = 'mm', res = 300)
ggarrange(nmds_habitat_Raritan, nmds_phy_Raritan, nmds_fish_abun_Raritan, 
          nmds_biomass_raritan, nmds_parasite_Raritan, common.legend = TRUE, ncol = 2, nrow = 3)
dev.off()

#############################################################################################
# pairwise procrustes and protests
###########################################################################################

# host abundance and habitat
pro_1 <- procrustes(ordO, ordO_hab)
#plot(pro_1, kind="1")
prot_ha_env <- protest(ordO,ordO_hab, permutations = 999)
a <- c(prot_ha_env$ss, prot_ha_env$signif)
# host biomass and habitat
pro_2 <- (procrustes(ordO_02, ordO_hab))
#plot(pro_2, kind="1")
prot_hb_env <- protest(ordO_02, ordO_hab, permutations = 999)
prot_hb_env
b<- c(prot_hb_env$ss, prot_hb_env$signif)
# host biomass and abundance
pro_3 <-(procrustes(ordO, ordO_02))
#jpeg(file="Raritan_procrustes_host_abundance_biomass.jpeg", width = 180, height = 180, units = 'mm', res = 300)
#plot(pro_3, kind="1")
#dev.off()
prot_hb_ha <- protest(ordO, ordO_02, permutations = 999)
prot_hb_ha
c<- c(prot_hb_ha$ss, prot_hb_ha$signif)
# host abundance and paraiste 
pro_4 <- procrustes(ordO, ordO_p)
#plot(pro_1, kind="1")
prot_ha_p <- protest(ordO, ordO_p, permutations = 999)
prot_ha_p
d<- c(prot_ha_p$ss, prot_ha_p$signif)
# host biomass and parasite
pro_5 <- (procrustes(ordO_02, ordO_p))
#plot(pro_5, kind="1")
prot_hb_p <- protest(ordO_02, ordO_p, permutations = 999)
prot_hb_p
e<- c(prot_hb_p$ss, prot_hb_p$signif)
# habitat and parasite 
pro_6 <- (procrustes(ordO_hab, ordO_p))
#plot(pro_5, kind="1")
prot_env_p <- protest(ordO_hab, ordO_p, permutations = 999)
prot_env_p
f<- c(prot_env_p$ss, prot_env_p$signif)

#############################################
# parasite and chemical 
prot_env_pchem <- protest(ordO_hab_chem, ordO_p, permutations = 999)
prot_env_pchem
g<- c(prot_env_pchem$ss, prot_env_pchem$signif)

#parasite and physical 
prot_env_pphy <- protest(ordO_hab_phy, ordO_p, permutations = 999)
prot_env_pphy
h<- c(prot_env_pphy$ss, prot_env_pphy$signif)

# host abundance and chemical
pro_12 <- procrustes(ordO, ordO_hab_chem)
#plot(pro_1, kind="1")
prot_ha_envchem <- protest(ordO,ordO_hab_chem, permutations = 999)
i <- c(prot_ha_envchem$ss, prot_ha_envchem$signif)

# host abundance and phys
pro_13 <- procrustes(ordO, ordO_hab_phy)
#plot(pro_1, kind="1")
prot_ha_envphy <- protest(ordO,ordO_hab_phy, permutations = 999)
j <- c(prot_ha_envphy$ss, prot_ha_envphy$signif)

# host biomass and chemical
pro_14 <- procrustes(ordO_02, ordO_hab_chem)
#plot(pro_1, kind="1")
prot_bio_envchem <- protest(ordO_02,ordO_hab_chem, permutations = 999)
k <- c(prot_bio_envchem$ss, prot_bio_envchem$signif)

# host biomass and phys
pro_15 <- procrustes(ordO_02, ordO_hab_phy)
#plot(pro_1, kind="1")
prot_bio_envphy <- protest(ordO_02,ordO_hab_phy, permutations = 999)
l <- c(prot_bio_envphy$ss, prot_bio_envphy$signif)

jpeg(file="Raritan_procrustes_errors.jpeg", width = 270, height = 270, units = 'mm', res = 300)
par(mfrow=c(3,3))
plot(pro_1, kind="1",  main="Fish-parasite")
plot(pro_5, kind="1",main="Fish biomass-parasite")
plot(prot_env_pchem , kind="1", main="Parasite-chemical")
plot(prot_env_pphy  , kind="1", main="Parasite-physical")
plot(pro_12  , kind="1",  main="Fish-chemical")
plot(pro_13  , kind="1", main="Fish-physical")
plot(pro_14  , kind="1", main="Fish biomass-chemical")
plot(pro_15  , kind="1", main="Fish biomass-physical")
plot(pro_3, kind="1", main="Fish biomass-abundance")
dev.off()
################################################################################################################
# put PROTEST CONCORDANCE in one df 
################################################################################################################

Raritan_protest <-data.frame(a,b,c,d,e,f,g,h,i,j,k,l)
blah <-c('fish_abundance-habitat', 'fish_biomass-habitat', 'fish_abundance-biomass', 
         'fish-abundance-parasite', 'fish_biomass-parasite', 'pararasite-habitat',
         'parasite-chemical', 'parasite-phyiscal', 'fish_abundance-chemical', 'fish_abundance-physical',
         'fish_biomass-chemical', 'fish_biomass-physical')
meh <-c('ms_square', "p_value")
colnames(Raritan_protest) <- blah
rownames(Raritan_protest) <- meh

###########################################################################################
# put PERMANOVA in one df 
################################################################################################################
parasite_permanova <- data.frame(p_para$aov.tab[1,1], p_para$aov.tab[3,1], p_para$aov.tab[1,6], p_para$aov.tab[1,4], p_para$aov.tab[1,5])
fish_numerical_permanova <- c(p_abun$aov.tab[1,1], p_abun$aov.tab[3,1], p_abun$aov.tab[1,6], p_abun$aov.tab[1,4], p_abun$aov.tab[1,5])
fish_biomass_permanova <- c(p_biomass$aov.tab[1,1], p_biomass$aov.tab[3,1], p_biomass$aov.tab[1,6], p_biomass$aov.tab[1,4], p_biomass$aov.tab[1,5])

Raritan_permanova <-rbind(parasite_permanova, fish_numerical_permanova, fish_biomass_permanova)
thing <-c('DF', 'DF2', 'P_value', 'F_model', 'R2')
yeah <-c('parasite', "fish_numerica", "fish_biomass")
colnames(Raritan_permanova) <- thing
rownames(Raritan_permanova) <- yeah

###########################################################################################
###########################################################################################
#############################################################################################
# export graphic of all ordinations 
#############################################################################################
#############################################################################################
ordination_agg <- cbind(parasite_all$NMDS1, parasite_all$NMDS2, 
                        fish_all$NMDS1, fish_all$NMDS2, 
                        fish_all_02$NMDS1, fish_all_02$NMDS2, 
                        habitat_chem$PC1, habitat_chem$PC2, 
                        habitat_phy$PC1, habitat_phy$PC2) 

parasite_o <-data.frame(t(c(parasite_all$NMDS1, parasite_all$NMDS2)))
fish_o <-data.frame(t(c(fish_all$NMDS1, fish_all$NMDS2)))
fish_bo <-data.frame(t(c(fish_all_02$NMDS1, fish_all_02$NMDS2)))
chem_o <-data.frame(t(c(habitat_chem$PC1, habitat_chem$PC2)))
phys_o <-data.frame(t(c(habitat_phy$PC1, habitat_phy$PC2)))

ord_agg <-rbind(parasite_o, fish_o, fish_bo, chem_o, phys_o)
#ord_agg_t <- t(ord_agg)
#ord_agg <-rbind(ordO_p$dist, ordO$dist, ordO_02$dist)
ord_ord<-rda(ordination_agg, scale=TRUE)# pca in vegan
ord_ord<-rda(ord_agg, scale=TRUE)# pca in vegan

plot(ord_ord)

ord_group <- data.frame(c("parasite", "fish", "fish biomass", "chemical", "physical"
))
colnames(ord_group) <- "Eco"

ord_labels <-(c("parasite", "fish", "fish biomass", "chem", "phys"))

sitecoordO_ord<-as.data.frame(scores(ord_ord, display="sites"))#extracts coordinates for plots
factcoordO<-as.data.frame(scores(ord_ord, display="species"))#extracts coordinates for parasite vectors
ordo_all <- bind_cols(ord_group, sitecoordO_ord)

topp<-max(ordo_all[,2:3]) 
bott<-min(ordo_all[,2:3]) 

###########################################################################################
###########################################################################################
pca_ordinations_raritan <-ggplot(ordo_all, aes(x= PC1, y=PC2))+
  geom_hline(yintercept = 0, lty=2, color="grey") +
  geom_vline(xintercept = 0, lty=2, color="grey") +
  geom_point( size =5, aes(fill=Eco), shape=21)+
  xlab("PC 1") +
  ylab("PC 2") +
  xlim(c( bott, topp)) +
  ylim(c( bott, topp)) +
  scale_shape_manual(values=c(21))+
  theme_pubr() +
  scale_fill_manual(values=c('#1b9e77','#d95f02','#7570b3'
                             ,'#e7298a','#66a61e'
  )) + 
  theme(legend.position = "top", legend.title = element_blank(), 
        legend.text = element_text(size=15), legend.box = "horizontal") +
  guides(colour = guide_legend(override.aes = list(size=70)))+
  theme(plot.title = element_text(hjust = 0.5, size=15, face="bold")) +
  theme(legend.title = element_text(size=8.5, face="bold"),
        legend.text = element_text(size=8))+ guides(shape=FALSE)+ 
  theme(legend.position=c(0.75,0.75))+
  guides(fill=guide_legend("",  keywidth=0.12,
                           keyheight=0.12, override.aes = list(shape=21),
                           default.unit="inch"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2)) # add black square panel around graphic


pca_ordinations_raritan <-ggplot(ordo_all, aes(x= PC1, y=PC2))+
  geom_blank()+
  xlab("PC 1") +
  ylab("PC 2") +
  xlim(c( bott, topp)) +
  ylim(c( bott, topp)) +
  geom_hline(yintercept = 0, lty=2, colour= '#bdbdbd') +
  geom_vline(xintercept = 0, lty=2, colour= '#bdbdbd') +
  theme_pubr() +
  geom_text(data=ordo_all, aes(x=PC1,y=PC2, label=ord_labels), fontface =2, colour = "black" , size=2.5) + # plant species label 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) # add black square panel around graphic


jpeg(file="Raritan_river_ordinations.jpeg", width = 100, height = 100, units = 'mm', res = 300)
pca_ordinations_raritan
dev.off()
