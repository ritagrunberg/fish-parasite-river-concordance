rm(list = ls())  # clear Rs brain 
###########################################################################################
library(readr)   # read csv file
library(tidyr)   # tidying data
library(dplyr)   # manipulating df
library(tidyverse)
library(ggplot2) # graphing
library(ggpubr)  # nice package for figures
library(vegan)   # nmds 
library(iNEXT)
###########################################################################################
# AIM TESTING FOR CONCORDANCE BETWEEN ALL PARASITES PRESENT IN ECOSYSTEM AND HOST + ENVIR VARIABLES 
# import dataseIt
river <- read_csv("C:/Users/grunberg/Dropbox/DISSERTATION_CHAPTERS/1_parasite_community_spatial/data/final_river_data.csv")
coordinates <- read_csv("Concordance_submitted_MS/Revision/data_files/river_coordinates.csv")
river <- river[,-c(38,57)] # remove notes that were in excel file... trash 
set.seed(1234567890)
river$plot <- as.factor(river$plot)
river$subplot <- as.factor(river$subplot)

envi <- river  %>% 
  group_by(river, plot) %>%
  # summarise_if(.predicate = function(x) is.numeric(x), .funs = funs(mean="mean"))
  summarise_all(funs(if(is.numeric(.)) mean(., na.rm = TRUE) else first(.))) %>%
  unite(river_plot, river, plot, season, sep = '_') #merge characters from columns into one thing


parasite <- river  %>% 
  group_by(river, plot) %>%
  # summarise_if(.predicate = function(x) is.numeric(x), .funs = funs(mean="mean"))
  summarise_all(funs(if(is.numeric(.)) sum(., na.rm = TRUE) else first(.))) %>%
  unite(river_plot, river, plot, season, sep = '_') #merge characters from columns into one thing


p_matrix<-parasite[,c(26:63)] # makes new matrix with parasite data (columns 28-65)
id <- unique(parasite$river_plot_seasaon)


rownames(p_matrix) <- c("passaic_1upstream", "passaic_2midstream", "passaic_3downstream",
                        "raritan_1upstream", "raritan_2midstream", "raritan_3downstream") 


#rownames(p_matrix) <- c("Raritan_1_fall"  , "Raritan_1_spring" ,"Raritan_1_summer" ,"Raritan_1_winter" 
 #                       ,"Raritan_2_fall"  , "Raritan_2_spring" ,"Raritan_2_summer", "Raritan_2_winter",
  #                      "Raritan_3_fall"  , "Raritan_3_spring" ,"Raritan_3_summer" ,"Raritan_3_winter")

p_matrix <- as.matrix(p_matrix)
#p_matrix<-p_matrix[rowSums(abs(p_matrix))!=0,]  # remove rows with no parasites
############################################################################################

plist = lapply(as.list(1:dim(p_matrix)[1]), function(x) p_matrix[x[1],])
plist = list(passaic1_upstream = p_matrix[1,],
             passaic2_midstream = p_matrix[2,],
             passaic3_downstream = p_matrix[3,]
             )


#out <- iNEXT(p_matrix, q=c(0), datatype="abundance")
m <- c(1, 100, 200, 500, 1000, 3000, 7500)
out <- iNEXT(plist, q=c(0), datatype="abundance", nboot=999, size = m)
out
out$iNextEst
out$AsyEst #lists the observed diversity, asymptotic estimates, estimated bootstrap s.e. and 95% confidence intervals for Hill numbers with q = 0, 1, and 2.
#################### Sample completeness curves #######################
passaic_rarefaction1 <-ggiNEXT(out, type=1, color.var="site", se=FALSE) +
  #ylim(c(0.9,1)) +
  theme_bw(base_size = 10) + 
  ggtitle("Passaic Paraste: sample based rarefaction")+
  theme(legend.position="bottom")+
  scale_fill_manual(values=c('#1b9e77','#d95f02','#7570b3'),
                    labels = c("upstream","midstream", "downstream"))+
  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'),
                     labels = c("upstream","midstream", "downstream"))+
  guides(shape=FALSE)

######################## Coverage-based R/E curves ##################

passaic_coverage1 <-ggiNEXT(out, type=3, color.var="site", facet.var="site", grey = TRUE) +
  #ylim(c(0.9,1)) +
  theme_bw(base_size = 10) + 
  ggtitle("Passaic Parasite: coverage based rarefaction")+
  theme(legend.position="none")

#################### Sample completeness curves #######################

passaic_coverage2 <-ggiNEXT(out, type=2, color.var ="site") + 
 # xlim(c(0,1)) +
  theme_bw(base_size = 10) +
  theme(legend.position="bottom",
        legend.title=element_blank())+
  ggtitle("Passaic River parasite coverage based rarifcation ")

jpeg(filename="passaic_parasite_coverage.jpeg", width=175, height=250,  units="mm", bg="white", res=300)
ggarrange(passaic_coverage1, passaic_coverage2, ncol=1, nrow=2)
dev.off()

############################################################################################
############################################################################################
############################################################################################

rlist = list(raritan1_upstream = p_matrix[4,],
             raritan2_midstream = p_matrix[5,],
             raritan3_downstream = p_matrix[6,]
)
 
m2 <- c(1, 100, 200, 500, 1000, 5000, 7000, 10000, 15000, 20000, 25000)
out2 <- iNEXT(rlist, q=c(0), datatype="abundance", size=m2, nboot=999)

###########rarefaction curve 
raritan_rarefaction1 <-ggiNEXT(out2, type=1, color.var="site", se=FALSE) +
  #ylim(c(0.9,1)) +
  theme_bw(base_size = 10) + 
  ggtitle("Raritan Parasite: sample based rarefaction")+
  theme(legend.position="bottom")+
  scale_fill_manual(values=c('#1b9e77','#d95f02','#7570b3'),
                    labels = c("upstream","midstream", "downstream"))+
  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'),
                     labels = c("upstream","midstream", "downstream"))+
  guides(shape=FALSE)

#################### Sample completeness curves ######################
raritan_coverage1 <-ggiNEXT(out2, type=3, color.var="site", facet.var="site", grey=TRUE) +
  #ylim(c(0.9,1)) +
  theme_bw(base_size = 10) + 
  theme(legend.position="none")+
  ggtitle("Raritan Paraste: coverage based rarefaction")
  
######################## Coverage-based R/E curves ##################
raritan_coverage2 <-ggiNEXT(out2, type=3, color.var ="site") + 
 # xlim(c(0.6,1)) +
  theme_bw(base_size = 10) +
  theme(legend.position="bottom",
        legend.title=element_blank())+
  ggtitle("Raritan River parasite coverage based rarifcation ")


jpeg(filename="raritan_parasite_coverage.jpeg", width=175, height=250, units="mm", bg="white", res=300)
ggarrange(raritan_coverage1, raritan_coverage2, ncol=1, nrow=2)
dev.off()

############################################################################################
##################### FISH FISH FISH #######################################################################
############################################################################################

# abundance based FISH MATRIX 
fish_01 <- river %>%  
  mutate(abun = rep(1)) %>% group_by(river, plot, season, subplot, host_species) %>%
  summarise_all(funs(if(is.numeric(.)) sum(., na.rm = TRUE) else first(.))) %>% spread(host_species, abun)
fish_01[,26:ncol((fish_01))][is.na(fish_01[,26:ncol((fish_01))])] <- 0
fish <- fish_01 %>%
  group_by(river, plot) %>%
  summarise_all(funs(if(is.numeric(.)) sum(., na.rm = TRUE) else first(.))) 
fish <- fish[,-c(1:65,91)] # remove parasite data

fish <- as.matrix(fish)
rownames(fish) <- c("passaic_upstream", "passaic_midstream", "passaic_downstream",
                    "raritan_upstream", "raritan_midstream", "raritan_downstream") 

pflist = list(passaic1_upstream = fish[1,],
             passaic2_midstream = fish[2,],
             passaic3_downstream = fish[3,]
)

m3 <- c(1, 10, 50, 100, 125, 150, 200)

out3 <- iNEXT(pflist, q=c(0), datatype ="abundance", nboot=999, size=m3)
###########rarefaction curve 
passaic_rarefaction1f <-ggiNEXT(out3, type=1, color.var="site", se=FALSE) +
  #ylim(c(0.9,1)) +
  theme_bw(base_size = 10) + 
  ggtitle("Passaic Fish: sampled based rarefaction")+
  theme(legend.position="bottom")+
  scale_fill_manual(values=c('#1b9e77','#d95f02','#7570b3'),
                    labels = c("upstream","midstream", "downstream"))+
  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'),
                     labels = c("upstream","midstream", "downstream"))+
  guides(shape=FALSE)

#################### Sample completeness curves #######################
passaic_coverage1f <-ggiNEXT(out3, type=3, color.var="site", facet.var="site", grey=TRUE) +
  #ylim(c(0.9,1)) +
  theme_bw(base_size = 10) + 
  theme(legend.position="none")+
  ggtitle("Passaic Fish: coverage based rarefaction")

######################## Coverage-based R/E curves ##################
passaic_coverage2f <-ggiNEXT(out3, type=3, color.var ="site") + 
  # xlim(c(0.9,1)) +
  theme_bw(base_size = 10) +
  theme(legend.position="bottom",
        legend.title=element_blank())+
  ggtitle("Passaic River fish coverage based rarifcation ")

jpeg(filename="passaic_fish_coverage.jpeg", width=175, height=250, units="mm", bg="white", res=300)
ggarrange(passaic_coverage1f, passaic_coverage2f, ncol=1, nrow=2)
dev.off()

############################################################################################
############################################################################################

rflist = list(raritan1_upstream = fish[4,],
              raritan2_midstream = fish[5,],
              raritan3_downstream = fish[6,]
)

m4 <- c(1, 10, 50, 100,  150, 200, 400, 600, 800)

out4 <- iNEXT(rflist, q=c(0), datatype="abundance", nboot=999, size=m4)
out4$DataInfo
###########rarefaction curve 
raritan_rarefaction1f <-ggiNEXT(out4, type=1, color.var="site", se=FALSE) +
  #ylim(c(0.9,1)) +
  theme_bw(base_size = 10) + 
  ggtitle("Raritan Fish: sampled based rarefaction")+
  theme(legend.position="bottom")+
  scale_fill_manual(values=c('#1b9e77','#d95f02','#7570b3'),
                    labels = c("upstream","midstream", "downstream"))+
  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'),
                     labels = c("upstream","midstream", "downstream"))+
  guides(shape=FALSE)

#################### Sample completeness curves #######################
raritan_coverage1f <-ggiNEXT(out4, type=3, color.var="site", facet.var="site", grey=TRUE) +
  #ylim(c(0.9,1)) +
  theme_bw(base_size = 10) + 
  theme(legend.position="none")  +
  ggtitle("Raritan Fish: coverage based rarefaction")


######################## Coverage-based R/E curves ##################
raritan_coverage2f <-ggiNEXT(out4, type=3, color.var ="site") + 
  # xlim(c(0.5,1)) +
  theme_bw(base_size = 10) +
  theme(legend.position="bottom",
        legend.title=element_blank())+
  ggtitle("Raritan River fish coverage based rarifcation ")

jpeg(filename="coverage_rarefaction.jpeg", width=175, height=250, units="mm", bg="white", res=300)
ggarrange(passaic_coverage1, passaic_coverage1f,
          raritan_coverage1, raritan_coverage1f, ncol=2, nrow=2)
dev.off()

jpeg(filename="coverage_rarifcation_rivers.jpeg", width=300, height=300, units="mm", bg="white", res=300)
ggarrange(passaic_coverage2, passaic_coverage2f, raritan_coverage2, raritan_coverage2f, ncol=2, nrow=2)
dev.off()

jpeg(filename="rarefaction_curves_river.jpeg", width=180, height=180, units="mm", bg="white", res=300)
ggarrange(raritan_rarefaction1, raritan_rarefaction1f, 
          passaic_rarefaction1, passaic_rarefaction1f,
          ncol=2, nrow=2, common.legend = TRUE)
dev.off()

jpeg(filename="rarefaction_curves_river.jpeg", width=200, height=300, units="mm", bg="white", res=300)
ggarrange(passaic_rarefaction1, passaic_coverage1, 
          passaic_rarefaction1f, passaic_coverage1f,
          raritan_rarefaction1, raritan_coverage1,
          raritan_rarefaction1f, raritan_coverage1f,
          widths  = c(0.9, 1), ncol=2, nrow=4, common.legend = TRUE)
dev.off()
############################################################################################
############################################################################################
############################################################################################

parasite.beta.raritan <- river  %>% 
  filter(river=="Raritan")%>%
  group_by(river, plot,season) %>%
  # summarise_if(.predicate = function(x) is.numeric(x), .funs = funs(mean="mean"))
  summarise_all(funs(if(is.numeric(.)) sum(., na.rm = TRUE) else first(.))) %>%
  unite(river_plot_season, river, plot, season, sep = '_') #merge characters from columns into one thing


pr_matrix<-parasite.beta.raritan[,c(26:63)] # makes new matrix with parasite data (columns 28-65)
id <- unique(parasite.beta.raritan$river_plot_season)
rownames(pr_matrix) <- id
pr_matrix<-apply(pr_matrix,2, function(x) {ifelse(x==0,0,1)})  #makes presence-absence matrix

library(betapart) # beta diversity analysis 

#computes qunatities needed for multiple site beta and pairwise dissimilarity
parasite.rar = betapart.core(pr_matrix)
parasite.rar.mult = beta.multi(pr_matrix , index.family="sor")
# distribution of values from resampled data below.. 
parasite.rar.samp = beta.sample(parasite.rar, sites=10, samples=100, index.family="sor")

# create matrix for beta diversity values 
parasite.rar.pair <- beta.pair(pr_matrix,index.family = 'sor')
mantel(parasite.rar.pair$beta.sim, parasite.rar.pair$beta.sne, method='pearson', permutations = 999, na.rm = TRUE)

boxplot((parasite.rar.samp$sampled.values$beta.SNE), (parasite.rar.samp$sampled.values$beta.SIM),  
        xaxt = "n", main='', ylab="beta diversity", lwd=1, ylim=c(0,1)) # nestedness
axis(1, at=1:2, labels=c( expression(beta['SNE']), expression(beta['SIM']))) 
mtext("A", side=3, line=-1.5, adj=0.05, cex=1.2)

############################################################################################

parasite.beta.passaic <- river  %>% 
  filter(river=="Passaic")%>%
  group_by(river, plot,season) %>%
  # summarise_if(.predicate = function(x) is.numeric(x), .funs = funs(mean="mean"))
  summarise_all(funs(if(is.numeric(.)) sum(., na.rm = TRUE) else first(.))) %>%
  unite(river_plot_season, river, plot, season, sep = '_') #merge characters from columns into one thing

pp_matrix<-parasite.beta.passaic[,c(26:63)] # makes new matrix with parasite data (columns 28-65)
id <- unique(parasite.beta.passaic$river_plot_season)
rownames(pp_matrix) <- id
pp_matrix<-apply(pp_matrix,2, function(x) {ifelse(x==0,0,1)})  #makes presence-absence matrix

#computes qunatities needed for multiple site beta and pairwise dissimilarity
parasite.pas = betapart.core(pp_matrix)
parasite.pas.mult = beta.multi(pp_matrix , index.family="sor")
# distribution of values from resampled data below.. 
parasite.pas.samp = beta.sample(parasite.pas, sites=10, samples=100, index.family="sor")

# create matrix for beta diversity values 
parasite.pas.pair <- beta.pair(pp_matrix,index.family = 'sor')
mantel(parasite.pas.pair$beta.sim, parasite.pas.pair$beta.sne, method='pearson', permutations = 999, na.rm = TRUE)

boxplot((parasite.pas.samp$sampled.values$beta.SNE), 
        (parasite.pas.samp$sampled.values$beta.SIM),  
        xaxt = "n", main='', ylab="beta diversity",
        lwd=1, 
        ylim=c(0,1)) # nestedness
axis(1, at=1:2, labels=c( expression(beta['SNE']),
                          expression(beta['SIM']))) 
mtext("B", side=3, line=-1.5, adj=0.05, cex=1.2)

############################################################################################
# abundance based FISH MATRIX 
fish_01.raritan <- river %>%  
  mutate(abun = rep(1)) %>% group_by(river, plot, season, subplot, host_species) %>%
  summarise_all(funs(if(is.numeric(.)) sum(., na.rm = TRUE) else first(.))) %>% spread(host_species, abun)
fish_01.raritan[,26:ncol((fish_01.raritan))][is.na(fish_01.raritan[,26:ncol((fish_01.raritan))])] <- 0
fish.raritan <- fish_01.raritan %>%
  group_by(river, plot, season) %>%
  summarise_all(funs(if(is.numeric(.)) sum(., na.rm = TRUE) else first(.))) %>%
  filter(river=="Raritan")
fish.raritan<- fish.raritan[,-c(1:65,91)] # remove parasite data

fr_matrix<-apply(fish.raritan,2, function(x) {ifelse(x==0,0,1)})  #makes presence-absence matrix

#computes qunatities needed for multiple site beta and pairwise dissimilarity
fish.rar = betapart.core(fr_matrix)
fish.rar.mult = beta.multi(fr_matrix , index.family="sor")
# distribution of values from resampled data below.. 
fish.rar.samp = beta.sample(fish.rar, sites=10, samples=100, index.family="sor")

# create matrix for beta diversity values 
fish.rar.pair <- beta.pair(fr_matrix,index.family = 'sor')
mantel(fish.rar.pair$beta.sim, fish.rar.pair$beta.sne, method='pearson', permutations = 999, na.rm = TRUE)

boxplot((fish.rar.samp$sampled.values$beta.SNE), 
        (fish.rar.samp$sampled.values$beta.SIM),  
        xaxt = "n", main='', ylab="beta diversity",
        lwd=1, 
        ylim=c(0,1)
        ) # nestedness
axis(1, at=1:2, labels=c( expression(beta['SNE']),
                          expression(beta['SIM']))) 
mtext("C", side=3, line=-1.5, adj=0.05, cex=1.2)

