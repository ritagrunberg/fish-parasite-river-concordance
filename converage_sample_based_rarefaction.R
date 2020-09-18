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
# import dataset 
river <- read_csv("C:/Users/grunberg/Dropbox/DISSERTATION_CHAPTERS/1_parasite_community_spatial/data/final_river_data.csv")
river <- river[,-c(38,57)] # remove notes that were in excel file... trash 
set.seed(1234567890)
river$plot <- as.factor(river$plot)
river$subplot <- as.factor(river$subplot)
river <- river %>% dplyr::select(-c(sex))
river <- river %>% drop_na()
river$host_species <-gsub("Etheostoma olmstedi_","Etheostoma olmstedi",  river$host_species)

# formatting data to get site level means 
envi <- river  %>% 
  group_by(river, plot) %>%
  # summarise_if(.predicate = function(x) is.numeric(x), .funs = funs(mean="mean"))
  summarise_all(funs(if(is.numeric(.)) mean(., na.rm = TRUE) else first(.))) %>%
  unite(river_plot, river, plot, season, sep = '_') #merge characters from columns into one thing
parasite <- river  %>% 
  group_by(river, plot, season) %>%
  # summarise_if(.predicate = function(x) is.numeric(x), .funs = funs(mean="mean"))
  summarise_all(funs(if(is.numeric(.)) sum(., na.rm = TRUE) else first(.))) %>%
  unite(river_plot_season, river, plot, season, season, sep = '_') #merge characters from columns into one thing

p_matrix<-parasite[,c(25:62)] # makes new matrix with parasite data (columns 28-65)
id <- unique(parasite$river_plot_season)


#rownames(p_matrix) <- c("passaic_1upstream", "passaic_2midstream", "passaic_3downstream",
 #                       "raritan_1upstream", "raritan_2midstream", "raritan_3downstream") 


#rownames(p_matrix) <- c("Raritan_1_fall"  , "Raritan_1_spring" ,"Raritan_1_summer" ,"Raritan_1_winter" 
 #                       ,"Raritan_2_fall"  , "Raritan_2_spring" ,"Raritan_2_summer", "Raritan_2_winter",
  #                      "Raritan_3_fall"  , "Raritan_3_spring" ,"Raritan_3_summer" ,"Raritan_3_winter")

p_matrix <- as.matrix(p_matrix)
#p_matrix<-p_matrix[rowSums(abs(p_matrix))!=0,]  # remove rows with no parasites
############################################################################################

plist = lapply(as.list(1:dim(p_matrix)[1]), function(x) p_matrix[x[1],])
#plist = list(passaic1_upstream = p_matrix[1,],
 #           passaic2_midstream = p_matrix[2,],
  #           passaic3_downstream = p_matrix[3,]
   #          )

plist = list(passaic1_upstream_fall = p_matrix[1,],
             passaic1_upstream_spring = p_matrix[2,],
             passaic1_upstream_summer = p_matrix[3,],
             passaic1_upstream_winter = p_matrix[4,],
             passaic2_midstream_fall = p_matrix[5,],
             passaic2_midstream_spring = p_matrix[6,],
             passaic2_midstream_summer = p_matrix[7,],
             passaic2_midstream_winter = p_matrix[8,],
             passaic3_downstream_fall = p_matrix[9,],
             passaic3_downstream_spring = p_matrix[10,],
             passaic3_downstream_summer = p_matrix[11,],
             passaic3_downstream_winter = p_matrix[12,]
             )


#out <- iNEXT(p_matrix, q=c(0), datatype="abundance")
#m <- c(1, 100, 200, 500, 1000, 3000, 7500)
m <- c(1, 100, 200, 500, 1000, 3000, 8000)
out <- iNEXT(plist, q=c(0), datatype="abundance", size =m)
out
out$iNextEst

#diversity estimates
passaic_richness  <- estimateD(plist,datatype="abundance", base="coverage", level=0.99, conf=0.95)

#### export info
passaic_diversity2 <-out$AsyEst #lists the observed diversity, asymptotic estimates, estimated bootstrap s.e. and 95% confidence intervals for Hill numbers with q = 0, 1, and 2.
passaic_richness2 <- passaic_diversity %>% filter(Diversity=='Species richness') 

#### graph richness 
pparasite <-
  ggplot()+
  geom_point(data= passaic_richness,
             aes(x=Estimator,y = Site),
             size=2)+
  geom_point(data = passaic_richness, 
             aes(x=Observed,y = Site),
             position = position_jitter( h = 0.25),
             pch=21, size=2)+
  labs(x="richness",y= "")+
  ggtitle("Parasite Passaic River")+
 # coord_flip()+
  theme_bw()

######### upstream site only 
plist.upstream = list(fall = p_matrix[1,],
                      spring = p_matrix[2,],
                      summer = p_matrix[3,],
                      winter = p_matrix[4,]
                      )
out.up.1 <- iNEXT(plist.upstream, q=c(0), datatype="abundance", size =m)

#################### Sample completeness curves #######################
passaic_rarefaction.upstream <-ggiNEXT(out.up.1, type=1, color.var="site", se=FALSE) +
  #ylim(c(0.9,1)) +
  theme_bw(base_size = 10) + 
  ggtitle("Upstream Passaic: Parasite richness sample rarefaction")+
  theme(legend.position="bottom")+
   scale_color_manual(values=c( '#a6611a','#dfc27d','#80cdc1','#018571')#,
                    #labels = c("upstream","midstream", "downstream")
                   )+
#  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'),
 #                    labels = c("upstream","midstream", "downstream"))+
  guides(shape=FALSE)

######################## Coverage-based R/E curves ##################
passaic_coverage.up.1 <-ggiNEXT(out.up.1, type=3, color.var="site", se=FALSE) +
  #ylim(c(0.9,1)) +
  theme_bw(base_size = 10) + 
  ggtitle("Upstream Passaic: Parasite richness coverage rarefaction")+
  theme(legend.position="bottom")+
  scale_color_manual(values=c( '#a6611a','#dfc27d','#80cdc1','#018571')#,
                     #labels = c("upstream","midstream", "downstream")
  )+
  #  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'),
  #                    labels = c("upstream","midstream", "downstream"))+
  guides(shape=FALSE)

##############################################################################################
##############################################################################################

plist.midstream = list(
             fall = p_matrix[5,],
             spring = p_matrix[6,],
             summer = p_matrix[7,],
             winter = p_matrix[8,]
             )

out.mid.1 <- iNEXT(plist.midstream, q=c(0), datatype="abundance", size =m)

#################### Sample completeness curves #######################
passaic_rarefaction.midstream <-ggiNEXT(out.mid.1, type=1, color.var="site", se=FALSE) +
  #ylim(c(0.9,1)) +
  theme_bw(base_size = 10) + 
  ggtitle("Midstream Passaic: Parasite richness sample rarefaction")+
  theme(legend.position="bottom")+
  scale_color_manual(values=c( '#a6611a','#dfc27d','#80cdc1','#018571')#,
                     #labels = c("upstream","midstream", "downstream")
  )+
  #  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'),
  #                    labels = c("upstream","midstream", "downstream"))+
  guides(shape=FALSE)

######################## Coverage-based R/E curves ##################
passaic_coverage.mid.1 <-ggiNEXT(out.mid.1, type=3, color.var="site", se=FALSE) +
  #ylim(c(0.9,1)) +
  theme_bw(base_size = 10) + 
  ggtitle("Midstream Passaic: Parasite richness coverage rarefaction")+
  theme(legend.position="bottom")+
  scale_color_manual(values=c( '#a6611a','#dfc27d','#80cdc1','#018571')#,
                     #labels = c("upstream","midstream", "downstream")
  )+
  #  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'),
  #                    labels = c("upstream","midstream", "downstream"))+
  guides(shape=FALSE)

##############################################################################################
##############################################################################################

plist.downstream = list(
                       fall = p_matrix[9,],
                       spring = p_matrix[10,],
                       summer = p_matrix[11,],
                       winter = p_matrix[12,]
)

out.down.1 <- iNEXT(plist.downstream, q=c(0), datatype="abundance", size =m)

#################### Sample completeness curves #######################
passaic_rarefaction.downstream <-ggiNEXT(out.down.1, type=1, color.var="site", se=FALSE) +
  #ylim(c(0.9,1)) +
  theme_bw(base_size = 10) + 
  ggtitle("Downstream Passaic: Parasite richness sample rarefaction")+
  theme(legend.position="bottom")+
  scale_color_manual(values=c( '#a6611a','#dfc27d','#80cdc1','#018571')#,
                     #labels = c("upstream","midstream", "downstream")
  )+
  #  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'),
  #                    labels = c("upstream","midstream", "downstream"))+
  guides(shape=FALSE)

######################## Coverage-based R/E curves ##################
passaic_coverage.down.1 <-ggiNEXT(out.down.1, type=3, color.var="site", se=FALSE) +
  #ylim(c(0.9,1)) +
  theme_bw(base_size = 10) + 
  ggtitle("Downstream Passaic: Parasite richness coverage rarefaction")+
  theme(legend.position="bottom")+
  scale_color_manual(values=c( '#a6611a','#dfc27d','#80cdc1','#018571')#,
                     #labels = c("upstream","midstream", "downstream")
  )+
  #  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'),
  #                    labels = c("upstream","midstream", "downstream"))+
  guides(shape=FALSE)

jpeg(filename="passaic_parasite_coverage.jpeg", width=250, height=300,  units="mm", bg="white", res=300)
ggarrange(passaic_rarefaction.upstream, passaic_coverage.up.1,
          passaic_rarefaction.midstream, passaic_coverage.mid.1,
          passaic_rarefaction.downstream, passaic_coverage.down.1,
          ncol=2, nrow=3, common.legend = TRUE)
dev.off()

############################################################################################
############### RARITAN PARASITES #############################################################################
############################################################################################

#rlist = list(raritan1_upstream = p_matrix[4,],
  #           raritan2_midstream = p_matrix[5,],
 #            raritan3_downstream = p_matrix[6,]
#)
 
rlist = list(raritan1_upstream_fall = p_matrix[13,],
             raritan1_upstream_spring = p_matrix[14,],
            # raritan1_upstream_summer = p_matrix[15,],
             raritan1_upstream_winter = p_matrix[16,],
             raritan2_midstream_fall = p_matrix[17,],
             raritan2_midstream_spring = p_matrix[18,],
             raritan2_midstream_summer = p_matrix[19,],
             raritan2_midstream_winter = p_matrix[20,],
             raritan3_downstream_fall = p_matrix[21,],
             raritan3_downstream_spring = p_matrix[22,],
             raritan3_downstream_summer = p_matrix[23,],
             raritan3_downstream_winter = p_matrix[24,]
)

m2 <- c(1, 100, 200, 500, 1000, 5000, 7000, 10000, 15000)
out2 <- iNEXT(rlist, q=c(0), datatype="abundance", size=m2, nboot=999)

#### export info
raritan_richness <- estimateD(rlist,datatype="abundance", base="coverage", level=0.99, conf=0.95)

raritan_diversity2 <-out2$AsyEst #lists the observed diversity, asymptotic estimates, estimated bootstrap s.e. and 95% confidence intervals for Hill numbers with q = 0, 1, and 2.
raritan_richness2 <- raritan_diversity %>% filter(Diversity=='Species richness') 

#### graph richness 
rparasite <-raritan_richness %>%
  ggplot()+
  geom_point(aes(x=Estimator,y = Site), size=2)+
  geom_point(aes(x=Observed,y = Site),
             position = position_jitter( h = 0.25),
             pch=21, size=2)+
  labs(x="richness",y= "")+
  ggtitle("Parasite Raritan River")+
  theme_bw()

######### upstream site only 
rlist.upstream = list(fall = p_matrix[13,],
                      spring = p_matrix[14,],
                     # summer = p_matrix[15,],
                      winter = p_matrix[16,]
)
out.up.2 <- iNEXT(rlist.upstream, q=c(0), datatype="abundance", size =m2)

#################### Sample completeness curves #######################
raritan_rarefaction.upstream <-ggiNEXT(out.up.2, type=1, color.var="site", se=FALSE) +
  #ylim(c(0.9,1)) +
  theme_bw(base_size = 10) + 
  ggtitle("Upstream Raritan: Parasite richness sample rarefaction")+
  theme(legend.position="bottom")+
  scale_color_manual(values=c( '#a6611a','#dfc27d','#018571')#,
                     #labels = c("upstream","midstream", "downstream")
  )+
  #  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'),
  #                    labels = c("upstream","midstream", "downstream"))+
  guides(shape=FALSE)

######################## Coverage-based R/E curves ##################
raritan_coverage.up.1 <-ggiNEXT(out.up.2, type=3, color.var="site", se=FALSE) +
  #ylim(c(0.9,1)) +
  theme_bw(base_size = 10) + 
  ggtitle("Upstream Raritan: Parasite richness coverage rarefaction")+
  theme(legend.position="bottom")+
  scale_color_manual(values=c( '#a6611a','#dfc27d','#018571')#,
                     #labels = c("upstream","midstream", "downstream")
  )+
  #  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'),
  #                    labels = c("upstream","midstream", "downstream"))+
  guides(shape=FALSE)

##############################################################################################
##############################################################################################
rlist.midstream = list(
  fall = p_matrix[17,],
  spring = p_matrix[18,],
  summer = p_matrix[19,],
  winter = p_matrix[20,]
)

out.mid.2 <- iNEXT(rlist.midstream, q=c(0), datatype="abundance", size =m2)

#################### Sample completeness curves #######################
raritan_rarefaction.midstream <-ggiNEXT(out.mid.2, type=1, color.var="site", se=FALSE) +
  #ylim(c(0.9,1)) +
  theme_bw(base_size = 10) + 
  ggtitle("Midstream Raritan: Parasite richness sample rarefaction")+
  theme(legend.position="bottom")+
  scale_color_manual(values=c( '#a6611a','#dfc27d','#80cdc1','#018571')#,
                     #labels = c("upstream","midstream", "downstream")
  )+
  #  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'),
  #                    labels = c("upstream","midstream", "downstream"))+
  guides(shape=FALSE)

######################## Coverage-based R/E curves ##################
raritan_coverage.mid.1 <-ggiNEXT(out.mid.2, type=3, color.var="site", se=FALSE) +
  #ylim(c(0.9,1)) +
  theme_bw(base_size = 10) + 
  ggtitle("Midstream Raritan: Parasite richness coverage rarefaction")+
  theme(legend.position="bottom")+
  scale_color_manual(values=c( '#a6611a','#dfc27d','#80cdc1','#018571')#,
                     #labels = c("upstream","midstream", "downstream")
  )+
  #  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'),
  #                    labels = c("upstream","midstream", "downstream"))+
  guides(shape=FALSE)

##############################################################################################
##############################################################################################

rlist.downstream = list(
  fall = p_matrix[21,],
  spring = p_matrix[22,],
  summer = p_matrix[23,],
  winter = p_matrix[24,]
)

out.down.2 <- iNEXT(rlist.downstream, q=c(0), datatype="abundance", size =m2)

#################### Sample completeness curves #######################
raritan_rarefaction.downstream <-ggiNEXT(out.down.2, type=1, color.var="site", se=FALSE) +
  #ylim(c(0.9,1)) +
  theme_bw(base_size = 10) + 
  ggtitle("Downstream Raritan: Parasite richness sample rarefaction")+
  theme(legend.position="bottom")+
  scale_color_manual(values=c( '#a6611a','#dfc27d','#80cdc1','#018571')#,
                     #labels = c("upstream","midstream", "downstream")
  )+
  #  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'),
  #                    labels = c("upstream","midstream", "downstream"))+
  guides(shape=FALSE)

######################## Coverage-based R/E curves ##################
raritan_coverage.down.1 <-ggiNEXT(out.down.1, type=3, color.var="site", se=FALSE) +
  #ylim(c(0.9,1)) +
  theme_bw(base_size = 10) + 
  ggtitle("Downstream raritan: Parasite richness coverage rarefaction")+
  theme(legend.position="bottom")+
  scale_color_manual(values=c( '#a6611a','#dfc27d','#80cdc1','#018571')#,
                     #labels = c("upstream","midstream", "downstream")
  )+
  #  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'),
  #                    labels = c("upstream","midstream", "downstream"))+
  guides(shape=FALSE)

jpeg(filename="raritan_parasite_coverage.jpeg", width=250, height=300,  units="mm", bg="white", res=300)
ggarrange(raritan_rarefaction.upstream, raritan_coverage.up.1,
          raritan_rarefaction.midstream, raritan_coverage.mid.1,
          raritan_rarefaction.downstream, raritan_coverage.down.1,
          ncol=2, nrow=3, common.legend = TRUE)
dev.off()

############################################################################################
##################### FISH FISH FISH #######################################################################
############################################################################################

# abundance based FISH MATRIX 
fish_01 <- river %>%  
  mutate(abun = rep(1)) %>% 
  group_by(river, plot, season, subplot, host_species) %>%
  summarise_all(funs(if(is.numeric(.)) sum(., na.rm = TRUE) else first(.))) %>% 
  spread(host_species, abun)

fish_01[,25:ncol((fish_01))][is.na(fish_01[,25:ncol((fish_01))])] <- 0
fish <- fish_01 %>%
  group_by(river, plot, season) %>%
  summarise_all(funs(if(is.numeric(.)) sum(., na.rm = TRUE) else first(.))) 
fish <- fish[,-c(1:64)] # remove parasite data

fish <- as.matrix(fish)
#rownames(fish) <- c("passaic_upstream", "passaic_midstream", "passaic_downstream",
 #                   "raritan_upstream", "raritan_midstream", "raritan_downstream") 

#pflist = list(passaic1_upstream = fish[1,],
 #            passaic2_midstream = fish[2,],
  #           passaic3_downstream = fish[3,]
#)

pflist = list(passaic1_upstream_fall = fish[1,],
             passaic1_upstream_spring = fish[2,],
             passaic1_upstream_summer = fish[3,],
             passaic1_upstream_winter = fish[4,],
             passaic2_midstream_fall = fish[5,],
             passaic2_midstream_spring = fish[6,],
             passaic2_midstream_summer = fish[7,],
             passaic2_midstream_winter = fish[8,],
             passaic3_downstream_fall = fish[9,],
             passaic3_downstream_spring = fish[10,],
             passaic3_downstream_summer = fish[11,],
             passaic3_downstream_winter = fish[12,]
)
m3 <- c(1, 10, 50, 100, 125, 150, 200)

out3 <- iNEXT(pflist, q=c(0), datatype ="abundance", nboot=999, size=m3)

passaicfish_richness <- estimateD(pflist,datatype="abundance", base="coverage", level=0.99, conf=0.95)

#### export info
passaicfish_diversity2 <-out3$AsyEst #lists the observed diversity, asymptotic estimates, estimated bootstrap s.e. and 95% confidence intervals for Hill numbers with q = 0, 1, and 2.
passaicfish_richness2 <- passaicfish_diversity %>% filter(Diversity=='Species richness') 


#### graph richness 
pfish <-passaicfish_richness %>%
  ggplot()+
  geom_point(aes(x=Estimator,y = Site),size=2)+
  geom_point(aes(x=Observed,y = Site),        
             position = position_jitter( h = 0.25),
             pch=21, size=2)+
  labs(x="richness",y= "")+
  ggtitle("Fish Passaic River")+
  theme_bw()

######### upstream site only 
pflist.upstream = list(fall = fish[1,],
                      spring = fish[2,],
                      summer = fish[3,],
                      winter = fish[4,]
)
out.up.1f <- iNEXT(pflist.upstream, q=c(0), datatype="abundance", size =m3)

#################### Sample completeness curves #######################
passaic_rarefaction.upstream.fish <-ggiNEXT(out.up.1f, type=1, color.var="site", se=FALSE) +
  #ylim(c(0.9,1)) +
  theme_bw(base_size = 10) + 
  ggtitle("Upstream Passaic: Fish richness sample rarefaction")+
  theme(legend.position="bottom")+
  scale_color_manual(values=c( '#a6611a','#dfc27d','#80cdc1','#018571')#,
                     #labels = c("upstream","midstream", "downstream")
  )+
  #  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'),
  #                    labels = c("upstream","midstream", "downstream"))+
  guides(shape=FALSE)

######################## Coverage-based R/E curves ##################
passaic_coverage.up.1.fish <-ggiNEXT(out.up.1f, type=3, color.var="site", se=FALSE) +
  #ylim(c(0.9,1)) +
  theme_bw(base_size = 10) + 
  ggtitle("Upstream Passaic: Fish richness coverage rarefaction")+
  theme(legend.position="bottom")+
  scale_color_manual(values=c( '#a6611a','#dfc27d','#80cdc1','#018571')#,
                     #labels = c("upstream","midstream", "downstream")
  )+
  #  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'),
  #                    labels = c("upstream","midstream", "downstream"))+
  guides(shape=FALSE)

##############################################################################################
##############################################################################################

pflist.midstream = list(
  fall = fish[5,],
  spring = fish[6,],
  summer = fish[7,],
  winter = fish[8,]
)

out.mid.1f <- iNEXT(pflist.midstream, q=c(0), datatype="abundance", size =m3)

#################### Sample completeness curves #######################
passaic_rarefaction.midstream.fish <-ggiNEXT(out.mid.1f, type=1, color.var="site", se=FALSE) +
  #ylim(c(0.9,1)) +
  theme_bw(base_size = 10) + 
  ggtitle("Midstream Passaic: Fish richness sample rarefaction")+
  theme(legend.position="bottom")+
  scale_color_manual(values=c( '#a6611a','#dfc27d','#80cdc1','#018571')#,
                     #labels = c("upstream","midstream", "downstream")
  )+
  #  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'),
  #                    labels = c("upstream","midstream", "downstream"))+
  guides(shape=FALSE)

######################## Coverage-based R/E curves ##################
passaic_coverage.mid.1.fish <-ggiNEXT(out.mid.1f, type=3, color.var="site", se=FALSE) +
  #ylim(c(0.9,1)) +
  theme_bw(base_size = 10) + 
  ggtitle("Midstream Passaic: Fish richness coverage rarefaction")+
  theme(legend.position="bottom")+
  scale_color_manual(values=c( '#a6611a','#dfc27d','#80cdc1','#018571')#,
                     #labels = c("upstream","midstream", "downstream")
  )+
  #  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'),
  #                    labels = c("upstream","midstream", "downstream"))+
  guides(shape=FALSE)

##############################################################################################
##############################################################################################

pflist.downstream = list(
  fall = fish[9,],
  spring = fish[10,],
  summer = fish[11,],
  winter = fish[12,]
)

out.down.1f<- iNEXT(pflist.downstream, q=c(0), datatype="abundance", size =m3)

#################### Sample completeness curves #######################
passaic_rarefaction.downstream.fish <-ggiNEXT(out.down.1f, type=1, color.var="site", se=FALSE) +
  #ylim(c(0.9,1)) +
  theme_bw(base_size = 10) + 
  ggtitle("Downstream Passaic: Fish richness sample rarefaction")+
  theme(legend.position="bottom")+
  scale_color_manual(values=c( '#a6611a','#dfc27d','#80cdc1','#018571')#,
                     #labels = c("upstream","midstream", "downstream")
  )+
  #  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'),
  #                    labels = c("upstream","midstream", "downstream"))+
  guides(shape=FALSE)

######################## Coverage-based R/E curves ##################
passaic_coverage.down.1.fish <-ggiNEXT(out.down.1f, type=3, color.var="site", se=FALSE) +
  #ylim(c(0.9,1)) +
  theme_bw(base_size = 10) + 
  ggtitle("Downstream Passaic: Fish richness coverage rarefaction")+
  theme(legend.position="bottom")+
  scale_color_manual(values=c( '#a6611a','#dfc27d','#80cdc1','#018571')#,
                     #labels = c("upstream","midstream", "downstream")
  )+
  #  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'),
  #                    labels = c("upstream","midstream", "downstream"))+
  guides(shape=FALSE)

jpeg(filename="passaic_fish_coverage.jpeg", width=250, height=300,  units="mm", bg="white", res=300)
ggarrange(passaic_rarefaction.upstream.fish, passaic_coverage.up.1.fish,
          passaic_rarefaction.midstream.fish, passaic_coverage.mid.1.fish,
          passaic_rarefaction.downstream.fish, passaic_coverage.down.1.fish,
          ncol=2, nrow=3, common.legend = TRUE)
dev.off()

############################################################################################
############################################################################################

#rflist = list(raritan1_upstream = fish[4,],
 #             raritan2_midstream = fish[5,],
  #            raritan3_downstream = fish[6,]
#)

rflist = list(raritan1_upstream_fall = fish[13,],
             raritan1_upstream_spring = fish[14,],
             # raritan1_upstream_summer = fish[15,],
             raritan1_upstream_winter = fish[16,],
             raritan2_midstream_fall = fish[17,],
             raritan2_midstream_spring = fish[18,],
             raritan2_midstream_summer = fish[19,],
             raritan2_midstream_winter = fish[20,],
             raritan3_downstream_fall = fish[21,],
             raritan3_downstream_spring = fish[22,],
             raritan3_downstream_summer = fish[23,],
             raritan3_downstream_winter = fish[24,]
)

m4 <- c(1, 10, 50, 100,  150, 200, 400, 600, 1200)

out4 <- iNEXT(rflist, q=c(0), datatype="abundance", nboot=999, size=m4)
out4$DataInfo

raritanfish_richness <- estimateD(rflist,datatype="abundance", base="coverage", level=0.985, conf=0.95)

#### export info
raritanfish_diversity2 <-out4$AsyEst #lists the observed diversity, asymptotic estimates, estimated bootstrap s.e. and 95% confidence intervals for Hill numbers with q = 0, 1, and 2.
raritanfish_richness2 <- raritanfish_diversity %>% filter(Diversity=='Species richness') 

#### graph richness 
rfish <-raritanfish_richness %>%
  ggplot()+
  geom_point(aes(x=Estimator,y = Site), size=2)+
  geom_point(aes(x=Observed,y = Site),         
             position = position_jitter( h = 0.25),
             pch=21, size=2)+
  labs(x="richness",y= "")+
  ggtitle("Fish Raritan River")+
  theme_bw()

############################################################################################
######### upstream site only 
rflist.upstream = list(fall = fish[13,],
                      spring = fish[14,],
                      # summer = fish[15,],
                      winter = fish[16,]
)
out.up.2f <- iNEXT(rflist.upstream, q=c(0), datatype="abundance", size =m4)

#################### Sample completeness curves #######################
raritan_rarefaction.upstream.fish <-ggiNEXT(out.up.2f, type=1, color.var="site", se=FALSE) +
  #ylim(c(0.9,1)) +
  theme_bw(base_size = 10) + 
  ggtitle("Upstream Raritan: Fish richness sample rarefaction")+
  theme(legend.position="bottom")+
  scale_color_manual(values=c( '#a6611a','#dfc27d','#018571')#,
                     #labels = c("upstream","midstream", "downstream")
  )+
  #  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'),
  #                    labels = c("upstream","midstream", "downstream"))+
  guides(shape=FALSE)

######################## Coverage-based R/E curves ##################
raritan_coverage.up.1.fish <-ggiNEXT(out.up.2f, type=3, color.var="site", se=FALSE) +
  #ylim(c(0.9,1)) +
  theme_bw(base_size = 10) + 
  ggtitle("Upstream Raritan: Fish richness coverage rarefaction")+
  theme(legend.position="bottom")+
  scale_color_manual(values=c( '#a6611a','#dfc27d','#018571')#,
                     #labels = c("upstream","midstream", "downstream")
  )+
  #  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'),
  #                    labels = c("upstream","midstream", "downstream"))+
  guides(shape=FALSE)

##############################################################################################
##############################################################################################
rflist.midstream = list(
  fall = fish[17,],
  spring = fish[18,],
  summer = fish[19,],
  winter = fish[20,]
)

out.mid.2f <- iNEXT(rflist.midstream, q=c(0), datatype="abundance", size =m4)

#################### Sample completeness curves #######################
raritan_rarefaction.midstream.fish <-ggiNEXT(out.mid.2f, type=1, color.var="site", se=FALSE) +
  #ylim(c(0.9,1)) +
  theme_bw(base_size = 10) + 
  ggtitle("Midstream Raritan: Fish richness sample rarefaction")+
  theme(legend.position="bottom")+
  scale_color_manual(values=c( '#a6611a','#dfc27d','#80cdc1','#018571')#,
                     #labels = c("upstream","midstream", "downstream")
  )+
  #  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'),
  #                    labels = c("upstream","midstream", "downstream"))+
  guides(shape=FALSE)

######################## Coverage-based R/E curves ##################
raritan_coverage.mid.1.fish <-ggiNEXT(out.mid.2f, type=3, color.var="site", se=FALSE) +
  #ylim(c(0.9,1)) +
  theme_bw(base_size = 10) + 
  ggtitle("Midstream Raritan: Fish richness coverage rarefaction")+
  theme(legend.position="bottom")+
  scale_color_manual(values=c( '#a6611a','#dfc27d','#80cdc1','#018571')#,
                     #labels = c("upstream","midstream", "downstream")
  )+
  #  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'),
  #                    labels = c("upstream","midstream", "downstream"))+
  guides(shape=FALSE)

##############################################################################################
##############################################################################################

rflist.downstream = list(
  fall = fish[21,],
  spring = fish[22,],
  summer = fish[23,],
  winter = fish[24,]
)

out.down.2f <- iNEXT(rflist.downstream, q=c(0), datatype="abundance", size =m4)

#################### Sample completeness curves #######################
raritan_rarefaction.downstream.fish <-ggiNEXT(out.down.2f, type=1, color.var="site", se=FALSE) +
  #ylim(c(0.9,1)) +
  theme_bw(base_size = 10) + 
  ggtitle("Downstream Raritan: Fish richness sample rarefaction")+
  theme(legend.position="bottom")+
  scale_color_manual(values=c( '#a6611a','#dfc27d','#80cdc1','#018571')#,
                     #labels = c("upstream","midstream", "downstream")
  )+
  #  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'),
  #                    labels = c("upstream","midstream", "downstream"))+
  guides(shape=FALSE)

######################## Coverage-based R/E curves ##################
raritan_coverage.down.1.fish <-ggiNEXT(out.down.2f, type=3, color.var="site", se=FALSE) +
  #ylim(c(0.9,1)) +
  theme_bw(base_size = 10) + 
  ggtitle("Downstream raritan: Fish richness coverage rarefaction")+
  theme(legend.position="bottom")+
  scale_color_manual(values=c( '#a6611a','#dfc27d','#80cdc1','#018571')#,
                     #labels = c("upstream","midstream", "downstream")
  )+
  #  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'),
  #                    labels = c("upstream","midstream", "downstream"))+
  guides(shape=FALSE)

jpeg(filename="raritan_fish_coverage.jpeg", width=250, height=300,  units="mm", bg="white", res=300)
ggarrange(raritan_rarefaction.upstream.fish, raritan_coverage.up.1.fish,
          raritan_rarefaction.midstream.fish, raritan_coverage.mid.1.fish,
          raritan_rarefaction.downstream.fish, raritan_coverage.down.1.fish,
          ncol=2, nrow=3, common.legend = TRUE)
dev.off()
##########################################################################

jpeg(filename="richness_estimator_values.jpeg", width=270, height=180, units="mm", bg="white", res=300)
ggarrange(pfish, pparasite,rfish, rparasite, ncol=2, nrow=2, labels='AUTO')
dev.off()

############################################################################################
# extract diversity estimates based on sample coverage and put into one dataframe
############################################################################################

# select q = 0; for richness only rarefaction and rename col headers 
#raritanfish.est <- raritanfish_richness %>% mutate(fish_richness = qD)%>% filter(order==0) %>% dplyr::select(site, fish_richness) 
#raritanparasiet.est <- raritan_richness %>% mutate(parasite_richness = qD)%>% filter(order==0) %>% dplyr::select(site, parasite_richness) 
#passaicfish.est <- passaicfish_richness %>% mutate(fish_richness = qD)%>% filter(order==0) %>% dplyr::select(site, fish_richness) 
#passaicparasiet.est <- passaic_richness %>% mutate(parasite_richness = qD)%>% filter(order==0) %>% dplyr::select(site, parasite_richness) 

raritan_rich_est_df <- merge(raritanfish_richnes, raritan_richness, by=c("Site")) %>% mutate(river ="Raritan")
passaic_rich_est_df <- merge(passaicfish_richness, passaic_richness, by=c("Site")) %>% mutate(river ="Passaic")

rich_est_all <- rbind(raritan_rich_est_df, passaic_rich_est_df) 
head(rich_est_all)

obs <-rich_est_all %>%
  ggplot()+
  geom_point(aes(x=Observed.x, y= Observed.y, fill=river), pch=21,size=3)+
  theme_bw(base_size = 15) +
  labs(x="fish richness (obs.)", y ="parasite richness (obs.)")+
  scale_fill_manual(values=c('white', 'black'),
                    guide=guide_legend(override.aes=list(shape=21))) +
  guides(fill=guide_legend("River",  keywidth=0.12,
                           keyheight=0.12, override.aes = list(shape=21),
                           default.unit="inch"), color= FALSE) +
  theme(legend.title = element_text(size=7, face="bold"),
        legend.text = element_text(size=6.5))+ 
  xlim(c(0,15))+
  ylim(c(0,27))+
  theme_bw()

obs.rltn <-lm(Observed.y~Observed.x, data= rich_est_all)
summary(obs.rltn)

est <-rich_est_all %>%
  ggplot()+
  geom_point(aes(x=Estimator.x, y= Estimator.ys, fill=river), size=3,pch=21)+
  theme_bw(base_size = 15) +
  labs(x="fish richness (est.)", y ="parasite richness (est.)")+
  scale_fill_manual(values=c('white', 'black'),
                    guide=guide_legend(override.aes=list(shape=21))) +
  guides(fill=guide_legend("River",  keywidth=0.12,
                           keyheight=0.12, override.aes = list(shape=21),
                           default.unit="inch"), color= FALSE) +
  theme(legend.title = element_text(size=7, face="bold"),
        legend.text = element_text(size=6.5))+ 
 # xlim(c(0,20))+
#  ylim(c(0,34))+
  theme_bw()

est.rltn <-lm(Estimator.y~Estimator.x, data= rich_est_all)
summary(est.rltn)

jpeg(filename="obs_est_richness_relationship.jpeg", width=180, height=90,units="mm", bg="white", res=300)
ggarrange(obs, est, labels='AUTO',common.legend = TRUE)
dev.off()
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

