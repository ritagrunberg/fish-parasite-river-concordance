rm(list = ls())  # clear Rs brain 
###########################################################################################
library(readr)   # read csv file
library(tidyr)   # tidying data
library(dplyr)   # manipulating df
library(ggplot2) # graphing
library(ggpubr)  # nice package for figures
library(vegan)   # nmds 
library(bipartite)
library(tidyverse)
library(circlize)#chord diagram
library(RColorBrewer) #color palette
###########################################################################################
# create host-parasite bipartite networks for each river
###########################################################################################

# import dataset
river <- read_csv("C:/Users/grunberg/Dropbox/DISSERTATION_CHAPTERS/1_parasite_community_spatial/data/final_river_data.csv")
river <- river[,-c(38,57)] # remove notes that were in excel file... trash 
set.seed(1234567890)
river$plot <- as.factor(river$plot)
river$subplot <- as.factor(river$subplot)
river <- river %>% dplyr::select(-c(sex))
river <- river %>% drop_na()
################################################################################################
#host and parasite bipartie networks 
##########################################################################################
parasite <- river  %>% 
  group_by(host_species, river) %>%
  summarise_all(funs(if(is.numeric(.)) sum(., na.rm = TRUE) else first(.))) %>%
  select(c(1,2, 27:64)) %>% 
  mutate_if(is.numeric, ~1 * (. != 0)) %>% 
  drop_na() %>% 
  gather(parasite, infection, 3:40)# %>%
  #filter(infection >0)
unique(parasite$parasite)
#changing names to look better for plotting 

parasite$parasite<-gsub("nematode_sp1","Larval nematode 1",  parasite$parasite)
parasite$parasite<-gsub("p_minimum","Posthodiplostomum minimum",  parasite$parasite)
parasite$parasite<-gsub("acanthocephalus_sp","Acanthocephalus sp.",  parasite$parasite)
parasite$parasite<-gsub("alloglossidium","Alloglossidium sp.",  parasite$parasite)
parasite$parasite<-gsub("bothriocephalus","Bothriocephalus sp.",  parasite$parasite)
parasite$parasite<-gsub("camallanus_adult","Camallanus sp.",  parasite$parasite)
parasite$parasite<-gsub("clinostomum_marginatum","Clinostomum marginatum",  parasite$parasite)
parasite$parasite<-gsub("crepidostomum_cooperi","Crepidostomum cooperi",  parasite$parasite)
parasite$parasite<-gsub("crepidostomum_isostomum","Crepidostomum isostomum",  parasite$parasite)
parasite$parasite<-gsub("echinorhynchus_adult","Echinorhynchus sp.",  parasite$parasite)
parasite$parasite<-gsub("eocolis_sp","Eocolis sp.",  parasite$parasite)
parasite$parasite<-gsub("eustronglyoides_sp","Eustronglyoides sp.",  parasite$parasite)
parasite$parasite<-gsub("Fessesentis_adult","Fessesentis sp.",  parasite$parasite)
parasite$parasite<-gsub("jelly_metacercariae","Metacercariae 1 (Heterophyidae)",  parasite$parasite)
parasite$parasite<-gsub("leptorhynchoides_larval","Leptorhynchoides sp. (la)",  parasite$parasite)
parasite$parasite<-gsub("leptorhynchoids_adult","Leptorhynchoides thecatus (ad)",  parasite$parasite)
parasite$parasite<-gsub("neoechinorhynchus_cristatus","Neoechinorhynchus cristatus",  parasite$parasite)
parasite$parasite<-gsub("neoechinorhynchus_cylindratus_ad","Neoechinorhynchus cylindratus (ad)",  parasite$parasite)
parasite$parasite<-gsub("neoechinorhynchus_larval","Neoechinorhynchus sp. (la)",  parasite$parasite)
parasite$parasite<-gsub("phyllodistomum","Phyllodistomum sp.",  parasite$parasite)
parasite$parasite<-gsub("plagioporus","Plagioporus sp.",  parasite$parasite)
parasite$parasite<-gsub("proteocephalus","Proteocephalus sp.",  parasite$parasite)
parasite$parasite<-gsub("rhabdochona","Rhabdochona sp.",  parasite$parasite)
parasite$parasite<-gsub("strigeoidae_larval","Metacercariae 2 (Strigeodiae)",  parasite$parasite)
parasite$parasite<-gsub("Trematode_sp1","Adult trematode 1",  parasite$parasite)
parasite$parasite<-gsub("triganodistomum","Triganodistomum sp.",  parasite$parasite)
parasite$parasite<-gsub("uvulifer_ambloplitis","Uvulifer ambloplitis",  parasite$parasite)
parasite$parasite<-gsub("Caecincola"," Caecincola sp.",  parasite$parasite)
parasite$parasite<-gsub("Philonema","Philonema sp.",  parasite$parasite)
parasite$parasite<-gsub("Raphidascaris","Raphidascaris sp.",  parasite$parasite)
parasite$parasite<-gsub("Isoglaridacris","Isoglaridacris sp.",  parasite$parasite)
parasite$parasite<-gsub("allocreadium_commune","Allocreadium commune",  parasite$parasite)
parasite$parasite<-gsub("red_intestinal_nematode","Adult red nematode 1",  parasite$parasite)
parasite$parasite<-gsub("pomporhynchus_adult","Pomporhynchus sp.",  parasite$parasite)
parasite$parasite<-gsub("pilum_adult","Pilum sp.",  parasite$parasite)
parasite$parasite<-gsub("Spinitectus","Spinitectus sp.",  parasite$parasite)
parasite$parasite<-gsub("tapeworm_sp2","Adult tapeworm 1",  parasite$parasite)
parasite$host_species<-gsub("Etheostoma olmstedi_","Etheostoma olmstedi",  parasite$host_species)

## chord diagram for network
parasite_passaic <- parasite %>%
  filter(river=="Passaic") %>% 
  select(-c(river)) %>%
  group_by(parasite, host_species) %>% 
  summarise_at(c("infection"), sum)

col_mat = rand_color(length(parasite_passaic), transparency = 0.5)


#### chord diagram instead of traditional diagram
jpeg(filename="chorddiagram_passaic_river.jpeg", width=180, height=180,units="mm", bg="white", res=300)
plot.new()
par(cex = 0.5, mar = c(0, 0, 0, 0))
chordDiagram(parasite_passaic, preAllocateTracks = 1, annotationTrack = "grid", directional = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.5, major.tick.percentage = 0.1, sector.index = sector.name, track.index = 2)
}, bg.border = NA)
dev.off()

#matrix for passaic river 
par_mat_passaic <- parasite %>% 
  filter(river=="Passaic") %>% 
  group_by(parasite) %>%
  filter(sum(infection)>0) %>%
  spread(parasite, infection)%>% remove_rownames() %>%
  column_to_rownames(var = 'host_species') %>% select(-c(river))

colSums(par_mat_passaic)
rowSums(par_mat_passaic)

#jpeg(filename="bipartie_passaic_river.jpeg", width=200, height=200,units="mm", bg="white", res=300)
plotweb(par_mat_passaic, y.width.low=0.05, y.width.high=0.05,  text.rot = 90,y.lim=c(-0.75,3.4),
        arrow="down",col.high="orange")
#dev.off()

#####################################################
#raritan
#####################################################

#matrix for raritan river 
## chord diagram for network
parasite_raritan <- parasite %>%
  filter(river=="Raritan") %>% 
  select(-c(river)) %>%
  group_by(parasite, host_species) %>% 
  summarise_at(c("infection"), sum)

#### chord diagram instead of traditional diagram

jpeg(filename="chorddiagram_raritan_river.jpeg", width=180, height=180,units="mm", bg="white", res=300)
plot.new()
par(cex = 0.5, mar = c(0, 0, 0, 0))
chordDiagram(parasite_raritan, preAllocateTracks = 1, annotationTrack = "grid", directional = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.5, major.tick.percentage = 0.1, sector.index = sector.name, track.index = 2)
}, bg.border = NA)
dev.off()

par_mat_raritan <- parasite %>% 
  filter(river=="Raritan") %>% 
  group_by(parasite) %>%
  filter(sum(infection)>0) %>%
  spread(parasite, infection)%>% remove_rownames() %>%
  column_to_rownames(var = 'host_species') %>% select(-c(river))

colSums(par_mat_raritan)
rowSums(par_mat_raritan)

#jpeg(filename="bipartie_raritan_river_sup.jpeg", width=200, height=200, units="mm", bg="white", res=300)
plotweb(par_mat_raritan, y.width.low=0.05, y.width.high=0.05, y.lim=c(-0.75,3.4),
        text.rot = 90, arrow="down", col.high="orange")
#dev.off()

outer = FALSE
line = 1.5
cex = 1.25
adj  = 0.5

#plot both networks in one graphic 

jpeg(filename="bipartie_networks.jpeg", width=180, height=300, units="mm", bg="white", res=300)
plot.new()
par(mfrow=c(2,1))
#layout(mat, c(3.5,1), c(1,3))
par(mar=c(0, 0, 0, 0))
plotweb(par_mat_passaic, y.width.low=0.05, y.width.high=0.05, text.rot = 90,y.lim=c(-0.3,2.85),
                  arrow="down",col.high="grey")
title(outer=outer,adj=adj,main="Passaic River",cex.main=cex,col="black",font=2,line=line)
plotweb(par_mat_raritan, y.width.low=0.05, y.width.high=0.05, y.lim=c(-0.3,2.85),
                  text.rot = 90, arrow="down", col.high="grey")
title(outer=outer,adj=adj,main="Raritan River",cex.main=cex,col="black",font=2,line=line)
dev.off()
#############################################################
### look over all parasite interactions 
par_mat_all <- parasite %>% 
 # filter(river=="Raritan") %>% 
  group_by(host_species, parasite) %>%
  summarise_all(funs(if(is.numeric(.)) sum(., na.rm = TRUE) else first(.))) %>%
  group_by(parasite, host_species) %>%
  filter(sum(infection)>0) %>%
  spread(parasite, infection)%>% remove_rownames() %>%
  column_to_rownames(var = 'host_species') %>% select(-c(river))

head(par_mat_all)
colSums(par_mat_all)
rowSums(par_mat_all)

plotweb(par_mat_all, y.width.low=0.05, y.width.high=0.05, text.rot = 90,y.lim=c(-0.3,2.85),
        arrow="down",col.high="grey")
