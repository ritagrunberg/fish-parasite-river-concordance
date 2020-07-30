rm(list = ls())  # clear Rs brain 
###########################################################################################
library(readr)   # read csv file
library(tidyr)   # tidying data
library(dplyr)   # manipulating df
library(ggplot2) # graphing
library(ggpubr)  # nice package for figures
library(vegan)   # nmds 
###########################################################################################
river <- read_csv("C:/Users/grunberg/Dropbox/DISSERTATION_CHAPTERS/1_parasite_community_spatial/data/final_river_data.csv")
river <- river[,-c(38,57)] # remove notes that were in excel file... trash 
set.seed(1234567890)
river$plot <- as.factor(river$plot)
river$subplot <- as.factor(river$subplot)
river <- river %>% dplyr::select(-c(sex))
river <- river %>% drop_na()
river$host_species <-gsub("Etheostoma olmstedi_","Etheostoma olmstedi",  river$host_species)
########################################################################################################
parasite_dist <- river  %>% 
  gather(parasite, abundance, 27:64)%>%
  filter(abundance >0) %>%
  group_by(river, plot, parasite) %>%
  summarise_at(c("abundance"), sum) %>%
  mutate(ranking = dense_rank(desc(abundance)),
         ranking = as.factor(ranking))

parasite_dist_relative <- river  %>% 
  gather(parasite, abundance, 27:64)%>%
  filter(abundance >0) %>%
  group_by(river, plot, season) %>%
  mutate(total_abundance = sum(abundance))%>% 
  ungroup()%>%
  group_by(river, plot, season, parasite) %>%
  mutate(sp_abundance = sum(abundance),
         relative_abundance = sp_abundance/total_abundance) %>%
  group_by(river, plot, season, parasite) %>%
  summarise_at(c("relative_abundance"), mean) %>%
  mutate(ranking = dense_rank(desc(relative_abundance)),
         ranking = as.factor(ranking))%>%
  unite(plot_season, plot, season,season, sep = '_') #merge characters from columns into one thing

###############################################################################################
parasite_dist_relative$parasite<-gsub("nematode_sp1","Larval nematode 1",  parasite_dist_relative$parasite)
parasite_dist_relative$parasite<-gsub("p_minimum","Posthodiplostomum minimum",  parasite_dist_relative$parasite)
parasite_dist_relative$parasite<-gsub("acanthocephalus_sp","Acanthocephalus sp.",  parasite_dist_relative$parasite)
parasite_dist_relative$parasite<-gsub("alloglossidium","Alloglossidium sp.",  parasite_dist_relative$parasite)
parasite_dist_relative$parasite<-gsub("bothriocephalus","Bothriocephalus sp.",  parasite_dist_relative$parasite)
parasite_dist_relative$parasite<-gsub("camallanus_adult","Camallanus sp.",  parasite_dist_relative$parasite)
parasite_dist_relative$parasite<-gsub("clinostomum_marginatum","Clinostomum marginatum",  parasite_dist_relative$parasite)
parasite_dist_relative$parasite<-gsub("crepidostomum_cooperi","Crepidostomum cooperi",  parasite_dist_relative$parasite)
parasite_dist_relative$parasite<-gsub("crepidostomum_isostomum","Crepidostomum isostomum",  parasite_dist_relative$parasite)
parasite_dist_relative$parasite<-gsub("echinorhynchus_adult","Echinorhynchus sp.",  parasite_dist_relative$parasite)
parasite_dist_relative$parasite<-gsub("eocolis_sp","Eocolis sp.",  parasite_dist_relative$parasite)
parasite_dist_relative$parasite<-gsub("eustronglyoides_sp","Eustronglyoides sp.",  parasite_dist_relative$parasite)
parasite_dist_relative$parasite<-gsub("Fessesentis_adult","Fessesentis sp.",  parasite_dist_relative$parasite)
parasite_dist_relative$parasite<-gsub("jelly_metacercariae","Metacercariae 1 (Heterophyidae)",  parasite_dist_relative$parasite)
parasite_dist_relative$parasite<-gsub("leptorhynchoides_larval","Leptorhynchoides (larval)",  parasite_dist_relative$parasite)
parasite_dist_relative$parasite<-gsub("leptorhynchoids_adult","Leptorhynchoides thecatus (adult)",  parasite_dist_relative$parasite)
parasite_dist_relative$parasite<-gsub("neoechinorhynchus_cristatus","Neoechinorhynchus cristatus",  parasite_dist_relative$parasite)
parasite_dist_relative$parasite<-gsub("neoechinorhynchus_cylindratus_ad","Neoechinorhynchus cylindratus (adult)",  parasite_dist_relative$parasite)
parasite_dist_relative$parasite<-gsub("neoechinorhynchus_larval","Neoechinorhynchus (larval)",  parasite_dist_relative$parasite)
parasite_dist_relative$parasite<-gsub("phyllodistomum","Phyllodistomum sp.",  parasite_dist_relative$parasite)
parasite_dist_relative$parasite<-gsub("plagioporus","Plagioporus sp.",  parasite_dist_relative$parasite)
parasite_dist_relative$parasite<-gsub("proteocephalus","Proteocephalus sp.",  parasite_dist_relative$parasite)
parasite_dist_relative$parasite<-gsub("rhabdochona","Rhabdochona sp.",  parasite_dist_relative$parasite)
parasite_dist_relative$parasite<-gsub("strigeoidae_larval","Metacercariae 2 (Strigeodiae)",  parasite_dist_relative$parasite)
parasite_dist_relative$parasite<-gsub("Trematode_sp1","Adult trematode 1",  parasite_dist_relative$parasite)
parasite_dist_relative$parasite<-gsub("triganodistomum","Triganodistomum sp.",  parasite_dist_relative$parasite)
parasite_dist_relative$parasite<-gsub("uvulifer_ambloplitis","Uvulifer ambloplitis",  parasite_dist_relative$parasite)
parasite_dist_relative$parasite<-gsub("allocreadium_commune","Allocreadium commune",  parasite_dist_relative$parasite)
parasite_dist_relative$parasite<-gsub("red_intestinal_nematode","Adult red nematode 1",  parasite_dist_relative$parasite)
parasite_dist_relative$parasite<-gsub("pomporhynchus_adult","Pomporhynchus sp.",  parasite_dist_relative$parasite)
parasite_dist_relative$parasite<-gsub("pilum_adult","Pilum sp.",  parasite_dist_relative$parasite)
parasite_dist_relative$parasite<-gsub("tapeworm_sp2","Adult tapeworm 1",  parasite_dist_relative$parasite)
##########################################################################################################
jpeg(filename="parasite_relative_cover_matrix.jpeg", width=180, height=140, units="mm", bg="white", res=300)
parasite_dist_relative %>% 
  ggplot(aes(x=plot_season, y = parasite, fill=(relative_abundance)))+
  geom_tile(color = "white")+
  geom_raster(aes(x=plot_season, y = parasite, fill=(relative_abundance)))+
  scale_fill_gradientn(colours = c("white", "#ccebc5", "#084081"), 
                       values = c(0,0.00000000001,1), 
                       name="relative abun.") +
  labs(x="", y="") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    #  legend.position = c(1.5, 0.8),
    legend.direction = "vertical"
  )+
  guides(fill = guide_colorbar(barwidth = 1, barheight = 7, title.position = "top", title.hjust = 0.5))+
  facet_wrap(~river)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  geom_vline(xintercept = 4.5)+
  geom_vline(xintercept = 8.5)+
  scale_x_discrete(name="",
                   labels=c("1_fall" = "upstream (fall)", 
                            "1_spring" = "upstream (spring)",
                            "1_summer"="upstream (summer)",
                            "1_winter"="upstream (winter)",
                            "2_fall" = "midstream (fall)", 
                            "2_spring" = "midstream (spring)",
                            "2_summer"="midstream (summer)",
                            "2_winter"="midstream (winter)",
                            "3_fall" = "downstream (fall)", 
                            "3_spring" = "downstream (spring)",
                            "3_summer"="downstream (summer)",
                            "3_winter"="downstream (winter)")
  )
dev.off()

###################################################################################3
fish_relative_abundance <-fish %>%
  gather(fish, abundance, 27:50) %>%
  # group_by(river, plot, season, fish) %>%
  # filter(abundance >0) %>%
  group_by(river, plot, season) %>%
  mutate(total_abundance = sum(abundance))%>% 
  ungroup()%>%
  group_by(river, plot, season, fish) %>%
  mutate(sp_abundance = sum(abundance),
         relative_abundance = sp_abundance/total_abundance) %>%
  #filter(sp_abundance >0) %>%
  group_by(river, plot, season, fish) %>%
  summarise_at(c("relative_abundance"), mean) %>%
  mutate(ranking = dense_rank(desc(relative_abundance)),
         ranking = as.factor(ranking))%>%
 unite(plot_season, plot, season,season, sep = '_') #merge characters from columns into one thing

jpeg(filename="fish_relative_cover_matrix.jpeg", width=180, height=120, units="mm", bg="white", res=300)
fish_relative_abundance %>% 
  ggplot(aes(x=plot_season, y = fish, fill=(relative_abundance)))+
  geom_tile(color = "white")+
  geom_raster(aes(x=plot_season, y = fish, fill=(relative_abundance)))+
  scale_fill_gradientn(colours = c("white", "#ccebc5", "#084081"), 
                      values = c(0,0.00000000001,1), 
                      name="relative abun.") +
  labs(x="", y="") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    #  legend.position = c(1.5, 0.8),
    legend.direction = "vertical"
  )+
  guides(fill = guide_colorbar(barwidth = 1, barheight = 7, title.position = "top", title.hjust = 0.5))+
  facet_wrap(~river)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  geom_vline(xintercept = 4.5)+
  geom_vline(xintercept = 8.5)+
  scale_x_discrete(name="",
                   labels=c("1_fall" = "upstream (fall)", 
                            "1_spring" = "upstream (spring)",
                            "1_summer"="upstream (summer)",
                            "1_winter"="upstream (winter)",
                            "2_fall" = "midstream (fall)", 
                            "2_spring" = "midstream (spring)",
                            "2_summer"="midstream (summer)",
                            "2_winter"="midstream (winter)",
                            "3_fall" = "downstream (fall)", 
                            "3_spring" = "downstream (spring)",
                            "3_summer"="downstream (summer)",
                            "3_winter"="downstream (winter)")
                   )
dev.off()
