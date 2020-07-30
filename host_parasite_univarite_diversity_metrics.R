rm(list = ls())  # clear Rs brain 
###########################################################################################
library(readr)   # read csv file
library(tidyr)   # tidying data
library(dplyr)   # manipulating df
library(ggplot2) # graphing
library(ggpubr)  # nice package for figures
library(vegan)   # nmds 
library(diverse)
library(hillR)
###########################################################################################
# here I analyze the relationship between host and parasite richness, shannon and simpson diversity for the RR and PR samples
###########################################################################################
river <- read_csv("C:/Users/grunberg/Dropbox/DISSERTATION_CHAPTERS/1_parasite_community_spatial/data/final_river_data.csv")
river <- river[,-c(38,57)] # remove notes that were in excel file... trash 
set.seed(1234567890)
river$plot <- as.factor(river$plot)
river$subplot <- as.factor(river$subplot)
river <- river %>% dplyr::select(-c(sex))
river <- river %>% drop_na()
river$host_species <-gsub("Etheostoma olmstedi_","Etheostoma olmstedi",  river$host_species)

##################################################################################################
# parasite
##################################################################################################
river.diversity <- river %>%
  group_by(river, season, plot) %>%
  summarise_all(funs(if(is.numeric(.)) sum(., na.rm = TRUE) else first(.)))

#calculate richness, hill shannon and hill simpson for parasite 
river.diversity$richness <- hill_taxa(river.diversity[27:64], q = 0, MARGIN = 1, base = exp(1) )
river.diversity$hill_shannon <- hill_taxa(river.diversity[27:64], q = 1, MARGIN = 1, base = exp(1))
river.diversity$hill_simpson <- hill_taxa(river.diversity[27:64], q = 2, MARGIN = 1, base = exp(1))

##################################################################################################
# fish 
##################################################################################################
# abundance based FISH MATRIX 
fish_01 <- river %>% 
  mutate(abun = rep(1)) %>% 
  group_by(river, plot, season, subplot, host_species) %>%
  summarise_all(funs(if(is.numeric(.)) sum(., na.rm = TRUE) else first(.))) %>% 
  spread(host_species, abun)

fish_01[,27:ncol((fish_01))][is.na(fish_01[,27:ncol((fish_01))])] <- 0
fish <- fish_01 %>%
  group_by(river, plot, season) %>%
  summarise_all(funs(if(is.numeric(.)) sum(., na.rm = TRUE) else first(.))) 

fish <- fish[,-c(26:64,89)] # remove parasite data

#calculate fish diversity metrics 
fish$fish_richness <- hill_taxa(fish[26:49], q = 0, MARGIN = 1, base = exp(1) )
fish$fish_hill_shannon <- hill_taxa(fish[26:49], q = 1, MARGIN = 1, base = exp(1))
fish$fish_hill_simpson <- hill_taxa(fish[26:49], q = 2, MARGIN = 1, base = exp(1))

#combine fish and parasite diversity calculations
diversity.data.full <- merge(river.diversity, fish, by =c("river", "plot", "season")) 
diversity.data <- diversity.data.full %>% 
  filter(richness> 1)

plot(fish_hill_shannon~fish_hill_simpson, data= diversity.data)

#################################################################################3
# host and parasite richness 
#################################################################################
rich.lm2 <- lm(richness~fish_richness*river, data= diversity.data.full)
summary(rich.lm2)
rich.lm <- lm(richness~fish_richness, data= diversity.data.full)
summary(rich.lm)
plot(rich.lm)
anova(rich.lm2, rich.lm)
## generate prediction frame
pframe <- with(diversity.data.full,
               expand.grid(fish_richness=seq(min(fish_richness),max(fish_richness),length=50) #,
                           #prog=levels(prog)
               )
)
## add predicted values (on response scale) to prediction frame
pframe$richness <- predict(rich.lm, newdata=pframe,type="response")

richn <-ggplot(diversity.data.full)+ 
  geom_point( aes(x=fish_richness, y = richness, fill=river),size=3, pch=21) +  
  theme_bw(base_size = 15) +
  labs(x="fish richness", y ="parasite richness")+
  scale_fill_manual(values=c('white', 'black'),
                    guide=guide_legend(override.aes=list(shape=21))) +
  guides(fill=guide_legend("River",  keywidth=0.12,
                           keyheight=0.12, override.aes = list(shape=21),
                           default.unit="inch"), color= FALSE) +
  theme(legend.title = element_text(size=7, face="bold"),
        legend.text = element_text(size=6.5))+ 
  # scale_shape_manual(values=(c(21,22,23,24)))+
  # geom_smooth(method = "lm", se = FALSE, linetype = "dashed")+
  geom_line(data=pframe, aes(x=fish_richness, y = richness), color="black", lty=2)+  ## use prediction data here
  xlim(c(0,13))+
  ylim(c(0,26))+
  theme(legend.position=c(0.15,0.85))

#################################################################################3
# host and parasite shannon diversity  
#################################################################################
shann <-ggplot(diversity.data, aes(x=fish_hill_shannon, y = hill_shannon, fill=river))+ 
  geom_point(pch=21, size=3) +  
  theme_bw(base_size = 15) +
  labs(x="fish shannon", y ="parasite shannon")+
  scale_fill_manual(values=c('white', 'black'),
                    guide=guide_legend(override.aes=list(shape=21))) +
  guides(fill=guide_legend("River",  keywidth=0.12,
                           keyheight=0.12, override.aes = list(shape=21),
                           default.unit="inch"), color= FALSE) +
  theme(legend.title = element_text(size=7, face="bold"),
        legend.text = element_text(size=6.5))+ 
  xlim(c(1,6.5))+
  ylim(c(1,6.5))+
  theme(legend.position=c(0.15,0.85))+
  guides(fill=FALSE)

shan.lm <- lm(hill_shannon~fish_hill_shannon, data= diversity.data)
summary(shan.lm)

#################################################################################3
# host and parasite simpson diversity 
#################################################################################
simps <-ggplot(diversity.data, aes(x=fish_hill_simpson, y = hill_simpson, fill=river))+ 
  geom_point(pch=21, size=3) +  
  theme_bw(base_size = 15) +
  labs(x="fish simpson", y ="parasite simpson")+
  scale_fill_manual(values=c('white', 'black'),
                    guide=guide_legend(override.aes=list(shape=21))) +
  guides(fill=guide_legend("River",  keywidth=0.12,
                           keyheight=0.12, override.aes = list(shape=21),
                           default.unit="inch"), color= FALSE) +
  theme(legend.title = element_text(size=7, face="bold"),
        legend.text = element_text(size=6.5))+ 
  #geom_smooth(method = "glm", , se = F, 
  #           method.args = list(family = "poisson"))+
  xlim(c(1,5.3))+
  ylim(c(1,5.3))+
  theme(legend.position=c(0.15,0.85))+
  guides(fill=FALSE)

simp.lm <- lm(hill_simpson~fish_hill_simpson, data= diversity.data)
summary(simp.lm)
plot(simp.lm)

#export graphics 
jpeg(filename="fish_parasite_diversity.jpeg", width=90, height=270,  units="mm", bg="white", res=300)
ggarrange(richn, shann, simps, ncol=1, nrow=3, labels = 'AUTO')
dev.off()