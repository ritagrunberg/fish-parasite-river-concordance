rm(list = ls())  # clear Rs brain 
library(readr)   # read csv file
library(tidyr)   # tidying data
library(dplyr)   # manipulating df
library(ggplot2) # graphing
library(ggpubr)  # nice package for figures
library(vegan)   # nmds 

residuals_passaic <- read_csv("C:/Users/grunberg/Documents/Concordance_submitted_MS/Revision/data_files/residuals_host_parasite_pa_passaic.csv")
residuals_raritan <- read_csv("C:/Users/grunberg/Documents/Concordance_submitted_MS/Revision/data_files/residuals_host_parasite_pa_raritan.csv")

residuals_host_passaic <- read_csv("C:/Users/grunberg/Documents/Concordance_submitted_MS/Revision/data_files/residuals_host_phys_pa_passaic.csv")
residuals_host_raritan <- read_csv("C:/Users/grunberg/Documents/Concordance_submitted_MS/Revision/data_files/residuals_host_phys_pa_raritan.csv")

passaic_river_protest <- read_csv("C:/Users/grunberg/Documents/Concordance_submitted_MS/Revision/data_files/passaic_river_protest.csv")
raritan_river_protest <- read_csv("C:/Users/grunberg/Documents/Concordance_submitted_MS/Revision/data_files/raritan_river_protest.csv")
###################################################

residuals_concordance <- rbind(residuals_passaic, residuals_raritan) 
residuals_host_concordance <- rbind(residuals_host_passaic, residuals_host_raritan) 

river_protest <- rbind(passaic_river_protest, raritan_river_protest)
river_protest <- river_protest %>%
  mutate(var1 = recode(var1, fish_abundance ='fish')) %>%
  mutate(var2 = recode(var2, fish_abundance ='fish')) 
  
pro1 <- river_protest %>% 
  group_by(var1, var2, river) %>% 
  summarise_at(c('ms_square'),mean) %>%
  ungroup()
pro1 <- data.frame(pro1)
pro2 <- data.frame(var1 = river_protest$var2, var2 = river_protest$var1, ms_square = river_protest$ms_square, river= river_protest$river)

pro3 <- data.frame(var1 = river_protest$var1, var2 = river_protest$var2, p_value = river_protest$p_value, river= river_protest$river)
pro4 <- data.frame(var1 = river_protest$var2, var2 = river_protest$var1, p_value = river_protest$p_value, river= river_protest$river)

protest_matrix <- (rbind(pro1, pro2))
pval_matrix <- (rbind(pro3, pro4))

#raritan matrix 
protest_matrix_raritan <- protest_matrix %>% 
  filter(river =="Raritan") %>%
  group_by(var1, var2) %>% 
  summarise_at(c('ms_square'),mean) %>% 
  ungroup()%>%
  spread(var2, ms_square)%>%
  arrange(factor(var1, levels = c("parasite",  "parasite_pa", "fish", "fish_biomass" ,"fish_pa", "physical" , "chemical"))) %>%
  select(var1, parasite,  parasite_pa, fish, fish_biomass ,fish_pa, physical , chemical)

pval_matrix_raritan <- pval_matrix %>% 
  filter(river =="Raritan") %>%
  group_by(var1, var2) %>% 
  summarise_at(c('p_value'),mean) %>% 
  ungroup()%>%
  spread(var2, p_value)%>%
  arrange(factor(var1, levels = c("parasite",  "parasite_pa", "fish", "fish_biomass" ,"fish_pa", "physical" , "chemical"))) %>%
  select(var1, parasite,  parasite_pa, fish, fish_biomass ,fish_pa, physical , chemical)

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
protest_matrix_raritan <- get_lower_tri(protest_matrix_raritan)
protest_matrix_raritan$river <- "Raritan"

pval_matrix_raritan <- get_lower_tri(pval_matrix_raritan)
pval_matrix_raritan$river <- "Raritan"

#passaic matrix 
protest_matrix_passaic <- protest_matrix %>%
  filter(river =="Passaic") %>%
  group_by(var1, var2) %>% 
  summarise_at(c('ms_square'),mean) %>% 
  ungroup()%>%
  spread(var2, ms_square)%>%
  arrange(factor(var1, levels = c("parasite",  "parasite_pa", "fish", "fish_biomass" ,"fish_pa", "physical" , "chemical"))) %>%
  select(var1, parasite,  parasite_pa, fish, fish_biomass ,fish_pa, physical , chemical)

pval_matrix_passaic <- pval_matrix %>% 
  filter(river =="Passaic") %>%
  group_by(var1, var2) %>% 
  summarise_at(c('p_value'),mean) %>% 
  ungroup()%>%
  spread(var2, p_value)%>%
  arrange(factor(var1, levels = c("parasite",  "parasite_pa", "fish", "fish_biomass" ,"fish_pa", "physical" , "chemical"))) %>%
  select(var1, parasite,  parasite_pa, fish, fish_biomass ,fish_pa, physical , chemical)


protest_matrix_passaic<- get_lower_tri(protest_matrix_passaic)
protest_matrix_passaic$river <- "Passaic"

pval_matrix_passaic<- get_lower_tri(pval_matrix_passaic)
pval_matrix_passaic$river <- "Passaic"

##############################################################
# create heatmap of procrustes analysis; shows pairwise comparisions of ordinations 
##############################################################
formatted_pro_mat <- rbind(protest_matrix_raritan, protest_matrix_passaic)
formatted_pro_mat <- formatted_pro_mat %>% gather(var2, ms_square, 2:8)

pval_mat <- rbind(pval_matrix_raritan, pval_matrix_passaic)
pval_mat <- pval_mat %>% gather(var2, p_value, 2:8)

#overide fucking ggplot ordering vars by alphabetical order 
formatted_pro_mat$var1 <- factor(formatted_pro_mat$var1, levels=c("parasite",  "parasite_pa", "fish", "fish_biomass" ,"fish_pa", "physical" , "chemical"))
formatted_pro_mat$var2 <- factor(formatted_pro_mat$var2, levels=c("parasite",  "parasite_pa", "fish", "fish_biomass" ,"fish_pa", "physical" , "chemical"))

pval_mat$var1 <- factor(pval_mat$var1, levels=c("parasite",  "parasite_pa", "fish", "fish_biomass" ,"fish_pa", "physical" , "chemical"))
pval_mat$var2 <- factor(pval_mat$var2, levels=c("parasite",  "parasite_pa", "fish", "fish_biomass" ,"fish_pa", "physical" , "chemical"))

formatted_pro_mat$p_value <- pval_mat$p_value
#protest_results$var1 <- factor(protest_results$var1, levels=c("parasite",  "parasite_pa", "fish", "fish_biomass" ,"fish_pa", "physical" , "chemical"))
#protest_results$var2 <- factor(protest_results$var2, levels=c("parasite",  "parasite_pa", "fish", "fish_biomass" ,"fish_pa", "physical" , "chemical"))

jpeg(filename="protest_matrix_results.jpeg",  width=180, height=120, units="mm", bg="white", res=300)
formatted_pro_mat %>% 
  ggplot(aes(x=var1, y = var2, fill=ms_square))+
 # geom_tile(color = "black")+
  geom_raster(aes(x=var1, y = var2, fill=ms_square))+
  scale_fill_gradient(high="#f1eef6", low="#045a8d", na.value = "white",
                      name="Tau", limits = c(0,1)) +
  labs(x="", y="", title="") +
  geom_text(aes(x=var1, y = var2, label=round(ms_square,2)), size=4.5, vjust=1)+ 
  geom_text(data=subset(formatted_pro_mat, p_value < 0.05), 
          aes(x=var1, y = var2, label="*"), size=7, color="white", vjust=0, hjust=-1)+ 
  # annotate ('text', , aes(x=Var1, y=Var2, label="*"), size=1) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.5, 0.8),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))+
  theme(axis.text.x = element_text( angle=90)) +guides(fill=FALSE)+ facet_wrap(~river)+
  theme(strip.background = element_rect(colour="black", fill="white"), strip.text.x = element_text(size=10)) + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))+ #add black square panel around graphic
  scale_x_discrete(labels=c("parasite_pa" = "parasite (PA)", "fish_biomass" = "fish (biomass)", "fish_pa" = "fish (PA)"))+
  scale_y_discrete(labels=c("parasite_pa" = "parasite (PA)", "fish_biomass" = "fish (biomass)", "fish_pa" = "fish (PA)"))
dev.off()

#################################################################################################################
#supplemental analysis; extracted procrustes residuals and asked whether they differ among seasons 
#################################################################################################################


### residuals of parasite and host presence absence ordination comparisons  
parasite_host <- residuals_concordance %>%  
  mutate_at(c('plot'), as.factor) %>%
  ggplot()+ 
  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = 0, ymax = 0.45,alpha = .15, fill="black")+
  annotate("rect", xmin = 3.5, xmax = 4.5, ymin = 0, ymax = 0.45,alpha = .15, fill="black")+
#  geom_bar(aes(x=(order_season), y = res_host_parasite_pa, fill=plot , group=river
 #              ),
  #         colour="black", stat = "identity", width= 0.1,   position = position_dodge(width = 0.5)  
   #        )+
  geom_point(aes(x=(order_season), y = res_host_parasite_pa, fill=plot, shape=river , group=plot
                 ), size=4 , 
             position = position_dodge(width = 0.5)
            # position = position_jitter(w = 0.25, h = 0)
             )+
  scale_fill_manual(values=c('#f7f7f7','#969696','#252525'), # breaks=c(1, 2,3),
                    labels = c("Upstream", "Midstream", "Downstream")) +
  scale_shape_manual(values = c(21,24))+
  theme_pubr() +labs(x="", y="Procrustes vector residuals") +
  guides(fill=guide_legend("Site",  keywidth=0.12,keyheight=0.12, override.aes = list(shape=21),default.unit="inch"))+
  guides(shape=guide_legend("River",  keywidth=0.12,keyheight=0.12, default.unit="inch"))+
  scale_x_discrete(name ="", limits=c("fall", "winter", "spring", "summer"))+
  ylim(c(0,0.45))+
  annotate("text", x = 1, y = 0.3, label = "a", size=5) + 
  annotate("text", x = 2, y = 0.42, label = "b", size=5) + 
  annotate("text", x = 3, y = 0.35, label = "ab", size=5) + 
  annotate("text", x = 4, y = 0.35, label = "ab", size=5) + 
  ggtitle("Parasite (PA) and fish (PA) comparision")

### residuals of host and physical envi ordination comparisons  
host_envi <- residuals_host_concordance %>%  mutate_at(c('plot'), as.factor) %>%
  ggplot()+ 
  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = 0, ymax = 0.45,alpha = .15, fill="black")+
  annotate("rect", xmin = 3.5, xmax = 4.5, ymin = 0, ymax = 0.45,alpha = .15, fill="black")+
  #  geom_bar(aes(x=(order_season), y = res_host_parasite_pa, fill=plot , group=river
  #              ),
  #         colour="black", stat = "identity", width= 0.1,   position = position_dodge(width = 0.5)  
  #        )+
  geom_point(aes(x=(order_season), y = res_host_phys_pa, fill=plot, shape=river, group=plot), size=4 , position = position_dodge(width = 0.5) )+
  scale_fill_manual(values=c('#f7f7f7','#969696','#252525'), # breaks=c(1, 2,3),
                    labels = c("Upstream", "Midstream", "Downstream")) +
  scale_shape_manual(values = c(21,24))+
  theme_pubr() +labs(x="", y="Procrustes vector residuals") +
  guides(fill=guide_legend("Site",  keywidth=0.12,keyheight=0.12, override.aes = list(shape=21),default.unit="inch"))+
  guides(shape=guide_legend("River",  keywidth=0.12,keyheight=0.12, default.unit="inch"))+
  scale_x_discrete(name ="", limits=c("fall", "winter", "spring", "summer"))+
  ylim(c(0,0.45))+
  ggtitle("Fish (PA) and physical habitat comparision")


#### export graphics 
jpeg(filename="residual_concordance_PA.jpeg",  width=270, height=120, units="mm", bg="white", res=300)
ggarrange(parasite_host, host_envi, common.legend = TRUE, nrow = 1, ncol=2, labels = "AUTO")
dev.off()


# parasite and host 
res.mod <- aov(res_host_parasite_pa ~ season + as.factor(plot)
               , data= residuals_concordance)
library(DescTools)
summary(res.mod)
res.mod
TukeyHSD(res.mod)
EtaSq(res.mod, type=1, anova=TRUE)

# host and envi
res.mod.2 <- aov(res_host_phys_pa ~ season + as.factor(plot)
               , data= residuals_host_concordance)
summary(res.mod.2)
EtaSq(res.mod.2, type=1, anova=TRUE)

#################################################################################################################
# pairwise protest graphic ; not used in the MS
#################################################################################################################

river_sub <-river_protest %>% 
  filter(var1=="parasite_pa") %>% 
  drop_na() %>% 
  mutate(significance = if_else(p_value >= 0.05, 'NS', 'Sign')) 
  #mutate_if(is.numeric, round, 3) 

pa_concor <- river_sub %>% ggplot()+
 # annotate("rect", xmin = 0, xmax = 6, ymin = 0.45, ymax = 0.65,alpha = .15)+
  geom_point(aes(x=order, y = ms_square , fill=significance, # fill= group,
                 shape=river, group=river), size=4 , position = position_dodge(width = 0))+ theme_pubr()+
  scale_shape_manual(values = c(21,24))+
  scale_x_discrete(name ="", limits=c("Fish PA" ,"Fish numerical", "Fish biomass", "Chemical", "Physical"))+
  labs(x="", y=expression(paste('',italic(m^2))))+
#  guides(fill=guide_legend("Group",  keywidth=0.12, keyheight=0.12, override.aes = list(shape=21),default.unit="inch"))+
  guides(shape=guide_legend("River",  keywidth=0.15,keyheight=0.15, default.unit="inch"))+
#  geom_text(data=subset(river_sub, p_value < 0.05), 
 #           aes(x=order, y = ms_square, label="*"), size=7, color="#636363", vjust=0.5, hjust=-1.0)+
  scale_fill_manual(values = c('white', 'black'))+
  guides(fill=FALSE, shape=FALSE)+
  ggtitle("Parasite (PA)")+
  theme(legend.box = "horizontal",
        legend.position=c(0.13,0.85)) +
  theme(axis.text.x = element_text( angle=90)) +
 # geom_hline(yintercept = 0.65, lty=2)+
  ylim(c(0.45,0.9))

river_sub2 <-river_protest %>% 
  filter(var1=="parasite") %>% 
  drop_na() %>% 
  mutate(significance = if_else(p_value >= 0.05, 'NS', 'Sign')) 

#mutate_if(is.numeric, round, 3) 
num_concor <- river_sub2 %>% ggplot()+
  # annotate("rect", xmin = 0, xmax = 6, ymin = 0.45, ymax = 0.65,alpha = .15)+
  geom_point(aes(x=order, y = ms_square , fill=significance, # fill= group,
                 shape=river, group=river), size=4 , position = position_dodge(width = 0))+ theme_pubr()+
  scale_shape_manual(values = c(21,24))+
  scale_x_discrete(name ="", limits=c("Fish PA" ,"Fish numerical", "Fish biomass", "Chemical", "Physical"))+
  labs(x="", y=expression(paste('',italic(m^2))))+
  #  guides(fill=guide_legend("Group",  keywidth=0.12, keyheight=0.12, override.aes = list(shape=21),default.unit="inch"))+
  guides(shape=guide_legend("River",  keywidth=0.15,keyheight=0.15, default.unit="inch"))+
  #  geom_text(data=subset(river_sub, p_value < 0.05), 
  #           aes(x=order, y = ms_square, label="*"), size=7, color="#636363", vjust=0.5, hjust=-1.0)+
  scale_fill_manual(values = c('white', 'black'))+
  guides(fill=FALSE)+
  ggtitle("Parasite (numerical density)")+
  theme(legend.box = "horizontal",
        legend.position=c(0.80,0.20)) +
  theme(axis.text.x = element_text( angle=90)) +
  # geom_hline(yintercept = 0.65, lty=2)+
  ylim(c(0.2,1))

jpeg(filename="concordance_ms_graphic.jpeg",  width=90, height=220, units="mm", bg="white", res=300)
ggarrange(num_concor, pa_concor, ncol=1, nrow = 2, labels = "AUTO")
dev.off()
########################################################################################################3
# calcualte fish richness and relationship with residuals from procrustes 
########################################################################################3
river <- read_csv("C:/Users/grunberg/Dropbox/DISSERTATION_CHAPTERS/1_parasite_community_spatial/data/final_river_data.csv")
river <- river[,-c(38,57)] # remove notes that were in excel file... trash 
river$plot <- as.factor(river$plot)
river$subplot <- as.factor(river$subplot)

fish <- river %>% #filter (river == "Raritan") %>% 
  mutate(abun = rep(1)) %>% 
  group_by(river, plot, season, subplot, host_species) %>%
  summarise_all(funs(if(is.numeric(.)) sum(., na.rm = TRUE) else first(.))) %>% spread(host_species, dry_weight_g)
fish[,20:ncol(fish)][is.na(fish[,20:ncol(fish)])] <- 0
fish_sp <- fish %>%
  group_by(river,season, plot, subplot) %>%
  summarise_all(funs(if(is.numeric(.)) sum(., na.rm = TRUE) else first(.))) %>%
  group_by(river ,plot, season) %>%  
  filter(!(river == "Raritan" & season == "summer" & plot ==1)) %>%  # remove summer plot 1 observation becaues no parasite data for that 
  summarise_all(funs(if(is.numeric(.)) mean(., na.rm = TRUE) else first(.))) 
fish_sp <- fish_sp[,-c(26:65, 91)]  # remove parasite data

fish_sp$sp_rich <- specnumber(fish_sp[,26:ncol(fish_sp)]) 

sp_residuals <- merge(fish_sp, residuals_concordance, by =c("season", "river", "plot")) %>%
  select(river, plot, season, sp_rich, res_host_parasite_pa )


# table of species richness 

fish_table <- fish_sp %>% select(river, plot, season, sp_rich) %>%
  spread(plot, sp_rich)

parasite <- river  %>% #filter(river == "Passaic") %>% # filter(season=='fall') %>%
  group_by(river, season, plot, subplot) %>%
  summarise_all(funs(if(is.numeric(.)) sum(., na.rm = TRUE) else first(.))) %>%
  group_by(river, plot, season) %>%
  # summarise_if(.predicate = function(x) is.numeric(x), .funs = funs(mean="mean"))
  summarise_all(funs(if(is.numeric(.)) mean(., na.rm = TRUE) else first(.))) # used to add observation without parasites 

parasite$sp_rich <- specnumber(parasite[,26:65])
parasite <- parasite %>%  select(river, plot, season, sp_rich) %>% drop_na()
parasite_table <- parasite %>%   spread(plot, sp_rich)

sp_residuals <- merge(sp_residuals, parasite, by =c("season", "river", "plot") )

sp_residuals <- sp_residuals %>%  mutate_at(c('plot'), as.factor) #%>% mutate_at(c('sp_rich'), as.numeric)

jpeg(filename="residuals_host_richness.jpeg",  width=120, height=120, units="mm", bg="white", res=300)
sp_residuals %>% ggplot()+
  geom_point(aes(x=sp_rich.x, y = res_host_parasite_pa, fill=plot, shape=river), size=3  )+
  scale_fill_manual(values=c('#f7f7f7','#969696','#252525'), # breaks=c(1, 2,3),
                    labels = c("Up", "Mid", "Down")) +
  scale_shape_manual(values = c(21,24))+
  theme_pubr() +labs(x="Host richness", y="Procrustes vector residuals") +
  guides(fill=guide_legend("Site",  keywidth=0.12,keyheight=0.12, override.aes = list(shape=21),default.unit="inch"))+
  guides(shape=guide_legend("River",  keywidth=0.12,keyheight=0.12, default.unit="inch"))+
  xlim(c(1,13))+ ylim(c(0,0.45))
dev.off()

sp_residuals %>% ggplot()+
  geom_point(aes(x=sp_rich.y, y = res_host_parasite_pa, fill=plot, shape=river), size=3  )+
  scale_fill_manual(values=c('#f7f7f7','#969696','#252525'), # breaks=c(1, 2,3),
                    labels = c("Up", "Mid", "Down")) +
  scale_shape_manual(values = c(21,24))+
  theme_pubr() +labs(x="Parasite richness", y="Procrustes vector residuals") +
  guides(fill=guide_legend("Site",  keywidth=0.12,keyheight=0.12, override.aes = list(shape=21),default.unit="inch"))+
  guides(shape=guide_legend("River",  keywidth=0.12,keyheight=0.12, default.unit="inch"))+
  xlim(c(1,27))+ facet_wrap(~season)+
  ylim(c(0,0.45))

sp.md <- aov(res_host_parasite_pa~sp_rich.x +season, data= sp_residuals)
summary(sp.md)
anova(sp.md)

sp.md.parasite <- aov(res_host_parasite_pa~ plot +season *sp_rich.y, data= sp_residuals)
plot(sp.md.parasite)

summary(sp.md.parasite)
anova(sp.md.parasite)
