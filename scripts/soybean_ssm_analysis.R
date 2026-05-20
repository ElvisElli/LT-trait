
## Script name: ====

## Script objective: ====

## Author: Elvis F. Elli ====

## Script was created on: 2025-11-03====

## Cleaning up environment ====
rm(list=ls())

## Libraries ====
library(tidyverse)
library(readxl)
library(cowplot)


## Set working directory ====
setwd("/home/user/LT-trait")

source("plot_theme.R")


ssm.model <- read_excel("data/SoybeanSSM_sims2.xlsx")

str(ssm.model)

ssm.model %>% 
  ggplot(aes(x=Crop,y=Yield_kgha)) +
  geom_boxplot()+
  geom_point()+
  facet_wrap(~Location)+
  temp+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##plot 1

plot1 <- ssm.model %>% 
  filter(Crop %in% c("SOYBEAN-MG4","SOYBEAN-LT1.5-MG4","SOYBEAN-LT1.5-MG3","SOYBEAN-MG3")) %>%
  mutate(MTMEAT = (MTMAXT+MTMINT)/2,
         Crop=ifelse(Crop=="SOYBEAN-MG4","Non-limited transpiration",
                     ifelse(Crop=="SOYBEAN-MG3","Non-limited transpiration","Limited transpiration")),
         Crop=as.factor(Crop),
         Crop=factor(Crop,levels = c("Non-limited transpiration","Limited transpiration"))) %>% 
  select(Crop,water_mgmt,Location,Planting_date,Yield_kgha,SRAINT,MTMEAT) %>% 
  group_by(Location,Planting_date) %>% 
  mutate(Rain.avg=mean(SRAINT),
         Rain.sd=sd(SRAINT),
         Rain.class=ifelse(SRAINT>Rain.avg+Rain.sd*0.5,"Wet years",
                           ifelse(SRAINT<Rain.avg-Rain.sd*0.5,"Dry years","Normal years"))) %>% 
  group_by(Crop,water_mgmt,Rain.class) %>%
  summarise(Yield_kgha=mean(Yield_kgha,na.rm=T)) %>% 
  filter(water_mgmt=="Rainfed")

plot1 %>%  
  ggplot(aes(x=Crop,y=Yield_kgha,fill=Crop)) +
  geom_col(position = position_dodge(),colour="black")+
  facet_wrap(~Rain.class)+
  temp+
  scale_y_continuous(limits = c(0,4100))+
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "top",
        plot.margin = margin(t = 10, r = 10, b = 10, l = 40))+
  scale_fill_manual(values = c("#31a354","#f7fcb9"))+
  labs(x=element_blank(),y="Seed yield (kg/ha)")

ggsave("figures/Fig_ssm_yield_rainclass.png",width=20,height=15,units ="cm",dpi=600,compression="lzw",bg="white")

##for USB, 12/12====
plot1 %>%  
  ggplot(aes(x=Crop,y=Yield_kgha,fill=Crop)) +
  geom_col(position = position_dodge(),colour="black")+
  facet_wrap(~Rain.class)+
  temp+
  scale_y_continuous(limits = c(0,4100))+
  theme(axis.text.x = element_blank(),
        legend.position = "top",
        plot.margin = margin(t = 10, r = 10, b = 10, l = 40))+
  scale_fill_manual(values = c("#31a354","#f7fcb9"))+
  labs(x=element_blank(),y="Seed yield (kg/ha)")

ggsave("figures/Fig_ssm_yield_USB.png",width=20,height=10,units ="cm",dpi=600,compression="lzw",bg="white")

##for discussion
plot1 %>% 
  pivot_wider(values_from = Yield_kgha,names_from = Crop) %>% 
  mutate(yield.change=(`Limited transpiration`-`Non-limited transpiration`)/`Non-limited transpiration`*100)

plot1.yv <- ssm.model %>% 
  filter(Crop %in% c("SOYBEAN-MG4","SOYBEAN-LT1.5-MG4","SOYBEAN-LT1.5-MG3","SOYBEAN-MG3")) %>%
  mutate(MTMEAT = (MTMAXT+MTMINT)/2,
         Crop=ifelse(Crop=="SOYBEAN-MG4","Non-limited transpiration",
                     ifelse(Crop=="SOYBEAN-MG3","Non-limited transpiration","Limited transpiration")),
         Crop=as.factor(Crop),
         Crop=factor(Crop,levels = c("Non-limited transpiration","Limited transpiration"))) %>% 
  group_by(Location,Crop,water_mgmt,Pyear) %>%
  summarise(Yield_kgha=mean(Yield_kgha,na.rm=T)) %>% 
  filter(water_mgmt=="Rainfed") %>% 
  pivot_wider(values_from = Yield_kgha,names_from = Crop) %>% 
  mutate(yield.change=(`Limited transpiration`-`Non-limited transpiration`)/`Non-limited transpiration`*100,
         yield.advantage = ifelse(yield.change>0,1,0)) %>% 
  group_by(Location) %>% 
  summarise(yield.advantage=sum(yield.advantage))





##plot 2

plot2 <- ssm.model %>% 
  filter(Crop %in% c("SOYBEAN-MG4","SOYBEAN-LT1.5-MG4")) %>%
  mutate(MTMEAT = (MTMAXT+MTMINT)/2,
         Crop=ifelse(Crop=="SOYBEAN-MG4","Non-limited transpiration",
                     ifelse(Crop=="SOYBEAN-MG3","Non-limited transpiration","Limited transpiration")),
         Crop=as.factor(Crop),
         Crop=factor(Crop,levels = c("Non-limited transpiration","Limited transpiration")),
         State = case_when(
           grepl("AR", Manag, ignore.case = TRUE) ~ "Arkansas",
           grepl("MO", Manag, ignore.case = TRUE) ~ "Missouri",
           grepl("NE", Manag, ignore.case = TRUE) ~ "Nebraska",
           TRUE ~ NA_character_
         )) %>% 
  select(Crop,water_mgmt,Location,State,Planting_date,Yield_kgha,SRAINT,MTMEAT,SUMETT) %>% 
  filter(water_mgmt=="Irrigated")

lm.lt <- lm(SUMETT ~ MTMEAT, data = plot2 %>% filter(Crop %in% c("Limited transpiration")))
slope.lm.lt <- coef(lm.lt)["MTMEAT"]
rel.slope.lm.lt <- slope.lm.lt/mean(plot2$SUMETT[plot2$Crop == "Limited transpiration"])*100

lm.nonlt <- lm(SUMETT ~ MTMEAT, data = plot2 %>% filter(Crop %in% c("Non-limited transpiration")))
slope.lm.nonlt <- coef(lm.nonlt)["MTMEAT"]
rel.slope.lm.nonlt <- slope.lm.lt/mean(plot2$SUMETT[plot2$Crop == "Non-limited transpiration"])*100


plot2 %>%  
  ggplot(aes(x=MTMEAT,y=SUMETT,fill=Crop)) +
  geom_point(colour="black",shape=21,size=4)+
  geom_smooth(method = "lm",se=F,size=1.5,aes(group=Crop),colour="black")+
  temp+
  scale_y_continuous(limits = c(0,2000))+
  scale_x_continuous(limits = c(13,28),breaks = c(12,16,20,24,28))+
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "top",
        plot.margin = margin(t = 10, r = 10, b = 10, l = 40))+
  scale_fill_manual(values = c("#31a354","#f7fcb9"))+
  scale_colour_manual(values = c("#31a354","#f7fcb9"))+
  labs(x='Season temperature (°C)',y="Season nonwater limited evapotranspiration (mm)")

ggsave("figures/Fig_ssm_et_temperature.png",width=20,height=20,units ="cm",dpi=600,compression="lzw",bg="white")

##plot 3

plot3 <- ssm.model %>% 
  filter(Crop %in% c("SOYBEAN-MG4","SOYBEAN-LT1.5-MG4","SOYBEAN-LT1.5-MG3","SOYBEAN-MG3")) %>%
  mutate(MTMEAT = (MTMAXT+MTMINT)/2,
         Crop=ifelse(Crop=="SOYBEAN-MG4","Non-limited transpiration",
                     ifelse(Crop=="SOYBEAN-MG3","Non-limited transpiration","Limited transpiration")),
         Crop=as.factor(Crop),
         Crop=factor(Crop,levels = c("Non-limited transpiration","Limited transpiration"))) %>% 
  select(Crop,water_mgmt,Location,Manag,Planting_date,SUMET2,SUMET3) %>% 
  filter(water_mgmt=="Rainfed") %>% 
  pivot_longer(!c(Crop:Planting_date)) %>% 
  mutate(name=ifelse(name=="SUMET2","Sowing to beginning seed","Beginning seed to harvest"),
         name=as.factor(name),
         name=factor(name,levels = c("Sowing to beginning seed","Beginning seed to harvest")),
         State = case_when(
           grepl("AR", Manag, ignore.case = TRUE) ~ "Arkansas",
           grepl("MO", Manag, ignore.case = TRUE) ~ "Missouri",
           grepl("NE", Manag, ignore.case = TRUE) ~ "Nebraska",
           TRUE ~ NA_character_
         )) %>% 
  group_by(Crop,name) %>% 
  summarise(value=mean(value,na.rm=T))


plot3 %>%  
  ggplot(aes(x=Crop,y=value,fill=Crop)) +
  geom_col(position = position_dodge(),colour="black")+
  facet_grid(~name)+
  temp+
  scale_y_continuous(limits = c(0,400))+
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "top",
        plot.margin = margin(t = 10, r = 10, b = 10, l = 40))+
  scale_fill_manual(values = c("#31a354","#f7fcb9"))+
  labs(x=element_blank(),y="Crop evapotranspiration (mm)")
  

ggsave("figures/Fig_ssm_et_stage.png",width=20,height=15,units ="cm",dpi=600,compression="lzw",bg="white")

##one more plot with year to year variability, and send presentation to co-authors

##for discussion
plot3 %>% 
  pivot_wider(values_from = value,names_from = Crop) %>% 
  mutate(percent.change=(`Limited transpiration`-`Non-limited transpiration`)/`Non-limited transpiration`*100,
         magnitude.change=(`Limited transpiration`-`Non-limited transpiration`))

##for discussion - across states

ssm.model %>% 
  filter(Crop %in% c("SOYBEAN-MG4","SOYBEAN-LT1.5-MG4","SOYBEAN-LT1.5-MG3","SOYBEAN-MG3")) %>%
  mutate(MTMEAT = (MTMAXT+MTMINT)/2,
         Crop=ifelse(Crop=="SOYBEAN-MG4","Non-limited transpiration",
                     ifelse(Crop=="SOYBEAN-MG3","Non-limited transpiration","Limited transpiration")),
         Crop=as.factor(Crop),
         Crop=factor(Crop,levels = c("Non-limited transpiration","Limited transpiration"))) %>% 
  select(Crop,water_mgmt,Location,Manag,Planting_date,SUMET2,SUMET3) %>% 
  filter(water_mgmt=="Rainfed") %>% 
  pivot_longer(!c(Crop:Planting_date)) %>% 
  mutate(name=ifelse(name=="SUMET2","Sowing to beginning seed","Beginning seed to harvest"),
         name=as.factor(name),
         name=factor(name,levels = c("Sowing to beginning seed","Beginning seed to harvest")),
         State = case_when(
           grepl("AR", Manag, ignore.case = TRUE) ~ "Arkansas",
           grepl("MO", Manag, ignore.case = TRUE) ~ "Missouri",
           grepl("NE", Manag, ignore.case = TRUE) ~ "Nebraska",
           TRUE ~ NA_character_
         )) %>% 
  group_by(Crop,name,State) %>% 
  summarise(value=mean(value,na.rm=T)) %>%  
  pivot_wider(values_from = value,names_from = Crop) %>% 
  mutate(percent.change=(`Limited transpiration`-`Non-limited transpiration`)/`Non-limited transpiration`*100,
         magnitude.change=(`Limited transpiration`-`Non-limited transpiration`))
          