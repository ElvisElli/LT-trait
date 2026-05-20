library(tidyverse)
library(ggplot2)
library(readxl)

#directory
data_all<-read_xlsx("../soybeanGC/AllDataGC.xlsx", sheet = "Final-mm") %>% 
  select(c(Genotype,Replication,Treatment,VPD,Transpiration,LeafArea_squarecm,StomatalConductance,TempLeaf,Day)) %>% 
  rename("LeafArea"=LeafArea_squarecm,
         "Conductance"=StomatalConductance) 


unique(data_all$VPD)
data_all$VPD <- as.character(data_all$VPD)
data_all$Genotype <- as.factor(data_all$Genotype)
str(data_all)
sum(is.na(dataDay1$Transpiration))

dataDay1<-data_all %>% 
filter(Day=="1") %>% 
  filter(!is.na(VPD) & !is.na(Transpiration))


ggplot(data = dataDay1,aes(x = VPD,y=Transpiration)) +
  geom_point() +
  geom_smooth(method="lm",se=F)+
  scale_y_continuous(breaks = seq(0,50,by=10))+
  facet_wrap(~Genotype)

ggplot(data = data_all,aes(x = VPD,y=Transpiration)) +
  geom_point() +
  geom_smooth(method="lm",se=F)+
  scale_y_continuous(breaks = seq(0,50,by=10))+
  facet_wrap(~Genotype)
