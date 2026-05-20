
## Script name: ====

## Script objective: ====

## Author: Elvis F. Elli ====

## Script was created on: 2025-02-07====

## Cleaning up environment ====
rm(list=ls())

## Libraries ====
library(rstudioapi)
library(tidyverse)
library(readxl)

## Set working directory ====
setwd(dirname(getActiveDocumentContext()$path))

sc <- read_excel("soybean_dataset_SC.xlsx",sheet = "dates_ee")

sc %>% pivot_longer(-c(Irrig:Phenologystage_8_30)) %>% 
  ggplot(aes(Genotype,y=value))+
  geom_boxplot(aes(colour=Irrig))+
  facet_wrap(Date~name,scales="free")

ggsave("p1.tiff",width=30,height=15,units ="cm",dpi=600,compression="lzw",bg="white")


sc %>% ggplot(aes(x=StomatalConductance,y=Yield))+
  geom_point(aes(colour=Irrig))+
  geom_smooth(method = "lm",se=F)+
  facet_wrap(~Date,scales="free")
ggsave("p2.tiff",width=30,height=15,units ="cm",dpi=600,compression="lzw",bg="white")


sc %>% ggplot(aes(x=LeafTemperature,y=Yield))+
  geom_point(aes(colour=Irrig))+
  geom_smooth(method = "lm",se=F)+
  geom_smooth(method = "lm",se=F,aes(colour = Irrig))+
  facet_wrap(~Date,scales="free")

ggsave("p3.tiff",width=30,height=15,units ="cm",dpi=600,compression="lzw",bg="white")


sc %>% filter(Irrig=="Irrigated") %>% 
  ggplot(aes(x=reorder(Genotype, - Yield),y=Yield))+
  geom_boxplot(aes(colour=Irrig))+
  facet_wrap(~Date,scales="free")+
  theme(axis.text.x = element_text(angle = 90))
ggsave("p4.tiff",width=30,height=15,units ="cm",dpi=600,compression="lzw",bg="white")


sc %>% filter( Irrig=="Irrigated") %>% 
  ggplot(aes(x=reorder(Genotype, - StomatalConductance),y=StomatalConductance))+
  geom_boxplot(aes(colour=Irrig))+
  facet_wrap(~Date,scales="free")+
  theme(axis.text.x = element_text(angle = 90))
ggsave("p5.tiff",width=30,height=15,units ="cm",dpi=600,compression="lzw",bg="white")


  
