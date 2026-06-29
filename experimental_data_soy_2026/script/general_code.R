library(tidyverse)
library(ggplot2)
library(dplyr)
library(readxl)


df_all <- read_excel("C:/Users/mmeira/Box/Experiments-2026/Soy-GC/Experiment 2026/data_all.xlsx")
View(data_all_rep1)

df_checks <- df_all %>%
  filter(Genotype %in% c("A5959","PI471938","PI416937")) %>% 
  mutate(day=as.factor(day)) %>% 
  mutate(Genotype=as.factor(Genotype)) %>% 
  mutate(avg_vpd = round(as.numeric(as.character(avg_vpd)), 1)) %>% 
  mutate(avg_vpd=as.factor(avg_vpd)) %>% 
  mutate(RepNo=as.factor(RepNo))

str(df_checks)


df_checks %>% 
  ggplot(aes(x=avg_vpd,y=TR_mg_sqm_seg, color = day)) +
  geom_point(size=2.5)+
  facet_grid(RepNo~Genotype)


ggsave(filename = "SoybeanAllReps.tiff", dpi=600,width = 9,height= 6, path = "../Desktop/")
