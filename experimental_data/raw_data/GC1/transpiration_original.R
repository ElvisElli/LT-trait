##Cleaning up environment====
rm(list=ls())

##Loading packages====
library(tidyverse)
library(raster)
library(readxl)
library(lubridate)
library(ggrepel)
library(cowplot)
library(rstudioapi)
library(dplyr)
library(tidyr)
library(lubridate)

##set working directory====
setwd(dirname(getActiveDocumentContext()$path))

source("plot_theme.R")

transpiration <- read_excel("Data analysis.xlsx") %>% 
  dplyr::select(-c(T1,RH1,T2,RH2,W7)) %>% 
  pivot_longer(names_to = "scale",-c(Hours,Minutes,Day,time2,time3,Date,Hour)) %>% 
  mutate(scale=as.factor(scale),
         scale=factor(scale, levels=c("W1","W2","W3","W4","W5","W6","W7","W8","W9","W10","W11","W12","W13","W14","W15","W16")))


##last 2 days
transpiration %>%
  #filter(Day>18) %>% 
  ggplot(aes(x=time3,y=value,colour=scale)) +
  geom_line()+
  temp

ggsave("fig1.tiff",width=30,height=15,units ="cm",dpi=600,compression="lzw",bg = "white")


#per day
codes <- read_excel("Data analysis.xlsx",sheet = "code")

# think more to make it relative

transpiration %>%
  left_join(codes) %>%
  filter(time3>8000) %>%
  #filter(value>1000,value<3000) %>%
  group_by(scale) %>% 
  mutate(rel.value=value/max(value)) %>% 
  group_by(Hours,scale,genotype) %>% 
  summarise(value=mean(value,na.rm=T),
            rel.value=mean(rel.value,na.rm=T)) %>%
  #group_by(Hours,genotype) %>% 
  #summarise(value=mean(value,na.rm=T),
  #          rel.value=mean(rel.value,na.rm=T)) %>%
  ggplot(aes(x=Hours,y=rel.value,colour=genotype)) +
  geom_line()+
  temp

ggsave("fig1.tiff",width=30,height=15,units ="cm",dpi=600,compression="lzw",bg = "white")


test <- transpiration %>%
  left_join(codes) %>% 
  mutate(current.pot.cap = value-pot_tare,
         current.pot.cap.rel=round(current.pot.cap/soil_full,2)) %>% 
  filter(current.pot.cap.rel<0.8,
        current.pot.cap.rel>0) %>% 
  group_by(scale,genotype,rep,Hours,Day) %>% 
  summarise(current.pot.cap.rel=mean(current.pot.cap.rel,na.rm=T))
  
test %>% 
  filter(Hours>320,
         Hours<520,
         !scale=="W9") %>% 
  ggplot(aes(x=Hours,y=current.pot.cap.rel,colour=scale)) +
  geom_line()+
  temp+
  labs(x="Time (hours)",y=element_blank())+
  theme(legend.position = "none",
        axis.text.x = element_blank())+
  scale_y_continuous(limits = c(0.2,1))+
  ggtitle("Water losses by transpiration (% of soil water capacity)")

ggsave("fig2.tiff",width=20,height=15,units ="cm",dpi=600,compression="lzw",bg = "white")


##gif
library(magick)

test.final <-
  test %>% 
  filter(Hours>320,
         Hours<520,#520
         !scale=="W9")

unique_hours <- unique(test.final$Hours)
test.list <- lapply(unique_hours, function(hour) {
  test.final %>% filter(Hours <= hour)
})

cm_to_px <- function(cm, dpi = 96) {
  (cm / 2.54) * dpi
}
cm_to_px(15)

# Define the desired dimensions in centimeters and DPI
width_cm <- 20
height_cm <- 15
dpi <- 150

# Convert the dimensions to pixels
width_px <- cm_to_px(width_cm, dpi)
height_px <- cm_to_px(height_cm, dpi)

# Open a graphics device with specified dimensions in pixels
img <- image_graph(width = width_px, height = height_px, res = dpi) 

# Apply the plotting function to each subset and print the plot
out <- lapply(test.list, function(data){
  
  p <-   ggplot(data, aes(x=Hours,y=current.pot.cap.rel,colour=scale)) +
    geom_line()+
    temp+
    labs(x="Time (hours)",y=element_blank())+
    theme(legend.position = "none",
          axis.text.x = element_blank())+
    scale_y_continuous(limits = c(0.2,1))+
    scale_x_continuous(limits = c(320,520))+
    ggtitle("Water losses by transpiration (% of soil water capacity)")
  
  print(p)
})

# Close the graphics device and create an animation
animation <- image_animate(img, fps = 10)

# Save the animated GIF
image_write(animation, "transpiration_animation.gif")

