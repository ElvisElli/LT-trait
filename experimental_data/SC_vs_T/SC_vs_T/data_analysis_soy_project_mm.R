## Libraries ====
library(rstudioapi)
library(tidyverse)
library(readxl)
library(ggplot2)
library(car)
library(multcomp)
install.packages("agricolae")
library(agricolae)
## Set working directory ====


sc <- read_excel("../QGIS/Elli_image_info/outputs/thermal_summary_attribute_table.xlsx",sheet = "mm")



#)

sc %>% #pivot_longer(-c(Irrig:phenology)) %>% 
  ggplot(aes(genotype,y=temp.mean))+
  geom_boxplot(aes(colour=Irrig))+
  facet_wrap(~date, scales = "free")

ggsave("temp.drone.tiff",width=30,height=15,units ="cm",dpi=600,compression="lzw",bg="white")


sc %>% ggplot(aes(x=conductance,y=yield))+
  geom_point(aes(colour=Irrig))+
  geom_smooth(method = "lm",se=F)+
  facet_wrap(~date,scales="free")
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

#####Only 08/09#####

sc_8_9 <- sc %>% 
  filter(date==as.Date("2024-08-09"))

sc_8_9 %>% #filter( Irrig=="Irrigated") %>% 
  ggplot(aes(x=reorder(genotype, - conductance),y=conductance))+
  geom_point(aes(colour=Irrig))+
  facet_wrap(~date,scales="free")+
  labs(y= expression(conductance~(mol~m^-2~s^-1)),x="genotype")+
  scale_colour_manual(values = c("Irrigated"="blue", "Rainfed"="red"))+
  theme(axis.text.x = element_text(angle = 90))
ggsave("conductanceVSgenotype.tiff",width=30,height=15,units ="cm",dpi=600,compression="lzw",bg="white")


###
sc_8_9 %>% filter( Irrig=="Irrigated") %>% 
  ggplot(aes(x=conductance,y=`yield(bu/ac)`))+
  geom_point(aes(colour=Irrig))+
  geom_text(data=df_label_irri,
            aes(label = genotype), size=5,color="black")+
  geom_smooth(method = "lm",se=F,aes(colour = Irrig))+
  labs(x= expression(conductance~(mol~m^-2~s^-1)),y=expression(yield~(bu~ac^-1)))+
  scale_colour_manual(values = c("Irrigated"="blue", "Rainfed"="red"))+
  theme_minimal()
ggsave("yieldVSconductance.tiff",width=30,height=15,units ="cm",dpi=600,compression="lzw",bg="white")

sc_8_9 %>%
  ggplot(aes(x=reorder(genotype, - diftemp),y=diftemp))+
  geom_boxplot(aes(colour=Irrig))+
  facet_wrap(~Irrig, scales = "free")+
  theme(axis.text.x = element_text(angle = 90))+
   scale_colour_manual(values = c("Irrigated"="blue", "Rainfed"="red"))+
  labs(y="diftemp(Tleaf-Timage)", x="Genotype")
ggsave("diftemp.tiff",width=30,height=15,units ="cm",dpi=600,compression="lzw",bg="white")

####YIELD vs TEMP####
genotype_high <-c("P42A84E","P48A14E","P52A14SE")
df_label<- sc_8_9 %>%  filter(genotype %in% genotype_high) 
df_label_irri<- df_label %>% 
  filter(Irrig=="Irrigated")


###CELICUS
sc_8_9 %>% 
  ggplot(aes(x=temp.mean, y =yield, label=genotype, color=Irrig))+
  geom_point(size= 3)+
  geom_text(data=df_label,
            aes(label = genotype), size=5,color="black")+
  scale_colour_manual(values = c("Irrigated"="blue", "Rainfed"="red"))+
  labs(x="CanopyTemperature(°C)",color="Irrig")+
  theme_minimal()
ggsave("CTvsyield.tiff",width=30,height=15,units ="cm",dpi=600,compression="lzw",bg="white")

###FAHRENHEIT
sc_8_9 %>% 
  ggplot(aes(x=`temp.mean(F)`, y =`yield(bu/ac)`, label=genotype, color=Irrig))+
  geom_point(size= 3)+
  geom_text(data=df_label,
            aes(label = genotype), size=5,color="black")+
  scale_colour_manual(values = c("Irrigated"="blue", "Rainfed"="red"))+
  labs(x="CanopyTemperature(°F)",y=expression(yield~(bu~ac^-1)),color="Irrig")+
  theme_minimal()
ggsave("CTvsyield_F.tiff",width=30,height=15,units ="cm",dpi=600,compression="lzw",bg="white")

###CELICUS
sc_8_9 %>% 
  ggplot(aes(x=leaftemp, y =yield, label=genotype, color=Irrig))+
  geom_point(size= 3)+
  geom_text(data=df_label,
            aes(label = genotype), size=5,color="black")+
  scale_colour_manual(values = c("Irrigated"="blue", "Rainfed"="red"))+
  labs(x="CanopyTemperature(°C)",color="Irrig")+
  theme_minimal()
ggsave("LeafTempvsyield.tiff",width=30,height=15,units ="cm",dpi=600,compression="lzw",bg="white")

###FAHRENHEIT
sc_8_9 %>% 
  ggplot(aes(x=`leaftemp(F)`, y =`yield(bu/ac)`, label=genotype, color=Irrig))+
  geom_point(size= 3)+
  geom_text(data=df_label,
            aes(label = genotype), size=5,color="black")+
  scale_colour_manual(values = c("Irrigated"="blue", "Rainfed"="red"))+
  labs(x="CanopyTemperature(°F)",y=expression(yield~(bu~ac^-1)),color="Irrig")+
  theme_minimal()
ggsave("LeafTempvsyieldF.tiff",width=30,height=15,units ="cm",dpi=600,compression="lzw",bg="white")

##vsT irrigated
sc_8_9_irri<-sc_8_9 %>% 
  filter(Irrig=="Irrigated")
lm_temp<-lm(temp.mean~leaftemp, data = sc_8_9_irri)
r2<-summary(lm_temp)$r.squared


sc_8_9 %>% filter( Irrig=="Irrigated") %>% 
  ggplot(aes(x = leaftemp,y = temp.mean))+
  geom_point(size=2, color="blue")+
  geom_smooth(method = "lm",se=F, color="black")+
  labs(x="LeafTemp porometer(°C)",y="Canopy temperature(°C)")+
  annotate("text", 
           x=max(sc_8_9$leaftemp)*0.97,
           y=max(sc_8_9$temp.mean)*0.74,
           label=paste("r2=", round(r2,3)),hjust=1, size=4)+
  theme_minimal()
ggsave("Temperatures_Irrig.tiff",width=30,height=15,units ="cm",dpi=600,compression="lzw",bg="white")

##vsT rainfed
sc_8_9_rainfed<-sc_8_9 %>% 
  filter(Irrig=="Rainfed")
lm_temp_rain<-lm(temp.mean~leaftemp, data = sc_8_9_rainfed)
r2_rf<-summary(lm_temp_rain)$r.squared


sc_8_9 %>% filter( Irrig=="Rainfed") %>% 
  ggplot(aes(x = leaftemp,y = temp.mean))+
  geom_point(size=2, color="red")+
  geom_smooth(method = "lm",se=F, color="black")+
  labs(x="LeafTemp porometer(°C)",y="Canopy temperature(°C)")+
  annotate("text", 
           x=max(sc_8_9$leaftemp)*1.02,
           y=max(sc_8_9$temp.mean)*0.9,
           label=paste("r2=", round(r2_rf,3)),hjust=1, size=4)+
  scale_y_continuous(limits = c(35,42))+
  scale_x_continuous(limits = c(30,36))+
  theme_minimal()
ggsave("Temperatures_Rainfed.tiff",width=30,height=15,units ="cm",dpi=600,compression="lzw",bg="white")



########
sc_mean<-df_label %>% 
  group_by(genotype,Irrig)

####
sc_mean_sd<-sc_8_9 %>% 
  group_by(genotype,Irrig) %>% 
  summarise(mean_Leaftemp=mean(leaftemp, na.rm=TRUE),
            sd_Leaftemp=sd(leaftemp,na.rm=TRUE),
            mean_yield=mean(yield,na.rm=TRUE))

sc_mean_sd %>% 
ggplot(aes(x=genotype, y =mean_Leaftemp, label=genotype, color=Irrig))+
  geom_point(size= 3)+
  #geom_text(data=sc_mean_sd,
            #aes(label = genotype), size=3,color="black")+
  geom_errorbar(aes(ymin = mean_Leaftemp - sd_Leaftemp, ymax = mean_Leaftemp + sd_Leaftemp), height = 0.2)+
  scale_colour_manual(values = c("Irrigated"="blue", "Rainfed"="red"))+
  labs(x="genotype",color="Irrig")+
  facet_wrap(~Irrig,scales="free")+
  #scale_y_continuous(limits = 25,40)+
  theme_minimal()

####YIELD####
sc_8_9.w <- sc_8_9 %>%
  filter(Irrig == "Irrigated")
sc_8_9.w$meanY <- mean(sc_8_9.w$yield) 
sc_8_9.w$sdY <- sd(sc_8_9.w$yield)

sc_8_9.r <- sc_8_9 %>%
  filter(Irrig == "Rainfed")
sc_8_9.r$meanY <- mean(sc_8_9.r$yield)
sc_8_9.r$sdY <- sd(sc_8_9.r$yield)

sc_total$Irrig<-as.factor(sc_total$Irrig)
str(sc_total)

sc_total<-rbind(sc_8_9.r,sc_8_9.w)


boxplot(yield~Irrig, data = sc_8_9)
t.test(yield~Irrig, data = sc_total)

ggplot(sc_8_9.r)+
  geom_point(aes(x = genotype,y = yield))


###ANOVA FOR RAINFED####
sc_8_9.r$genotype<-as.factor(sc_8_9.r$genotype)
str(sc_8_9.w)
yield_aov <- aov(yield~genotype, data = sc_8_9.r)
hist(yield_aov$residuals)

plot(yield_aov, which = 3)

boxplot(yield~genotype ,sc_8_9.w)

summary(yield_aov)

TukeyHSD(yield_aov)
lsd.t4est<-LSD.test(yield_aov, "genotype")
print(lsd.t4est)
summary(lsd.t4est)


###ANOVA FOR RAINFED####

sc_8_9.w$genotype<-as.factor(sc_8_9.w$genotype)
str(sc_8_9.w)
yield_aov_w <- aov(yield~genotype, data = sc_8_9.w)
hist(yield_aov_w$residuals)

plot(yield_aov_w, which = 3)

boxplot(yield~genotype ,sc_8_9.w)

summary(yield_aov_w)

TukeyHSD(yield_aov_w)
lsd.t4est.w<-LSD.test(yield_aov_w, "genotype")
print(lsd.t4est.w)
summary(lsd.t4est.w)



#####mn
skimr::skim(sc_8_9)


mod1 <- aov(yield ~ genotype * Irrig,data = sc_8_9)
residuals.mod1 <- residuals(mod1)
qqnorm(residuals.mod1)
qqline(residuals.mod1)
lsd.test2 <- LSD.test(mod1, "Irrig","genotype")
print(lsd.test2)
TukeyHSD(mod1)
anova(mod1)
po <- SNK.test(mod1,"genotype*Irrig")
print(po)

sc_8_9 %>%
  ggplot(aes(x= genotype, y= yield, color= Irrig))+
  geom_point()
mod3 <- aov(yield~ genotype, data= sc_8_9)
r.mod3 <- residuals(mod3)
qqnorm(r.mod3)
qqline(r.mod3)

lsd.mod3 <- LSD.test(mod3, "genotype")
print(lsd.mod3)

p <- SNK.test(mod3,"genotype", alpha= 0.05)

print(p)

sc_8_9.r$sdT <- sd(sc_8_9.r$leaftemp)

 #sc_8_9.r %>% 
 # group_by(genotype)%>%
 sc_8_9.r$sdyield <- sc_8_9.r$yield- sc_8_9.r$meanyield  
