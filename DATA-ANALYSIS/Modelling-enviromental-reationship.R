library(lme4)
library(MuMIn)
library(blmeco)
library(piecewiseSEM)
library(insight)
library(DHARMa)
library(performance)
library(glmmTMB)
library(interactions)
library(vegan)
library(dplyr)


attach (RA_PROSTREDI_kopie_2)
summary (RA_PROSTREDI_kopie_2)

chrastal_data <- openxlsx::read.xlsx("data/chrast_data.xlsx")


model<-(glmmTMB(RA_Celkem~pocet+(1|Oblast),RA_PROSTREDI_kopie_2,family=nbinom1))
summary(model) #druha hypotéza - jestli je zde vliv prostředí na toho chřástala 
plot(RA_Celkem, pocet, pch = 16, cex = 1,abline (lm(RA_Celkem~pocet)), col = "blue", main = "Početnost chřástala v rámci počtu různých typů stanovišť", xlab = "Počet chřástalů", ylab = "počet stanovišť") #grafickéznázornění
check_overdispersion (model)

RA_PROSTREDI_kopie_2 %>% 
  mutate(shannon = diversity(dplyr::select(., rakos, orobinec, zblochan, ostrice, dreviny, ostatni), index = "shannon"),
         shannon_perc = diversity(dplyr::select(., 'rakos%', 'orobinec%', 'zblochan%', 'ostrice%', 'dreviny%', 'ostatni%'), index = "shannon"))


model<-(glmmTMB(RA_Celkem~ zarust+(1|Oblast),RA_PROSTREDI_kopie_2,family=nbinom1))
summary(model) #vliv zárůstu (teoreticky vydělit 10)
plot(RA_Celkem, zarust, pch = 16, cex = 1,abline (lm(RA_Celkem~pocet)), col = "blue", main = "Početnost chřástala v rámci počtu různých typů stanovišť", xlab = "Počet chřástalů", ylab = "počet stanovišť") #grafickéznázornění
check_overdispersion (model)

model<-(glmmTMB(RA_Celkem~ `rakos%`+`orobinec%`+`zblochan%`+`ostrice%`+`dreviny%`+`ostatni%`+(1|Oblast),RA_PROSTREDI_kopie_2,family=nbinom1))
summary(model) #model na zjištění, zda má nějaký prostředí na toho chřástala vliv a který to je 

model<-(glmmTMB(RA_Celkem~`zblochan%`+`ostrice%`+`dreviny%`+`ostatni%`+(1|Oblast),RA_PROSTREDI_kopie_2,family=nbinom1))
summary(model) #model z kolika procent mi vysvetluje rokas a orobinec, jak jsou důležitý 

r.squaredGLMM(model) #kolik procent variability mi ukazuj ty data

#orobinec a rákos jsou pro něj nejdůležitější, ale můžou být pro něj důležité ostatní prostředí kvůli hnízdění 

detach(RA_PROSTREDI_kopie_2)

chrastal_data <- openxlsx::read.xlsx("S:/Složky uživatelů/Gaigr/chrast_data.xlsx") %>%
  rename(pocet_celkem = `RA-Celkem`) %>%
  mutate(shannon = diversity(dplyr::select(., rakos, orobinec, zblochan, ostrice, dreviny, ostatni), index = "shannon"),
         shannon_perc = diversity(dplyr::select(., 'rakos%', 'orobinec%', 'zblochan%', 'ostrice%', 'dreviny%', 'ostatni%'), index = "shannon"))

model<-(glmmTMB(pocet_celkem~X19+(1|Oblast),chrastal_data,family=nbinom1))
summary(model) #druha hypotéza - jestli je zde vliv prostředí na toho chřástala 
model1<-(glmmTMB(pocet_celkem~shannon+(1|Oblast),chrastal_data,family=nbinom1))
summary(model1) #druha hypotéza - jestli je zde vliv prostředí na toho chřástala 
plot(chrastal_data, X19, pch = 16, cex = 1,abline (lm(pocet_celkem~X19)), col = "blue", main = "Početnost chřástala v rámci počtu různých typů stanovišť", xlab = "Počet chřástalů", ylab = "počet stanovišť") #grafickéznázornění
check_overdispersion(model)

chrast_sum <- chrastal_data %>%
  dplyr::group_by(Oblast) %>%
  dplyr::summarise(pocet_pozorovani = n(),
                   pocet_celkem = mean(pocet_celkem, na.rm = TRUE),
                   pocet_biotopu = mean(X19, na.rm = TRUE),
                   shannon = mean(shannon, na.rm = TRUE)) %>%
  dplyr::ungroup()
ggplot(data = chrast_sum, aes(x = pocet_biotopu, y = pocet_celkem, size = pocet_pozorovani)) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("\nprůměrný počet biotopů na lokalitě") +
  ylab("průměrný počet jedinců na lokalitě\n") +
  theme_classic()
