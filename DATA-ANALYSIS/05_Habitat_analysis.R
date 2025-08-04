libs <- c(
  "lme4", "MuMIn", "blmeco", "piecewiseSEM", "insight",
  "DHARMa", "performance", "glmmTMB", "interactions",
  "vegan", "dplyr", "ggplot2", "ggeffects", "gt", "officer",
  "flextable"
)

for (lib in libs) {
  if (!require(lib, character.only = TRUE, quietly = TRUE)) {
    install.packages(lib, dependencies = TRUE)
    library(lib, character.only = TRUE)
  }
}

data_habitat <- openxlsx::read.xlsx(
  "DATA-ANALYSIS/data/data_habitat.xlsx",
  sheet = 3
  )
summary(data_habitat)

model<-(glmmTMB(RA_Celkem~pocet+(1|Oblast),data_habitat,family=nbinom1))
summary(model) #druha hypotéza - jestli je zde vliv prostředí na toho chřástala 
plot(data = data_habitat, 
     RA_Celkem ~ pocet, 
     pch = 16, cex = 1,
     abline (lm(RA_Celkem~pocet)), col = "blue", 
     main = "Početnost chřástala v rámci počtu různých typů stanovišť", 
     xlab = "Počet chřástalů", ylab = "počet stanovišť") #grafickéznázornění
check_overdispersion (model)

model<-(glmmTMB(RA_Celkem~ zarust+(1|Oblast),data_habitat,family=nbinom1))
summary(model) #vliv zárůstu (teoreticky vydělit 10)
library(ggplot2)

ggplot(data_habitat, aes(x = pocet, y = RA_Celkem)) +
  geom_jitter(color = "black", size = 2) +                         # body
  geom_smooth(method = "lm", color = "blue", se = FALSE) +       # regresní přímka
  labs(
    title = "Početnost chřástala v rámci počtu různých typů stanovišť",
    x = "Počet stanovišť",
    y = "Počet chřástalů"
  ) +
  theme_minimal()

check_overdispersion (model)

# Model se všemi proměnnými
full_model <- glmmTMB(RA_Celkem ~ `rakos%` + `orobinec%` + `zblochan%` + `ostrice%` + `dreviny%` + `ostatni%` + (1 | Oblast),
                      data = data_habitat, family = nbinom1)
r.squaredGLMM(full_model) #kolik procent variability mi ukazuj ty data
# Shrnutí fixních efektů
model_df <- tidy(full_model, effects = "fixed") |> 
  dplyr::select(term, estimate, std.error, statistic, p.value) |> 
  dplyr::mutate(
    estimate = round(estimate, 3),
    std.error = round(std.error, 3),
    statistic = round(statistic, 2),
    p.value = round(p.value, 4)
  )

# Model bez rakosu a orobince
reduced_model <- glmmTMB(RA_Celkem ~ `zblochan%` + `ostrice%` + `dreviny%` + `ostatni%` + (1 | Oblast),
                         data = data_habitat, family = nbinom1)
r.squaredGLMM(reduced_model) #kolik procent variability mi ukazuj ty data

#orobinec a rákos jsou pro něj nejdůležitější, ale můžou být pro něj důležité ostatní prostředí kvůli hnízdění 

library(broom.mixed)  # pokud není nainstalováno: install.packages("broom.mixed")

full_model_summary <- broom.mixed::tidy(full_model, effects = "fixed")
reduced_model_summary <- broom.mixed::tidy(reduced_model, effects = "fixed")

# Výpis tabulek
full_model_summary
reduced_model_summary

library(performance)

compare_models <- performance::compare_performance(full_model, reduced_model)
print(compare_models)

## To MS Word -----
# Převod na flextable
ft <- flextable(model_df)
ft <- set_header_labels(ft,
                        term = "Proměnná",
                        estimate = "Odhad",
                        std.error = "Směrodatná chyba",
                        statistic = "Z-statistika",
                        p.value = "P-hodnota")
ft <- autofit(ft)
ft <- add_header_lines(ft, values = "Modelový odhad vlivu typů prostředí na početnost chřástala")

# Vytvoření Word dokumentu
doc <- read_docx()
doc <- body_add_flextable(doc, ft)
print(doc, target = "DATA-ANALYSIS/data/model_vystup.docx")


library(knitr)
kable(full_model_summary, digits = 3, caption = "Koeficienty full modelu")

gt(model_df) |>
  tab_header(
    title = "Modelový odhad vlivu typů prostředí na početnost chřástala"
  ) |>
  cols_label(
    term = "Proměnná",
    estimate = "Odhad",
    std.error = "Směrodatná chyba",
    statistic = "Z-statistika",
    p.value = "P-hodnota"
  ) |>
  fmt_number(columns = where(is.numeric), decimals = 3) |>
  tab_options(
    table.font.names = "Arial",
    heading.title.font.size = 14,
    column_labels.font.weight = "bold"
  )


library(interactions)
library(ggeffects)

# Vizuální odhad efektu jednotlivých prediktorů
plot(ggpredict(full_model, terms = c("rakos%", "orobinec%", "zblochan%")))


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
