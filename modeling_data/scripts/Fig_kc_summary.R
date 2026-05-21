## Fig 6 — Kc summary: LT vs Non-LT vs FAO-56 standard, by state

## ── Path configuration ────────────────────────────────────────────────────────
## Windows (local): edit and uncomment the line below
# project_root <- "C:/Users/YourName/Projects/LT-trait"
if (.Platform$OS.type == "windows") setwd(project_root) else setwd("/home/user/LT-trait")
## ─────────────────────────────────────────────────────────────────────────────
library(dplyr); library(ggplot2); library(readxl)

pdoy_eto <- read.csv("modeling_data/data/nasapower_ETo_pdoy.csv")   # cols: State, Year, Planting_date, ETo_pdoy_mm
ssm      <- read_excel("modeling_data/data/SoybeanSSM_sims2.xlsx")

d <- ssm %>%
  filter(Crop %in% c("SOYBEAN-MG4","SOYBEAN-LT1.5-MG4",
                     "SOYBEAN-LT1.5-MG3","SOYBEAN-MG3"),
         water_mgmt == "Irrigated") %>%
  mutate(
    Genotype = factor(
      ifelse(Crop %in% c("SOYBEAN-MG4","SOYBEAN-MG3"),
             "Non-limited transpiration","Limited transpiration"),
      levels = c("Non-limited transpiration","Limited transpiration")),
    State = case_when(
      grepl("AR", Manag) ~ "AR",
      grepl("MO", Manag) ~ "MO",
      grepl("NE", Manag) ~ "NE"),
    State_label = factor(
      case_when(grepl("AR",Manag) ~ "Arkansas",
                grepl("MO",Manag) ~ "Missouri",
                grepl("NE",Manag) ~ "Nebraska"),
      levels = c("Arkansas","Missouri","Nebraska"))
  ) %>%
  rename(Year = Pyear) %>%
  left_join(pdoy_eto, by = c("State","Year","Planting_date")) %>%
  filter(!is.na(ETo_pdoy_mm), !is.na(SUMETT)) %>%
  mutate(Kc = SUMETT / ETo_pdoy_mm)

## Summary: mean ± SD per genotype × state
kc_summary <- d %>%
  group_by(Genotype, State_label) %>%
  summarise(Kc_mean = mean(Kc, na.rm = TRUE),
            Kc_sd   = sd(Kc,   na.rm = TRUE),
            .groups = "drop")

## Pooled annotation
kc_pooled <- d %>%
  group_by(Genotype) %>%
  summarise(Kc_mean = mean(Kc, na.rm = TRUE),
            Kc_sd   = sd(Kc,   na.rm = TRUE),
            .groups = "drop") %>%
  mutate(label = paste0("Overall Kc = ", round(Kc_mean, 2),
                        " ± ", round(Kc_sd, 2)))

fao_kc <- 1.15

p <- ggplot(kc_summary,
            aes(x = State_label, y = Kc_mean,
                colour = Genotype, fill = Genotype)) +

  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = 1.05, ymax = 1.25, fill = "grey70", alpha = 0.25) +
  geom_hline(yintercept = fao_kc, linetype = "dashed",
             colour = "grey40", linewidth = 0.8) +
  annotate("text", x = 3.45, y = fao_kc + 0.05,
           label = "FAO-56  Kc = 1.15", size = 3,
           colour = "grey35", hjust = 1, fontface = "italic") +

  geom_errorbar(aes(ymin = Kc_mean - Kc_sd, ymax = Kc_mean + Kc_sd),
                width = 0.18, linewidth = 0.7,
                position = position_dodge(width = 0.5)) +
  geom_point(size = 5, shape = 21, stroke = 0.7,
             position = position_dodge(width = 0.5)) +

  geom_text(data = kc_pooled,
            aes(x = 3.52, y = Kc_mean, label = label, colour = Genotype),
            size = 2.9, hjust = 0, fontface = "italic",
            inherit.aes = FALSE, show.legend = FALSE) +

  scale_colour_manual(values = c("Non-limited transpiration" = "#31a354",
                                 "Limited transpiration"      = "#d95f02")) +
  scale_fill_manual(  values = c("Non-limited transpiration" = "#31a354",
                                 "Limited transpiration"      = "#d95f02")) +
  scale_y_continuous(limits = c(0.75, 2.6),
                     breaks = seq(0.8, 2.4, 0.2),
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_x_discrete() +
  coord_cartesian(clip = "off") +
  labs(x        = NULL,
       y        = "Seasonal crop coefficient  (Kc = SUMETT / ETo)",
       colour   = NULL, fill = NULL,
       title    = "Genotype-dependent seasonal crop coefficient (irrigated)",
       subtitle = "Points = mean ± 1 SD  ·  Planting-date-specific ETo windows  ·  Shaded band = FAO-56 range") +
  theme_bw(base_size = 12) +
  theme(legend.position    = "bottom",
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title         = element_text(face = "bold", size = 12),
        plot.subtitle      = element_text(colour = "grey40", size = 9),
        plot.margin        = margin(5, 110, 5, 5))

ggsave("modeling_data/figures/Fig_kc_summary.png",  plot = p, width = 18, height = 14,
       units = "cm", dpi = 150, bg = "white")
ggsave("modeling_data/figures/Fig_kc_summary.png", plot = p, width = 18, height = 14,
message("Saved: Fig_kc_summary")
