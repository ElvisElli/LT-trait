## Fig 3 — ETo vs SUMETT (irrigated only) — clean Kc regression by genotype

## ── Path configuration ────────────────────────────────────────────────────────
## Windows (local): edit and uncomment the line below
# project_root <- "C:/Users/YourName/Projects/LT-trait"
if (.Platform$OS.type == "windows") setwd(project_root) else setwd("/home/user/LT-trait")
## ─────────────────────────────────────────────────────────────────────────────
library(dplyr); library(tidyr); library(ggplot2); library(readxl); library(scales)

pdoy_eto <- read.csv("data/nasapower_ETo_pdoy.csv")
ssm      <- read_excel("data/SoybeanSSM_sims2.xlsx")

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
      grepl("AR",Manag) ~ "AR", grepl("MO",Manag) ~ "MO",
      grepl("NE",Manag) ~ "NE")
  ) %>%
  rename(Year = Pyear) %>%
  left_join(pdoy_eto, by = c("State","Year","Planting_date")) %>%
  filter(!is.na(ETo_pdoy_mm), !is.na(SUMETT)) %>%
  mutate(Kc = SUMETT / ETo_pdoy_mm) %>%
  rename(ETo_season_mm = ETo_pdoy_mm)

## Kc annotations per genotype
kc_labels <- d %>%
  group_by(Genotype) %>%
  summarise(
    Kc_mean = mean(Kc, na.rm = TRUE),
    Kc_sd   = sd(Kc,   na.rm = TRUE),
    x = quantile(ETo_season_mm, 0.08),
    y = max(SUMETT) * c(0.97, 0.87)[cur_group_id()],
    .groups = "drop"
  ) %>%
  mutate(label = paste0("Kc = ", round(Kc_mean, 2), " ± ", round(Kc_sd, 2)))

## 1:1 line range
lim_x <- range(d$ETo_season_mm)
ref_line <- tibble(x = lim_x, y = lim_x)   # Kc = 1 reference

p <- ggplot(d, aes(x = ETo_season_mm, y = SUMETT,
                   colour = Genotype, fill = Genotype)) +
  ## Kc = 1 reference line
  geom_line(data = ref_line, aes(x = x, y = y),
            inherit.aes = FALSE, linetype = "dotted",
            colour = "grey50", linewidth = 0.7) +
  annotate("text", x = max(lim_x) - 10, y = max(lim_x) + 15,
           label = "Kc = 1", size = 3, colour = "grey50", hjust = 1) +
  geom_point(alpha = 0.4, shape = 21, size = 2.0,
             stroke = 0.3, colour = "grey20") +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.2, alpha = 0.15) +
  geom_label(data = kc_labels,
             aes(x = x, y = y, label = label, colour = Genotype),
             fill = "white", size = 3.2, fontface = "bold",
             label.size = 0.3, show.legend = FALSE) +
  facet_wrap(~ State, nrow = 1,
             labeller = labeller(State = c(AR = "Arkansas",
                                           MO = "Missouri",
                                           NE = "Nebraska"))) +
  scale_colour_manual(values = c("Non-limited transpiration" = "#31a354",
                                 "Limited transpiration"      = "#d95f02")) +
  scale_fill_manual(  values = c("Non-limited transpiration" = "#31a354",
                                 "Limited transpiration"      = "#d95f02")) +
  scale_x_continuous(breaks = breaks_pretty(4)) +
  scale_y_continuous(breaks = breaks_pretty(4)) +
  labs(x        = "Seasonal ETo (mm, DOY 120–300)",
       y        = "Seasonal crop ET — SUMETT (mm)",
       colour   = NULL, fill = NULL,
       title    = "Genotype-dependent crop coefficient under irrigated conditions",
       subtitle = "Planting-date-specific ETo (Pdoy → R7 physiological maturity)  ·  Dotted line = Kc 1.0") +
  theme_bw(base_size = 12) +
  theme(legend.position  = "bottom",
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text       = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        plot.title       = element_text(face = "bold", size = 12),
        plot.subtitle    = element_text(colour = "grey40", size = 9))

ggsave("figures/Fig_kc_scatter.png",  plot = p, width = 24, height = 13,
       units = "cm", dpi = 150, bg = "white")
ggsave("figures/Fig_kc_scatter.png", plot = p, width = 24, height = 13,
message("Saved: Fig_kc_scatter")
