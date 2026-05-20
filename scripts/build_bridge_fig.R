## Build two-panel bridge figure:
## Panel A — Sinclair et al. (2008) leaf transpiration vs VPD curves
## Panel B — Our seasonal Kc scatter (SUMETT vs ETo) from SSM simulations

setwd("/home/user/LT-trait")

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(readxl)

## ── Panel A: Sinclair transpiration-VPD curves ────────────────────────────────

vpd_seq <- seq(0, 4, by = 0.01)

sinclair <- bind_rows(
  tibble(VPD = vpd_seq, Genotype = "A5959 (non-LT)",
         Tr = -14.2 + 26.0 * VPD),
  tibble(VPD = vpd_seq, Genotype = "Graham (non-LT)",
         Tr = 5.7 + 22.5 * VPD),
  tibble(VPD = vpd_seq, Genotype = "PI 471938 (LT)",
         Tr = -1.8 + 21.0 * VPD),
  tibble(VPD = vpd_seq, Genotype = "PI 416937 (LT)",
         Tr = ifelse(VPD < 2.13,
                     -21.0 + 30.8 * VPD,
                      56.6 -  5.6 * VPD))
) %>%
  mutate(Tr = pmax(0, Tr),
         Type = ifelse(grepl("non-LT", Genotype),
                       "Non-limited transpiration", "Limited transpiration"),
         Type = factor(Type,
                       levels = c("Non-limited transpiration",
                                  "Limited transpiration")),
         Genotype = factor(Genotype,
                           levels = c("A5959 (non-LT)", "Graham (non-LT)",
                                      "PI 471938 (LT)", "PI 416937 (LT)")))

## threshold annotation for PI 416937
threshold <- tibble(x = 2.13, xend = 2.13, y = 0,
                    yend = -21.0 + 30.8 * 2.13)

pA <- ggplot(sinclair, aes(x = VPD, y = Tr,
                           colour = Type, linetype = Genotype)) +
  geom_line(linewidth = 0.9, na.rm = TRUE) +
  geom_segment(data = threshold,
               aes(x = x, xend = xend, y = y, yend = yend),
               inherit.aes = FALSE,
               linetype = "dotted", colour = "#d95f02", linewidth = 0.6) +
  annotate("text", x = 2.13, y = -4, label = "VPD threshold\n(2.13 kPa)",
           size = 2.8, colour = "#d95f02", hjust = 0.5) +
  scale_colour_manual(values = c("Non-limited transpiration" = "#31a354",
                                 "Limited transpiration"      = "#d95f02")) +
  scale_linetype_manual(values = c("A5959 (non-LT)"   = "solid",
                                   "Graham (non-LT)"   = "dashed",
                                   "PI 471938 (LT)"    = "solid",
                                   "PI 416937 (LT)"    = "dashed")) +
  scale_x_continuous(limits = c(0, 4), breaks = 0:4) +
  scale_y_continuous(limits = c(-5, 100), breaks = seq(0, 100, 25)) +
  labs(x = "VPD (kPa)",
       y = expression(Transpiration~(mg~H[2]*O~m^{-2}~s^{-1})),
       colour = NULL, linetype = "Genotype",
       tag = "(a)",
       title = "Leaf-level: transpiration response to VPD",
       subtitle = "Sinclair et al. (2008)") +
  theme_bw(base_size = 11) +
  theme(legend.position   = "bottom",
        legend.box        = "vertical",
        legend.key.width  = unit(1.2, "cm"),
        panel.grid.minor  = element_blank(),
        plot.tag          = element_text(face = "bold"),
        plot.title        = element_text(face = "bold", size = 10),
        plot.subtitle     = element_text(colour = "grey40", size = 8)) +
  guides(colour   = guide_legend(order = 1, nrow = 2),
         linetype = guide_legend(order = 2, nrow = 2))

## ── Panel B: Seasonal Kc scatter from SSM simulations ────────────────────────

eto_seasonal <- read.csv("data/nasapower_ETo_seasonal.csv")
ssm.model    <- read_excel("data/SoybeanSSM_sims2.xlsx")

ssm_proc <- ssm.model %>%
  filter(Crop %in% c("SOYBEAN-MG4", "SOYBEAN-LT1.5-MG4",
                     "SOYBEAN-LT1.5-MG3", "SOYBEAN-MG3")) %>%
  mutate(
    Genotype = ifelse(Crop %in% c("SOYBEAN-MG4", "SOYBEAN-MG3"),
                      "Non-limited transpiration", "Limited transpiration"),
    Genotype = factor(Genotype,
                      levels = c("Non-limited transpiration",
                                 "Limited transpiration")),
    State = case_when(
      grepl("AR", Manag, ignore.case = TRUE) ~ "AR",
      grepl("MO", Manag, ignore.case = TRUE) ~ "MO",
      grepl("NE", Manag, ignore.case = TRUE) ~ "NE",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(State)) %>%
  rename(Year = Pyear)

merged <- ssm_proc %>%
  left_join(eto_seasonal, by = c("State", "Year")) %>%
  filter(!is.na(ETo_season_mm), !is.na(SUMETT))

kc_slopes <- merged %>%
  group_by(Genotype) %>%
  summarise(
    slope    = coef(lm(SUMETT ~ ETo_season_mm))[["ETo_season_mm"]],
    intercept = coef(lm(SUMETT ~ ETo_season_mm))[["(Intercept)"]],
    x_pos    = quantile(ETo_season_mm, 0.15, na.rm = TRUE),
    y_pos    = max(SUMETT, na.rm = TRUE) * c(0.92, 0.82)[cur_group_id()],
    .groups  = "drop"
  ) %>%
  mutate(label = paste0("Kc = ", round(slope, 2)))

pB <- ggplot(merged, aes(x = ETo_season_mm, y = SUMETT,
                         colour = Genotype, fill = Genotype)) +
  geom_point(alpha = 0.35, shape = 21, size = 1.8,
             stroke = 0.2, colour = "grey30") +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.1, alpha = 0.15) +
  geom_text(data = kc_slopes,
            aes(x = x_pos, y = y_pos, label = label, colour = Genotype),
            size = 3.2, fontface = "bold.italic", show.legend = FALSE,
            hjust = 0) +
  scale_colour_manual(values = c("Non-limited transpiration" = "#31a354",
                                 "Limited transpiration"      = "#d95f02")) +
  scale_fill_manual(values   = c("Non-limited transpiration" = "#31a354",
                                 "Limited transpiration"      = "#d95f02")) +
  scale_x_continuous(breaks = scales::breaks_pretty(4)) +
  scale_y_continuous(breaks = scales::breaks_pretty(4)) +
  labs(x = "Seasonal ETo (mm)",
       y = "Seasonal crop ET — SUMETT (mm)",
       colour = NULL, fill = NULL,
       tag = "(b)",
       title = "Canopy-season: crop ET response to ETo",
       subtitle = "DSSAT SSM simulations — AR, MO, NE (1990–2023)") +
  theme_bw(base_size = 11) +
  theme(legend.position  = "bottom",
        legend.box       = "vertical",
        panel.grid.minor = element_blank(),
        plot.tag         = element_text(face = "bold"),
        plot.title       = element_text(face = "bold", size = 10),
        plot.subtitle    = element_text(colour = "grey40", size = 8)) +
  guides(colour = guide_legend(nrow = 2),
         fill   = guide_legend(nrow = 2))

## ── Combine panels ────────────────────────────────────────────────────────────

p_bridge <- plot_grid(pA, pB, ncol = 2, align = "h", axis = "tb",
                      rel_widths = c(1, 1))

ggsave("figures/Fig_bridge_sinclair_vs_kc.png",  plot = p_bridge,
       width = 24, height = 14, units = "cm", dpi = 150, bg = "white")
ggsave("figures/Fig_bridge_sinclair_vs_kc.png", plot = p_bridge,
       width = 24, height = 14, units = "cm", dpi = 300,

message("Saved: figures/Fig_bridge_sinclair_vs_kc.png")
