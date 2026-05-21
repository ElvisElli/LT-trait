## Fig 4 — Water productivity (yield / SUMETT) by genotype × water_mgmt × state

## ── Path configuration ────────────────────────────────────────────────────────
## Windows (local): edit and uncomment the line below
# project_root <- "C:/Users/YourName/Projects/LT-trait"
if (.Platform$OS.type == "windows") setwd(project_root) else setwd("/home/user/LT-trait")
## ─────────────────────────────────────────────────────────────────────────────
library(dplyr); library(ggplot2); library(readxl)

ssm <- read_excel("modeling_data/data/SoybeanSSM_sims2.xlsx")

d <- ssm %>%
  filter(Crop %in% c("SOYBEAN-MG4","SOYBEAN-LT1.5-MG4",
                     "SOYBEAN-LT1.5-MG3","SOYBEAN-MG3")) %>%
  mutate(
    Genotype = factor(
      ifelse(Crop %in% c("SOYBEAN-MG4","SOYBEAN-MG3"),
             "Non-limited transpiration","Limited transpiration"),
      levels = c("Non-limited transpiration","Limited transpiration")),
    WP      = Yield_kgha / SUMETT,
    water_mgmt = factor(water_mgmt,
                        levels = c("Rainfed","Irrigated"),
                        labels = c("Rainfed\n(water-limited)",
                                   "Irrigated\n(potential)"))
  )

## summary for significance brackets (t-test LT vs Non-LT per panel)
sig <- d %>%
  group_by(water_mgmt) %>%
  summarise(
    p     = t.test(WP ~ Genotype)$p.value,
    y_max = max(WP, na.rm = TRUE) * 1.05,
    .groups = "drop"
  ) %>%
  mutate(label = ifelse(p < 0.001, "***",
                 ifelse(p < 0.01,  "**",
                 ifelse(p < 0.05,  "*", "ns"))))

p <- ggplot(d, aes(x = Genotype, y = WP, fill = Genotype)) +
  geom_violin(trim = TRUE, alpha = 0.35, colour = NA) +
  geom_boxplot(width = 0.25, outlier.shape = 21, outlier.size = 1.2,
               outlier.alpha = 0.4, colour = "grey25", linewidth = 0.5) +
  geom_text(data = sig, aes(x = 1.5, y = y_max, label = label),
            inherit.aes = FALSE, size = 5, fontface = "bold") +
  geom_segment(data = sig, aes(x = 1, xend = 2, y = y_max * 0.985, yend = y_max * 0.985),
               inherit.aes = FALSE, linewidth = 0.4, colour = "grey30") +
  facet_wrap(~ water_mgmt, nrow = 1) +
  scale_fill_manual(values = c("Non-limited transpiration" = "#31a354",
                                "Limited transpiration"     = "#d95f02")) +
  scale_x_discrete(labels = c("Non-limited\ntranspiration",
                               "Limited\ntranspiration")) +
  labs(x        = NULL,
       y        = expression(Water~productivity~(kg~mm^{-1})),
       fill     = NULL,
       title    = "Water productivity: yield per unit seasonal crop ET",
       subtitle = "Violin + boxplot across 6 locations × 30 years  ·  *** p < 0.001") +
  theme_bw(base_size = 12) +
  theme(legend.position  = "none",
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text       = element_text(face = "bold", size = 11),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title       = element_text(face = "bold", size = 12),
        plot.subtitle    = element_text(colour = "grey40", size = 9),
        axis.text.x      = element_text(size = 9))

ggsave("modeling_data/figures/Fig_water_productivity.png",  plot = p, width = 18, height = 14,
       units = "cm", dpi = 150, bg = "white")
ggsave("modeling_data/figures/Fig_water_productivity.png", plot = p, width = 18, height = 14,
message("Saved: Fig_water_productivity")
