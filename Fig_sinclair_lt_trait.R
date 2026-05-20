## Standalone figure: Sinclair et al. (2008) leaf transpiration vs VPD

setwd("/home/user/LT-trait")

library(dplyr)
library(ggplot2)

vpd_seq <- seq(0, 4, by = 0.01)

sinclair <- bind_rows(
  tibble(VPD = vpd_seq, Genotype = "A5959 (non-LT)",
         Tr = -14.2 + 26.0 * VPD),
  tibble(VPD = vpd_seq, Genotype = "Graham (non-LT)",
         Tr =   5.7 + 22.5 * VPD),
  tibble(VPD = vpd_seq, Genotype = "PI 471938 (LT)",
         Tr =  -1.8 + 21.0 * VPD),
  tibble(VPD = vpd_seq, Genotype = "PI 416937 (LT)",
         Tr = ifelse(VPD < 2.13,
                     -21.0 + 30.8 * VPD,
                      56.6 -  5.6 * VPD))
) %>%
  mutate(
    Tr   = pmax(0, Tr),
    Type = ifelse(grepl("non-LT", Genotype),
                  "Non-limited transpiration", "Limited transpiration"),
    Type = factor(Type, levels = c("Non-limited transpiration",
                                   "Limited transpiration")),
    Genotype = factor(Genotype,
                      levels = c("A5959 (non-LT)", "Graham (non-LT)",
                                 "PI 471938 (LT)", "PI 416937 (LT)"))
  )

threshold_pt <- tibble(x = 2.13, y = -21.0 + 30.8 * 2.13)

p <- ggplot(sinclair, aes(x = VPD, y = Tr,
                          colour = Type, linetype = Genotype)) +
  geom_line(linewidth = 1.0) +
  geom_segment(data = threshold_pt,
               aes(x = x, xend = x, y = 0, yend = y),
               inherit.aes = FALSE,
               linetype = "dotted", colour = "#d95f02", linewidth = 0.7) +
  annotate("text", x = 2.18, y = 2,
           label = "VPD threshold (2.13 kPa)",
           size = 3, colour = "#d95f02", hjust = 0, fontface = "italic") +
  scale_colour_manual(values = c("Non-limited transpiration" = "#31a354",
                                 "Limited transpiration"      = "#d95f02")) +
  scale_linetype_manual(values = c("A5959 (non-LT)"  = "solid",
                                   "Graham (non-LT)"  = "longdash",
                                   "PI 471938 (LT)"   = "solid",
                                   "PI 416937 (LT)"   = "longdash")) +
  scale_x_continuous(limits = c(0, 4), breaks = 0:4,
                     expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x     = "VPD (kPa)",
       y     = expression(Transpiration ~ (mg ~ H[2]*O ~ m^{-2} ~ s^{-1})),
       colour   = NULL,
       linetype = "Genotype",
       title    = "Leaf transpiration response to vapour pressure deficit",
       subtitle = "Sinclair et al. (2008) — four soybean genotypes") +
  theme_bw(base_size = 12) +
  theme(legend.position  = "bottom",
        legend.box       = "vertical",
        legend.key.width = unit(1.4, "cm"),
        panel.grid.minor = element_blank(),
        plot.title       = element_text(face = "bold", size = 12),
        plot.subtitle    = element_text(colour = "grey40", size = 9)) +
  guides(colour   = guide_legend(order = 1, override.aes = list(linewidth = 1.2)),
         linetype = guide_legend(order = 2))

ggsave("Fig_sinclair_lt_trait.png",  plot = p,
       width = 16, height = 13, units = "cm", dpi = 150, bg = "white")
ggsave("Fig_sinclair_lt_trait.tiff", plot = p,
       width = 16, height = 13, units = "cm", dpi = 300,
       compression = "lzw", bg = "white")
message("Saved: Fig_sinclair_lt_trait (.png + .tiff)")
