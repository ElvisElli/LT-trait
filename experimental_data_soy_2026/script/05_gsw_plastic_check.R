## 05_gsw_plastic_check.R
## Quick check: stomatal conductance (gsw) at ~1 kPa (treatment 1) and ~3 kPa
## (treatment 5), grouped as plastic / high WUE vs nonplastic / low WUE
## (PI398198, the MAGIC parent, is grouped with the nonplastic lines as in the
## paper figures) and the wilting checks.
##
## Uses Day 1 and Day 2 measurements because gsw at treatment 1 was only taken
## in replication 5 on Day 1. Input: data/data_clean.csv (from 01_data_inspection_QC.R).
## Output: figures/gsw_plastic_vs_nonplastic_check.png
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(readr); library(patchwork) })

d <- read_csv("data/data_clean.csv", show_col_types = FALSE) |>
  filter(!is.na(gsw), Treatment %in% c(1, 5)) |>
  mutate(
    VPD_level = factor(ifelse(Treatment == 1, "~1 kPa", "~3 kPa"), levels = c("~1 kPa", "~3 kPa")),
    Group = case_when(
      WUE_class == "Plastic, High WUE"                          ~ "Plastic / high WUE",
      WUE_class %in% c("Nonplastic, Low WUE", "MAGIC parent") ~ "Nonplastic / low WUE",
      TRUE                                                     ~ "Wilting checks"),
    Group = factor(Group, levels = c("Plastic / high WUE", "Nonplastic / low WUE", "Wilting checks")))

COLS <- c("Plastic / high WUE" = "#4393c3", "Nonplastic / low WUE" = "#d4824a", "Wilting checks" = "#9d2235")

grp <- d |> group_by(Group, VPD_level) |>
  summarise(n = n(), n_geno = n_distinct(Genotype), mean = mean(gsw), se = sd(gsw) / sqrt(n()), .groups = "drop")
geno <- d |> group_by(Group, Genotype, VPD_level) |>
  summarise(n = n(), mean = mean(gsw), .groups = "drop")

cat("\nGroup means (gsw, mol m-2 s-1):\n"); print(grp)
wide <- grp |> select(Group, VPD_level, mean) |> pivot_wider(names_from = VPD_level, values_from = mean) |>
  mutate(decline_pct = 100 * (`~3 kPa` - `~1 kPa`) / `~1 kPa`)
cat("\nDecline from ~1 to ~3 kPa:\n"); print(wide)
cat("\nGenotype means:\n"); print(geno |> select(-n) |> pivot_wider(names_from = VPD_level, values_from = mean) |> arrange(Group), n = Inf)

thm <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "#e6e6e6"),
        panel.border = element_rect(colour = "black", linewidth = 0.7),
        axis.text = element_text(colour = "black"), legend.position = "top",
        plot.background = element_rect(fill = "white", colour = NA))

pd <- position_dodge(width = 0.5)
pA <- ggplot() +
  geom_point(data = d, aes(VPD_level, gsw, colour = Group),
             position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.5, seed = 1),
             alpha = 0.35, size = 1.8) +
  geom_line(data = grp, aes(VPD_level, mean, colour = Group, group = Group), position = pd, linewidth = 0.9) +
  geom_errorbar(data = grp, aes(VPD_level, ymin = mean - se, ymax = mean + se, colour = Group),
                position = pd, width = 0, linewidth = 0.9) +
  geom_point(data = grp, aes(VPD_level, mean, fill = Group), position = pd, shape = 21,
             colour = "white", size = 4, stroke = 0.8) +
  geom_text(data = grp, aes(VPD_level, mean, label = sprintf("%.2f", mean), group = Group),
            position = pd, vjust = -1.4, size = 3.2, colour = "#333333") +
  scale_colour_manual(values = COLS, name = NULL) + scale_fill_manual(values = COLS, name = NULL) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
  labs(x = "VPD treatment", y = expression(g[sw]~(mol~m^{-2}~s^{-1})),
       title = "A  Plants (Day 1 and Day 2), group mean +/- SE") + thm

pB <- ggplot(geno, aes(VPD_level, mean, group = Genotype, colour = Group)) +
  geom_line(linewidth = 0.7, alpha = 0.9) +
  geom_point(size = 2.2) +
  ggrepel::geom_text_repel(data = filter(geno, VPD_level == "~3 kPa"), aes(label = Genotype),
            direction = "y", hjust = 0, nudge_x = 0.12, size = 2.7, segment.size = 0.25,
            segment.colour = "grey60", box.padding = 0.12, min.segment.length = 0, seed = 1,
            show.legend = FALSE) +
  scale_colour_manual(values = COLS, name = NULL) +
  scale_x_discrete(expand = expansion(mult = c(0.1, 0.6))) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
  labs(x = "VPD treatment", y = NULL, title = "B  Genotype means") + thm + theme(legend.position = "none")

p <- pA + pB + plot_layout(widths = c(1.1, 1))
ggsave("figures/gsw_plastic_vs_nonplastic_check.png", p, width = 11, height = 5.5, dpi = 200, bg = "white")
cat("\nSaved figures/gsw_plastic_vs_nonplastic_check.png\n")
