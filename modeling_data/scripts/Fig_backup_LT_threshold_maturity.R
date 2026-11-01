## Fig_backup_LT_threshold_maturity.R
## Two backup figures for the ASA 2026 talk, from the 10-site SSM run:
##   A. LT threshold sweep: mean rainfed and irrigated yield for Non-LT, LT2.5,
##      LT2 and LT1.5 (mean of the 10 site means, as on the trade-off slide),
##      labelled with the change vs Non-LT.
##   B. Maturity check: share of seasons reaching physiological maturity
##      (MATYP == 1) and mean days from sowing to R7, by state and water regime
##      (Non-LT), next to the Arkansas-only trade-off (sites where the crop matures).
##
## Input : modeling_data/data/updated-outputs/all_yearly.csv
## Output: modeling_data/figures/Fig_backup_LT_threshold.png, Fig_backup_maturity.png
## Run from the repo root:  Rscript modeling_data/scripts/Fig_backup_LT_threshold_maturity.R
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(readr); library(patchwork) })

yr <- read_csv("modeling_data/data/updated-outputs/all_yearly.csv", show_col_types = FALSE) |>
  mutate(Location = trimws(gsub("­", "", Location)),
         Yield_kgha = as.numeric(Yield_kgha), R7PM = as.numeric(`R7.PM`),
         Genotype = factor(Genotype_class, levels = c("Non-LT", "LT2.5", "LT2", "LT1.5"),
                           labels = c("Non-LT", "LT 2.5 kPa", "LT 2.0 kPa", "LT 1.5 kPa")),
         water_mgmt = factor(water_mgmt, levels = c("Rainfed", "Irrigated")))

WATER_COL <- c(Rainfed = "#31a354", Irrigated = "#3182bd")
theme_bk <- theme_bw(base_size = 15) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(colour = "black"), legend.position = "top",
        legend.title = element_blank(),
        strip.background = element_rect(fill = "white", colour = "black"),
        plot.background = element_rect(fill = "white", colour = NA))

## ── A. Threshold sweep ─────────────────────────────────────────────────────────
site_mean <- function(d) d |>
  group_by(Genotype, water_mgmt, Location) |> summarise(y = mean(Yield_kgha, na.rm = TRUE), .groups = "drop") |>
  group_by(Genotype, water_mgmt) |> summarise(Yield = mean(y), se = sd(y) / sqrt(n()), .groups = "drop") |>
  group_by(water_mgmt) |> mutate(pct = 100 * (Yield / Yield[Genotype == "Non-LT"] - 1)) |> ungroup() |>
  mutate(lab = ifelse(Genotype == "Non-LT", "", sprintf("%+.0f%%", pct)))

sweep <- site_mean(yr)
print(sweep)
pA <- ggplot(sweep, aes(Genotype, Yield, fill = water_mgmt)) +
  geom_col(position = position_dodge(0.8), width = 0.75, colour = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = Yield - se, ymax = Yield + se), position = position_dodge(0.8), width = 0.2) +
  geom_text(aes(y = Yield + se + 250, label = lab), position = position_dodge(0.8), size = 5, fontface = "bold") +
  scale_fill_manual(values = WATER_COL) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Simulated genotype (VPD threshold)", y = expression(Yield~(kg~ha^{-1}))) +
  theme_bk
ggsave("modeling_data/figures/Fig_backup_LT_threshold.png", pA, width = 9, height = 5.2, dpi = 300)

## ── B. Maturity check + Arkansas-only trade-off ─────────────────────────────────
mat <- yr |> filter(Genotype == "Non-LT") |>
  group_by(State, water_mgmt) |>
  summarise(matured = 100 * mean(MATYP == 1), r7 = mean(R7PM, na.rm = TRUE), .groups = "drop")
print(mat)
pB1 <- ggplot(mat, aes(State, matured, fill = water_mgmt)) +
  geom_col(position = position_dodge(0.8), width = 0.75, colour = "black", linewidth = 0.3) +
  geom_text(aes(y = matured + 6, label = sprintf("%.0f d", r7)), position = position_dodge(0.8), size = 4.3) +
  scale_fill_manual(values = WATER_COL) +
  scale_y_continuous(limits = c(0, 110), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  labs(x = NULL, y = "Seasons reaching maturity (%)",
       title = "Non-LT: seasons reaching R7 (labels: mean days to R7)") +
  theme_bk + theme(plot.title = element_text(size = 13))

ar <- site_mean(yr |> filter(State == "Arkansas", Genotype %in% c("Non-LT", "LT 1.5 kPa"))) |>
  mutate(Genotype = droplevels(Genotype))
print(ar)
pB2 <- ggplot(ar, aes(Genotype, Yield, fill = water_mgmt)) +
  geom_col(position = position_dodge(0.8), width = 0.75, colour = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = Yield - se, ymax = Yield + se), position = position_dodge(0.8), width = 0.2) +
  geom_text(aes(y = Yield + se + 300, label = lab), position = position_dodge(0.8), size = 5, fontface = "bold") +
  scale_fill_manual(values = WATER_COL) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = expression(Yield~(kg~ha^{-1})), title = "Arkansas sites only (4 sites)") +
  theme_bk + theme(plot.title = element_text(size = 13))

pB <- pB1 + pB2 + plot_layout(widths = c(1.2, 1), guides = "collect") & theme(legend.position = "top")
ggsave("modeling_data/figures/Fig_backup_maturity.png", pB, width = 11, height = 5.2, dpi = 300)
cat("Saved Fig_backup_LT_threshold.png, Fig_backup_maturity.png\n")
