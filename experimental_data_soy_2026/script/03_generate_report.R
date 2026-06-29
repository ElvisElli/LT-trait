# =============================================================================
# 03_generate_report.R
# Growth Chamber 2026 - Soybean Limited-Transpiration Experiment
#
# Generates a self-contained HTML report (GC2026_analysis_report.html)
# embedding all figures and tables. Figures must already exist in figures/.
#
# Run from project root  (experimental_data_soy_2026/):
#   source("script/03_generate_report.R")
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(base64enc)
})

# --------------------------------------------------------------------------- #
# HELPERS
# --------------------------------------------------------------------------- #
UA_RED  <- "#9D2235"
UA_WINE <- "#6E1424"

# Embed a PNG as a base-64 <img> tag
img_tag <- function(fname, alt = "", width = "100%") {
  path <- file.path("figures", fname)
  if (!file.exists(path))
    return(sprintf('<p style="color:#c00">Missing: %s</p>', fname))
  b64 <- base64enc::base64encode(path)
  sprintf('<img src="data:image/png;base64,%s" alt="%s" style="width:%s;max-width:100%%;">',
          b64, alt, width)
}

# Plain HTML table (no kableExtra needed)
ua_table <- function(df, caption = NULL) {
  cap <- if (!is.null(caption))
    sprintf('<caption style="font-weight:bold;text-align:left;padding:4px 0 6px">%s</caption>', caption)
  else ""
  th <- paste(sprintf("<th>%s</th>", names(df)), collapse = "")
  rows <- apply(df, 1, function(r) {
    cells <- paste(sprintf("<td>%s</td>", r), collapse = "")
    paste0("<tr>", cells, "</tr>")
  })
  paste0('<div style="overflow-x:auto"><table>', cap,
         "<thead><tr>", th, "</tr></thead>",
         "<tbody>", paste(rows, collapse = ""), "</tbody></table></div>")
}

# Colour-code the Best_model column
colour_model <- function(x) {
  dplyr::case_when(
    x == "Bilinear"     ~ sprintf('<span style="color:%s;font-weight:bold">Bilinear</span>',     UA_RED),
    x == "Lin. plateau" ~ sprintf('<span style="color:%s;font-weight:bold">Lin. plateau</span>', UA_WINE),
    x == "Linear"       ~ '<span style="color:#555">Linear</span>',
    x == "FAILED"       ~ '<span style="color:#aaa;font-style:italic">FAILED</span>',
    TRUE                ~ x
  )
}

# --------------------------------------------------------------------------- #
# LOAD RESULTS
# --------------------------------------------------------------------------- #
safe_csv <- function(f) {
  p <- file.path("data", f)
  if (file.exists(p)) read_csv(p, show_col_types = FALSE) else tibble()
}

tr_best   <- safe_csv("TR_best_model_summary.csv")
gsw_best  <- safe_csv("gsw_best_model_summary.csv")
tr_params <- safe_csv("TR_model_parameters.csv")
gsw_params<- safe_csv("gsw_model_parameters.csv")

df_clean  <- safe_csv("data_clean.csv")
day1      <- filter(df_clean, Day_label == "Day 1")

n_geno    <- n_distinct(day1$Genotype)
n_rep     <- n_distinct(day1$Rep)
vpd_range <- sprintf("%.2f-%.2f", min(day1$VPD, na.rm=TRUE), max(day1$VPD, na.rm=TRUE))
tr_range  <- sprintf("%.1f-%.1f", min(day1$TR, na.rm=TRUE), max(day1$TR, na.rm=TRUE))

# --------------------------------------------------------------------------- #
# BUILD TABLE HTML FRAGMENTS
# --------------------------------------------------------------------------- #
# VPD design table
vpd_tbl <- tibble(
  Treatment          = 1:6,
  `Temp (degC)`        = 32,
  `Target RH (%)`    = c(79,68,57,47,37,24),
  `Target VPD (kPa)` = c(1.00,1.52,2.04,2.52,3.00,3.61),
  `Mean VPD (kPa)`   = c(1.003,1.517,2.056,2.502,3.024,3.666)
)

# TR summary table with coloured model column
tr_disp <- tr_best |>
  mutate(
    Best_model = colour_model(Best_model),
    Slope1     = ifelse(is.na(Slope1), "-", sprintf("%.3f", Slope1)),
    BP_kPa     = ifelse(is.na(BP_kPa), "-", sprintf("%.3f", BP_kPa)),
    Slope2     = ifelse(is.na(Slope2), "-", sprintf("%.3f", Slope2))
  ) |>
  rename(
    `WUE class`         = WUE_class,
    `Expected phenotype`= Expected_behavior,
    `Best model`        = Best_model,
    `Slope1`            = Slope1,
    `BP (kPa)`          = BP_kPa,
    `Slope2`            = Slope2
  ) |>
  select(Genotype, `WUE class`, `Expected phenotype`,
         `Best model`, `Slope1`, `BP (kPa)`, `Slope2`)

gsw_disp <- gsw_best |>
  mutate(
    Best_model = colour_model(Best_model),
    Slope1     = ifelse(is.na(Slope1), "-", sprintf("%.4f", Slope1)),
    BP_kPa     = ifelse(is.na(BP_kPa), "-", sprintf("%.3f", BP_kPa)),
    Slope2     = ifelse(is.na(Slope2), "-", sprintf("%.4f", Slope2))
  ) |>
  rename(
    `WUE class`         = WUE_class,
    `Expected phenotype`= Expected_behavior,
    `Best model`        = Best_model,
    `Slope1`            = Slope1,
    `BP (kPa)`          = BP_kPa,
    `Slope2`            = Slope2
  ) |>
  select(Genotype, `WUE class`, `Expected phenotype`,
         `Best model`, `Slope1`, `BP (kPa)`, `Slope2`)

# AIC tables
tr_aic_disp <- tr_params |>
  select(Genotype, Model, AIC, Delta_AIC) |>
  rename(DeltaAIC = Delta_AIC) |>
  mutate(across(c(AIC, DeltaAIC), ~round(., 1)))

gsw_aic_disp <- gsw_params |>
  select(Genotype, Model, AIC, Delta_AIC) |>
  rename(DeltaAIC = Delta_AIC) |>
  mutate(across(c(AIC, DeltaAIC), ~round(., 1)))

# --------------------------------------------------------------------------- #
# CSS
# --------------------------------------------------------------------------- #
CSS <- sprintf('
* { box-sizing: border-box; margin: 0; padding: 0; }
body { font-family: Georgia, serif; font-size: 14px; color: #1A1A1A;
       background: #FAFAFA; line-height: 1.65; }
a { color: %s; }
header { background: %s; color: white; padding: 28px 40px 20px; }
header h1 { font-size: 1.8em; font-weight: bold; letter-spacing: 0.02em; }
header .subtitle { font-size: 0.95em; opacity: 0.88; margin-top: 4px; }
nav { background: %s; padding: 8px 40px; position: sticky; top: 0; z-index: 100; }
nav a { color: white; text-decoration: none; margin-right: 20px;
        font-size: 0.88em; font-family: sans-serif; opacity: 0.9; }
nav a:hover { opacity: 1; text-decoration: underline; }
main { max-width: 1100px; margin: 0 auto; padding: 30px 40px 60px; }
section { margin-bottom: 48px; }
h2 { font-size: 1.25em; color: %s; border-bottom: 2px solid %s;
     padding-bottom: 5px; margin-bottom: 18px; font-family: sans-serif; }
h3 { font-size: 1.05em; color: %s; margin: 20px 0 8px; font-family: sans-serif; }
p  { margin-bottom: 10px; }
ul { margin: 8px 0 14px 24px; }
.stats-grid { display: grid; grid-template-columns: repeat(auto-fill,minmax(180px,1fr));
              gap: 14px; margin: 18px 0; }
.stat-card { background: #FAF0F2; border-left: 4px solid %s; border-radius: 4px;
             padding: 12px 16px; }
.stat-card .value { font-size: 1.5em; font-weight: bold; color: %s;
                    font-family: sans-serif; }
.stat-card .label { font-size: 0.82em; color: #555; font-family: sans-serif; margin-top: 2px; }
.fig-wrap { background: white; border: 1px solid #ddd; border-radius: 6px;
            padding: 14px; margin: 14px 0;
            box-shadow: 0 1px 4px rgba(0,0,0,0.06); }
.fig-caption { font-size: 0.85em; color: #555; margin-top: 8px;
               font-style: italic; font-family: sans-serif; }
.qc-warn { background: #FFF3CD; border-left: 4px solid #E6A817; border-radius: 4px;
           padding: 12px 16px; margin: 12px 0; font-size: 0.9em; font-family: sans-serif; }
.qc-ok   { background: #D4EDDA; border-left: 4px solid #28A745; border-radius: 4px;
           padding: 12px 16px; margin: 12px 0; font-size: 0.9em; font-family: sans-serif; }
footer { background: %s; color: rgba(255,255,255,0.75); text-align: center;
         padding: 18px; font-size: 0.82em; font-family: sans-serif; margin-top: 40px; }
',
UA_RED, UA_RED, UA_WINE,
UA_RED, UA_RED, UA_WINE,
UA_RED, UA_RED, UA_WINE)

# --------------------------------------------------------------------------- #
# ASSEMBLE HTML
# --------------------------------------------------------------------------- #
today <- format(Sys.Date(), "%Y-%m-%d")

html_body <- paste0('
<header>
  <h1>Growth Chamber 2026 - Soybean Limited-Transpiration Experiment</h1>
  <div class="subtitle">
    Transpiration and stomatal conductance response to VPD &nbsp;|&nbsp;
    Elli Lab, University of Arkansas &nbsp;|&nbsp;
    Report generated: ', today, '
  </div>
</header>

<nav>
  <a href="#overview">Overview</a>
  <a href="#design">Design</a>
  <a href="#qc">QC</a>
  <a href="#tr-models">TR models</a>
  <a href="#gsw-models">gsw models</a>
  <a href="#tables">Tables</a>
  <a href="#methods">Methods</a>
</nav>

<main>

<!-- OVERVIEW -->
<section id="overview">
<h2>Overview</h2>
<div class="stats-grid">
  <div class="stat-card"><div class="value">', n_geno, '</div><div class="label">Genotypes</div></div>
  <div class="stat-card"><div class="value">', n_rep, '</div><div class="label">Replicates</div></div>
  <div class="stat-card"><div class="value">6</div><div class="label">VPD treatments</div></div>
  <div class="stat-card"><div class="value">', vpd_range, '</div><div class="label">VPD range (kPa)</div></div>
  <div class="stat-card"><div class="value">', tr_range, '</div><div class="label">TR range (mg m&#x207B;&#xB2; s&#x207B;&#xB9;)</div></div>
  <div class="stat-card"><div class="value">', nrow(day1), '</div><div class="label">Day-1 observations</div></div>
</div>
<p>Transpiration rate (TR) and stomatal conductance (g<sub>sw</sub>) responses
to stepped VPD levels in 14 soybean genotypes across five WUE phenotype classes.
Three model types fitted to Day-1 data and compared by AIC.</p>
</section>

<!-- DESIGN -->
<section id="design">
<h2>Experimental design</h2>
<h3>VPD treatments</h3>',
ua_table(vpd_tbl, "Table 1. VPD treatments - targets and achieved values."),
'<h3>Genotypes</h3>',
ua_table(
  day1 |> distinct(Genotype, WUE_class) |> arrange(WUE_class, Genotype) |>
    rename(`WUE class` = WUE_class),
  "Table 2. Genotypes and WUE phenotype class."
),

'</section>

<!-- QC -->
<section id="qc">
<h2>Quality control</h2>
<div class="qc-ok">&#x2713; No TR outliers detected (3&#xD7;IQR rule within Genotype &#xD7; Day &#xD7; Treatment).</div>
<div class="qc-warn">&#x26A0; LICOR g<sub>sw</sub> data present for only 41 of 253 Day-1 rows (16%) and 9 of 14 genotypes.</div>

<h3>QC1 - TR distribution</h3>
<div class="fig-wrap">',
img_tag("QC1_TR_distribution.png"),
'<div class="fig-caption">Fig. QC1. TR boxplots by genotype and day. No extreme outliers.</div></div>

<h3>QC2 - TR vs VPD (both days)</h3>
<div class="fig-wrap">',
img_tag("QC2_TR_vs_VPD_by_genotype.png"),
'<div class="fig-caption">Fig. QC2. Raw TR vs VPD per genotype, OLS trend per day.</div></div>

<h3>QC4 - g<sub>sw</sub> distribution</h3>
<div class="fig-wrap">',
img_tag("QC4_gsw_distribution.png"),
'<div class="fig-caption">Fig. QC4. g<sub>sw</sub> distribution (LICOR data only).</div></div>

<h3>QC5 - g<sub>sw</sub> vs VPD, Day 1</h3>
<div class="fig-wrap">',
img_tag("QC5_gsw_vs_VPD_Day1.png"),
'<div class="fig-caption">Fig. QC5. Stomatal conductance vs VPD for genotypes with Day-1 LICOR data.</div></div>

<h3>QC6 - Leaf area by replicate</h3>
<div class="fig-wrap">',
img_tag("QC6_leaf_area_by_rep_genotype.png"),
'<div class="fig-caption">Fig. QC6. Leaf area (cm&#xB2;) per replicate and genotype.</div></div>

<h3>QC7 - Chamber conditions</h3>
<div class="fig-wrap">',
img_tag("QC7_chamber_conditions.png"),
'<div class="fig-caption">Fig. QC7. Temperature, RH, and VPD achieved per treatment.</div></div>

<h3>QC8 - Day 1 vs Day 2 TR</h3>
<div class="fig-wrap">',
img_tag("QC8_TR_day1_vs_day2.png"),
'<div class="fig-caption">Fig. QC8. Day 1 vs Day 2 TR. Points near the diagonal indicate day-to-day consistency.</div></div>
</section>

<!-- TR MODELS -->
<section id="tr-models">
<h2>Transpiration rate (TR) models</h2>
<p>Models fitted per genotype (Day-1 data, Rep as random effect via nlme):</p>
<ul>
  <li><strong>Linear</strong>: TR = a + b&#xB7;VPD</li>
  <li><strong>Bilinear</strong>: two slopes joined at a breakpoint BP</li>
  <li><strong>Linear plateau</strong>: linear increase then flat at TR&#x2098;&#x2090;&#x2093; (classic LT response)</li>
</ul>

<h3>Fig 1 - All models per genotype</h3>
<div class="fig-wrap">',
img_tag("Fig1_TR_models_per_genotype.png"),
'<div class="fig-caption"><strong>Fig. 1.</strong> TR vs VPD per genotype.
Points coloured by WUE class. All three model curves shown; bold = best model by AIC.
Inverted triangles mark estimated breakpoints.</div></div>

<h3>Fig 2 - Pooled by WUE class</h3>
<div class="fig-wrap">',
img_tag("Fig2_TR_by_WUEclass.png"),
'<div class="fig-caption"><strong>Fig. 2.</strong> TR vs VPD by WUE class.
Thin = per-genotype best model; thick = pooled OLS &#xB1;95% CI.</div></div>

<h3>Fig 3 - Breakpoint VPD</h3>
<div class="fig-wrap">',
img_tag("Fig3_TR_breakpoints.png"),
'<div class="fig-caption"><strong>Fig. 3.</strong> Estimated breakpoint VPD for genotypes with non-linear
best model. Dashed line at 1.73 kPa (literature LT breakpoint).</div></div>

<h3>Fig 4 - Transpiration sensitivity (Slope&#x2081;)</h3>
<div class="fig-wrap">',
img_tag("Fig4_TR_slopes.png"),
'<div class="fig-caption"><strong>Fig. 4.</strong> Initial slope of TR vs VPD per genotype.
Higher = greater sensitivity to VPD in the linear range.</div></div>

<h3>Fig 6 - AIC model comparison</h3>
<div class="fig-wrap">',
img_tag("Fig6_TR_AIC_heatmap.png"),
'<div class="fig-caption"><strong>Fig. 6.</strong> &#x394;AIC heatmap. Red = best model (&#x394;AIC = 0).</div></div>

<h3>Fig 7 - TR mean &#xB1; SD per genotype</h3>
<div class="fig-wrap">',
img_tag("Fig7_TR_mean_SD.png"),
'<div class="fig-caption"><strong>Fig. 7.</strong> Transpiration rate: genotype mean +/- SD across replicates per VPD treatment.
Each symbol = one genotype; error bars = 1 SD across replicates.</div></div>
</section>

<!-- GSW MODELS -->
<section id="gsw-models">
<h2>Stomatal conductance (g<sub>sw</sub>) models</h2>
<p>g<sub>sw</sub> (mol m&#x207B;&#xB2; s&#x207B;&#xB9;) measured with LICOR for 9 of 14 genotypes
on Day 1 (41 observations). Mean of three LICOR readings used.</p>

<h3>Fig 5 - g<sub>sw</sub> per genotype</h3>
<div class="fig-wrap">',
img_tag("Fig5_gsw_models_per_genotype.png"),
'<div class="fig-caption"><strong>Fig. 5.</strong> Stomatal conductance vs VPD per genotype.
Points show the three individual LICOR readings per plant; model curve fitted to the per-plant mean.
Negative slopes indicate stomatal closure in response to increasing VPD.</div></div>

<h3>Fig 8 - g<sub>sw</sub> mean &#xB1; SD per genotype</h3>
<div class="fig-wrap">',
img_tag("Fig8_gsw_mean_SD.png"),
'<div class="fig-caption"><strong>Fig. 8.</strong> Stomatal conductance: genotype mean +/- SD across the 3 LICOR readings per VPD treatment.
Each symbol = one genotype; error bars = 1 SD across the three LICOR readings.</div></div>
</section>

<!-- TABLES -->
<section id="tables">
<h2>Summary tables</h2>
<h3>Table 3 - TR best-model parameters</h3>
<p style="font-size:0.88em;color:#555">Slope&#x2081; = initial slope; BP = breakpoint VPD (kPa);
Slope&#x2082; = second-segment slope (bilinear only). &#x2014; = not applicable.</p>',
ua_table(tr_disp, "Table 3. Transpiration rate - best-model parameters per genotype."),
'<h3>Table 4 - g<sub>sw</sub> best-model parameters</h3>
<p style="font-size:0.88em;color:#555">Negative Slope&#x2081; = stomatal closure with increasing VPD. FAILED = too few observations.</p>',
ua_table(gsw_disp, "Table 4. Stomatal conductance - best-model parameters per genotype."),
'<h3>Table 5 - AIC comparison (TR)</h3>',
ua_table(tr_aic_disp, "Table 5. AIC and DeltaAIC for all three TR model types per genotype."),
'<h3>Table 6 - AIC comparison (g<sub>sw</sub>)</h3>',
ua_table(gsw_aic_disp, "Table 6. AIC and DeltaAIC for all three gsw model types per genotype."),
'</section>

<!-- METHODS -->
<section id="methods">
<h2>Methods summary</h2>
<h3>Data</h3>
<p>Source: <code>data/data_all.xlsx</code> (sheet "dataAll"). ', nrow(df_clean),
' observations, 14 genotypes, 4 replicates, 2 days, 6 VPD treatments.
TR (mg m&#x207B;&#xB2; s&#x207B;&#xB9;) = (Initial&#x2009;&#x2212;&#x2009;Final weight) &#xD7; 1000 / (leaf area m&#xB2; &#xD7; duration s).</p>

<h3>Model fitting</h3>
<p>Day-1 data only. Random effect: Rep (replicate) fitted via <code>nlme::lme()</code>
for the linear model and <code>nlme::nlme()</code> for bilinear and linear-plateau models
(nls fallback when nlme fails to converge). Model selection by AIC.</p>

<h3>Scripts</h3>
<ul>
  <li><code>script/01_data_inspection_QC.R</code> &mdash; data loading, QC, diagnostic figures</li>
  <li><code>script/02_model_fitting_TR_gs_VPD.R</code> &mdash; model fitting, publication figures, parameter tables</li>
  <li><code>script/03_generate_report.R</code> &mdash; this HTML report</li>
</ul>
</section>

</main>

<footer>
  Growth Chamber 2026 Soybean LT Experiment &nbsp;|&nbsp;
  Elli Lab - University of Arkansas &nbsp;|&nbsp;
  Generated ', today, '
</footer>
')

full_html <- paste0(
  '<!DOCTYPE html>\n<html lang="en">\n<head>\n',
  '<meta charset="UTF-8">\n',
  '<meta name="viewport" content="width=device-width, initial-scale=1.0">\n',
  '<title>GC2026 Soybean LT - Analysis Report</title>\n',
  '<style>', CSS, '</style>\n',
  '</head>\n<body>\n',
  html_body,
  '\n</body>\n</html>'
)

writeLines(full_html, "GC2026_analysis_report.html")
cat("Report written: GC2026_analysis_report.html (",
    round(file.size("GC2026_analysis_report.html") / 1e6, 1), "MB)\n")
