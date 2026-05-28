"""
update_meta_analysis.py — Add new entries to limited_transpiration_soybean_meta_analysis.xlsx
based on PDFs present in the meta_analysis/ folder that are not yet in the table.
Run from /home/user/LT-trait:
    python3 meta_analysis/update_meta_analysis.py
"""

import os, shutil, copy
from datetime import date
import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

XLSX = "meta_analysis/limited_transpiration_soybean_meta_analysis.xlsx"
BACKUP = f"meta_analysis/limited_transpiration_soybean_meta_analysis_BACKUP_{date.today()}.xlsx"

# ── New entries identified from PDF review ────────────────────────────────────
# Columns (row 2 headers):
# 0  Study ID         1  Authors           2  Year              3  Journal
# 4  DOI/Link         5  Open Access       6  Study Type        7  Genotype(s) Studied
# 8  Key Genotype Class  9  Population/Material  10 Env. Conditions
# 11 VPD Range Tested(kPa)  12 VPD Threshold(kPa)  13 Transpiration Measure
# 14 Phenotyping Method  15 Aquaporin Inhibitor?  16 Stomatal Conductance Data?
# 17 Leaf Hydraulic Conductance?  18 QTL/Genomic Region(s)  19 Chromosomes
# 20 Marker Platform  21 Yield Data?  22 WUE Data?  23 Simulation/Model?
# 24 Temperature Effect Noted?  25 Slow-Wilting Phenotype?
# 26 Main Finding – TRlim Expression  27 Effect on Drought Tolerance  28 Notes/Limitations

NEW_ENTRIES = [
    {   # S16 — Bunce 1984 — earliest known genotypic variation in stomatal VPD sensitivity
        0:  "S16",
        1:  "Bunce",
        2:  1984,
        3:  "Annals of Applied Biology",
        4:  "10.1111/j.1744-7348.1984.tb03017.x",
        5:  "No",
        6:  "Greenhouse + field – gas exchange",
        7:  "Kent, Williams, Clark, Dorman, Lincoln; Kent intra-cultivar lines",
        8:  "VPD-sensitive vs VPD-insensitive",
        9:  "Selected cultivars + intra-cultivar selfed progeny",
        10: "Controlled environment (28/22°C, 14 h light) + field (Beltsville MD)",
        11: "1.0 – 3.0",
        12: "~2.0 (inferred from gₛ reduction between 1.0 and 3.0 kPa)",
        13: "Stomatal conductance + net photosynthesis (steady-state gas exchange)",
        14: "Open gas exchange system; single attached leaf",
        15: "No",
        16: "Yes",
        17: "Indirect (root:leaf ratio used as hydraulic proxy)",
        18: "None",
        19: "N/A",
        20: "N/A",
        21: "No",
        22: "No",
        23: "No",
        24: "No",
        25: "Not assessed",
        26: (
            "First demonstration of genotypic variation in stomatal and photosynthetic "
            "sensitivity to VPD within soybean. Plants with low root:leaf area ratio had "
            "greater reductions in leaf conductance between VPD 1.0 and 3.0 kPa. "
            "Rankings were consistent in controlled and field conditions. Precursor to "
            "the LT trait concept."
        ),
        27: "Suggests that selection for high root:leaf ratio could reduce stomatal sensitivity to VPD, conserving water under high atmospheric demand.",
        28: "Only two VPD levels tested (1.0 and 3.0 kPa); no breakpoint characterisation. Study predates the 'limited transpiration' terminology.",
    },
    {   # S17 — Sinclair, Zwieniecki, Holbrook 2008 Physiologia Plantarum — hydraulic mechanism
        0:  "S17",
        1:  "Sinclair, Zwieniecki & Holbrook",
        2:  2008,
        3:  "Physiologia Plantarum",
        4:  "10.1111/j.1399-3054.2007.01028.x",
        5:  "No",
        6:  "Greenhouse – hydraulic physiology",
        7:  "PI 416937; two comparison genotypes",
        8:  "Slow-wilting (LT) vs commercial",
        9:  "Contrasting pair",
        10: "Controlled environment; varied xylem water potential",
        11: "Not reported (mechanism study, not VPD-threshold screening)",
        12: "~2.0 (VPD response confirmed in PI 416937)",
        13: "Transpiration rate; leaf rehydration; carbon assimilation after petiole cutting",
        14: "Gravimetric + gas exchange; leaf hydraulic conductance measurement",
        15: "No",
        16: "Yes",
        17: "Yes — core result: PI 416937 has low leaf hydraulic conductance",
        18: "None",
        19: "N/A",
        20: "N/A",
        21: "No",
        22: "No",
        23: "No",
        24: "No",
        25: "Yes",
        26: (
            "Confirms VPD-threshold response in PI 416937. Three experiments demonstrate that "
            "low leaf hydraulic conductance is the underlying mechanism: (1) VPD response persists "
            "when xylem water potential is maintained at zero; (2) slower leaf rehydration in "
            "PI 416937; (3) more rapid reduction in carbon assimilation after petiole cutting. "
            "Provides mechanistic basis for the LT trait: restricted leaf hydraulic conductance "
            "limits transpiration at high VPD while allowing normal gas exchange at moderate VPD."
        ),
        27: "Low leaf hydraulic conductance conserves water under high evaporative demand; allows normal gas exchange under moderate VPD, preserving carbon gain.",
        28: (
            "Note: this Physiologia Plantarum 132:446-451 paper is distinct from the same "
            "group's Plant Cell & Environment 31:796-803 (2008) paper (S02). Both describe "
            "hydraulic conductance mechanisms but use different experimental approaches."
        ),
    },
    {   # S18 — Seversike et al. 2012/2013 Physiologia Plantarum — temperature × VPD
        0:  "S18",
        1:  "Seversike, Sermons, Sinclair, Carter & Rufty",
        2:  2013,
        3:  "Physiologia Plantarum",
        4:  "10.1111/j.1399-3054.2012.01693.x",
        5:  "No",
        6:  "Greenhouse – temperature × VPD interaction",
        7:  "PI 416937; Hutcheson; 4 Glycine soja genotypes",
        8:  "Slow-wilting G. max + wild soybean G. soja vs commercial G. max",
        9:  "Two G. max cultivars + four wild soybean (G. soja) accessions",
        10: "Controlled environment; 25°C or 30°C growth temperature with abrupt 5°C increase",
        11: "1.0 – 4.0+ (wide VPD range)",
        12: (
            "~2.0 kPa in G. max at 30°C (growth optimum); G. soja threshold lower (~1.5 kPa) "
            "at 25°C; threshold lost or shifted after temperature increase to 35°C (G. max) "
            "or 30°C (G. soja)"
        ),
        13: "Transpiration rate (gravimetric, whole-plant)",
        14: "Whole-plant gravimetric; growth chamber temperature manipulation",
        15: "No",
        16: "No",
        17: "No",
        18: "None",
        19: "N/A",
        20: "N/A",
        21: "No",
        22: "No",
        23: "No",
        24: "Yes — central finding: LT response is temperature-sensitive",
        25: "Yes (PI 416937)",
        26: (
            "TR control in G. max (PI 416937) is evident above ~2.0 kPa at 30°C (growth optimum) "
            "but is lost when temperature shifts to 35°C. Three of four G. soja genotypes show TR "
            "control at 25°C (their growth optimum) at lower VPD threshold than G. max, but lose "
            "control at 30°C. LT response is kinetically definable and temperature-sensitive, "
            "suggesting a temperature-regulated hydraulic or biochemical mechanism. Only PI 416937 "
            "retained TR control after shift to 35°C among G. max genotypes."
        ),
        27: "TR control under high VPD is temperature-dependent; beneficial only near growth optimum. G. soja offers genetic diversity for the LT trait but requires decoupling from low-temperature dependency for use in warm climates.",
        28: (
            "Distinct from S14 (Seversike et al. 2015 Crop Science). This is Physiologia Plantarum "
            "148:62-73 (2013). Highlights that LT trait expression in warm production environments "
            "(e.g., US Mid-South summer) may be compromised if temperatures regularly exceed 33°C."
        ),
    },
    {   # S19 — Sadok & Sinclair 2009 Crop Science — genetic variability across broader germplasm
        0:  "S19",
        1:  "Sadok & Sinclair",
        2:  2009,
        3:  "Crop Science",
        4:  "10.2135/cropsci2008.09.0535",
        5:  "No",
        6:  "Greenhouse – transpiration phenotyping",
        7:  "PI 416937; 3 broad-background cultivars; 4 PI 416937 progeny lines",
        8:  "LT (PI 416937 and progenies) vs commercial cultivars",
        9:  "7 genotypes: 3 commercial cultivars + 4 progeny lines",
        10: "Controlled environment; full VPD range",
        11: "~1.0 – 3.5",
        12: "~2.0 kPa (two-segment response in selected genotypes)",
        13: "Whole-plant transpiration rate (gravimetric)",
        14: "Gravimetric whole-plant; VPD ramp protocol",
        15: "No",
        16: "No",
        17: "No",
        18: "None",
        19: "N/A",
        20: "N/A",
        21: "No",
        22: "No",
        23: "No",
        24: "No",
        25: "Yes (PI 416937 and some progeny)",
        26: (
            "Extended the LT trait screen to a wider genetic base (7 genotypes). Two cultivars "
            "without PI 416937 in their pedigree showed two-segment TR response to VPD "
            "(breakpoint ~2.0 kPa), demonstrating that the LT trait has more than one genetic "
            "source. Three PI 416937 progeny lines expressing slow-wilting in the field did NOT "
            "show limited TR, indicating slow-wilting can arise from mechanisms other than "
            "VPD-threshold transpiration limitation."
        ),
        27: "Multiple genetic sources for VPD-threshold TR limitation identified. Slow-wilting field phenotype and controlled-environment LT expression are not perfectly correlated, implying complex genetic architecture.",
        28: (
            "Note: this is a DIFFERENT paper from S03 (Sadok & Sinclair 2009, JXB erp350) and S04 "
            "(Sadok & Sinclair 2010, Crop Science). This Crop Science 2009 paper (published May) "
            "doi 10.2135/cropsci2008.09.0535 focuses on broader germplasm screening."
        ),
    },
]

# ── Update workbook ───────────────────────────────────────────────────────────
shutil.copy2(XLSX, BACKUP)
print(f"Backup saved: {BACKUP}")

wb = openpyxl.load_workbook(XLSX)
ws = wb["Meta-Analysis Table"]

# Find the last row with data
last_row = 2  # header is row 2
for row in ws.iter_rows(min_row=3):
    if any(c.value is not None for c in row):
        last_row = row[0].row

print(f"Last data row: {last_row}  →  appending from row {last_row+1}")

# Get header row to map column names to indices
headers = {ws.cell(2, j+1).value: j+1 for j in range(ws.max_column)}

# Style from last data row
def copy_style(src_row_idx, dst_row_idx):
    for col in range(1, ws.max_column + 1):
        src = ws.cell(src_row_idx, col)
        dst = ws.cell(dst_row_idx, col)
        if src.has_style:
            dst.font      = copy.copy(src.font)
            dst.fill      = copy.copy(src.fill)
            dst.border    = copy.copy(src.border)
            dst.alignment = copy.copy(src.alignment)

for i, entry in enumerate(NEW_ENTRIES):
    new_row = last_row + 1 + i
    copy_style(last_row, new_row)
    for col_idx, value in entry.items():
        ws.cell(new_row, col_idx + 1).value = value
    print(f"  Added {entry[0]}: {entry[1]} ({entry[2]})")

# Update Summary & Notes sheet
ws2 = wb["Summary & Notes"]
n_entries = last_row - 2 + len(NEW_ENTRIES)  # minus 2 header rows
for row in ws2.iter_rows():
    for cell in row:
        if isinstance(cell.value, str) and "studies included" in cell.value.lower():
            cell.value = f"Studies included (total): {n_entries}"
            break

wb.save(XLSX)
print(f"\nSaved updated workbook: {XLSX}")
print(f"Total entries now: {n_entries}")
