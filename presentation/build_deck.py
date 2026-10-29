#!/usr/bin/env python3
"""Build the ASA 2026 deck (ASA-2026-Elli-LT-trait.pptx at the repo root).

Starts from last year's ASA deck (presentation/template/ASA-2025-Elli-template.pptx),
keeps its title, intro, LT-trait, team and thank-you slides, and adds new slides
built from the figures in the repo (growth chamber 2026 + SSM simulations) and the
static images in presentation/img/ (taken from the Corteva poster and the 2025 deck).

Needs the pptx skill's helper scripts (add_slide.py, clean.py). Point
PPTX_SKILL_DIR at the skill folder, otherwise it is searched under ~/.claude/skills.

Run from anywhere:  python3 presentation/build_deck.py
Then check:         python3 $PPTX_SKILL_DIR/scripts/office/validate.py ASA-2026-Elli-LT-trait.pptx \
                        --original presentation/template/ASA-2025-Elli-template.pptx
"""
import os, re, shutil, subprocess, sys, html, glob, tempfile, zipfile
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
TEMPLATE = os.path.join(HERE, "template", "ASA-2025-Elli-template.pptx")
IMG = os.path.join(HERE, "img")
PF = f"{REPO}/experimental_data_soy_2026/figures/paper_figures"
MF = f"{REPO}/modeling_data/figures"
OUT = os.path.join(REPO, "ASA-2026-Elli-LT-trait.pptx")

def _find_skill():
    env = os.environ.get("PPTX_SKILL_DIR")
    if env and os.path.exists(os.path.join(env, "scripts", "add_slide.py")):
        return env
    hits = glob.glob(os.path.expanduser("~/.claude/skills/**/pptx/scripts/add_slide.py"), recursive=True)
    if not hits:
        sys.exit("pptx skill not found: set PPTX_SKILL_DIR to the folder that contains scripts/add_slide.py")
    return os.path.dirname(os.path.dirname(hits[0]))
SK = _find_skill()

## the template is unpacked into a temporary folder on each run
_TMP = tempfile.mkdtemp(prefix="asa_deck_")
SRC = os.path.join(_TMP, "template")
WORK = os.path.join(_TMP, "work")
zipfile.ZipFile(TEMPLATE).extractall(SRC)

EMU = 914400
RED = "A61E36"
DARK = "333333"
GREY = "BFBFBF"
MUTED = "6B6B6B"
TINT = "F7EEF0"     # very light red tint for cards
LIGHT = "F2F2F2"
GREEN = "31A354"

def E(inches):
    return int(round(inches * EMU))

def esc(s):
    return html.escape(s, quote=False)

# ---------------------------------------------------------------- XML pieces
CAL = '<a:latin typeface="Calibri" panose="020F0502020204030204" pitchFamily="34" charset="0"/><a:cs typeface="Calibri" panose="020F0502020204030204" pitchFamily="34" charset="0"/>'

def rpr(sz, b=False, i=False, color=DARK, u=False):
    return (f'<a:rPr lang="en-US" sz="{sz}"{" b=\"1\"" if b else ""}{" i=\"1\"" if i else ""}{" u=\"sng\"" if u else ""} dirty="0">'
            f'<a:solidFill><a:srgbClr val="{color}"/></a:solidFill>{CAL}</a:rPr>')

def run(text, sz=1600, b=False, i=False, color=DARK, u=False):
    return f'<a:r>{rpr(sz,b,i,color,u)}<a:t>{esc(text)}</a:t></a:r>'

def para(runs, algn="l", bullet=False, spc_after=600, spc_before=0, line=None, lvl_indent=None):
    """runs: str or list of run() strings. bullet -> red round bullet."""
    if isinstance(runs, str):
        runs = run(runs)
    elif isinstance(runs, list):
        runs = "".join(runs)
    ppr = f'<a:pPr algn="{algn}"'
    if bullet:
        ind = lvl_indent or 285750
        ppr += f' marL="{ind}" indent="-{ind}"'
    ppr += '>'
    if line:
        ppr += f'<a:lnSpc><a:spcPct val="{line}"/></a:lnSpc>'
    if spc_before:
        ppr += f'<a:spcBef><a:spcPts val="{spc_before}"/></a:spcBef>'
    ppr += f'<a:spcAft><a:spcPts val="{spc_after}"/></a:spcAft>'
    if bullet:
        ppr += f'<a:buClr><a:srgbClr val="{RED}"/></a:buClr><a:buFont typeface="Arial" panose="020B0604020202020204" pitchFamily="34" charset="0"/><a:buChar char="&#8226;"/>'
    else:
        ppr += '<a:buNone/>'
    ppr += '</a:pPr>'
    return f'<a:p>{ppr}{runs}</a:p>'


class Slide:
    def __init__(self, name):
        self.name = name
        self.shapes = []
        self.images = []   # (rId, media filename)
        self._id = 2
        self.notes = None

    def nid(self):
        self._id += 1
        return self._id

    # --- primitives -------------------------------------------------------
    def text(self, x, y, w, h, paras, anchor="t", fill=None, line=None, ins=(0.1, 0.05, 0.1, 0.05),
             autofit=False, rot=0, radius=None, name="TextBox", shadow=False):
        i = self.nid()
        geom = 'roundRect' if radius is not None else 'rect'
        av = f'<a:gd name="adj" fmla="val {radius}"/>' if radius is not None else ''
        sppr = f'<a:xfrm{" rot=\"%d\"" % rot if rot else ""}><a:off x="{E(x)}" y="{E(y)}"/><a:ext cx="{E(w)}" cy="{E(h)}"/></a:xfrm><a:prstGeom prst="{geom}"><a:avLst>{av}</a:avLst></a:prstGeom>'
        sppr += f'<a:solidFill><a:srgbClr val="{fill}"/></a:solidFill>' if fill else '<a:noFill/>'
        sppr += f'<a:ln w="12700"><a:solidFill><a:srgbClr val="{line}"/></a:solidFill></a:ln>' if line else '<a:ln><a:noFill/></a:ln>'
        if shadow:
            sppr += '<a:effectLst><a:outerShdw blurRad="50800" dist="25400" dir="5400000" algn="t" rotWithShape="0"><a:srgbClr val="000000"><a:alpha val="20000"/></a:srgbClr></a:outerShdw></a:effectLst>'
        l, t, r, b = [E(v) for v in ins]
        body = f'<a:bodyPr wrap="square" lIns="{l}" tIns="{t}" rIns="{r}" bIns="{b}" rtlCol="0" anchor="{anchor}">{"<a:normAutofit/>" if autofit else "<a:noAutofit/>"}</a:bodyPr><a:lstStyle/>'
        if isinstance(paras, str):
            paras = [para(paras)]
        self.shapes.append(
            f'<p:sp><p:nvSpPr><p:cNvPr id="{i}" name="{name} {i}"/><p:cNvSpPr txBox="1"/><p:nvPr/></p:nvSpPr>'
            f'<p:spPr>{sppr}</p:spPr><p:txBody>{body}{"".join(paras)}</p:txBody></p:sp>')
        return i

    def rect(self, x, y, w, h, fill, line=None, radius=None, name="Shape", shadow=False):
        return self.text(x, y, w, h, [para([run("", 800)])], fill=fill, line=line, radius=radius, name=name, shadow=shadow)

    def hline(self, x, y, w, color=RED, wt=28575):
        i = self.nid()
        self.shapes.append(
            f'<p:cxnSp><p:nvCxnSpPr><p:cNvPr id="{i}" name="Line {i}"/><p:cNvCxnSpPr/><p:nvPr/></p:nvCxnSpPr>'
            f'<p:spPr><a:xfrm><a:off x="{E(x)}" y="{E(y)}"/><a:ext cx="{E(w)}" cy="0"/></a:xfrm><a:prstGeom prst="line"><a:avLst/></a:prstGeom>'
            f'<a:ln w="{wt}"><a:solidFill><a:srgbClr val="{color}"/></a:solidFill></a:ln></p:spPr></p:cxnSp>')

    def image(self, path, x, y, w, h, align="c", valign="c", crop=None):
        """Fit image in (x,y,w,h) preserving aspect ratio."""
        im = Image.open(path)
        iw, ih = im.size
        if crop:   # (l,t,r,b) fractions
            iw = iw * (1 - crop[0] - crop[2]); ih = ih * (1 - crop[1] - crop[3])
        s = min(w / iw, h / ih)
        fw, fh = iw * s, ih * s
        fx = x + (w - fw) / 2 if align == "c" else (x if align == "l" else x + w - fw)
        fy = y + (h - fh) / 2 if valign == "c" else (y if valign == "t" else y + h - fh)
        rid = f"rId{len(self.images) + 2}"
        self.images.append((rid, path))
        i = self.nid()
        src = ''
        if crop:
            src = f'<a:srcRect l="{int(crop[0]*100000)}" t="{int(crop[1]*100000)}" r="{int(crop[2]*100000)}" b="{int(crop[3]*100000)}"/>'
        self.shapes.append(
            f'<p:pic><p:nvPicPr><p:cNvPr id="{i}" name="Picture {i}"/><p:cNvPicPr><a:picLocks noChangeAspect="1"/></p:cNvPicPr><p:nvPr/></p:nvPicPr>'
            f'<p:blipFill><a:blip r:embed="{rid}"/>{src}<a:stretch><a:fillRect/></a:stretch></p:blipFill>'
            f'<p:spPr><a:xfrm><a:off x="{E(fx)}" y="{E(fy)}"/><a:ext cx="{E(fw)}" cy="{E(fh)}"/></a:xfrm><a:prstGeom prst="rect"><a:avLst/></a:prstGeom></p:spPr></p:pic>')
        return fx, fy, fw, fh

    def table(self, x, y, colw, rows, header_fill=RED, row_h=0.42, sz=1400, band=True):
        """rows[0] is header. colw list of inches."""
        i = self.nid()
        grid = "".join(f'<a:gridCol w="{E(w)}"/>' for w in colw)
        trs = []
        for r, row in enumerate(rows):
            tcs = []
            for c, cell in enumerate(row):
                hdr = (r == 0)
                txt = run(cell, sz, b=hdr, color=("FFFFFF" if hdr else DARK))
                fill = header_fill if hdr else (LIGHT if (band and r % 2 == 0) else "FFFFFF")
                tcs.append(
                    f'<a:tc><a:txBody><a:bodyPr/><a:lstStyle/><a:p><a:pPr algn="{"l" if c == 0 else "ctr"}"/>{txt}</a:p></a:txBody>'
                    f'<a:tcPr marL="91440" marR="91440" marT="36000" marB="36000" anchor="ctr">'
                    f'<a:lnL><a:noFill/></a:lnL><a:lnR><a:noFill/></a:lnR><a:lnT><a:noFill/></a:lnT><a:lnB w="6350"><a:solidFill><a:srgbClr val="D9D9D9"/></a:solidFill></a:lnB>'
                    f'<a:solidFill><a:srgbClr val="{fill}"/></a:solidFill></a:tcPr></a:tc>')
            trs.append(f'<a:tr h="{E(row_h)}">{"".join(tcs)}</a:tr>')
        self.shapes.append(
            f'<p:graphicFrame><p:nvGraphicFramePr><p:cNvPr id="{i}" name="Table {i}"/><p:cNvGraphicFramePr><a:graphicFrameLocks noGrp="1"/></p:cNvGraphicFramePr><p:nvPr/></p:nvGraphicFramePr>'
            f'<p:xfrm><a:off x="{E(x)}" y="{E(y)}"/><a:ext cx="{E(sum(colw))}" cy="{E(row_h*len(rows))}"/></p:xfrm>'
            f'<a:graphic><a:graphicData uri="http://schemas.openxmlformats.org/drawingml/2006/table"><a:tbl><a:tblPr firstRow="1" bandRow="1"/>'
            f'<a:tblGrid>{grid}</a:tblGrid>{"".join(trs)}</a:tbl></a:graphicData></a:graphic></p:graphicFrame>')

    # --- composites -------------------------------------------------------
    def title(self, t, sz=3200):
        # same position/style as the template's content titles
        self.text(0.49, 0.14, 12.35, 0.64, [para([run(t, sz, b=True, color=RED)], spc_after=0)], anchor="ctr", name="Title")

    def card(self, x, y, w, h, head, body, head_sz=1800, body_sz=1500, fill=TINT, head_color=RED):
        ps = [para([run(head, head_sz, b=True, color=head_color)], spc_after=300)]
        if isinstance(body, str):
            body = [body]
        for b in body:
            ps.append(para([run(b, body_sz, color=DARK)], spc_after=200, line=100000))
        self.text(x, y, w, h, ps, anchor="ctr", fill=fill, radius=6000, ins=(0.18, 0.12, 0.18, 0.12), name="Card")

    def stat(self, x, y, w, h, number, label, num_sz=4000, lab_sz=1400, fill=LIGHT, num_color=RED):
        ps = [para([run(number, num_sz, b=True, color=num_color)], algn="ctr", spc_after=200),
              para([run(label, lab_sz, color=DARK)], algn="ctr", spc_after=0, line=100000)]
        self.text(x, y, w, h, ps, anchor="ctr", fill=fill, radius=6000, ins=(0.12, 0.1, 0.12, 0.1), name="Stat")

    def caption(self, x, y, w, t, algn="l"):
        self.text(x, y, w, 0.5, [para([run(t, 1200, i=True, color=MUTED)], algn=algn, spc_after=0)], anchor="t", name="Caption")

    def objective_rows(self, items, active=None, icon=f"{SRC}/ppt/media/image8.png", top=1.55, sz=3200, gap=1.75):
        """Template-style objective list: target icon + text; non-active greyed."""
        for k, (txt) in enumerate(items):
            on = (active is None) or (active == k)
            col = DARK if on else GREY
            y = top + k * gap
            fx, fy, fw, fh = self.image(icon, 0.6, y, 1.25, 1.25)
            self.text(2.1, y - 0.05, 10.6, 1.35, [para([run(txt, sz, color=col)], spc_after=0, line=95000)], anchor="ctr")

    # --- write ------------------------------------------------------------
    def xml(self):
        return ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
                '<p:sld xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships" xmlns:p="http://schemas.openxmlformats.org/presentationml/2006/main">'
                '<p:cSld><p:spTree><p:nvGrpSpPr><p:cNvPr id="1" name=""/><p:cNvGrpSpPr/><p:nvPr/></p:nvGrpSpPr>'
                '<p:grpSpPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="0" cy="0"/><a:chOff x="0" y="0"/><a:chExt cx="0" cy="0"/></a:xfrm></p:grpSpPr>'
                + "".join(self.shapes) +
                '</p:spTree></p:cSld><p:clrMapOvr><a:masterClrMapping/></p:clrMapOvr></p:sld>')

    def rels(self, media_names):
        r = ['<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/slideLayout" Target="../slideLayouts/slideLayout7.xml"/>']
        for (rid, _), mn in zip(self.images, media_names):
            r.append(f'<Relationship Id="{rid}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/image" Target="../media/{mn}"/>')
        return ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
                '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">' + "".join(r) + '</Relationships>')


# ====================================================================== content
def build_slides():
    S = []

    # --- 3. Rising atmospheric demand ------------------------------------
    s = Slide("vpd_rising"); S.append(s)
    s.title("Atmospheric demand on soybean is rising")
    s.image(f"{IMG}/vpd_projection.png", 0.5, 1.05, 6.3, 3.5)
    s.caption(0.5, 4.6, 6.3, "July mean VPD in Arkansas, UKESM1-0-LL model, SSP3-7.0 scenario")
    s.image(f"{IMG}/vpd_heatmap.png", 7.1, 1.25, 5.8, 3.0)
    s.caption(7.1, 4.4, 5.8, "Hours per day above 2 kPa, 10 sites in Arkansas, Missouri and Nebraska")
    s.card(0.5, 5.25, 6.3, 1.65, "+0.8% per year",
           "Projected increase in July VPD through 2100: more hours of high evaporative demand each season")
    s.card(7.1, 5.25, 5.8, 1.65, "Why it matters",
           "How genotypes respond to VPD is the basis for strategies that optimize seasonal water use and raise water use efficiency (WUE)")

    # --- 6. Hypothesis -----------------------------------------------------
    s = Slide("hypothesis"); S.append(s)
    s.title("Hypothesis: WUE plasticity is a signature of the LT trait")
    s.text(0.5, 1.15, 6.4, 1.9, [
        para([run("Chamarthi et al. (2022)", 1800, b=True, color=RED), run(" classified soybean accessions by the plasticity of intrinsic WUE across environments, inferred from carbon isotope discrimination", 1800)], spc_after=600, line=100000),
    ], anchor="t")
    s.card(0.5, 3.15, 6.4, 1.55, "Plastic, high WUE",
           "WUE rises with evaporative demand: expected to limit transpiration under high VPD", fill=TINT)
    s.card(0.5, 4.9, 6.4, 1.55, "Nonplastic, low WUE",
           "WUE stays flat across environments: expected to keep transpiring as VPD increases", fill=LIGHT, head_color=DARK)
    s.image(f"{IMG}/hypothesis.png", 7.3, 1.1, 5.6, 4.9)
    s.caption(7.3, 6.1, 5.6, "Expected relationship between WUE plasticity and the limited transpiration (LT) response", algn="ctr")

    # --- 7. Objectives -----------------------------------------------------
    s = Slide("objectives"); S.append(s)
    s.title("Objectives")
    s.objective_rows(["Identify genotypes with the LT trait among accessions differing in WUE plasticity",
                      "Quantify the yield impact and trade-offs of the LT trait across the US Mid-South and Mid-West"], top=1.7, gap=2.3)

    # --- 8. Materials and methods ------------------------------------------
    s = Slide("methods"); S.append(s)
    s.title("Material and Methods")
    # row 1: controlled environment
    s.rect(0.5, 1.0, 12.35, 2.75, "EAF4F7", radius=5000)
    s.text(0.65, 1.05, 2.0, 0.45, [para([run("Objective 1", 1600, b=True, color=RED)], spc_after=0)], anchor="ctr")
    s.image(f"{SRC}/ppt/media/image11.jpg", 0.7, 1.5, 1.75, 2.15)
    s.image(f"{IMG}/chamber.png", 2.6, 1.3, 3.1, 2.35)
    s.text(5.9, 1.1, 6.85, 2.6, [
        para([run("Controlled-environment experiment", 1800, b=True, color=DARK)], spc_after=400),
        para([run("14 soybean MG IV accessions: 5 plastic / high WUE, 4 nonplastic / low WUE, slow- and high-wilting checks, 1 MAGIC parent", 1500)], bullet=True, spc_after=200, line=100000),
        para([run("Stepwise VPD from 1.0 to 3.5 kPa at constant temperature, 5 replications (Dec 2025 to Aug 2026)", 1500)], bullet=True, spc_after=200),
        para([run("Gravimetric transpiration per leaf area (TR); stomatal conductance (gsw) on 13 accessions", 1500)], bullet=True, spc_after=200),
        para([run("Linear, linear-plateau and bilinear mixed models (replicate as random effect, R nlme) compared by AIC to locate VPD breakpoints", 1500)], bullet=True, spc_after=0, line=100000),
    ], anchor="ctr")
    # row 2: modeling
    s.rect(0.5, 3.95, 12.35, 2.95, "EAF4F7", radius=5000)
    s.text(0.65, 4.0, 2.0, 0.45, [para([run("Objective 2", 1600, b=True, color=RED)], spc_after=0)], anchor="ctr")
    s.text(0.7, 4.5, 5.0, 2.3, [
        para([run("Systems-modeling experiment", 1800, b=True, color=DARK)], spc_after=400),
        para([run("SSM soybean model (Soltani and Sinclair, 2012)", 1500)], bullet=True, spc_after=200),
        para([run("Biomass: radiation use efficiency; water use: hourly VPD-driven transpiration efficiency", 1500)], bullet=True, spc_after=200, line=100000),
        para([run("LT vs non-LT genotype, rainfed vs irrigated", 1500)], bullet=True, spc_after=200),
        para([run("10 locations (AR, MO, NE) x 30 weather-years x 3 sowing dates", 1500)], bullet=True, spc_after=0),
    ], anchor="ctr")
    s.image(f"{IMG}/acreage_map.png", 5.9, 4.05, 2.9, 2.8)
    s.image(f"{IMG}/vpd_heatmap.png", 8.9, 4.3, 3.85, 2.0)
    s.caption(8.9, 6.3, 3.85, "Hours above 2 kPa, study sites", algn="ctr")

    # --- 9. Objective 1 marker ----------------------------------------------
    s = Slide("obj1"); S.append(s)
    s.title("Objectives")
    s.objective_rows(["Identify genotypes with the LT trait among accessions differing in WUE plasticity",
                      "Quantify the yield impact and trade-offs of the LT trait across the US Mid-South and Mid-West"], active=0, top=1.7, gap=2.3)

    # --- 10. TR response 9 panels -------------------------------------------
    s = Slide("tr_response"); S.append(s)
    s.title("Diverging transpiration responses to VPD across accessions")
    s.image(f"{REPO}/experimental_data_soy_2026/figures/presentation/Fig_TR_response_VPD_presentation.png", 0.4, 0.9, 7.4, 5.4)
    s.caption(0.4, 6.35, 7.4, "Plant-level transpiration under increasing VPD, Day 1 of each run, accessions with at least 3 replications; lines are the AIC-selected mixed model")
    s.card(8.1, 1.0, 4.8, 1.75, "4 of 9 limit transpiration",
           "A linear-plateau or bilinear model beats the linear model above a VPD breakpoint")
    s.card(8.1, 2.95, 4.8, 1.75, "Checks mostly behave as expected",
           "Slow-wilting PI471938 plateaus at 2.6 kPa; A5959 keeps increasing; PI416937 shows no breakpoint", fill=LIGHT, head_color=DARK)
    s.card(8.1, 4.9, 4.8, 1.75, "Plasticity is an imperfect guide",
           "3 of 4 high-plasticity accessions respond (PI398492, PI398617, PI603418D); PI424231 does not")

    # --- 11. Breakpoints ----------------------------------------------------
    s = Slide("breakpoints"); S.append(s)
    s.title("Limited transpiration emerges between 1.7 and 3.4 kPa")
    s.image(f"{PF}/Figure2_TR_breakpoints.png", 0.4, 1.0, 7.2, 5.6, crop=(0.0, 0.0, 0.0, 0.11))
    s.caption(0.4, 6.6, 7.2, "VPD breakpoint of the eight accessions with a linear-plateau or bilinear response")
    s.stat(7.9, 1.1, 2.4, 1.6, "8 of 14", "accessions showed a nonlinear TR response (57%)")
    s.stat(10.5, 1.1, 2.4, 1.6, "4 of 5", "of the plastic / high WUE accessions are among them")
    s.stat(7.9, 2.9, 5.0, 1.45, "2.76 kPa", "mean breakpoint across the eight accessions")
    s.card(7.9, 4.55, 5.0, 2.1, "Early closers",
           ["PI398617 (1.69 kPa) restricts transpiration before typical midday field VPD of 2.0 to 2.5 kPa; PI377574, PI471938 and PI404191 follow between 2.3 and 2.7 kPa",
            "PI603418D only declines at the last step (3.4 kPa)"])

    # --- 12. Quadrants -----------------------------------------------------
    s = Slide("quadrants"); S.append(s)
    s.title("Two routes to low water use under high VPD")
    s.image(f"{PF}/Figure3_intercept_slope_quadrants.png", 0.4, 1.0, 7.4, 5.6)
    s.caption(0.4, 6.6, 7.4, "Initial slope vs intercept of the TR response; dashed lines are the medians")
    s.card(8.1, 1.1, 4.8, 1.65, "Low baseline, strong response",
           "PI398617 and PI377574: intercept below 12 and slope above 25 mg m\u207b\u00b2 s\u207b\u00b9 kPa\u207b\u00b9, with early breakpoints")
    s.card(8.1, 2.95, 4.8, 1.65, "Slow-wilting check PI471938",
           "Sits in the low-intercept / high-sensitivity quadrant and plateaus at 2.6 kPa")
    s.card(8.1, 4.8, 4.8, 1.65, "Plasticity class is not enough",
           "The five plastic / high WUE accessions spread across three quadrants", fill=LIGHT, head_color=DARK)

    # --- 13. Stomatal conductance --------------------------------------------
    s = Slide("gsw"); S.append(s)
    s.title("Stomata close in every genotype, but transpiration does not")
    s.image(f"{PF}/Figure5_gsw_predicted.png", 0.4, 0.95, 12.5, 3.9)
    s.caption(0.4, 4.85, 12.5, "Predicted stomatal conductance (gsw) at 1.5 and 3.5 kPa for the 13 accessions measured with a porometer", algn="ctr")
    s.stat(0.5, 5.3, 3.9, 1.6, "\u221247% to \u2212100%", "decline in gsw from 1.5 to 3.5 kPa, in all 13 accessions")
    s.stat(4.7, 5.3, 3.9, 1.6, "6 of 14", "accessions still increased transpiration linearly up to 3.5 kPa")
    s.card(8.9, 5.3, 4.0, 1.6, "Beyond the stomata",
           "At high demand, transpiration follows the vapor pressure gradient and boundary layer, not aperture alone", body_sz=1400)

    # --- 14. Objective 1 summary ----------------------------------------------
    s = Slide("obj1_summary"); S.append(s)
    s.title("Candidate LT accessions identified")
    rows = [["Accession", "WUE class (Chamarthi et al., 2022)", "Breakpoint (kPa)", "Response"],
            ["PI398617", "Plastic / high WUE", "1.69", "Bilinear"],
            ["PI377574", "Nonplastic / low WUE", "2.34", "Plateau"],
            ["PI471938", "Slow-wilting check", "2.60", "Plateau"],
            ["PI404191", "Nonplastic / low WUE", "2.68", "Plateau"],
            ["PI424231", "Plastic / high WUE", "3.06", "Plateau"],
            ["AG43XF2", "High-wilting check", "3.13", "Plateau"],
            ["PI398492", "Plastic / high WUE", "3.14", "Plateau"],
            ["PI603418D", "Plastic / high WUE", "3.42", "Decline at 3.4 kPa"]]
    s.table(0.5, 1.1, [1.6, 3.2, 1.6, 2.0], rows, row_h=0.44, sz=1300)
    s.caption(0.5, 5.15, 8.4, "Linear up to 3.5 kPa: A5959, PI398198, PI398272, PI416937, PI442012A, PI603171")
    s.stat(9.2, 1.1, 3.7, 1.45, "8 of 14", "accessions limit transpiration under high VPD")
    s.stat(9.2, 2.75, 3.7, 1.45, "4 of 5", "plastic / high WUE accessions carry the LT response")
    s.stat(9.2, 4.4, 3.7, 1.45, "1.7 to 3.4 kPa", "breakpoint range: a tunable trait, not a switch")
    s.card(0.5, 5.75, 8.4, 1.15, "Take-home",
           "WUE plasticity flags LT candidates, but phenotyping the VPD response is needed to confirm the trait and its threshold", body_sz=1400)

    # --- 15. Objective 2 marker ---------------------------------------------
    s = Slide("obj2"); S.append(s)
    s.title("Objectives")
    s.objective_rows(["Identify genotypes with the LT trait among accessions differing in WUE plasticity",
                      "Quantify the yield impact and trade-offs of the LT trait across the US Mid-South and Mid-West"], active=1, top=1.7, gap=2.3)

    # --- 16. Yield gap maps ----------------------------------------------------
    s = Slide("yield_maps"); S.append(s)
    s.title("The LT trait narrows the yield gap caused by water deficit")
    s.image(f"{MF}/Fig_map_yields_10loc_A.png", 0.3, 0.9, 9.3, 5.1, align="l")
    s.caption(0.4, 6.1, 9.0, "Water-limited yield (green) and yield gap by water deficit (yellow) as shares of potential yield at 10 sites, 30 weather-years x 3 sowing dates each; labels give the water-limited share")
    s.stat(9.7, 1.1, 3.2, 2.15, "61%", "of potential yield lost to water deficit without the LT trait", num_sz=4400)
    s.stat(9.7, 3.5, 3.2, 2.15, "37%", "lost with the LT trait", num_sz=4400, fill=TINT)
    s.text(9.7, 5.9, 3.2, 1.0, [para([run("Largest gains in Missouri and Nebraska; Arkansas sites keep the widest gap even with the LT trait", 1300, color=DARK)], spc_after=0, line=100000)], anchor="t")

    # --- 17. Trade-off -----------------------------------------------------------
    s = Slide("tradeoff"); S.append(s)
    s.title("The trade-off: yield potential for water-limited yield")
    s.image(f"{MF}/Fig_map_yields_10loc_B.png", 0.5, 1.0, 5.0, 5.0)
    s.caption(0.5, 6.05, 5.0, "Mean attainable yield (green) plus yield gap by water deficit (yellow) across 10 locations; error bars are SE among sites", algn="ctr")
    s.stat(5.9, 1.1, 3.35, 2.1, "+8%", "water-limited yield with the LT trait", num_sz=4400, fill=TINT)
    s.stat(9.55, 1.1, 3.35, 2.1, "+50%", "potential yield without the LT trait", num_sz=4400)
    s.card(5.9, 3.45, 7.0, 1.55, "Where it pays",
           "Rainfed fields and dry years: conserved soil water carries the crop through seed filling")
    s.card(5.9, 5.2, 7.0, 1.55, "Where it costs",
           "Irrigated or wet seasons: restricted transpiration caps carbon gain and yield potential", fill=LIGHT, head_color=DARK)

    # --- 18. Conclusions -------------------------------------------------------
    s = Slide("conclusions"); S.append(s)
    s.title("Conclusions")
    items = ["Transpiration restriction under high VPD partly explains the WUE plasticity inferred from carbon isotopes (Chamarthi et al., 2022)",
             "Most high-plasticity accessions limited transpiration, with breakpoints from 1.7 to 3.4 kPa",
             "Not all plastic accessions showed the response: other water-conservation mechanisms also shape WUE plasticity",
             "The LT trait trades yield potential for water-limited yield: +8% under deficit and a yield gap cut from 61% to 37% across 10 sites"]
    for k, t in enumerate(items):
        y = 1.1 + k * 1.45
        s.image(f"{SRC}/ppt/media/image8.png", 0.6, y + 0.1, 1.0, 1.0)
        s.text(1.85, y, 10.9, 1.25, [para([run(t, 2000, color=DARK)], spc_after=0, line=95000)], anchor="ctr")

    # --- 19. Future directions ---------------------------------------------
    s = Slide("future"); S.append(s)
    s.title("Next steps")
    s.image(f"{SRC}/ppt/media/image2.jpg", 0.5, 1.15, 4.2, 5.2, crop=(0.068, 0.10, 0.068, 0.10))
    s.card(5.0, 1.15, 7.9, 1.6, "Field testing",
           "Diurnal stomatal conductance of LT candidates under field temperature and humidity, Arkansas 2026-2027")
    s.card(5.0, 2.95, 7.9, 1.6, "Breeding donors",
           "PI398617 and PI377574 as sources of an early, strong LT response for MG IV backgrounds")
    s.card(5.0, 4.75, 7.9, 1.6, "Model the threshold",
           "Genotype-specific breakpoints in SSM to map where, and in which seasons, the LT trait pays off", fill=LIGHT, head_color=DARK)

    return S


# ====================================================================== assembly
def main():
    if os.path.exists(WORK):
        shutil.rmtree(WORK)
    shutil.copytree(SRC, WORK)
    # remove stray converted previews in media (png copies of tiff/emf we made while inspecting)
    for f in os.listdir(f"{WORK}/ppt/media"):
        base, ext = os.path.splitext(f)
        if ext == ".png" and (os.path.exists(f"{WORK}/ppt/media/{base}.tiff") or os.path.exists(f"{WORK}/ppt/media/{base}.emf")):
            os.remove(f"{WORK}/ppt/media/{f}")

    slides = build_slides()

    # 1) structural: duplicate slide11 once per new slide (appended at end)
    created = []
    for s in slides:
        out = subprocess.run([sys.executable, f"{SK}/scripts/add_slide.py", WORK, "slide11.xml"],
                             capture_output=True, text=True, check=True).stdout
        m = re.search(r"(ppt/slides/slide\d+\.xml)", out)
        created.append(m.group(1))
        print(s.name, "->", m.group(1))

    # 2) write content + rels + media
    media_counter = 100
    for s, path in zip(slides, created):
        names = []
        for rid, src in s.images:
            media_counter += 1
            ext = os.path.splitext(src)[1].lower()
            mn = f"image{media_counter}{ext}"
            shutil.copy(src, f"{WORK}/ppt/media/{mn}")
            names.append(mn)
        with open(f"{WORK}/{path}", "w") as f:
            f.write(s.xml())
        relp = f"{WORK}/ppt/slides/_rels/{os.path.basename(path)}.rels"
        with open(relp, "w") as f:
            f.write(s.rels(names))

    # 3) edit kept template slides
    def sub(path, old, new, count=1):
        p = f"{WORK}/{path}"
        x = open(p).read()
        assert old in x, (path, old)
        x = x.replace(old, new, count)
        open(p, "w").write(x)

    NEW_TITLE = "Impact of the limited transpiration trait on soybean potential and water-limited yield"
    OLD_TITLE = "Modeling crop physiological strategies to alleviate soybean drought stress in the US-Midsouth"
    for sl in ("ppt/slides/slide1.xml", "ppt/slides/slide17.xml"):
        sub(sl, f"<a:t>{OLD_TITLE}</a:t>", f"<a:t>{NEW_TITLE}</a:t>")
    # co-author line on the title slide (after "Elvis Elli"), with a taller author box
    coauth = ('<a:p><a:r><a:rPr lang="en-US" sz="1500" dirty="0"><a:solidFill><a:srgbClr val="A51E36"/></a:solidFill>' + CAL + '</a:rPr>'
              '<a:t>Lucas Romano, Martina Meira, Felix Fritschi, Afsin Soltani, Caio Canella Vieira</a:t></a:r></a:p>')
    x = open(f"{WORK}/ppt/slides/slide1.xml").read()
    i = x.find("<a:t>Elvis Elli</a:t>")
    j = x.find("</a:p>", i) + len("</a:p>")
    x = x[:j] + coauth + x[j:]
    x = x.replace('<a:off x="136776" y="3429000"/><a:ext cx="6877110" cy="1698169"/>',
                  '<a:off x="136776" y="3380000"/><a:ext cx="7100000" cy="2250000"/>')
    open(f"{WORK}/ppt/slides/slide1.xml", "w").write(x)
    # thank-you slide: ASA 2026 line instead of nothing? keep as is.

    # 4) final order
    order = (["ppt/slides/slide1.xml", "ppt/slides/slide2.xml", created[0],
              "ppt/slides/slide3.xml", "ppt/slides/slide4.xml"]
             + created[1:]
             + ["ppt/slides/slide15.xml", "ppt/slides/slide17.xml"])
    pres = f"{WORK}/ppt/presentation.xml"
    rels = f"{WORK}/ppt/_rels/presentation.xml.rels"
    rx = open(rels).read()
    px = open(pres).read()
    target2rid = {m.group(2): m.group(1) for m in re.finditer(r'<Relationship Id="(rId\d+)" Type="[^"]*/slide" Target="(slides/slide\d+\.xml)"/>', rx)}
    ids = re.findall(r'<p:sldId id="(\d+)" r:id="(rId\d+)"/>', px)
    rid2id = {r: i for i, r in ids}
    nid = max(int(i) for i, _ in ids) + 1
    new_lst = []
    for p in order:
        rid = target2rid[p.replace("ppt/", "")]
        sid = rid2id.get(rid)
        if sid is None:
            sid = str(nid); nid += 1
        new_lst.append(f'<p:sldId id="{sid}" r:id="{rid}"/>')
    px = re.sub(r'<p:sldIdLst>.*?</p:sldIdLst>', '<p:sldIdLst>' + "".join(new_lst) + '</p:sldIdLst>', px, flags=re.S)
    open(pres, "w").write(px)

    # 5) clean orphans, zip
    subprocess.run([sys.executable, f"{SK}/scripts/clean.py", WORK], check=True)
    if os.path.exists(OUT):
        os.remove(OUT)
    subprocess.run(["zip", "-Xr", "-q", OUT, "."], cwd=WORK, check=True)
    print("wrote", OUT, os.path.getsize(OUT) // 1024, "KB")


if __name__ == "__main__":
    main()
