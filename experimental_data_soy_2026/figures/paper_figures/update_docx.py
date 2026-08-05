"""
update_docx.py
Insert 3 key-result bullet points after each figure caption in the paper
Word document.  Bullets are read from figure_bullets.json, which is produced
by figure_bullets.R.

Run:
    python3 update_docx.py
"""
import json
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

BASE   = Path(__file__).parent
SRC    = Path("/root/.claude/uploads/df313029-2814-5867-990d-e761f419a647/06ebe891-paper0628.docx")
OUT    = BASE / "paper_figures_annotated.docx"
JSON   = BASE / "figure_bullets.json"

# Paragraph IDs of the five figure captions (from word/document.xml)
PARA_IDS = {
    "F1": "1840C233",   # Figure 1 — TR response
    "F2": "0B6914F2",   # Figure 2 — breakpoints
    "F3": "2D77048E",   # Figure 3 — intercept/slope quadrant
    "F4": "25BD8838",   # Figure 4 — gsw response
    "F5": "375D88BA",   # Figure 5 — predicted bars
}

# --------------------------------------------------------------------------- #

def bullet_xml(text: str) -> str:
    """Return a Word paragraph XML node formatted as a bullet."""
    for ch, esc in [("&", "&amp;"), ("<", "&lt;"), (">", "&gt;"), ('"', "&quot;")]:
        text = text.replace(ch, esc)
    return (
        '<w:p>'
          '<w:pPr>'
            '<w:pStyle w:val="ListParagraph"/>'
            '<w:ind w:left="720" w:hanging="360"/>'
          '</w:pPr>'
          '<w:r>'
            '<w:rPr><w:sz w:val="22"/></w:rPr>'
            f'<w:t xml:space="preserve">•  {text}</w:t>'
          '</w:r>'
        '</w:p>'
    )


def insert_after_para(xml: str, para_id: str, bullets: list[str]) -> str:
    """Find the paragraph with para_id and insert bullet paragraphs after it."""
    anchor = f'w14:paraId="{para_id}"'
    pos = xml.find(anchor)
    if pos == -1:
        print(f"  WARNING: paraId {para_id} not found")
        return xml
    close_tag = "</w:p>"
    after = xml.find(close_tag, pos)
    if after == -1:
        print(f"  WARNING: closing </w:p> not found after {para_id}")
        return xml
    insert_at = after + len(close_tag)
    block = "\n" + "\n".join(bullet_xml(b) for b in bullets) + "\n"
    return xml[:insert_at] + block + xml[insert_at:]


# --------------------------------------------------------------------------- #

print("Reading bullets from:", JSON)
with open(JSON, encoding="utf-8") as f:
    bullets_by_fig = json.load(f)

tmp = Path(tempfile.mkdtemp())
shutil.copy2(SRC, tmp / "source.docx")
subprocess.run(["unzip", "-q", str(tmp / "source.docx"), "-d", str(tmp / "unpacked")], check=True)

xml_path = tmp / "unpacked" / "word" / "document.xml"
xml = xml_path.read_text(encoding="utf-8")

for fig_key, para_id in PARA_IDS.items():
    blist = bullets_by_fig.get(fig_key, [])
    xml = insert_after_para(xml, para_id, blist)
    print(f"  Inserted {len(blist)} bullets for {fig_key}")

xml_path.write_text(xml, encoding="utf-8")

out_tmp = tmp / "output.docx"
subprocess.run(
    ["zip", "-qr", str(out_tmp), "."],
    cwd=str(tmp / "unpacked"),
    check=True
)
shutil.copy2(out_tmp, OUT)
shutil.rmtree(tmp)

print("Saved:", OUT)
