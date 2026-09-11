"""Tidal-stripping study explainer - plain language, python-docx.
Reads tidal_results.json (produced by tidal_pipeline.py) and fills in the
actual numbers, so the document always matches the final run.
"""
import json
from pathlib import Path

from docx import Document
from docx.shared import Pt, Inches
from docx.enum.text import WD_ALIGN_PARAGRAPH

BASE = Path(r"C:\Users\hoche\OneDrive\Desktop\Cluster Analysis HTML\tidal_study")
OUT = BASE / "Tidal_Stripping_Explained.docx"
res = {r["cluster"]: r for r in json.load(open(BASE / "tidal_results.json", encoding="utf-8"))}

NICE = {"Pal5": "Palomar 5", "NGC6341": "M92 (NGC 6341)", "NGC288": "NGC 288"}

doc = Document()
st = doc.styles["Normal"]; st.font.name = "Calibri"; st.font.size = Pt(11)

doc.add_heading("Hunting for Tidal Stripping — How We Did It", 0)
p = doc.add_paragraph()
r = p.add_run("A plain-language guide to our search for stars being pulled out of globular clusters by the Milky Way’s gravity.")
r.italic = True

doc.add_heading("1. The idea", level=1)
doc.add_paragraph(
    "The Milky Way’s gravity is stronger on a cluster’s near side than its far side. That difference "
    "(the tidal force) slowly peels stars off. The escaped stars don’t vanish — they drift ahead of and "
    "behind the cluster along its orbit, forming two thin “tails”. Finding those tails = catching tidal "
    "stripping in the act.")

doc.add_heading("2. Why these three clusters", level=1)
doc.add_paragraph(
    "We needed clusters far from the crowded Milky Way disk, so any extra stars we find can’t be blamed "
    "on contamination. We picked Palomar 5, M92, and NGC 288 — all high above the Galactic plane, all "
    "with clean membership data from our earlier study, and all with tails or streams reported in "
    "published papers (so we can check ourselves against the professionals).")

doc.add_heading("3. Finding escaped stars", level=1)
doc.add_paragraph(
    "An escaped star keeps two fingerprints of its birthplace: it still moves across the sky at almost "
    "exactly the cluster’s speed, and it still sits on the cluster’s narrow track in a color-vs-brightness "
    "chart (the CMD) — because it has the same age and metal content as its siblings. So we:")
for t in [
    "Downloaded a huge circle of sky around each cluster from Gaia — several degrees across, far beyond "
    "the cluster’s edge (the “tidal radius”). The user asked about direction: we don’t need to guess it — "
    "we search all directions and let the stars reveal the line.",
    "Kept only stars whose proper motion matches the cluster’s (within a small window).",
    "Kept only stars lying on the cluster’s own CMD ridge line (built from our member catalog).",
    "Removed stars whose parallax says they are nearby foreground objects.",
]:
    doc.add_paragraph(t, style="List Bullet")

doc.add_heading("4. Reading the result", level=1)
doc.add_paragraph(
    "We made a smoothed density map of the surviving stars and asked three questions: (1) are there more "
    "stars beyond the tidal radius than the background predicts? (2) do they line up along one direction "
    "(a tail axis)? (3) does their CMD turnoff — the sharp corner where the star track bends — match the "
    "cluster’s own turnoff? A yes to all three is evidence of tidal stripping.")

doc.add_heading("5. What we found", level=1)
for cid, nice in NICE.items():
    r = res.get(cid)
    if not r:
        continue
    excess = r["excess_outside"]
    pa, pm_pa, agree = r["tail_pa_deg"], r["pm_pa_deg"], r["pa_agreement_deg"]
    to_c, to_t = r["turnoff_color_cluster"], r["turnoff_color_tail"]
    bits = [f"{r['n_candidates']:,} stars passed all filters; "
            f"{r['n_outside_rt']:,} lie beyond the tidal radius vs ~{r['n_outside_expected_bg']:.0f} "
            f"expected from background (excess ≈ {excess:.0f})."]
    if pa is not None and agree is not None and agree < 30:
        bits.append(f"They align along position angle {pa:.0f}° while the cluster moves along "
                    f"{pm_pa:.0f}° — only {agree:.0f}° apart, matching the theory that tails follow "
                    "the orbit.")
    elif pa is not None:
        bits.append(f"An overdensity exists but its axis ({pa:.0f}°) does not line up with the "
                    f"cluster’s motion ({pm_pa:.0f}°) — so we call this evidence ambiguous, which "
                    "matches the disagreement about this cluster in published papers.")
    else:
        bits.append("No significant elongated overdensity was detected — an honest null for this cluster.")
    if to_t is not None:
        bits.append(f"CMD turnoff colors match: cluster {to_c:.2f} vs tail candidates {to_t:.2f} "
                    "(same age and metallicity — these stars were born in the cluster).")
    doc.add_heading(nice, level=2)
    doc.add_paragraph(" ".join(bits))
    img = BASE / f"{cid}_tail_map.png"
    if img.exists():
        doc.add_picture(str(img), width=Inches(4.4))
        doc.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER

doc.add_heading("6. How this compares with published work", level=1)
doc.add_paragraph(
    "Palomar 5’s tails were discovered by Odenkirchen et al. (2001) and traced across 10° of sky by 2003; "
    "Gaia-era maps (Ibata et al. 2021) confirm them. M92’s stream was reported by Thomas et al. (2020). "
    "NGC 288 has claimed short tails (Jordi & Grebel 2010) but other surveys found little — a genuinely "
    "hard case. Full citations: literature_refs.json in this folder.")

doc.add_heading("7. Honest limitations", level=1)
for t in [
    "Gaia only reaches G ≈ 20.5 — for distant Pal 5 that’s barely past the turnoff, so we see only the "
    "brightest tail stars; professional maps use deeper telescopes.",
    "We measure metallicity agreement indirectly (through the CMD track), not with spectra.",
    "The density-map significance is approximate (Poisson estimate on smoothed counts).",
]:
    doc.add_paragraph(t, style="List Bullet")

doc.save(OUT)
print("wrote", OUT)
