"""
Venn diagram — Figma-ready SVG.
Z-Hunter (green) ∩ Z-DNABERT (blue), black labels.
"""

import svgwrite
from pathlib import Path

DATA_DIR = Path("/Users/dassagaripova/Downloads/project")

# ── Data ──────────────────────────────────────────────────────
ZHUNT_N      = "1,781,183"
ZDNABERT_N   = "1,588,238"
ZHUNT_ONLY   = "985,176"
INTERSECT    = "681,135"
ZDNABERT_ONLY= "1,021,974"

# ── Geometry ──────────────────────────────────────────────────
W, H   = 900, 700
R      = 210          # circle radius
CY     = H // 2 + 20  # circle centres y
CX_L   = W // 2 - 110 # Z-Hunter centre x
CX_R   = W // 2 + 110 # Z-DNABERT centre x

dwg = svgwrite.Drawing(
    str(DATA_DIR / "venn_figma.svg"),
    size=(f"{W}px", f"{H}px"),
    profile="full"
)
dwg.attribs["font-family"] = "Arial, Helvetica, sans-serif"
dwg.add(dwg.rect(insert=(0,0), size=("100%","100%"), fill="white"))

# ── Defs: clip paths for intersection fill ────────────────────
defs = dwg.defs

# clip to left circle
clip_l = defs.add(dwg.clipPath(id="clip-left"))
clip_l.add(dwg.circle(center=(CX_L, CY), r=R))

# clip to right circle
clip_r = defs.add(dwg.clipPath(id="clip-right"))
clip_r.add(dwg.circle(center=(CX_R, CY), r=R))

# ── Left circle fill (Z-Hunter green) ────────────────────────
dwg.add(dwg.circle(
    center=(CX_L, CY), r=R,
    fill="#52B788", fill_opacity=0.55,
    stroke="#238B45", stroke_width=2.5
))

# ── Right circle fill (Z-DNABERT blue) ───────────────────────
dwg.add(dwg.circle(
    center=(CX_R, CY), r=R,
    fill="#74A9D8", fill_opacity=0.45,
    stroke="#2171B5", stroke_width=2.5
))

# ── Intersection highlight (cyan-ish) ────────────────────────
# Draw right circle clipped to left → gives intersection area
inter_grp = dwg.add(dwg.g(clip_path="url(#clip-left)"))
inter_grp.add(dwg.circle(
    center=(CX_R, CY), r=R,
    fill="#7ECECA", fill_opacity=0.70,
    stroke="none"
))

# ── Labels: titles ────────────────────────────────────────────
# Z-Hunter title (green)
dwg.add(dwg.text("Z-Hunter",
    insert=(CX_L, 55),
    text_anchor="middle",
    font_size="32px", font_weight="bold", fill="#238B45"
))
# n =
dwg.add(dwg.text(f"n = {ZHUNT_N}",
    insert=(CX_L, 90),
    text_anchor="middle",
    font_size="16px", fill="#333"
))

# Z-DNABERT title (blue)
dwg.add(dwg.text("Z-DNABERT",
    insert=(CX_R, 55),
    text_anchor="middle",
    font_size="32px", font_weight="bold", fill="#2171B5"
))
dwg.add(dwg.text(f"n = {ZDNABERT_N}",
    insert=(CX_R, 90),
    text_anchor="middle",
    font_size="16px", fill="#333"
))

# ── Numbers inside circles ─────────────────────────────────────
# Z-Hunter only (left)
dwg.add(dwg.text(ZHUNT_ONLY,
    insert=(CX_L - 80, CY + 8),
    text_anchor="middle",
    font_size="26px", font_weight="bold", fill="#111"
))

# Intersection (centre)
ix = (CX_L + CX_R) // 2
dwg.add(dwg.text(INTERSECT,
    insert=(ix, CY + 8),
    text_anchor="middle",
    font_size="26px", font_weight="bold", fill="#111"
))

# Z-DNABERT only (right)
dwg.add(dwg.text(ZDNABERT_ONLY,
    insert=(CX_R + 80, CY + 8),
    text_anchor="middle",
    font_size="26px", font_weight="bold", fill="#111"
))

dwg.save()
print(f"✅ {DATA_DIR}/venn_figma.svg")
