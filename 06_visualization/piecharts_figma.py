"""
Three pie charts — Figma-ready SVG.
Style: clean callout labels with leader lines, Arial font.
Same color palette as the chromosome ideogram figure.
"""

import numpy as np
import math
import svgwrite
from pathlib import Path

DATA_DIR = Path("/Users/dassagaripova/Downloads/project")

# ── Colors (same as figure) ───────────────────────────────────
COLORS = {
    "Intergenic":  "#A5A5A5",
    "Introns":     "#4472C4",
    "Promoters":   "#FFC000",
    "Exons":       "#FF0000",
    "UTR":         "#70AD47",
    "CDS":         "#ED7D31",
    "Downstream":  "#7030A0",
}
# Order matches reference: small slices grouped top-right and bottom-left,
# large slices (Introns, Intergenic) fill right and left halves
CAT_ORDER = ["Promoters","Exons","Introns","Intergenic","Downstream","CDS","UTR"]

# ── Data ──────────────────────────────────────────────────────
PIES = [
    {
        "letter": "B",
        "title":  "Genome",
        "subtitle": "",
        "data": {
            "Intergenic": 54.50, "Introns": 40.63, "Promoters": 2.06,
            "Exons": 2.27, "UTR": 1.12, "CDS": 1.14, "Downstream": 1.76,
        }
    },
    {
        "letter": "C",
        "title":  "Z-DNABERT",
        "subtitle": "N = 1,588,238",
        "data": {
            "Intergenic": 61.26, "Introns": 29.27, "Promoters": 1.97,
            "Exons": 0.47, "UTR": 0.40, "CDS": 0.07, "Downstream": 1.29,
        }
    },
    {
        "letter": "D",
        "title":  "Z-Hunter",
        "subtitle": "N = 1,781,183",
        "data": {
            "Intergenic": 58.84, "Introns": 38.09, "Promoters": 1.93,
            "Exons": 0.84, "UTR": 0.82, "CDS": 0.02, "Downstream": 1.24,
        }
    },
]

# ── SVG geometry ──────────────────────────────────────────────
PIE_R    = 100        # pie radius px
CX_GAP   = 420        # horizontal distance between pie centres
MARGIN_X = 220        # left margin
MARGIN_Y = 160        # top margin (room for title + letter)
SVG_W    = MARGIN_X * 2 + CX_GAP * 2
SVG_H    = 420

dwg = svgwrite.Drawing(
    str(DATA_DIR / "piecharts_figma.svg"),
    size=(f"{SVG_W}px", f"{SVG_H}px"),
    profile="full"
)
dwg.attribs["font-family"] = "Arial, Helvetica, sans-serif"
dwg.add(dwg.rect(insert=(0,0), size=("100%","100%"), fill="white"))

# ── Helper: arc path ──────────────────────────────────────────
def arc_path(cx, cy, r, a_start_deg, a_end_deg):
    """SVG arc path for a pie wedge."""
    a0 = math.radians(a_start_deg)
    a1 = math.radians(a_end_deg)
    x0 = cx + r * math.cos(a0)
    y0 = cy + r * math.sin(a0)
    x1 = cx + r * math.cos(a1)
    y1 = cy + r * math.sin(a1)
    large = 1 if (a_end_deg - a_start_deg) > 180 else 0
    return f"M {cx} {cy} L {x0:.3f} {y0:.3f} A {r} {r} 0 {large} 1 {x1:.3f} {y1:.3f} Z"

# ── Helper: push labels apart ─────────────────────────────────
def push_apart(ys, gap=16, lo=-95, hi=95):
    ys = list(ys)
    n = len(ys)
    for _ in range(800):
        moved = False
        for i in range(n-1):
            d = ys[i] - ys[i+1]
            if d < gap:
                push = (gap - d) / 2
                ys[i]   += push
                ys[i+1] -= push
                moved = True
        if not moved:
            break
    if ys[0]  > hi: delta = ys[0]-hi;   ys=[y-delta for y in ys]
    if ys[-1] < lo: delta = lo-ys[-1];  ys=[y+delta for y in ys]
    return ys

# ── Draw each pie ─────────────────────────────────────────────
for pi, pie in enumerate(PIES):
    cx = MARGIN_X + pi * CX_GAP
    cy = MARGIN_Y

    data  = pie["data"]
    total = sum(data[k] for k in CAT_ORDER)
    pcts  = [data[k] / total * 100 for k in CAT_ORDER]

    grp = dwg.add(dwg.g(id=f"pie-{pie['letter']}"))

    # ── Panel letter ─────────────────────────────────────────
    grp.add(dwg.text(
        pie["letter"],
        insert=(cx - PIE_R - 80, cy - PIE_R - 30),
        font_size="16px", font_weight="bold", fill="#111"
    ))

    # ── Title ─────────────────────────────────────────────────
    grp.add(dwg.text(
        pie["title"],
        insert=(cx, cy - PIE_R - 30),
        text_anchor="middle",
        font_size="13px", font_weight="bold", fill="#111"
    ))
    if pie["subtitle"]:
        grp.add(dwg.text(
            pie["subtitle"],
            insert=(cx, cy - PIE_R - 14),
            text_anchor="middle",
            font_size="12px", font_weight="bold", fill="#111"
        ))

    # ── Wedges ────────────────────────────────────────────────
    # Start at 90° (top), go clockwise → SVG y-axis is flipped,
    # so clockwise = increasing angle in SVG coords
    angle = -90.0   # start at top
    wedge_angles = []   # (cat, pct, mid_angle)
    wedges_grp = grp.add(dwg.g(id=f"wedges-{pie['letter']}"))

    for cat, pct in zip(CAT_ORDER, pcts):
        sweep    = pct / 100 * 360
        a_start  = angle
        a_end    = angle + sweep
        a_mid    = (a_start + a_end) / 2
        color    = COLORS[cat]
        d        = arc_path(cx, cy, PIE_R, a_start, a_end)

        wedges_grp.add(dwg.path(
            d=d, fill=color,
            stroke="white", stroke_width=0.8,
            id=f"wedge-{pie['letter']}-{cat}"
        ))
        wedge_angles.append((cat, pct, a_mid))
        angle = a_end

    # ── Callout labels ────────────────────────────────────────
    R_EDGE  = PIE_R * 1.05
    R_KINK  = PIE_R * 1.25
    X_R     = cx + PIE_R * 1.78
    X_L     = cx - PIE_R * 1.78
    FS      = "9px"

    labels_grp = grp.add(dwg.g(id=f"labels-{pie['letter']}"))

    right = [(cat, pct, a) for cat, pct, a in wedge_angles
             if math.cos(math.radians(a)) >= 0]
    left  = [(cat, pct, a) for cat, pct, a in wedge_angles
             if math.cos(math.radians(a)) <  0]

    right = sorted(right, key=lambda e: -math.sin(math.radians(e[2])))
    left  = sorted(left,  key=lambda e: -math.sin(math.radians(e[2])))

    for group, x_txt, anchor in [(right, X_R, "start"), (left, X_L, "end")]:
        if not group:
            continue
        ys_init = [cy + math.sin(math.radians(a)) * PIE_R for _,_,a in group]
        ys      = push_apart(ys_init, gap=15, lo=cy-PIE_R*1.45, hi=cy+PIE_R*1.45)

        for (cat, pct, ang_deg), y_lbl in zip(group, ys):
            if pct < 0.01:
                continue
            rad  = math.radians(ang_deg)
            x0   = cx + R_EDGE * math.cos(rad)
            y0   = cy + R_EDGE * math.sin(rad)
            x1   = cx + R_KINK * math.cos(rad)
            y1   = cy + R_KINK * math.sin(rad)
            x2   = x_txt * 0.965 + cx * 0.035  # slightly inward of anchor

            # leader line: edge → kink → text
            labels_grp.add(dwg.polyline(
                points=[(x0,y0),(x1,y1),(x2,y_lbl)],
                fill="none", stroke="#666",
                stroke_width=0.7,
                stroke_linejoin="round"
            ))

            txt = f"{cat}, {pct:.2f}%"
            labels_grp.add(dwg.text(
                txt,
                insert=(x_txt, y_lbl + 3.5),
                text_anchor=anchor,
                font_size=FS, fill="#111"
            ))

# ── Shared legend at bottom ───────────────────────────────────
leg_y   = SVG_H - 40
leg_grp = dwg.add(dwg.g(id="legend"))
n_cats  = len(CAT_ORDER)
box_w   = 12
spacing = SVG_W / (n_cats + 1)

for i, cat in enumerate(CAT_ORDER):
    x = spacing * (i + 1) - box_w / 2
    leg_grp.add(dwg.rect(insert=(x, leg_y), size=(box_w, box_w),
                         fill=COLORS[cat], rx=2))
    leg_grp.add(dwg.text(cat,
                          insert=(x + box_w + 4, leg_y + 9),
                          font_size="9px", fill="#222"))

dwg.save()
print(f"✅ SVG: {DATA_DIR}/piecharts_figma.svg")
print("   → Перетащи в Figma — каждый пирог отдельная группа (pie-B, pie-C, pie-D)")
