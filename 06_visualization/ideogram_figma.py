"""
Generates a clean, Figma-ready SVG of the Z-DNA ideogram.
Structure:
  <g id="chromosome-1"> ... <g id="zdnabert-1"/><g id="zhunter-1"/> </g>
  ...
Each chromosome = named group → fully editable in Figma.
"""

import numpy as np
import pandas as pd
import svgwrite
from svgwrite import cm, mm
from pathlib import Path

DATA_DIR = Path("/Users/dassagaripova/Downloads/project")

# ── Karyotype ─────────────────────────────────────────────────
CHR_LENGTHS = [
    199866856,192480835,168055261,147789705,127475568,
    117078234,110520290, 97793126, 96880575, 94596256,
     80862590, 68457376, 65055609, 60575670, 58152616,
     57488271, 57420321, 55529144, 55115435, 47890911,
     40058617, 37711260, 37670884, 36152331, 35487451,
     26644192, 22551867, 18633596, 11972102, 10081836,
]
N_CHR     = 30
CHR_NAMES = [f"Chr{i}" for i in range(1, N_CHR+1)]
MAX_LEN   = max(CHR_LENGTHS)
NC_IDS    = ([f"NC_06{n}.1" for n in range(8981,9000)] +
             [f"NC_06{n}.1" for n in range(9000,9011)])
CHR_MAP   = dict(zip(NC_IDS, CHR_NAMES))

# ── Read bedgraph ─────────────────────────────────────────────
def read_bg(path):
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith(("track","browser","#")) or not line.strip(): continue
            p = line.split()
            if len(p) < 4 or p[0] not in CHR_MAP: continue
            try: rows.append((CHR_MAP[p[0]], int(p[1]), int(p[2]), float(p[3])))
            except: continue
    return pd.DataFrame(rows, columns=["chr","start","end","score"])

print("Reading bedgraphs...")
bert  = read_bg(DATA_DIR/"GCF_001194135.2_ASM119413v2_zdna_thr025.bedgraph.gz")
zhunt = read_bg(DATA_DIR/"zhunter_z-dna.bedgraph")

# ── Density 1 Mb ──────────────────────────────────────────────
WINDOW = 1_000_000
def density(df):
    res = {}
    for name, length in zip(CHR_NAMES, CHR_LENGTHS):
        sub = df[df["chr"]==name]
        brk = np.arange(0, length+WINDOW, WINDOW)
        vals = []
        for s,e in zip(brk[:-1], brk[1:]):
            e = min(e, length)
            ov = (np.minimum(sub["end"],e)-np.maximum(sub["start"],s)).clip(lower=0).sum() if len(sub) else 0
            vals.append(ov/(e-s)*100)
        res[name] = np.array(vals)
    return res

print("Computing densities...")
bd = density(bert);  vb = np.percentile(np.concatenate(list(bd.values())), 99)
zd = density(zhunt); vz = np.percentile(np.concatenate(list(zd.values())), 99)

# ── Color interpolation ───────────────────────────────────────
def lerp_color(t, stops):
    """Interpolate hex color from list of (t, r, g, b) stops."""
    t = max(0.0, min(1.0, t))
    for i in range(len(stops)-1):
        t0, r0, g0, b0 = stops[i]
        t1, r1, g1, b1 = stops[i+1]
        if t0 <= t <= t1:
            f = (t-t0)/(t1-t0)
            r = int(r0 + f*(r1-r0))
            g = int(g0 + f*(g1-g0))
            b = int(b0 + f*(b1-b0))
            return f"#{r:02X}{g:02X}{b:02X}"
    return f"#{int(stops[-1][1]):02X}{int(stops[-1][2]):02X}{int(stops[-1][3]):02X}"

# white→blue  (Z-DNABERT)
BERT_STOPS = [
    (0.0,  255,255,255),
    (0.25, 198,219,239),
    (0.5,  107,174,214),
    (0.75,  33,113,181),
    (1.0,    8, 48,107),
]
# white→green (Z-Hunter)
ZHUNT_STOPS = [
    (0.0,  255,255,255),
    (0.25, 199,233,192),
    (0.5,  116,196,118),
    (0.75,  35,139, 69),
    (1.0,    0, 68, 27),
]

# ── SVG layout ────────────────────────────────────────────────
# Each chromosome column = COL_PX wide; height = MAX_LEN scaled to PAGE_H
PAGE_H   = 900      # px total height
PAGE_PAD = 50       # top/bottom padding
CHR_AREA = PAGE_H - 2*PAGE_PAD
SCALE    = CHR_AREA / MAX_LEN  # px per bp

COL_W    = 22       # total column width per chromosome (px)
CHR_W    = 14       # chromosome body width (px) — split in two halves
HALF_W   = CHR_W / 2  # each track occupies exactly half
GAP      = 0        # no gap
COL_SEP  = COL_W   # x step between chromosomes

TOTAL_W  = N_CHR * COL_SEP + 120   # + right margin for colorbar

dwg = svgwrite.Drawing(
    str(DATA_DIR / "octopus_ideogram_figma.svg"),
    size=(f"{TOTAL_W}px", f"{PAGE_H + 80}px"),
    profile="full"
)
dwg.attribs["font-family"] = "Arial, Helvetica, sans-serif"

# White background
dwg.add(dwg.rect(insert=(0,0), size=("100%","100%"), fill="white"))

# ── Title ─────────────────────────────────────────────────────
dwg.add(dwg.text(
    "Z-DNA distribution across Octopus bimaculoides chromosomes",
    insert=(TOTAL_W/2, 22),
    text_anchor="middle", font_size="11px", font_weight="bold", fill="#111"
))

# ── Scale bar 100 Mb ──────────────────────────────────────────
sb_x   = 8
sb_y0  = PAGE_PAD
sb_y1  = PAGE_PAD + 100e6 * SCALE
sb_grp = dwg.add(dwg.g(id="scale-bar"))
sb_grp.add(dwg.line(start=(sb_x, sb_y0), end=(sb_x, sb_y1),
                    stroke="#333", stroke_width=1.2))
for yy in [sb_y0, sb_y1]:
    sb_grp.add(dwg.line(start=(sb_x-3,yy), end=(sb_x+3,yy),
                        stroke="#333", stroke_width=0.8))
sb_grp.add(dwg.text("100 Mb", insert=(sb_x-4, (sb_y0+sb_y1)/2),
                    text_anchor="end", font_size="7px", fill="#444",
                    transform=f"rotate(-90,{sb_x-4},{(sb_y0+sb_y1)/2})"))

print("Drawing chromosomes...")
for idx, (name, length) in enumerate(zip(CHR_NAMES, CHR_LENGTHS)):
    x_center = 35 + idx * COL_SEP + COL_SEP / 2
    chr_h    = length * SCALE
    y_bot    = PAGE_PAD + CHR_AREA - chr_h    # top-aligned (tallest at top)
    y_top    = PAGE_PAD

    # chromosome group
    grp = dwg.add(dwg.g(id=f"chromosome-{idx+1}"))

    nw   = len(bd[name])
    wh   = chr_h / nw     # window height in px

    x_left  = x_center - HALF_W   # left edge of chromosome
    x_mid   = x_center            # divider (centre)
    x_right = x_center + HALF_W  # right edge

    # ── Z-DNABERT — LEFT half (blue) ──────────────────────────
    bert_grp = grp.add(dwg.g(id=f"zdnabert-{idx+1}"))
    for wi, val in enumerate(bd[name]):
        t     = min(val / vb, 1.0)
        color = lerp_color(t, BERT_STOPS)
        yw    = y_bot + wi * wh
        bert_grp.add(dwg.rect(
            insert=(x_left, yw),
            size=(HALF_W, wh + 0.5),
            fill=color, stroke="none"
        ))

    # ── Z-Hunter — RIGHT half (green) ─────────────────────────
    zhunt_grp = grp.add(dwg.g(id=f"zhunter-{idx+1}"))
    for wi, val in enumerate(zd[name]):
        t     = min(val / vz, 1.0)
        color = lerp_color(t, ZHUNT_STOPS)
        yw    = y_bot + wi * wh
        zhunt_grp.add(dwg.rect(
            insert=(x_mid, yw),
            size=(HALF_W, wh + 0.5),
            fill=color, stroke="none"
        ))

    # ── Chromosome outline (on top, rounded border) ───────────
    body_grp = grp.add(dwg.g(id=f"body-{idx+1}"))
    rx = HALF_W * 0.55
    body_grp.add(dwg.rect(
        insert=(x_left, y_bot),
        size=(CHR_W, chr_h),
        rx=rx, ry=rx,
        fill="none",
        stroke="#2a5080", stroke_width=0.6
    ))

    # ── Chromosome number BOTTOM ───────────────────────────────
    grp.add(dwg.text(
        str(idx+1),
        insert=(x_center, PAGE_PAD + CHR_AREA + 12),
        text_anchor="middle",
        font_size="6.5px", font_weight="bold", fill="#111"
    ))

    # ── Length label below number ──────────────────────────────
    grp.add(dwg.text(
        f"{length/1e6:.0f}",
        insert=(x_center, PAGE_PAD + CHR_AREA + 22),
        text_anchor="middle",
        font_size="5px", fill="#888"
    ))

# ── Colorbars ─────────────────────────────────────────────────
cb_x     = TOTAL_W - 90
cb_h     = 160
cb_w     = 10

for stops, label, y0 in [
    (BERT_STOPS,  "Z-DNABERT\n(% bp / 1 Mb)", PAGE_PAD + 10),
    (ZHUNT_STOPS, "Z-Hunter\n(% bp / 1 Mb)",  PAGE_PAD + 200),
]:
    cb_grp = dwg.add(dwg.g(id=f"colorbar-{label[:6]}"))
    n_steps = 50
    for i in range(n_steps):
        t  = 1.0 - i / n_steps
        yy = y0 + i * (cb_h / n_steps)
        color = lerp_color(t, stops)
        cb_grp.add(dwg.rect(
            insert=(cb_x, yy),
            size=(cb_w, cb_h/n_steps + 0.5),
            fill=color, stroke="none"
        ))
    # border
    cb_grp.add(dwg.rect(
        insert=(cb_x, y0), size=(cb_w, cb_h),
        fill="none", stroke="#888", stroke_width=0.4
    ))
    cb_grp.add(dwg.text("High", insert=(cb_x+cb_w+3, y0+5),
                         font_size="7px", fill="#555"))
    cb_grp.add(dwg.text("Low",  insert=(cb_x+cb_w+3, y0+cb_h),
                         font_size="7px", fill="#555"))
    for line in label.split("\n"):
        cb_grp.add(dwg.text(line,
                             insert=(cb_x + cb_w/2, y0 + cb_h + 14 + label.split("\n").index(line)*10),
                             text_anchor="middle", font_size="7px", fill="#333"))

# ── Legend ────────────────────────────────────────────────────
leg_x, leg_y = TOTAL_W - 90, PAGE_PAD + 380
leg = dwg.add(dwg.g(id="legend"))
leg.add(dwg.rect(insert=(leg_x-4,leg_y-10), size=(88,52),
                 fill="white", stroke="#ccc", stroke_width=0.5, rx=3))
leg.add(dwg.text("Z-DNA predictor", insert=(leg_x, leg_y),
                  font_size="7.5px", font_weight="bold", fill="#222"))
for i,(color,lbl) in enumerate([("#2171B5","Z-DNABERT (left)"),
                                  ("#238B45","Z-Hunter (right)")]):
    yy = leg_y + 14 + i*14
    leg.add(dwg.rect(insert=(leg_x,yy-7), size=(10,8),
                     fill=color, rx=1))
    leg.add(dwg.text(lbl, insert=(leg_x+13, yy),
                      font_size="7.5px", fill="#222"))

dwg.save()
print(f"\n✅ Figma SVG: {DATA_DIR}/octopus_ideogram_figma.svg")
print("   → Перетащи файл в Figma (File → Place image или просто drag & drop)")
