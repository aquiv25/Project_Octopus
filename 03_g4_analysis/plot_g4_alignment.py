"""
Publication-ready figure: G4 structural alignment in promoters of epigenetic genes
across 5 cephalopod species.

Shows per-gene panels with:
  - Colored sequence blocks (G-runs = dark green, loops = orange gradient by length)
  - Loop length annotations
  - G-run count annotations
  - Variation markers (bulge, non-uniform G)
"""
import re, csv
from pathlib import Path
from collections import defaultdict
import matplotlib
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["font.family"] = "Arial"
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import numpy as np

BASE = Path("/Users/dassagaripova/Documents/Claude/Projects/Article/g4_promoter/results")

SP_COLORS = {
    "O_bimaculoides": "#D6604D",
    "O_vulgaris":     "#4393C3",
    "O_sinensis":     "#2CA02C",
    "A_fangsiao":     "#FF7F0E",
    "A_hians":        "#8C564B",
}
SP_LABELS = {
    "O_bimaculoides": "O. bimaculoides",
    "O_vulgaris":     "O. vulgaris",
    "O_sinensis":     "O. sinensis",
    "A_fangsiao":     "A. fangsiao",
    "A_hians":        "A. hians",
}
SP_ORDER = list(SP_LABELS.keys())

G_COLOR   = "#1B5E20"   # dark green: G-runs
LOOP_COLORS = ["#FFB300","#F57C00","#E64A19"]  # short→long loops: amber→deep orange→red-orange

def parse_g4_structure(seq):
    m = re.match(
        r'(G{2,7})([ACGTN]{1,15})(G{2,7})([ACGTN]{1,15})(G{2,7})([ACGTN]{1,15})(G{2,7})',
        seq, re.IGNORECASE)
    if not m: return None
    g1, l1, g2, l2, g3, l3, g4 = m.groups()
    return {"G1":g1,"L1":l1,"G2":g2,"L2":l2,"G3":g3,"L3":l3,"G4":g4,
            "g_runs":[len(g1),len(g2),len(g3),len(g4)],
            "loop_lens":[len(l1),len(l2),len(l3)]}

def loop_color(loop_len):
    if loop_len <= 2:   return LOOP_COLORS[0]
    elif loop_len <= 6: return LOOP_COLORS[1]
    else:               return LOOP_COLORS[2]

def draw_gene_panel(ax, rows_by_sp, gene_name):
    """Draw one gene panel with all species' G4 sequences aligned."""
    ax.set_facecolor("#FAFAFA")

    # Keep only the BEST hit per species (first/strict preferred)
    best = {}
    for sp in SP_ORDER:
        sp_rows = [r for r in rows_by_sp if r["species"] == sp]
        if sp_rows:
            # prefer strict
            strict = [r for r in sp_rows if r["type"] == "strict"]
            best[sp] = (strict or sp_rows)[0]

    if not best:
        ax.text(0.5,0.5,f"{gene_name}\n—",ha="center",va="center",
                transform=ax.transAxes,color="#aaa",fontsize=12); ax.axis("off"); return

    y_tick_pos = []
    y_tick_lab = []

    for row_idx, sp in enumerate(reversed(SP_ORDER)):
        if sp not in best: continue
        row  = best[sp]
        y    = row_idx
        struct = parse_g4_structure(row["G4_seq"])
        if not struct: continue

        sp_color = SP_COLORS[sp]

        # Draw colored blocks + characters
        x = 0
        parts = [
            (struct["G1"], G_COLOR, None),
            (struct["L1"], loop_color(len(struct["L1"])), len(struct["L1"])),
            (struct["G2"], G_COLOR, None),
            (struct["L2"], loop_color(len(struct["L2"])), len(struct["L2"])),
            (struct["G3"], G_COLOR, None),
            (struct["L3"], loop_color(len(struct["L3"])), len(struct["L3"])),
            (struct["G4"], G_COLOR, None),
        ]
        for seq_part, bg, loop_len in parts:
            w = len(seq_part)
            rect = FancyBboxPatch(
                (x, y - 0.42), w, 0.84,
                boxstyle="square,pad=0",
                linewidth=0.3, edgecolor="#fff", facecolor=bg, alpha=0.92
            )
            ax.add_patch(rect)
            # characters
            for ci, ch in enumerate(seq_part):
                ax.text(x + ci + 0.5, y, ch, ha="center", va="center",
                        fontsize=7.5, fontweight="bold", color="white",
                        fontfamily="monospace")
            # loop length annotation above block
            if loop_len is not None:
                ax.text(x + w/2, y + 0.52, str(loop_len),
                        ha="center", va="bottom", fontsize=6.5,
                        color=bg, fontweight="bold")
            x += w

        # G-run annotation at right
        g_str = "-".join(str(g) for g in struct["g_runs"])
        ax.text(x + 0.3, y, f"G({g_str})", va="center", ha="left",
                fontsize=7, color=G_COLOR, fontweight="bold")

        # Type badge
        badge_col = "#1565C0" if row["type"] == "strict" else "#888"
        ax.text(x + 3.5, y, row["type"][0].upper(), va="center", ha="left",
                fontsize=6.5, color=badge_col, fontweight="bold")

        y_tick_pos.append(y)
        y_tick_lab.append(SP_LABELS[sp])

    max_len = max(len(r["G4_seq"]) for r in best.values()) + 6
    ax.set_xlim(-0.5, max_len)
    ax.set_ylim(-0.7, len(best) - 0.1)
    ax.set_yticks(y_tick_pos)
    ax.set_yticklabels(y_tick_lab, fontsize=9, fontstyle="italic", color="#222")
    ax.set_title(gene_name, fontsize=12, fontweight="bold", color="#111", pad=5)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", bottom=False, labelbottom=False)
    for spine in ["top","right","bottom"]: ax.spines[spine].set_visible(False)
    ax.spines["left"].set_color("#ccc")

# ── Load table ────────────────────────────────────────────────────────────
table_path = BASE / "g4_promoters_table.tsv"
rows = []
with open(table_path) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        rows.append(row)

by_gene = defaultdict(list)
for r in rows:
    by_gene[r["gene"]].append(r)

GENE_ORDER = ["DNMT1","TET3","EZH2","KMT2A","CHD1","BRD3","ADAR1","UHRF1"]
genes_present = [g for g in GENE_ORDER if g in by_gene]

# ── Layout: 2 rows × 4 cols ───────────────────────────────────────────────
ncols, nrows = 4, 2
fig, axes = plt.subplots(nrows, ncols, figsize=(26, 9), facecolor="white",
                         gridspec_kw={"hspace": 0.55, "wspace": 0.35})

for idx, gene in enumerate(genes_present):
    ax = axes[idx // ncols][idx % ncols]
    draw_gene_panel(ax, by_gene[gene], gene)

# Hide empty
for idx in range(len(genes_present), nrows * ncols):
    axes[idx // ncols][idx % ncols].axis("off")

# ── Legend ─────────────────────────────────────────────────────────────────
leg_patches = [
    mpatches.Patch(fc=G_COLOR,         label="G-run (tetrad)",    ec="none"),
    mpatches.Patch(fc=LOOP_COLORS[0],  label="Short loop (1–2 nt)",  ec="none"),
    mpatches.Patch(fc=LOOP_COLORS[1],  label="Medium loop (3–6 nt)", ec="none"),
    mpatches.Patch(fc=LOOP_COLORS[2],  label="Long loop (>6 nt)",   ec="none"),
    mpatches.Patch(fc="#1565C0",       label="Strict G4 (G≥3)",     ec="none", alpha=0.7),
    mpatches.Patch(fc="#888",          label="Relaxed G4 (G≥2)",    ec="none", alpha=0.7),
]
fig.legend(handles=leg_patches, loc="lower center", ncol=6, fontsize=10,
           frameon=False, bbox_to_anchor=(0.5, -0.025),
           handlelength=1.2, handleheight=0.9)

# Species color dots
for i, sp in enumerate(SP_ORDER):
    fig.text(0.07 + i*0.17, -0.045, f"● {SP_LABELS[sp]}",
             color=SP_COLORS[sp], fontsize=9.5, fontstyle="italic",
             transform=fig.transFigure, ha="left", va="bottom")

fig.suptitle(
    "G-quadruplex motifs in promoters of epigenetic regulatory genes across Cephalopoda\n"
    "Numbers above loops indicate loop length (nt)  ·  G(x-x-x-x) = G-run lengths  "
    "·  S = strict (G≥3)  ·  R = relaxed (G≥2)",
    fontsize=12, fontweight="bold", y=1.015, color="#111"
)

out_svg = BASE / "g4_promoters_alignment_v2.svg"
out_png = BASE / "g4_promoters_alignment_v2.png"
fig.savefig(out_svg, format="svg", bbox_inches="tight", facecolor="white")
fig.savefig(out_png, format="png", bbox_inches="tight", dpi=150, facecolor="white")
print(f"✅ {out_svg}")
print(f"✅ {out_png}")
plt.close()
