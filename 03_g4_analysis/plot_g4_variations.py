"""
G4 structural alignment with variation annotations.
For each gene, compares G4 structures across species and marks:
  - Loop length variation   ([v] above loop if differs between species)
  - G-run length variation  (★ if G-run count differs)
  - Bulge                   ([B] if any G-run > 4)
  - Non-uniform G-runs      ([!] if G-runs within one motif differ)
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

G_COLOR = "#1B5E20"
LOOP_COLORS = ["#FFB300", "#F57C00", "#BF360C"]  # short / medium / long

VAR_COLORS = {
    "loop_var":   "#C62828",   # red    — loop length varies across species
    "grun_var":   "#1565C0",   # blue   — G-run length varies
    "bulge":      "#6A1B9A",   # purple — G-run ≥ 5
    "nonuniform": "#00838F",   # teal   — non-uniform G-runs within motif
}

def parse_struct(seq):
    m = re.match(
        r'(G{2,7})([ACGTN]{1,15})(G{2,7})([ACGTN]{1,15})(G{2,7})([ACGTN]{1,15})(G{2,7})',
        seq, re.IGNORECASE)
    if not m: return None
    g1,l1,g2,l2,g3,l3,g4 = m.groups()
    return dict(G1=g1,L1=l1,G2=g2,L2=l2,G3=g3,L3=l3,G4=g4,
                g_runs=[len(g1),len(g2),len(g3),len(g4)],
                loop_lens=[len(l1),len(l2),len(l3)])

def loop_color(n):
    if n <= 2:   return LOOP_COLORS[0]
    elif n <= 6: return LOOP_COLORS[1]
    else:        return LOOP_COLORS[2]

def draw_gene_panel(ax, all_rows, gene_name):
    ax.set_facecolor("white")

    # Best hit per species
    best = {}
    for sp in SP_ORDER:
        rows = [r for r in all_rows if r["species"] == sp]
        if rows:
            strict = [r for r in rows if r["type"] == "strict"]
            best[sp] = (strict or rows)[0]

    if not best:
        ax.text(0.5, 0.5, f"{gene_name}\n—", ha="center", va="center",
                transform=ax.transAxes, color="#bbb", fontsize=11)
        ax.axis("off"); return

    structs = {sp: parse_struct(r["G4_seq"]) for sp, r in best.items()}
    structs = {sp: s for sp, s in structs.items() if s}

    # ── Compute inter-species variation ──────────────────────────────────
    # Loop length variation per loop position
    loop_vars = []
    for li in range(3):
        lens = [structs[sp]["loop_lens"][li] for sp in structs]
        loop_vars.append(len(set(lens)) > 1)   # True if varies

    # G-run length variation per G position
    grun_vars = []
    for gi in range(4):
        runs = [structs[sp]["g_runs"][gi] for sp in structs]
        grun_vars.append(len(set(runs)) > 1)

    # ── Draw rows ─────────────────────────────────────────────────────────
    y_tick_pos, y_tick_lab = [], []

    for row_idx, sp in enumerate(reversed(SP_ORDER)):
        if sp not in structs: continue
        s = structs[sp]
        y = row_idx
        sp_color = SP_COLORS[sp]

        # Variation flags for this species
        has_bulge     = any(g >= 5 for g in s["g_runs"])
        non_uniform   = len(set(s["g_runs"])) > 1

        parts = [
            (s["G1"], G_COLOR,              None,           0, "g"),
            (s["L1"], loop_color(len(s["L1"])), len(s["L1"]), 0, "l"),
            (s["G2"], G_COLOR,              None,           1, "g"),
            (s["L2"], loop_color(len(s["L2"])), len(s["L2"]), 1, "l"),
            (s["G3"], G_COLOR,              None,           2, "g"),
            (s["L3"], loop_color(len(s["L3"])), len(s["L3"]), 2, "l"),
            (s["G4"], G_COLOR,              None,           3, "g"),
        ]

        x = 0
        g_idx = 0
        l_idx = 0
        for seq_part, bg, loop_len, idx, kind in parts:
            w = len(seq_part)

            # Check if this position varies
            varies = (grun_vars[idx] if kind == "g" else loop_vars[idx])
            edge_color = VAR_COLORS["loop_var"] if (kind == "l" and varies) else \
                         VAR_COLORS["grun_var"] if (kind == "g" and varies) else "none"
            lw = 2.0 if edge_color != "none" else 0

            rect = FancyBboxPatch(
                (x, y - 0.41), w, 0.82,
                boxstyle="square,pad=0",
                linewidth=lw, edgecolor=edge_color,
                facecolor=bg, alpha=0.90, zorder=2
            )
            ax.add_patch(rect)

            # Characters
            for ci, ch in enumerate(seq_part):
                ax.text(x + ci + 0.5, y, ch, ha="center", va="center",
                        fontsize=7, fontweight="bold", color="white",
                        fontfamily="monospace", zorder=3)

            # Loop length annotation (number above block)
            if loop_len is not None:
                txt_color = VAR_COLORS["loop_var"] if loop_vars[idx] else loop_color(loop_len)
                ax.text(x + w/2, y + 0.52, str(loop_len),
                        ha="center", va="bottom", fontsize=6.5,
                        color=txt_color, fontweight="bold", zorder=4)

                # Triangle marker if loop varies
                if loop_vars[idx]:
                    ax.text(x + w/2, y + 0.88, "^",
                            ha="center", va="bottom", fontsize=6,
                            color=VAR_COLORS["loop_var"], zorder=4)

            # G-run annotation
            if kind == "g":
                g_len = s["g_runs"][idx]
                # Mark G-run variation
                if grun_vars[idx]:
                    ax.text(x + w/2, y - 0.62, f"G{g_len}",
                            ha="center", va="top", fontsize=5.5,
                            color=VAR_COLORS["grun_var"], fontweight="bold", zorder=4)
                # Mark bulge
                if g_len >= 5:
                    ax.text(x + w/2, y + 0.55, "+",
                            ha="center", va="bottom", fontsize=5.5,
                            color=VAR_COLORS["bulge"], zorder=4)

            x += w

        # G-run summary at right
        g_str = "-".join(str(g) for g in s["g_runs"])
        g_col = VAR_COLORS["nonuniform"] if non_uniform else "#2E7D32"
        ax.text(x + 0.4, y, f"G({g_str})", va="center", ha="left",
                fontsize=7, color=g_col, fontweight="bold")

        # Non-uniform marker
        if non_uniform:
            ax.text(x + 4.2, y, "!=", va="center", ha="left",
                    fontsize=8, color=VAR_COLORS["nonuniform"], fontweight="bold")

        # Bulge marker
        if has_bulge:
            ax.text(x + 5.0, y, "+", va="center", ha="left",
                    fontsize=7, color=VAR_COLORS["bulge"])

        # Type (S/R)
        type_col = "#1565C0" if best[sp]["type"] == "strict" else "#999"
        ax.text(x + 6.2, y, best[sp]["type"][0].upper(), va="center", ha="left",
                fontsize=6.5, color=type_col, fontweight="bold")

        y_tick_pos.append(y)
        y_tick_lab.append(SP_LABELS[sp])

    if not y_tick_pos: ax.axis("off"); return

    max_len = max(len(r["G4_seq"]) for r in best.values()) + 10
    ax.set_xlim(-0.5, max_len)
    ax.set_ylim(-0.9, len(y_tick_pos) - 0.1)
    ax.set_yticks(y_tick_pos)
    ax.set_yticklabels(y_tick_lab, fontsize=9, fontstyle="italic")
    ax.set_title(gene_name, fontsize=12, fontweight="bold", pad=5)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", bottom=False, labelbottom=False)
    for sp in ["top", "right", "bottom"]: ax.spines[sp].set_visible(False)
    ax.spines["left"].set_color("#ddd")

    # ── Variation summary box ─────────────────────────────────────────────
    var_notes = []
    if any(loop_vars):
        lpos = [f"L{i+1}" for i, v in enumerate(loop_vars) if v]
        var_notes.append(f"Loop var: {', '.join(lpos)}")
    if any(grun_vars):
        gpos = [f"G{i+1}" for i, v in enumerate(grun_vars) if v]
        var_notes.append(f"G-run var: {', '.join(gpos)}")
    bulge_sp = [SP_LABELS[sp][:3] for sp in structs if any(g>=5 for g in structs[sp]["g_runs"])]
    if bulge_sp:
        var_notes.append(f"Bulge (G≥5): {', '.join(bulge_sp)}")

    if var_notes:
        ax.text(0.01, -0.12, "  ".join(var_notes),
                transform=ax.transAxes, fontsize=7, color="#555",
                va="top", ha="left")

# ── Load data ─────────────────────────────────────────────────────────────
rows = []
with open(BASE / "g4_promoters_table.tsv") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        rows.append(r)

by_gene = defaultdict(list)
for r in rows:
    by_gene[r["gene"]].append(r)

GENE_ORDER = ["DNMT1", "TET3", "EZH2", "CHD1", "BRD3", "ADAR1", "UHRF1"]
genes = [g for g in GENE_ORDER if g in by_gene]

# ── Layout ────────────────────────────────────────────────────────────────
ncols, nrows = 4, 2
fig, axes = plt.subplots(nrows, ncols,
                         figsize=(26, 10),
                         facecolor="white",
                         gridspec_kw={"hspace": 0.7, "wspace": 0.3})

for idx, gene in enumerate(genes):
    draw_gene_panel(axes[idx // ncols][idx % ncols], by_gene[gene], gene)

for idx in range(len(genes), nrows * ncols):
    axes[idx // ncols][idx % ncols].axis("off")

# ── Legend ────────────────────────────────────────────────────────────────
struct_leg = [
    mpatches.Patch(fc=G_COLOR,         label="G-run (tetrad)",       ec="none"),
    mpatches.Patch(fc=LOOP_COLORS[0],  label="Short loop (1–2 nt)",  ec="none"),
    mpatches.Patch(fc=LOOP_COLORS[1],  label="Medium loop (3–6 nt)", ec="none"),
    mpatches.Patch(fc=LOOP_COLORS[2],  label="Long loop (>6 nt)",    ec="none"),
]
var_leg = [
    mpatches.Patch(fc="white", ec=VAR_COLORS["loop_var"],  lw=2,
                   label="[v] Loop length varies across species"),
    mpatches.Patch(fc="white", ec=VAR_COLORS["grun_var"],  lw=2,
                   label="★ G-run length varies across species"),
    mpatches.Patch(fc="white", ec="none",
                   label="[B] Bulge (G-run ≥ 5)"),
    mpatches.Patch(fc="white", ec="none",
                   label="[!] Non-uniform G-runs within motif"),
    mpatches.Patch(fc="white", ec="none",
                   label="S = strict G≥3 · R = relaxed G≥2"),
]
# Color dots for species
sp_leg = [mpatches.Patch(fc=SP_COLORS[sp], label=SP_LABELS[sp], ec="none")
          for sp in SP_ORDER]

fig.legend(handles=struct_leg + var_leg, loc="lower left",
           ncol=4, fontsize=9, frameon=True,
           fancybox=False, edgecolor="#ccc",
           bbox_to_anchor=(0.01, -0.07),
           title="Structure", title_fontsize=9)

fig.legend(handles=sp_leg, loc="lower right",
           ncol=5, fontsize=9, frameon=True,
           fancybox=False, edgecolor="#ccc",
           bbox_to_anchor=(0.99, -0.07),
           title="Species", title_fontsize=9)

fig.suptitle(
    "G-quadruplex motifs in promoters of epigenetic regulatory genes — Cephalopoda\n"
    "Variation markers: [v] loop length  ·  colored border = varies across species  "
    "·  G(x-x-x-x) = G-run lengths  ·  [!] non-uniform  ·  [B] bulge (G≥5)",
    fontsize=12, fontweight="bold", y=1.02, color="#111"
)

out_svg = BASE / "g4_alignment_variations.svg"
out_png = BASE / "g4_alignment_variations.png"
fig.savefig(out_svg, format="svg", bbox_inches="tight", facecolor="white")
fig.savefig(out_png, format="png", bbox_inches="tight", dpi=150, facecolor="white")
print(f"✅ {out_svg}")
print(f"✅ {out_png}")
plt.close()
