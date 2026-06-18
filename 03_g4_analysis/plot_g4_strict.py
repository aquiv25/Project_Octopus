"""
G4 structural alignment — STRICT pattern only:
G{3,6}[ACGT]{1,7}G{3,6}[ACGT]{1,7}G{3,6}[ACGT]{1,7}G{3,6}
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

G_COLOR      = "#1B5E20"
LOOP_COLORS  = ["#FFB300", "#F57C00", "#BF360C"]   # short / medium / long
VAR_LOOP_COL = "#C62828"
VAR_GRUN_COL = "#1565C0"
BULGE_COL    = "#6A1B9A"

def parse_struct(seq):
    m = re.match(
        r'(G{3,6})([ACGT]{1,7})(G{3,6})([ACGT]{1,7})(G{3,6})([ACGT]{1,7})(G{3,6})',
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

def draw_panel(ax, rows, gene_name):
    ax.set_facecolor("white")

    # one row per species (keep first hit)
    best = {}
    for sp in SP_ORDER:
        sp_rows = [r for r in rows if r["species"] == sp]
        if sp_rows:
            best[sp] = sp_rows[0]

    if not best:
        ax.text(0.5, 0.5, f"{gene_name}\n—", ha="center", va="center",
                transform=ax.transAxes, color="#bbb", fontsize=11)
        ax.axis("off"); return

    structs = {sp: parse_struct(r["G4_seq"]) for sp, r in best.items()}
    structs  = {sp: s for sp, s in structs.items() if s}

    # variation flags
    loop_vars = [len({structs[sp]["loop_lens"][i] for sp in structs}) > 1 for i in range(3)]
    grun_vars = [len({structs[sp]["g_runs"][i]    for sp in structs}) > 1 for i in range(4)]

    y_tick_pos, y_tick_lab = [], []

    for row_idx, sp in enumerate(reversed(SP_ORDER)):
        if sp not in structs: continue
        s = structs[sp]
        y = row_idx

        has_bulge   = any(g >= 5 for g in s["g_runs"])
        non_uniform = len(set(s["g_runs"])) > 1

        parts = [
            (s["G1"], G_COLOR,              None,         0, "g"),
            (s["L1"], loop_color(len(s["L1"])), len(s["L1"]), 0, "l"),
            (s["G2"], G_COLOR,              None,         1, "g"),
            (s["L2"], loop_color(len(s["L2"])), len(s["L2"]), 1, "l"),
            (s["G3"], G_COLOR,              None,         2, "g"),
            (s["L3"], loop_color(len(s["L3"])), len(s["L3"]), 2, "l"),
            (s["G4"], G_COLOR,              None,         3, "g"),
        ]

        x = 0
        for seq_part, bg, loop_len, idx, kind in parts:
            w = len(seq_part)
            varies    = grun_vars[idx] if kind == "g" else loop_vars[idx]
            edge_col  = (VAR_LOOP_COL if kind == "l" else VAR_GRUN_COL) if varies else "none"
            lw        = 2.0 if edge_col != "none" else 0

            ax.add_patch(FancyBboxPatch(
                (x, y - 0.42), w, 0.84,
                boxstyle="square,pad=0",
                linewidth=lw, edgecolor=edge_col,
                facecolor=bg, alpha=0.92, zorder=2))

            for ci, ch in enumerate(seq_part):
                ax.text(x + ci + 0.5, y, ch, ha="center", va="center",
                        fontsize=8.5, fontweight="bold", color="white",
                        fontfamily="monospace", zorder=3)

            if loop_len is not None:
                txt_col = VAR_LOOP_COL if loop_vars[idx] else loop_color(loop_len)
                ax.text(x + w/2, y + 0.53, str(loop_len),
                        ha="center", va="bottom", fontsize=8,
                        color=txt_col, fontweight="bold", zorder=4)
                if loop_vars[idx]:
                    ax.text(x + w/2, y + 0.92, "^",
                            ha="center", va="bottom", fontsize=7,
                            color=VAR_LOOP_COL, zorder=4)

            if kind == "g":
                g_len = s["g_runs"][idx]
                if grun_vars[idx]:
                    ax.text(x + w/2, y - 0.65, f"G{g_len}",
                            ha="center", va="top", fontsize=6.5,
                            color=VAR_GRUN_COL, fontweight="bold", zorder=4)
                if g_len >= 5:
                    ax.text(x + w/2, y + 0.56, "+",
                            ha="center", va="bottom", fontsize=6,
                            color=BULGE_COL, zorder=4)
            x += w

        # G-run summary
        g_str = "-".join(str(g) for g in s["g_runs"])
        g_col = "#00838F" if non_uniform else "#2E7D32"
        ax.text(x + 0.5, y, f"G({g_str})", va="center", ha="left",
                fontsize=8, color=g_col, fontweight="bold")
        if non_uniform:
            ax.text(x + 5.5, y, "!=", va="center", ha="left",
                    fontsize=8, color="#00838F", fontweight="bold")
        if has_bulge:
            ax.text(x + 7.0, y, "[B]", va="center", ha="left",
                    fontsize=7, color=BULGE_COL)

        y_tick_pos.append(y)
        y_tick_lab.append(SP_LABELS[sp])

    if not y_tick_pos:
        ax.axis("off"); return

    max_len = max(len(r["G4_seq"]) for r in best.values()) + 12
    ax.set_xlim(-0.5, max_len)
    ax.set_ylim(-0.9, len(y_tick_pos) - 0.1)
    ax.set_yticks(y_tick_pos)
    ax.set_yticklabels(y_tick_lab, fontsize=10, fontstyle="italic")
    ax.set_title(gene_name, fontsize=13, fontweight="bold", pad=6)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", bottom=False, labelbottom=False)
    for sp in ["top","right","bottom"]: ax.spines[sp].set_visible(False)
    ax.spines["left"].set_color("#ddd")

    # variation summary below panel
    notes = []
    if any(loop_vars):
        notes.append("Loop var: " + ", ".join(f"L{i+1}" for i,v in enumerate(loop_vars) if v))
    if any(grun_vars):
        notes.append("G-run var: " + ", ".join(f"G{i+1}" for i,v in enumerate(grun_vars) if v))
    bulge_sp = [SP_LABELS[sp][:3] for sp in structs if any(g>=5 for g in structs[sp]["g_runs"])]
    if bulge_sp:
        notes.append(f"Bulge (G≥5): {', '.join(bulge_sp)}")
    if notes:
        ax.text(0.01, -0.14, "  |  ".join(notes),
                transform=ax.transAxes, fontsize=7.5, color="#555",
                va="top", ha="left")

# ── Load data ──────────────────────────────────────────────────────────────
rows = []
with open(BASE / "g4_strict_only.tsv") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        rows.append(r)

by_gene = defaultdict(list)
for r in rows:
    by_gene[r["gene"]].append(r)

# genes present, sorted by number of species (desc) then alphabetically
genes = sorted(by_gene.keys(),
               key=lambda g: (-len({r["species"] for r in by_gene[g]}), g))
print(f"Genes with strict G4: {genes}")
for g in genes:
    sps = [r["species"] for r in by_gene[g]]
    print(f"  {g}: {sps}")

# ── Layout: 1 row if ≤4 genes, else 2 rows ─────────────────────────────────
n = len(genes)
ncols = min(n, 4)
nrows = (n + ncols - 1) // ncols
fig_w = ncols * 6.5
fig_h = nrows * 3.8

fig, axes = plt.subplots(nrows, ncols, figsize=(fig_w, fig_h),
                         facecolor="white",
                         gridspec_kw={"hspace": 0.75, "wspace": 0.3})

# Normalise axes to 2D list
if nrows == 1 and ncols == 1:
    axes = [[axes]]
elif nrows == 1:
    axes = [list(axes)]
elif ncols == 1:
    axes = [[a] for a in axes]
else:
    axes = [list(row) for row in axes]

for idx, gene in enumerate(genes):
    draw_panel(axes[idx // ncols][idx % ncols], by_gene[gene], gene)

for idx in range(n, nrows * ncols):
    axes[idx // ncols][idx % ncols].axis("off")

# ── Legend ─────────────────────────────────────────────────────────────────
struct_leg = [
    mpatches.Patch(fc=G_COLOR,        label="G-run (G≥3)",          ec="none"),
    mpatches.Patch(fc=LOOP_COLORS[0], label="Short loop  1–2 nt",   ec="none"),
    mpatches.Patch(fc=LOOP_COLORS[1], label="Medium loop 3–6 nt",   ec="none"),
    mpatches.Patch(fc=LOOP_COLORS[2], label="Long loop   >6 nt",    ec="none"),
    mpatches.Patch(fc="white", ec=VAR_LOOP_COL, lw=2,
                   label="^ loop varies between species"),
    mpatches.Patch(fc="white", ec=VAR_GRUN_COL, lw=2,
                   label="G-run varies between species"),
    mpatches.Patch(fc="white", ec="none",
                   label="[B] = bulge (G≥5)  |  != = non-uniform G-runs"),
]
fig.legend(handles=struct_leg, loc="lower center", ncol=4,
           fontsize=9.5, frameon=True, fancybox=False,
           edgecolor="#ccc", bbox_to_anchor=(0.5, -0.06),
           title="Legend", title_fontsize=9.5)

fig.suptitle(
    "Canonical G-quadruplex motifs (G{3,6}·N{1,7}·G{3,6}·N{1,7}·G{3,6}·N{1,7}·G{3,6})"
    "  in promoters of epigenetic genes — Cephalopoda\n"
    "Numbers above loops = loop length (nt)  ·  G(x-x-x-x) = G-run lengths  "
    "·  coloured border = varies across species",
    fontsize=11.5, fontweight="bold", y=1.03, color="#111"
)

out_svg = BASE / "g4_strict_alignment.svg"
out_png = BASE / "g4_strict_alignment.png"
fig.savefig(out_svg, format="svg", bbox_inches="tight", facecolor="white")
fig.savefig(out_png, format="png", bbox_inches="tight", dpi=150, facecolor="white")
print(f"\n✅ {out_svg}")
print(f"✅ {out_png}")
plt.close()
