"""
G4 structural alignment for 5 target genes (TPM, MED28, SIK1, ABDA, ZNF271)
across cephalopod species.
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

SP_LABELS = {
    "O_bimaculoides": "O. bimaculoides",
    "O_vulgaris":     "O. vulgaris",
    "O_sinensis":     "O. sinensis",
    "A_fangsiao":     "A. fangsiao",
    "A_hians":        "A. hians",
}
SP_ORDER = list(SP_LABELS.keys())

G_COLOR     = "#1B5E20"
LOOP_COLORS = ["#FFB300", "#F57C00", "#BF360C"]
VAR_LOOP    = "#C62828"
VAR_GRUN    = "#1565C0"
BULGE_COL   = "#6A1B9A"
STRICT_COL  = "#1565C0"
RELAXED_COL = "#9E9E9E"

GENE_LABELS = {
    "TPM":    "TPM (Tropomyosin)",
    "MED28":  "MED28 (Mediator 28)",
    "SIK1":   "SIK1 (Salt-inducible kinase)",
    "ABDA":   "ABDA (Homeobox abd-A)",
    "ZNF271": "ZNF271 (Zinc finger 271)",
}

def parse_struct(seq):
    """Try strict first, then relaxed."""
    m = re.match(
        r'(G{3,6})([ACGT]{1,7})(G{3,6})([ACGT]{1,7})(G{3,6})([ACGT]{1,7})(G{3,6})',
        seq, re.IGNORECASE)
    if m:
        g1,l1,g2,l2,g3,l3,g4 = m.groups()
        return dict(G1=g1,L1=l1,G2=g2,L2=l2,G3=g3,L3=l3,G4=g4,
                    g_runs=[len(g1),len(g2),len(g3),len(g4)],
                    loop_lens=[len(l1),len(l2),len(l3)], stype="strict")
    m = re.match(
        r'(G{2,6})([ACGT]{1,12})(G{2,6})([ACGT]{1,12})(G{2,6})([ACGT]{1,12})(G{2,6})',
        seq, re.IGNORECASE)
    if m:
        g1,l1,g2,l2,g3,l3,g4 = m.groups()
        return dict(G1=g1,L1=l1,G2=g2,L2=l2,G3=g3,L3=l3,G4=g4,
                    g_runs=[len(g1),len(g2),len(g3),len(g4)],
                    loop_lens=[len(l1),len(l2),len(l3)], stype="relaxed")
    return None

def loop_color(n):
    if n <= 2:   return LOOP_COLORS[0]
    elif n <= 7: return LOOP_COLORS[1]
    else:        return LOOP_COLORS[2]

def draw_panel(ax, rows, gene_name):
    ax.set_facecolor("white")

    best = {}
    for sp in SP_ORDER:
        sp_rows = [r for r in rows if r["species"] == sp]
        if sp_rows:
            best[sp] = sp_rows[0]

    if not best:
        ax.text(0.5, 0.5, f"{gene_name}\n—", ha="center", va="center",
                transform=ax.transAxes, color="#bbb", fontsize=11)
        ax.axis("off"); return

    structs = {}
    for sp, r in best.items():
        s = parse_struct(r["G4_seq"])
        if s: structs[sp] = s

    if not structs:
        ax.axis("off"); return

    loop_vars = [len({structs[sp]["loop_lens"][i] for sp in structs}) > 1 for i in range(3)]
    grun_vars = [len({structs[sp]["g_runs"][i]    for sp in structs}) > 1 for i in range(4)]

    y_tick_pos, y_tick_lab = [], []

    for row_idx, sp in enumerate(reversed(SP_ORDER)):
        if sp not in structs: continue
        s = structs[sp]
        y = row_idx
        is_strict    = (s["stype"] == "strict")
        has_bulge    = any(g >= 5 for g in s["g_runs"])
        non_uniform  = len(set(s["g_runs"])) > 1

        parts = [
            (s["G1"], G_COLOR, None, 0, "g"),
            (s["L1"], loop_color(len(s["L1"])), len(s["L1"]), 0, "l"),
            (s["G2"], G_COLOR, None, 1, "g"),
            (s["L2"], loop_color(len(s["L2"])), len(s["L2"]), 1, "l"),
            (s["G3"], G_COLOR, None, 2, "g"),
            (s["L3"], loop_color(len(s["L3"])), len(s["L3"]), 2, "l"),
            (s["G4"], G_COLOR, None, 3, "g"),
        ]

        x = 0
        for seq_part, bg, loop_len, idx, kind in parts:
            w = len(seq_part)
            varies   = grun_vars[idx] if kind == "g" else loop_vars[idx]
            edge_col = (VAR_LOOP if kind == "l" else VAR_GRUN) if varies else "none"
            lw = 2.0 if edge_col != "none" else 0
            alpha = 0.95 if is_strict else 0.60

            ax.add_patch(FancyBboxPatch(
                (x, y - 0.42), w, 0.84,
                boxstyle="square,pad=0",
                linewidth=lw, edgecolor=edge_col,
                facecolor=bg, alpha=alpha, zorder=2))

            for ci, ch in enumerate(seq_part):
                ax.text(x + ci + 0.5, y, ch, ha="center", va="center",
                        fontsize=8, fontweight="bold", color="white",
                        fontfamily="monospace", zorder=3)

            if loop_len is not None:
                txt_col = VAR_LOOP if loop_vars[idx] else loop_color(loop_len)
                ax.text(x + w/2, y + 0.53, str(loop_len),
                        ha="center", va="bottom", fontsize=7.5,
                        color=txt_col, fontweight="bold", zorder=4)
                if loop_vars[idx]:
                    ax.text(x + w/2, y + 0.92, "^",
                            ha="center", va="bottom", fontsize=7,
                            color=VAR_LOOP, zorder=4)

            if kind == "g":
                g_len = s["g_runs"][idx]
                if grun_vars[idx]:
                    ax.text(x + w/2, y - 0.65, f"G{g_len}",
                            ha="center", va="top", fontsize=6,
                            color=VAR_GRUN, fontweight="bold", zorder=4)
                if g_len >= 5:
                    ax.text(x + w/2, y + 0.56, "+",
                            ha="center", va="bottom", fontsize=6,
                            color=BULGE_COL, zorder=4)
            x += w

        # G-run summary
        g_str = "-".join(str(g) for g in s["g_runs"])
        g_col = "#00838F" if non_uniform else "#2E7D32"
        ax.text(x + 0.5, y, f"G({g_str})", va="center", ha="left",
                fontsize=7.5, color=g_col, fontweight="bold")
        if non_uniform:
            ax.text(x + 5.2, y, "!=", va="center", ha="left",
                    fontsize=7.5, color="#00838F", fontweight="bold")
        if has_bulge:
            ax.text(x + 6.8, y, "[B]", va="center", ha="left",
                    fontsize=7, color=BULGE_COL)

        # strict/relaxed badge
        badge_col = STRICT_COL if is_strict else RELAXED_COL
        ax.text(x + 9.0, y, "S" if is_strict else "R",
                va="center", ha="left", fontsize=7,
                color=badge_col, fontweight="bold")

        y_tick_pos.append(y)
        y_tick_lab.append(SP_LABELS[sp])

    if not y_tick_pos:
        ax.axis("off"); return

    max_len = max(len(r["G4_seq"]) for r in best.values()) + 13
    ax.set_xlim(-0.5, max_len)
    ax.set_ylim(-0.9, len(y_tick_pos) - 0.1)
    ax.set_yticks(y_tick_pos)
    ax.set_yticklabels(y_tick_lab, fontsize=10, fontstyle="italic")
    ax.set_title(GENE_LABELS.get(gene_name, gene_name), fontsize=12,
                 fontweight="bold", pad=6)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", bottom=False, labelbottom=False)
    for sp in ["top","right","bottom"]: ax.spines[sp].set_visible(False)
    ax.spines["left"].set_color("#ddd")


# ── Load data ──────────────────────────────────────────────────────────────
rows = []
with open(BASE.parent / "results/target_g4_v2.tsv") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        rows.append(r)

# Filter polyG runs (no real loops)
def is_valid(r):
    s = parse_struct(r["G4_seq"])
    if s is None: return False
    # at least one loop must be >= 1
    return all(ll >= 1 for ll in s["loop_lens"])

rows = [r for r in rows if is_valid(r)]

by_gene = defaultdict(list)
for r in rows:
    by_gene[r["gene"]].append(r)

GENE_ORDER = ["SIK1","ZNF271","MED28"]
genes = [g for g in GENE_ORDER if g in by_gene]
print(f"Genes: {genes}")
for g in genes:
    sps = [r['species'] for r in by_gene[g]]
    print(f"  {g}: {sps}")

# ── Layout ─────────────────────────────────────────────────────────────────
n = len(genes)
ncols = 1
nrows = n
fig, axes = plt.subplots(nrows, ncols,
                         figsize=(11, nrows * 3.6),
                         facecolor="white",
                         gridspec_kw={"hspace": 0.70})

axes = [[ax] for ax in axes] if nrows > 1 else [[axes]]

for idx, gene in enumerate(genes):
    draw_panel(axes[idx][0], by_gene[gene], gene)
for idx in range(n, nrows):
    axes[idx][0].axis("off")

# ── Legend ─────────────────────────────────────────────────────────────────
from matplotlib.lines import Line2D

col1 = [
    mpatches.Patch(fc=G_COLOR,        ec="none", label="G-run (tetrad)"),
    mpatches.Patch(fc=LOOP_COLORS[0], ec="none", label="Short loop (1–2 nt)"),
    mpatches.Patch(fc=LOOP_COLORS[1], ec="none", label="Medium loop (3–7 nt)"),
    mpatches.Patch(fc=LOOP_COLORS[2], ec="none", label="Long loop (>7 nt)"),
]
col2 = [
    mpatches.Patch(fc="#e8f0fe", ec=VAR_LOOP, lw=2,
                   label="Loop length varies across species"),
    mpatches.Patch(fc="#e8f0fe", ec=VAR_GRUN, lw=2,
                   label="G-run length varies across species"),
    mpatches.Patch(fc="white", ec="none",
                   label="[B] — G-run ≥ 5 (bulge potential)"),
    mpatches.Patch(fc="white", ec="none",
                   label="!= — non-uniform G-runs within motif"),
]
col3 = [
    mpatches.Patch(fc=STRICT_COL,  ec="none", alpha=0.85,
                   label="S — strict G4  (G≥3, loop ≤7 nt)"),
    mpatches.Patch(fc=RELAXED_COL, ec="none", alpha=0.85,
                   label="R — relaxed G4 (G≥2, loop ≤12 nt)"),
    mpatches.Patch(fc="white", ec="none",
                   label="Numbers above blocks = loop length (nt)"),
    mpatches.Patch(fc="white", ec="none",
                   label="G(x-x-x-x) = lengths of the four G-runs"),
]

# Draw three side-by-side legend boxes
leg1 = fig.legend(handles=col1, loc="lower left",   bbox_to_anchor=(0.01, -0.13),
                  fontsize=9, frameon=True, fancybox=False, edgecolor="#bbb",
                  title="Sequence blocks", title_fontsize=9.5, ncol=1)
leg2 = fig.legend(handles=col2, loc="lower center", bbox_to_anchor=(0.50, -0.13),
                  fontsize=9, frameon=True, fancybox=False, edgecolor="#bbb",
                  title="Variation markers", title_fontsize=9.5, ncol=1)
leg3 = fig.legend(handles=col3, loc="lower right",  bbox_to_anchor=(0.99, -0.13),
                  fontsize=9, frameon=True, fancybox=False, edgecolor="#bbb",
                  title="Motif type & annotations", title_fontsize=9.5, ncol=1)
fig.add_artist(leg1)
fig.add_artist(leg2)


out_png = BASE / "target_g4_alignment.png"
out_svg = BASE / "target_g4_alignment.svg"
fig.savefig(out_png, format="png", bbox_inches="tight", dpi=150, facecolor="white")
fig.savefig(out_svg, format="svg", bbox_inches="tight", facecolor="white")
print(f"\n✅ {out_png}")
print(f"✅ {out_svg}")
plt.close()
