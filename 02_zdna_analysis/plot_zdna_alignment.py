"""
Z-DNA forming sequence alignment in promoters across 3 octopus species.
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

BASE = Path("/Users/dassagaripova/Documents/Claude/Projects/Article/zdna_promoter")

SP_ORDER  = ["O_bimaculoides", "O_vulgaris", "A_fangsiao"]
SP_LABELS = {
    "O_bimaculoides": "O. bimaculoides",
    "O_vulgaris":     "O. vulgaris",
    "A_fangsiao":     "A. fangsiao",
}

# Color dinucleotide repeats
DINU_COLORS = {
    "CG": "#B71C1C",   # dark red
    "GC": "#C62828",   # red
    "CA": "#1565C0",   # blue
    "TG": "#1976D2",   # light blue
    "AT": "#2E7D32",   # green
    "TA": "#388E3C",   # light green
    "AC": "#6A1B9A",   # purple
    "GT": "#7B1FA2",   # light purple
}
DEFAULT_COL = "#546E7A"   # blue-grey for other

GENE_LABELS = {
    "guanylate kinase":                        "Guanylate kinase",
    "metalloendopeptidase oma1, mitochondrial":"Metalloendopeptidase OMA1",
    "histone h2a":                             "Histone H2A",
    "60s ribosomal protein l3":                "60S Ribosomal protein L3",
    "neurotrypsin":                            "Neurotrypsin",
    "adar1":                                    "ADAR1 (dsRNA adenosine deaminase)",
}

def classify_dinu(seq):
    """Find dominant dinucleotide repeat in sequence."""
    seq = seq.upper()
    best = ("XX", 0)
    for i in range(0, len(seq)-1, 2):
        di = seq[i:i+2]
        if len(di) == 2:
            count = seq.count(di)
            if count > best[1]:
                best = (di, count)
    return best[0]

def nt_color(ch, pos, seq):
    """Color nucleotide based on dinucleotide context."""
    if pos + 1 < len(seq):
        di = seq[pos:pos+2].upper()
        if di in DINU_COLORS: return DINU_COLORS[di]
    if pos > 0:
        di = seq[pos-1:pos+1].upper()
        if di in DINU_COLORS: return DINU_COLORS[di]
    return DEFAULT_COL

def draw_panel(ax, rows, gene_name):
    ax.set_facecolor("white")
    best = {r["species"]: r for r in rows}
    if not best: ax.axis("off"); return

    # Variation: positions where nucleotide differs across species
    seqs = {sp: best[sp]["zdna_seq"].upper() for sp in SP_ORDER if sp in best}
    max_len = max(len(s) for s in seqs.values())

    # Pad sequences to same length for comparison
    padded = {sp: s.ljust(max_len, "-") for sp, s in seqs.items()}

    y_tick_pos, y_tick_lab = [], []

    for row_idx, sp in enumerate(reversed(SP_ORDER)):
        if sp not in seqs: continue
        seq = seqs[sp]
        y = row_idx

        # Classify dominant dinucleotide
        dominant = classify_dinu(seq)
        bg_col = DINU_COLORS.get(dominant, DEFAULT_COL)

        # Check position variation
        pos_varies = []
        for pi in range(len(seq)):
            nts = set()
            for s2 in padded.values():
                if pi < len(s2): nts.add(s2[pi])
            pos_varies.append(len(nts) > 1)

        # Draw block
        w = len(seq)
        ax.add_patch(FancyBboxPatch(
            (0, y - 0.42), w, 0.84,
            boxstyle="square,pad=0",
            linewidth=0, facecolor=bg_col, alpha=0.25, zorder=1))

        # Draw each nucleotide
        for ci, ch in enumerate(seq):
            varies = pos_varies[ci] if ci < len(pos_varies) else False
            col = nt_color(ch, ci, seq)
            # Highlight varying positions
            if varies:
                ax.add_patch(FancyBboxPatch(
                    (ci, y - 0.40), 1, 0.80,
                    boxstyle="square,pad=0",
                    linewidth=0, facecolor="#FFD600", alpha=0.5, zorder=2))
            ax.text(ci + 0.5, y, ch, ha="center", va="center",
                    fontsize=8, fontweight="bold", color=col,
                    fontfamily="monospace", zorder=3)

        # Length annotation
        ax.text(w + 0.5, y, f"{w} nt  [{dominant}]n",
                va="center", ha="left", fontsize=8,
                color=DINU_COLORS.get(dominant, DEFAULT_COL), fontweight="bold")

        y_tick_pos.append(y)
        y_tick_lab.append(SP_LABELS[sp])

    if not y_tick_pos: ax.axis("off"); return

    ax.set_xlim(-0.5, max_len + 14)
    ax.set_ylim(-0.8, len(y_tick_pos) - 0.1)
    ax.set_yticks(y_tick_pos)
    ax.set_yticklabels(y_tick_lab, fontsize=10, fontstyle="italic")
    ax.set_title(GENE_LABELS.get(gene_name, gene_name.title()),
                 fontsize=12, fontweight="bold", pad=6)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", bottom=False, labelbottom=False)
    for sp in ["top","right","bottom"]: ax.spines[sp].set_visible(False)
    ax.spines["left"].set_color("#ddd")

# ── Load data ─────────────────────────────────────────────────────────────────
rows = []
with open(BASE / "zdna_orthologs.tsv") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        rows.append(r)

by_gene = defaultdict(list)
for r in rows: by_gene[r["gene"]].append(r)
genes = list(by_gene.keys())
print(f"Genes: {genes}")

# ── Layout: vertical stack ────────────────────────────────────────────────────
n = len(genes)
fig, axes = plt.subplots(n, 1, figsize=(12, n * 2.8),
                          facecolor="white",
                          gridspec_kw={"hspace": 0.65})
if n == 1: axes = [axes]

for idx, gene in enumerate(genes):
    draw_panel(axes[idx], by_gene[gene], gene)

# ── Legend ────────────────────────────────────────────────────────────────────
dinu_patches = [
    mpatches.Patch(fc=DINU_COLORS["CG"], ec="none", label="CG repeats"),
    mpatches.Patch(fc=DINU_COLORS["CA"], ec="none", label="CA/TG repeats"),
    mpatches.Patch(fc=DINU_COLORS["AT"], ec="none", label="AT/TA repeats"),
    mpatches.Patch(fc=DINU_COLORS["AC"], ec="none", label="AC/GT repeats"),
    mpatches.Patch(fc=DEFAULT_COL,       ec="none", label="Other"),
    mpatches.Patch(fc="#FFD600", ec="none", alpha=0.8,
                   label="Position varies between species"),
]
fig.legend(handles=dinu_patches, loc="lower center", ncol=3,
           fontsize=9.5, frameon=True, fancybox=False, edgecolor="#bbb",
           bbox_to_anchor=(0.5, -0.05),
           title="Z-DNA forming dinucleotide repeats  ·  [XY]n = dominant repeat unit  ·  length in nt",
           title_fontsize=9.5)

out_png = BASE / "zdna_alignment.png"
out_svg = BASE / "zdna_alignment.svg"
fig.savefig(out_png, format="png", bbox_inches="tight", dpi=150, facecolor="white")
fig.savefig(out_svg, format="svg", bbox_inches="tight", facecolor="white")
print(f"✅ {out_png}")
print(f"✅ {out_svg}")
plt.close()
