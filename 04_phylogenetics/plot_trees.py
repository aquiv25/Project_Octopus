"""
Visualize 8 phylogenetic trees in one figure — publication style
"""
from Bio import Phylo
from io import StringIO
import matplotlib
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["font.family"] = "Arial"
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

TREES_DIR = Path("/Users/dassagaripova/Documents/Claude/Projects/Article/trees")
OUT_DIR   = Path("/Users/dassagaripova/Documents/Claude/Projects/Article")

GENES = ["DNMT1","TET3","EZH2","KMT2A","CHD1","BRD3","ADAR1","UHRF1"]

# Species colors
SP_COLORS = {
    "O_bimaculoides": "#D6604D",
    "O_vulgaris":     "#4393C3",
    "O_sinensis":     "#2CA02C",
    "A_fangsiao":     "#FF7F0E",
    "E_cirrhosa":     "#9467BD",
    "O_americanus":   "#8C564B",
    "O_rubescens":    "#E377C2",
    "O_mimus":        "#7F7F7F",
}
SP_LABELS = {
    "O_bimaculoides": "O. bimaculoides",
    "O_vulgaris":     "O. vulgaris",
    "O_sinensis":     "O. sinensis",
    "A_fangsiao":     "A. fangsiao",
    "E_cirrhosa":     "E. cirrhosa",
    "O_americanus":   "O. americanus",
    "O_rubescens":    "O. rubescens",
    "O_mimus":        "O. mimus",
}

fig, axes = plt.subplots(2, 4, figsize=(22, 11), facecolor="white")
axes = axes.flatten()

for ax, gene in zip(axes, GENES):
    treefile = TREES_DIR / f"{gene}_tree.treefile"
    if not treefile.exists():
        ax.text(0.5, 0.5, f"{gene}\n(no tree)", ha="center", va="center",
                fontsize=12, color="#999", transform=ax.transAxes)
        ax.axis("off")
        continue

    tree = Phylo.read(str(treefile), "newick")
    tree.root_at_midpoint()
    tree.ladderize()

    # Draw tree
    Phylo.draw(tree, axes=ax, do_show=False,
               label_func=lambda c: "", show_confidence=False)

    # Color leaf labels manually
    ax.set_facecolor("white")
    for text_obj in ax.texts:
        raw = text_obj.get_text().strip()
        sp_key = raw.split("_", 1)[1] if "_" in raw else raw
        col = SP_COLORS.get(sp_key, "#333")
        label = SP_LABELS.get(sp_key, sp_key.replace("_", " "))
        text_obj.set_text(label)
        text_obj.set_color(col)
        text_obj.set_fontsize(10.5)
        text_obj.set_fontweight("bold")

    # Add bootstrap values to internal nodes
    for clade in tree.find_clades(order="level"):
        if clade.confidence is not None and not clade.is_terminal():
            # get position from axes
            pass  # Phylo.draw doesn't expose positions easily; skip for now

    # Title
    ax.set_title(gene, fontsize=14, fontweight="bold", color="#111", pad=6)
    ax.set_xlabel("")
    ax.set_ylabel("")
    # Clean up spines
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)

# Legend
handles = [mpatches.Patch(fc=SP_COLORS[k], label=SP_LABELS[k], ec="white")
           for k in SP_COLORS if SP_LABELS[k] in
           [t.get_text() for ax in axes for t in ax.texts]]
# Simpler: show all 5 species we have
used = ["O_bimaculoides","O_vulgaris","O_sinensis","A_fangsiao","E_cirrhosa"]
handles = [mpatches.Patch(fc=SP_COLORS[k], label=SP_LABELS[k], ec="white", linewidth=0)
           for k in used]
fig.legend(handles=handles, loc="lower center", ncol=5, fontsize=11,
           frameon=False, handlelength=1.4, bbox_to_anchor=(0.5, -0.01))

fig.suptitle("Phylogenetic trees of epigenetic genes — Cephalopoda\n"
             "(Maximum Likelihood, LG+G4, 1000 bootstraps | MAFFT + IQ-TREE)",
             fontsize=13, fontweight="bold", y=1.01)

plt.tight_layout(rect=[0, 0.05, 1, 1])

out_svg = OUT_DIR / "epigenetic_trees.svg"
out_png = OUT_DIR / "epigenetic_trees.png"
fig.savefig(out_svg, format="svg", bbox_inches="tight", facecolor="white")
fig.savefig(out_png, format="png", bbox_inches="tight", dpi=150, facecolor="white")
print(f"✅ {out_svg}")
print(f"✅ {out_png}")
plt.close()
