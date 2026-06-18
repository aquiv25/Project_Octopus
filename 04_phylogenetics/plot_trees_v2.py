"""
Manual phylogenetic tree drawing — publication style with colored labels
"""
import numpy as np
import matplotlib
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["font.family"] = "Arial"
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from Bio import Phylo
from pathlib import Path

TREES_DIR = Path("/Users/dassagaripova/Documents/Claude/Projects/Article/trees")
OUT_DIR   = Path("/Users/dassagaripova/Documents/Claude/Projects/Article")

GENES = ["DNMT1","TET3","EZH2","KMT2A","CHD1","BRD3","ADAR1","UHRF1"]

SP_COLORS = {
    "O_bimaculoides": "#D6604D",
    "O_vulgaris":     "#4393C3",
    "O_sinensis":     "#2CA02C",
    "A_fangsiao":     "#FF7F0E",
    "E_cirrhosa":     "#9467BD",
    "A_hians":        "#8C564B",
}
SP_LABELS = {
    "O_bimaculoides": "O. bimaculoides",
    "O_vulgaris":     "O. vulgaris",
    "O_sinensis":     "O. sinensis",
    "A_fangsiao":     "A. fangsiao",
    "E_cirrhosa":     "E. cirrhosa",
    "A_hians":        "A. hians",
}

def get_sp_key(clade_name):
    """Extract species key from clade name like GENE_Species_name"""
    if clade_name and "_" in clade_name:
        parts = clade_name.split("_", 1)
        return parts[1] if len(parts) > 1 else clade_name
    return clade_name

def assign_y(tree):
    """Assign y positions to leaves (0,1,2,...) and compute internal y as mean of children."""
    leaves = tree.get_terminals()
    y_pos = {}
    for i, leaf in enumerate(leaves):
        y_pos[id(leaf)] = i
    def set_internal(clade):
        if clade.is_terminal():
            return y_pos[id(clade)]
        child_ys = [set_internal(c) for c in clade.clades]
        y_pos[id(clade)] = np.mean(child_ys)
        return y_pos[id(clade)]
    set_internal(tree.root)
    return y_pos

def assign_x(tree):
    """Assign x positions based on branch lengths from root."""
    x_pos = {}
    def set_x(clade, cumlen):
        bl = clade.branch_length if clade.branch_length else 0
        x_pos[id(clade)] = cumlen + bl
        for child in clade.clades:
            set_x(child, cumlen + bl)
    # root has no branch length
    x_pos[id(tree.root)] = 0
    for child in tree.root.clades:
        set_x(child, 0)
    return x_pos

def draw_tree(ax, tree, title):
    tree.root_at_midpoint()
    tree.ladderize()

    y_pos = assign_y(tree)
    x_pos = assign_x(tree)

    # Normalize x to [0, 1]
    max_x = max(x_pos.values()) or 1
    x_norm = {k: v / max_x for k, v in x_pos.items()}

    def draw_clade(clade):
        x = x_norm[id(clade)]
        y = y_pos[id(clade)]
        for child in clade.clades:
            xc = x_norm[id(child)]
            yc = y_pos[id(child)]
            # Horizontal line to child
            ax.plot([x, xc], [yc, yc], color="#333", lw=1.3, solid_capstyle="round")
            # Vertical connector
        if clade.clades:
            ys = [y_pos[id(c)] for c in clade.clades]
            ax.plot([x, x], [min(ys), max(ys)], color="#333", lw=1.3)
        # Bootstrap on internal nodes
        if not clade.is_terminal() and clade.confidence is not None and clade != tree.root:
            bs = int(clade.confidence)
            if bs >= 50:
                ax.text(x - 0.03, y + 0.05, str(bs),
                        fontsize=7, color="#666", ha="right", va="bottom")
        for child in clade.clades:
            draw_clade(child)

    draw_clade(tree.root)

    # Leaf labels
    max_label_x = 0
    for leaf in tree.get_terminals():
        x = x_norm[id(leaf)]
        y = y_pos[id(leaf)]
        sp_key = get_sp_key(leaf.name)
        color = SP_COLORS.get(sp_key, "#333")
        label = SP_LABELS.get(sp_key, sp_key.replace("_", " "))
        ax.text(x + 0.025, y, label, fontsize=10, color=color,
                fontweight="bold", va="center", ha="left")
        max_label_x = max(max_label_x, x)

    # Scale bar (10% of tree)
    n_leaves = len(tree.get_terminals())
    sb_y = -0.7
    ax.plot([0, 0.1], [sb_y, sb_y], color="#555", lw=1.5)
    ax.plot([0, 0], [sb_y - 0.08, sb_y + 0.08], color="#555", lw=1.0)
    ax.plot([0.1, 0.1], [sb_y - 0.08, sb_y + 0.08], color="#555", lw=1.0)
    ax.text(0.05, sb_y - 0.22, f"{max_x*0.1:.2f}", ha="center",
            fontsize=7.5, color="#555")

    ax.set_xlim(-0.08, 1.75)
    ax.set_ylim(-1.0, n_leaves - 0.2)
    ax.set_title(title, fontsize=13, fontweight="bold", color="#111", pad=5)
    ax.axis("off")

# ── Figure ────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 4, figsize=(24, 12), facecolor="white")
axes = axes.flatten()

for ax, gene in zip(axes, GENES):
    treefile = TREES_DIR / f"{gene}_tree.treefile"
    ax.set_facecolor("white")
    if not treefile.exists():
        ax.text(0.5, 0.5, f"{gene}\n(no tree)", ha="center", va="center",
                fontsize=12, color="#bbb", transform=ax.transAxes)
        ax.axis("off")
        continue
    tree = Phylo.read(str(treefile), "newick")
    draw_tree(ax, tree, gene)

# Legend
handles = [mpatches.Patch(fc=SP_COLORS[k], label=SP_LABELS[k],
                           ec="none", linewidth=0)
           for k in ["O_bimaculoides","O_vulgaris","O_sinensis","A_fangsiao","E_cirrhosa","A_hians"]]
fig.legend(handles=handles, loc="lower center", ncol=6, fontsize=11.5,
           frameon=False, handlelength=1.5, handleheight=1.0,
           bbox_to_anchor=(0.5, -0.005))

fig.suptitle(
    "Phylogenetic trees of epigenetic regulatory genes across Cephalopoda\n"
    "Maximum Likelihood · LG+G4 model · 1000 ultrafast bootstraps · MAFFT + IQ-TREE3",
    fontsize=13, fontweight="bold", y=1.015, color="#111")

plt.tight_layout(rect=[0, 0.05, 1, 1])

out_svg = OUT_DIR / "epigenetic_trees_v2.svg"
out_png = OUT_DIR / "epigenetic_trees_v2.png"
fig.savefig(out_svg, format="svg", bbox_inches="tight", facecolor="white")
fig.savefig(out_png, format="png", bbox_inches="tight", dpi=150, facecolor="white")
print(f"✅ {out_svg}")
print(f"✅ {out_png}")
plt.close()
