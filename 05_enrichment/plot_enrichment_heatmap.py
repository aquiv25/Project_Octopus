import pandas as pd
import numpy as np
import matplotlib
matplotlib.rcParams["font.family"] = "Arial"
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap

df = pd.read_csv("/Users/dassagaripova/Downloads/enrichment.all.tsv", sep="\t", comment="#",
                 names=["category","term_id","description","obs","bg","strength","signal","fdr","ids","labels"])

# Focus on GO Process, GO Component, KEGG — most biologically informative
cats_use = ["GO Process", "GO Component", "KEGG"]
cat_colors = {
    "GO Process":   "#1565C0",
    "GO Component": "#6A1B9A",
    "KEGG":         "#E65100",
}

TOP = 8
rows = []
for cat in cats_use:
    sub = df[df["category"] == cat].sort_values("fdr").head(TOP)
    rows.append(sub)
plot_df = pd.concat(rows, ignore_index=True)
plot_df["-log10fdr"] = -np.log10(plot_df["fdr"].astype(float))
plot_df["signal"] = plot_df["signal"].astype(float)

# Short descriptions
plot_df["desc_short"] = plot_df["description"].apply(
    lambda x: x[:50] + "…" if len(x) > 50 else x)

# Build gene × term matrix
# Collect all unique genes across selected terms
all_genes = set()
for labels in plot_df["labels"]:
    all_genes.update(str(labels).split(","))
all_genes = sorted(all_genes)

# For heatmap: value = -log10(fdr) if gene is in term, else 0
matrix = pd.DataFrame(0.0, index=plot_df["desc_short"], columns=all_genes)
for _, row in plot_df.iterrows():
    genes = str(row["labels"]).split(",")
    for g in genes:
        if g in matrix.columns:
            matrix.loc[row["desc_short"], g] = row["-log10fdr"]

# Remove genes that appear in fewer than 2 terms
matrix = matrix.loc[:, (matrix > 0).sum() >= 2]

# --- Plot ---
n_terms, n_genes = matrix.shape
fig_w = min(max(n_genes * 0.18 + 3, 12), 22)
fig_h = n_terms * 0.45 + 2.5

fig, ax = plt.subplots(figsize=(fig_w, fig_h), facecolor="white")
ax.set_facecolor("white")

# Custom colormap: white → light blue → dark blue
cmap = LinearSegmentedColormap.from_list(
    "custom", ["#ffffff", "#bbdefb", "#1565C0"], N=256)

im = ax.imshow(matrix.values, aspect="auto", cmap=cmap,
               interpolation="nearest", vmin=0)

# Axes
ax.set_xticks(np.arange(n_genes))
ax.set_xticklabels(matrix.columns, fontsize=6.5, rotation=90)
ax.set_yticks(np.arange(n_terms))
ax.set_yticklabels(matrix.index, fontsize=9)

# Colorbar
cbar = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01)
cbar.set_label("-log10(FDR)", fontsize=9)
cbar.ax.tick_params(labelsize=8)

# Category color bar on left
cat_ax = fig.add_axes([0.01, ax.get_position().y0,
                       0.012, ax.get_position().height])
cat_ax.set_xlim(0, 1)
cat_ax.set_ylim(0, n_terms)
cat_ax.axis("off")

# Draw category blocks
prev_cat = None
cat_start = 0
cat_blocks = []
for i, row in plot_df.iterrows():
    cat = row["category"]
    if cat != prev_cat:
        if prev_cat is not None:
            cat_blocks.append((prev_cat, cat_start, i))
        cat_start = i
        prev_cat = cat
cat_blocks.append((prev_cat, cat_start, len(plot_df)))

for cat, start, end in cat_blocks:
    color = cat_colors.get(cat, "#888")
    y0 = n_terms - end
    height = end - start
    cat_ax.add_patch(mpatches.Rectangle(
        (0, y0), 1, height,
        facecolor=color, edgecolor="white", linewidth=0.5))

# Separator lines between categories
for cat, start, end in cat_blocks[1:]:
    ax.axhline(end - 0.5, color="white", lw=1.5)

# Legend
legend_patches = [mpatches.Patch(fc=cat_colors[c], label=c) for c in cats_use]
ax.legend(handles=legend_patches, loc="upper right",
          bbox_to_anchor=(1.18, 1.02), fontsize=8.5,
          title="Category", title_fontsize=9,
          frameon=True, edgecolor="#ccc")

plt.tight_layout(rect=[0.03, 0, 1, 1])

out = "/Users/dassagaripova/Documents/Claude/Projects/Article/g4_promoter/results/enrichment_heatmap_g4.png"
out_svg = out.replace(".png", ".svg")
fig.savefig(out, dpi=180, bbox_inches="tight", facecolor="white")
fig.savefig(out_svg, format="svg", bbox_inches="tight", facecolor="white")
print(f"✅ {out}")
plt.close()
