import pandas as pd
import numpy as np
import matplotlib
matplotlib.rcParams["font.family"] = "Arial"
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt

df = pd.read_csv("/Users/dassagaripova/Downloads/enrichment.all.tsv", sep="\t", comment="#",
                 names=["category","term_id","description","obs","bg","strength","signal","fdr","ids","labels"])

# Select categories to show
cats = {
    "GO Process":    ("Biological Process", "#1565C0"),
    "GO Function":   ("Molecular Function", "#2E7D32"),
    "GO Component":  ("Cellular Component", "#6A1B9A"),
    "KEGG":          ("KEGG Pathways",       "#E65100"),
    "Reactome":      ("Reactome",            "#AD1457"),
}

# Keep top N per category by FDR
TOP = 5
rows = []
for cat_key, (cat_label, color) in cats.items():
    sub = df[df["category"] == cat_key].copy()
    sub = sub.sort_values("fdr").head(TOP)
    sub["cat_label"] = cat_label
    sub["color"] = color
    rows.append(sub)

plot_df = pd.concat(rows, ignore_index=True)
plot_df["-log10fdr"] = -np.log10(plot_df["fdr"].astype(float))
plot_df["signal"] = plot_df["signal"].astype(float)

# Truncate long descriptions
plot_df["desc_short"] = plot_df["description"].apply(lambda x: x[:55] + "…" if len(x) > 55 else x)

# Sort: by category then by fdr within category
cat_order = [v[0] for v in cats.values()]
plot_df["cat_order"] = plot_df["cat_label"].map({v: i for i, v in enumerate(cat_order)})
plot_df = plot_df.sort_values(["cat_order", "fdr"], ascending=[True, False]).reset_index(drop=True)

fig, ax = plt.subplots(figsize=(10, len(plot_df) * 0.42 + 1.5), facecolor="white")
ax.set_facecolor("white")

y_pos = np.arange(len(plot_df))

scatter = ax.scatter(
    plot_df["-log10fdr"],
    y_pos,
    s=plot_df["obs"] * 6,
    c=plot_df["color"],
    alpha=0.85,
    edgecolors="white",
    linewidths=0.5,
    zorder=3
)

ax.set_yticks(y_pos)
ax.set_yticklabels(plot_df["desc_short"], fontsize=9.5)
ax.set_xlabel("-log10(FDR)", fontsize=11)
ax.axvline(x=-np.log10(0.05), color="#aaa", lw=0.8, ls="--", zorder=1)
ax.text(-np.log10(0.05) + 0.05, len(plot_df) - 0.5, "FDR=0.05", fontsize=8, color="#aaa")

# Category separators
cat_borders = []
prev = None
for i, row in plot_df.iterrows():
    if row["cat_label"] != prev:
        cat_borders.append(i)
        prev = row["cat_label"]

for b in cat_borders[1:]:
    ax.axhline(b - 0.5, color="#ddd", lw=0.8, zorder=1)

# Category labels on right
for b_idx, b in enumerate(cat_borders):
    end = cat_borders[b_idx + 1] if b_idx + 1 < len(cat_borders) else len(plot_df)
    mid = (b + end - 1) / 2
    cat_name = plot_df.iloc[b]["cat_label"]
    color = plot_df.iloc[b]["color"]
    ax.text(ax.get_xlim()[1] if ax.get_xlim()[1] > 1 else 10,
            mid, cat_name, fontsize=8.5, color=color,
            fontweight="bold", va="center", ha="left")

# Legend for bubble size
for size, label in [(10, "10"), (30, "30"), (60, "60")]:
    ax.scatter([], [], s=size * 6, c="#888", alpha=0.7, label=f"{label} genes")
ax.legend(title="Gene count", loc="lower right", fontsize=8.5,
          title_fontsize=9, frameon=True, edgecolor="#ccc")

for sp in ["top", "right"]:
    ax.spines[sp].set_visible(False)
ax.spines["left"].set_color("#ccc")
ax.spines["bottom"].set_color("#ccc")
ax.grid(axis="x", color="#f0f0f0", zorder=0)

plt.tight_layout()
out = "/Users/dassagaripova/Documents/Claude/Projects/Article/g4_promoter/results/enrichment_bubble_g4.png"
out_svg = out.replace(".png", ".svg")
fig.savefig(out, dpi=180, bbox_inches="tight", facecolor="white")
fig.savefig(out_svg, format="svg", bbox_inches="tight", facecolor="white")
print(f"✅ {out}")
plt.close()
