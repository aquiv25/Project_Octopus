import re, numpy as np
import matplotlib
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["font.family"] = "Arial"
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

txt = open("meme_output/meme.txt").read()
NT_COLORS = {"A": "#2CA02C", "C": "#1565C0", "G": "#FFB300", "T": "#D62728"}
NTS = ["A", "C", "G", "T"]

# Parse Motif 1 only
sec = re.split(r'(?=MOTIF\s+\S+\s+MEME-)', txt)[1]
consensus = re.search(r'MOTIF\s+(\S+)', sec).group(1)
evalue    = re.search(r'E-value\s*=\s*([\d.e+\-]+)', sec).group(1)
sites     = re.search(r'sites\s*=\s*(\d+)', sec).group(1)
lp        = re.search(r'letter-probability matrix[^\n]*\n((?:[\s\d.]+\n)+)', sec)
matrix = []
for row in lp.group(1).strip().split("\n"):
    vals = list(map(float, row.split()))
    if len(vals) == 4:
        matrix.append(vals)

fig, ax = plt.subplots(figsize=(8, 4), facecolor="white")
ax.set_facecolor("white")

for pos, freqs in enumerate(matrix):
    ic = float(2 + sum(float(f) * np.log2(float(f) + 1e-10) for f in freqs))
    ic = max(ic, 0.0)
    order = sorted(range(4), key=lambda i: float(freqs[i]))
    y = 0.0
    for nt_idx in order:
        h = float(freqs[nt_idx]) * ic
        if h < 0.005: continue
        ax.add_patch(mpatches.FancyBboxPatch(
            (pos + 0.06, y), 0.88, h,
            boxstyle="round,pad=0.03",
            facecolor=NT_COLORS[NTS[nt_idx]], edgecolor="white",
            linewidth=0.5, zorder=2))
        if h > 0.22:
            ax.text(pos + 0.5, y + h / 2, NTS[nt_idx],
                    ha="center", va="center",
                    fontsize=max(9, int(h * 13)),
                    fontweight="bold", color="white", zorder=3)
        y += h

n = len(matrix)
ax.set_xlim(-0.3, n + 0.3)
ax.set_ylim(0, 2.25)
ax.axhline(2, color="#ccc", lw=0.8, ls="--", zorder=1)
ax.set_xticks(np.arange(n) + 0.5)
ax.set_xticklabels(range(1, n + 1), fontsize=11)
ax.set_xlabel("Position", fontsize=12)
ax.set_ylabel("Information content (bits)", fontsize=12)
ax.set_yticks([0, 0.5, 1.0, 1.5, 2.0])
for sp in ["top", "right"]: ax.spines[sp].set_visible(False)
ax.spines["left"].set_color("#aaa")
ax.spines["bottom"].set_color("#aaa")

# Title + stats box
stats_text = f"E-value = {evalue}   ·   found in {sites} / 18 sequences   ·   width = {n} nt"
ax.set_title(f"Conserved Z-DNA motif in octopus promoters\n{stats_text}",
             fontsize=14, fontweight="bold", pad=12, linespacing=1.6)

# Legend
leg = [mpatches.Patch(fc=c, ec="none", label=nt) for nt, c in NT_COLORS.items()]
ax.legend(handles=leg, loc="upper left", ncol=4, fontsize=10.5,
          frameon=True, fancybox=False, edgecolor="#ccc",
          title="Nucleotide", title_fontsize=10)

# Annotation: motif type
ax.text(n - 0.2, 2.12,
        "Purine–pyrimidine alternation (R·Y)n — Z-DNA forming",
        ha="right", va="bottom", fontsize=9.5, color="#555", style="italic")

out_png = "meme_output/motif1_logo.png"
out_svg = "meme_output/motif1_logo.svg"
fig.savefig(out_png, dpi=180, bbox_inches="tight", facecolor="white")
fig.savefig(out_svg, format="svg", bbox_inches="tight", facecolor="white")
print(f"✅ {out_png}")
plt.close()
