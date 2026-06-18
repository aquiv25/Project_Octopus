import re, numpy as np
import matplotlib
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["font.family"] = "Arial"
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

txt = open("meme_output/meme.txt").read()
NT_COLORS = {"A": "#2CA02C", "C": "#1565C0", "G": "#FFB300", "T": "#D62728"}
NTS = ["A", "C", "G", "T"]

# Parse motifs
motif_blocks = []
sections = re.split(r'(?=MOTIF\s+\S+\s+MEME-)', txt)
for sec in sections:
    m  = re.search(r'MOTIF\s+(\S+)\s+MEME-(\d+)', sec)
    w  = re.search(r'width\s*=\s*(\d+)', sec)
    s  = re.search(r'sites\s*=\s*(\d+)', sec)
    e  = re.search(r'E-value\s*=\s*([\d.e+\-]+)', sec)
    lp = re.search(r'letter-probability matrix[^\n]*\n((?:[\s\d.]+\n)+)', sec)
    if not (m and w and s and e and lp): continue
    matrix = []
    for row in lp.group(1).strip().split("\n"):
        vals = list(map(float, row.split()))
        if len(vals) == 4:
            matrix.append(vals)
    if matrix:
        motif_blocks.append({
            "consensus": m.group(1), "num": int(m.group(2)),
            "width": int(w.group(1)), "sites": int(s.group(1)),
            "evalue": e.group(1), "matrix": matrix
        })

print(f"Parsed {len(motif_blocks)} motifs")

def draw_logo(ax, mb):
    matrix = mb["matrix"]
    for pos, freqs in enumerate(matrix):
        ic = float(2 + sum(float(f) * np.log2(float(f) + 1e-10) for f in freqs))
        ic = max(ic, 0.0)
        order = sorted(range(4), key=lambda i: float(freqs[i]))
        y = 0.0
        for nt_idx in order:
            h = float(freqs[nt_idx]) * ic
            if h < 0.005: continue
            ax.add_patch(mpatches.FancyBboxPatch(
                (pos + 0.05, y), 0.9, h,
                boxstyle="square,pad=0",
                facecolor=NT_COLORS[NTS[nt_idx]], edgecolor="none", zorder=2))
            if h > 0.25:
                ax.text(pos + 0.5, y + h / 2, NTS[nt_idx],
                        ha="center", va="center",
                        fontsize=max(6, int(h * 10)), fontweight="bold",
                        color="white", zorder=3)
            y += h

    n = len(matrix)
    ax.set_xlim(0, n)
    ax.set_ylim(0, 2.15)
    ax.set_xticks(np.arange(n) + 0.5)
    ax.set_xticklabels(range(1, n + 1), fontsize=8)
    ax.set_ylabel("bits", fontsize=9)
    ev = mb["evalue"]
    ax.set_title(f"Motif {mb['num']}  [{mb['consensus']}]\nE-value = {ev}   sites = {mb['sites']}",
                 fontsize=9.5, fontweight="bold", pad=5)
    ax.axhline(2, color="#ddd", lw=0.7, ls="--")
    for sp in ["top", "right"]: ax.spines[sp].set_visible(False)

n = len(motif_blocks)
fig, axes = plt.subplots(1, n, figsize=(n * 3.8, 3.8),
                          facecolor="white",
                          gridspec_kw={"wspace": 0.45})
if n == 1: axes = [axes]

for i, mb in enumerate(motif_blocks):
    draw_logo(axes[i], mb)

fig.legend(
    handles=[mpatches.Patch(fc=c, ec="none", label=nt) for nt, c in NT_COLORS.items()],
    loc="lower center", ncol=4, fontsize=10,
    bbox_to_anchor=(0.5, -0.10), frameon=False)

fig.suptitle("MEME: conserved motifs in Z-DNA forming promoter sequences — Cephalopoda",
             fontsize=11, fontweight="bold", y=1.03)

out = "meme_output/meme_logos.png"
fig.savefig(out, dpi=150, bbox_inches="tight", facecolor="white")
print(f"✅ {out}")
plt.close()
